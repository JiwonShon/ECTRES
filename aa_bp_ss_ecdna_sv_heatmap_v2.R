# ecdna_sv_heatmap.R
#
# R 콘솔 / RStudio에서 위에서부터 그냥 실행하는 스크립트 (Rscript CLI 아님).
# 아래 "설정" 부분만 고쳐서 source() 하거나 한 줄씩 돌리면 됨.
#
# sample_mapping csv에서 sample_id == "parental"인 row가 parental 샘플임.
#
# 1. ecDNA_SVABA_filtered_combined.csv (annotate_svaba_vcf.py + combine_filtered_vcf_to_csv.py 산출물)
# 2. sample_mapping csv (aliquot_barcode, source_barcode, sample_id)
# 를 불러와서:
#   - source_barcode 하나(예: H2170)로 그룹핑
#   - BND mate pair(ID "숫자:1"/"숫자:2")를 SV junction 단위로 합침
#   - 각 샘플이 독립적으로 AA를 돌렸다고 가정, parental 좌표를 기준으로
#     각 clone의 junction을 window(bp) 내 최근접 매칭 (orientation 무관)
#   - VAF = AD / (AD + DP)
#   - AA_CYCLECLASS에 "ecDNA-like"가 있으면 그 SV는 ecDNA-defining으로 표시(red box)
#   - presence/absence clustermap + VAF annotated heatmap 두 장 저장
# ECGI1
# EFM19
# H2170
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(data.table)
})

# =============================== 설정 ======================================
combined_csv    <- "/mnt/NAS3/home/jiwon/ECTRES/summary/aa_bp_seq_seek/ecDNA_SVABA_filtered_combined.csv"
mapping_csv     <- "/mnt/NAS3/home/jiwon/ECTRES/manifest/sample_mapping_20260507.csv"
source_barcode  <- "EFM19"     # 이 source(환자)에 속한 샘플들만
window_bp       <- 500         # clone 자체 breakpoint를 parental 것과 같은 SV로 볼 최대 거리(bp)
out_dir         <- "/mnt/NAS3/home/jiwon/ECTRES/summary/aa_bp_seq_seek/plots/"
out_prefix      <- "H2170_ecDNA"

# gene 매핑 (BED 형식: chrom, start, end, gene_name, 헤더 없음, tab-separated).
# 없으면 NULL로 두면 gene 주석 없이 그대로 진행됨.
# gene 매핑: BED(chrom,start,end,gene_name, 헤더 없음) 또는 GTF 파일 다 지원 (확장자로 자동 판별).
# CGC(Cancer Gene Census) gtf처럼 암 관련 유전자만 있는 파일이면 결과 해석에 더 유용함.
# 없으면 NULL/빈 문자열로 두면 gene 주석 없이 그대로 진행됨.
gene_annotation_path <- "/mnt/NAS3/home/jiwon/ECTRES/data/CGC.gtf"
gene_flank_bp   <- 10000       # breakpoint 좌우로 이만큼(bp) 안에 있으면 "근처 gene"으로 인정
# ============================================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(suffix) file.path(out_dir, paste0(out_prefix, suffix))

# ---------------------------------------------------------------------------
# 1. load + restrict to one source_barcode, parental 자동 판별 (sample_id == "parental")
# ---------------------------------------------------------------------------
combined <- read.csv(combined_csv, stringsAsFactors = FALSE)
mapping  <- read.csv(mapping_csv, stringsAsFactors = FALSE)

mapping_sub <- mapping %>% filter(source_barcode == !!source_barcode)
if (nrow(mapping_sub) == 0) {
  stop(sprintf("source_barcode '%s' not found in mapping csv", source_barcode))
}

n_parental <- sum(mapping_sub$sample_id == "parental")
if (n_parental != 1) {
  stop(sprintf("source_barcode '%s' should have exactly 1 sample_id=='parental' row, found %d",
               source_barcode, n_parental))
}

dat <- combined %>%
  filter(aliquot_barcode %in% mapping_sub$aliquot_barcode) %>%
  left_join(mapping_sub %>% select(aliquot_barcode, sample_id), by = "aliquot_barcode") %>%
  mutate(AD = as.numeric(AD), DP = as.numeric(DP))

if (nrow(dat) == 0) {
  stop("No rows in combined csv matched this source_barcode's aliquot_barcodes -- check aliquot_barcode spelling / that annotate_svaba_vcf.py + combine_filtered_vcf_to_csv.py were run for these samples.")
}

# ---------------------------------------------------------------------------
# 2. BND mate pair(ID = "<base_id>:1" / "<base_id>:2")를 SV junction 단위로 합침
# ---------------------------------------------------------------------------
dat <- dat %>%
  mutate(base_id = sub(":[12]$", "", ID),
         mate_num = sub("^.*:", "", ID))

junctions <- dat %>%
  group_by(sample_id, base_id) %>%
  filter(n() == 2) %>%           # 짝 안 맞는 BND row는 제외
  summarise(
    chrom1 = CHROM[mate_num == "1"][1],
    pos1   = POS[mate_num == "1"][1],
    chrom2 = CHROM[mate_num == "2"][1],
    pos2   = POS[mate_num == "2"][1],
    AD = AD[1], DP = DP[1],
    FILTER = FILTER[1],
    AA_CYCLECLASS = AA_CYCLECLASS[1],
    AA_AMPLICON = AA_AMPLICON[1],
    AA_CYCLE = AA_CYCLE[1],
    .groups = "drop"
  ) %>%
  mutate(
    VAF = AD / (AD + DP),
    is_ecdna = grepl("ecDNA-like", AA_CYCLECLASS, fixed = TRUE)
  )

# ---------------------------------------------------------------------------
# 3. breakpoint reconciliation
#    reference_mode == "union": 전체 샘플(parental+모든 clone)의 junction을 다 모아서
#      서로 window_bp 이내면 같은 SV로 묶음 (union-find). parental에 없던 clone-only
#      SV도 행(row)으로 잡히고, 그 경우 parental 칸은 자연히 NA가 됨.
#    reference_mode == "parental": (예전 방식) parental 좌표만 기준 row로 삼고,
#      parental에 없는 clone-only SV는 제외.
# ---------------------------------------------------------------------------
reference_mode <- "union"   # "union" 또는 "parental"

# 두 breakend를 chrom:pos 문자열로 정렬해 orientation을 고정 (같은 SV면 항상 같은 순서로 저장되게)
canon_orient <- function(chrom1, pos1, chrom2, pos2) {
  a <- sprintf("%s:%015d", chrom1, pos1)
  b <- sprintf("%s:%015d", chrom2, pos2)
  swap <- a > b
  list(
    chrom1 = ifelse(swap, chrom2, chrom1), pos1 = ifelse(swap, pos2, pos1),
    chrom2 = ifelse(swap, chrom1, chrom2), pos2 = ifelse(swap, pos1, pos2)
  )
}
co <- canon_orient(junctions$chrom1, junctions$pos1, junctions$chrom2, junctions$pos2)
junctions$chrom1 <- co$chrom1; junctions$pos1 <- co$pos1
junctions$chrom2 <- co$chrom2; junctions$pos2 <- co$pos2

if (reference_mode == "union") {
  N <- nrow(junctions)
  uf <- seq_len(N)
  find_root <- function(x) { while (uf[x] != x) { uf[x] <<- uf[uf[x]]; x <- uf[x] }; x }
  union_set <- function(a, b) { ra <- find_root(a); rb <- find_root(b); if (ra != rb) uf[ra] <<- rb }

  chrom_key <- paste(junctions$chrom1, junctions$chrom2, sep = "__")
  for (k in unique(chrom_key)) {
    idx <- which(chrom_key == k)
    if (length(idx) < 2) next
    for (a in seq_len(length(idx) - 1)) {
      for (b in (a + 1):length(idx)) {
        i <- idx[a]; j <- idx[b]
        d <- abs(junctions$pos1[i] - junctions$pos1[j]) + abs(junctions$pos2[i] - junctions$pos2[j])
        if (d <= window_bp) union_set(i, j)
      }
    }
  }
  junctions$ref_junction_id <- paste0("cluster_", sapply(seq_len(N), find_root))

  rep_tbl <- junctions %>%
    group_by(ref_junction_id) %>%
    summarise(
      has_parental = any(sample_id == "parental"),
      rep_chrom1 = if (any(sample_id == "parental")) chrom1[sample_id == "parental"][1] else chrom1[which.min(pos1)],
      rep_pos1   = if (any(sample_id == "parental")) pos1[sample_id == "parental"][1]   else pos1[which.min(pos1)],
      rep_chrom2 = if (any(sample_id == "parental")) chrom2[sample_id == "parental"][1] else chrom2[which.min(pos1)],
      rep_pos2   = if (any(sample_id == "parental")) pos2[sample_id == "parental"][1]   else pos2[which.min(pos1)],
      .groups = "drop"
    )

  n_novel <- sum(!rep_tbl$has_parental)
  if (n_novel > 0) {
    message(sprintf("[note] parental에는 없고 clone에만 있는 breakpoint %d개가 union 기준으로 포함됨 (parental 칸은 NA로 표시됨).",
                     n_novel))
  }

} else {  # reference_mode == "parental" (예전 방식)
  ref <- junctions %>% filter(sample_id == "parental")
  ref <- ref %>% mutate(ref_junction_id = sprintf("%s:%s-%s:%s", chrom1, pos1, chrom2, pos2))

  match_to_ref <- function(row, ref, window) {
    cand <- ref[ref$chrom1 == row$chrom1 & ref$chrom2 == row$chrom2, ]
    if (nrow(cand) == 0) return(NA_character_)
    d <- abs(cand$pos1 - row$pos1) + abs(cand$pos2 - row$pos2)
    best <- which.min(d)
    if (d[best] <= window) cand$ref_junction_id[best] else NA_character_
  }

  junctions$ref_junction_id <- NA_character_
  junctions$ref_junction_id[junctions$sample_id == "parental"] <- ref$ref_junction_id
  for (i in which(junctions$sample_id != "parental")) {
    junctions$ref_junction_id[i] <- match_to_ref(junctions[i, ], ref, window_bp)
  }

  n_novel <- sum(is.na(junctions$ref_junction_id) & junctions$sample_id != "parental")
  if (n_novel > 0) {
    message(sprintf("[note] %d개 clone junction이 parental 기준 %dbp 이내에서 매칭 안 됨 -- 그림에서 제외.",
                     n_novel, window_bp))
  }

  rep_tbl <- ref %>%
    transmute(ref_junction_id, has_parental = TRUE,
              rep_chrom1 = chrom1, rep_pos1 = pos1, rep_chrom2 = chrom2, rep_pos2 = pos2)
}

# ---------------------------------------------------------------------------
# 3b. intra/inter-chromosomal 분류 + gene 주석 (선택) + 정렬
#     정렬 순서: intra-chromosomal 먼저, 그 안에서 chrom(1..22,X,Y) -> pos 순.
#     이후 inter-chromosomal, chrom1 -> pos1 순.
# ---------------------------------------------------------------------------
rep_tbl <- rep_tbl %>%
  mutate(sv_class = ifelse(rep_chrom1 == rep_chrom2, "intra-chromosomal", "inter-chromosomal"))

chrom_levels <- c(as.character(1:22), "X", "Y")
chrom_rank <- function(chrom) {
  r <- match(as.character(chrom), chrom_levels)
  ifelse(is.na(r), length(chrom_levels) + 1L, r)
}

annotate_genes <- function(chrom_vec, pos_vec, gene_dt, flank) {
  if (is.null(gene_dt) || length(chrom_vec) == 0) return(rep(NA_character_, length(chrom_vec)))
  q <- data.table(chrom = as.character(chrom_vec),
                   start = pmax(0, pos_vec - flank), end = pos_vec + flank,
                   qid = seq_along(chrom_vec))
  setkey(q, chrom, start, end)
  hits <- foverlaps(q, gene_dt, type = "any", nomatch = NULL)
  if (nrow(hits) == 0) return(rep(NA_character_, length(chrom_vec)))
  agg <- hits[, .(genes = paste(sort(unique(gene_name)), collapse = ",")), by = qid]
  out <- rep(NA_character_, length(chrom_vec))
  out[agg$qid] <- agg$genes
  out
}

load_gene_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("gtf", "gff", "gff3")) {
    gtf <- fread(path, header = FALSE, sep = "\t", quote = "", skip = "#",
                 col.names = c("chrom", "source", "feature", "start", "end",
                               "score", "strand", "frame", "attribute"))
    gtf <- gtf[feature == "gene"]
    if (nrow(gtf) == 0) {
      # 일부 gtf는 "gene" feature가 없고 transcript/exon만 있을 수 있음 -> 그 경우 전체에서 gene_name만 뽑아 first-occurrence로 축약
      gtf <- fread(path, header = FALSE, sep = "\t", quote = "", skip = "#",
                   col.names = c("chrom", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attribute"))
    }
    gtf[, gene_name := sub('.*gene_name[ =]"?([^";]+)"?.*', "\\1", attribute)]
    gtf <- gtf[!is.na(gene_name) & gene_name != attribute]
    gtf[, chrom := sub("^chr", "", chrom)]
    out <- gtf[, .(start = min(start), end = max(end)), by = .(chrom, gene_name)][, .(chrom, start, end, gene_name)]
  } else {
    out <- fread(path, header = FALSE, col.names = c("chrom", "start", "end", "gene_name"))
    out[, chrom := sub("^chr", "", as.character(chrom))]
  }
  setkey(out, chrom, start, end)
  out
}

gene_dt <- NULL
if (!is.null(gene_annotation_path) && nzchar(gene_annotation_path)) {
  if (file.exists(gene_annotation_path)) {
    gene_dt <- load_gene_table(gene_annotation_path)
    message(sprintf("[info] gene 주석 파일에서 %d개 gene 로드: %s", nrow(gene_dt), gene_annotation_path))
  } else {
    message(sprintf("[warn] gene_annotation_path '%s' 파일이 없어서 gene 주석 없이 진행합니다.", gene_annotation_path))
  }
}

rep_tbl$gene1 <- annotate_genes(rep_tbl$rep_chrom1, rep_tbl$rep_pos1, gene_dt, gene_flank_bp)
rep_tbl$gene2 <- annotate_genes(rep_tbl$rep_chrom2, rep_tbl$rep_pos2, gene_dt, gene_flank_bp)
rep_tbl$gene_label <- mapply(function(g1, g2) {
  g <- unique(stats::na.omit(c(g1, g2)))
  if (length(g) == 0) "" else paste0(" [", paste(g, collapse = "/"), "]")
}, rep_tbl$gene1, rep_tbl$gene2)

rep_tbl <- rep_tbl %>%
  mutate(
    bp_label = sprintf("%s:%s-%s:%s%s", rep_chrom1, format(rep_pos1, big.mark = ","),
                        rep_chrom2, format(rep_pos2, big.mark = ","), gene_label),
    chrom_rank1 = chrom_rank(rep_chrom1)
  ) %>%
  arrange(factor(sv_class, levels = c("inter-chromosomal", "intra-chromosomal")), chrom_rank1, rep_pos1)

matched <- junctions %>%
  filter(!is.na(ref_junction_id)) %>%
  left_join(rep_tbl %>% select(ref_junction_id, bp_label), by = "ref_junction_id")
bp_labels <- rep_tbl$bp_label

# ---------------------------------------------------------------------------
# 4. matrix 구성: VAF, ecDNA flag (row=bp_label, col=sample)
# ---------------------------------------------------------------------------
samples <- unique(mapping_sub$sample_id)
samples <- c("parental", setdiff(samples, "parental"))

to_matrix <- function(value_col) {
  m <- matched %>%
    select(bp_label, sample_id, all_of(value_col)) %>%
    pivot_wider(names_from = sample_id, values_from = all_of(value_col)) %>%
    as.data.frame()
  rownames(m) <- m$bp_label
  m$bp_label <- NULL
  as.matrix(m[bp_labels, intersect(samples, colnames(m)), drop = FALSE])
}

vaf_mat   <- to_matrix("VAF")
ecdna_mat <- to_matrix("is_ecdna"); ecdna_mat[is.na(ecdna_mat)] <- FALSE


row_sv_class <- setNames(rep_tbl$sv_class, rep_tbl$bp_label)[bp_labels]  # Heatmap row_split용
row_sv_class <- factor(row_sv_class, levels = c("inter-chromosomal", "intra-chromosomal"))

# ---------------------------------------------------------------------------
# 5. presence/absence clustermap (clone만)
# ---------------------------------------------------------------------------
clones <- setdiff(colnames(vaf_mat), "parental")
if (length(clones) >= 2 && nrow(vaf_mat) >= 2) {
  pres_abs <- ifelse(is.na(vaf_mat[, clones, drop = FALSE]), 0, 1)
  d <- dist(t(pres_abs), method = "binary")
  hc <- hclust(d, method = "average")
  clone_order <- clones[hc$order]

  ht <- Heatmap(pres_abs,
                name = "detected",
                col = c("0" = "#e8e8e8", "1" = "#2b6cb0"),
                cluster_rows = FALSE,
                row_split = row_sv_class,
                cluster_row_slices = FALSE,
                row_title_rot = 0,
                cluster_columns = hc,
                column_title = sprintf("%s: presence/absence of ecDNA-linked SV breakpoints across clones\n(clustered by Jaccard-style distance; rows grouped intra- vs inter-chromosomal)", source_barcode),
                row_names_side = "left",
                rect_gp = gpar(col = "white", lwd = 0.5))

  draw(ht)   # <- 현재 R plot 창(RStudio Plots pane 등)에 바로 표시

  png(out_path("_clustermap.png"),
      width = max(800, 90 * length(clones)), height = max(500, 45 * nrow(pres_abs)), res = 120)
  draw(ht)
  dev.off()
  message(sprintf("[done] %s", out_path("_clustermap.png")))
} else {
  clone_order <- clones
  message("[warn] 클론 또는 breakpoint가 2개 미만이라 클러스터링 스킵")
}

# ---------------------------------------------------------------------------
# 6. VAF annotated heatmap: parental + clustered clone 순서, red box = ecDNA-defining SV
# ---------------------------------------------------------------------------
sample_order <- c("parental", clone_order)
vaf_mat   <- vaf_mat[, sample_order, drop = FALSE]
ecdna_mat <- ecdna_mat[, sample_order, drop = FALSE]


col_fun <- colorRamp2(c(0, 1), c("#e8e8e8", "#2c5282"))

cell_fun <- function(j, i, x, y, width, height, fill) {
  v <- vaf_mat[i, j]
  if (is.na(v)) {
    grid.rect(x, y, width, height, gp = gpar(fill = "#f0f0f0", col = "white"))
    grid.text("NA", x, y, gp = gpar(fontsize = 8, col = "grey40"))
  } else {
    grid.rect(x, y, width, height, gp = gpar(fill = col_fun(v), col = "white"))
    txt_col <- if (v > 0.6) "white" else "black"
    label <- sprintf("%.2f", v)
    grid.text(label, x, y, gp = gpar(fontsize = 8, col = txt_col))
  }
  if (isTRUE(ecdna_mat[i, j])) {
    grid.rect(x, y, width, height, gp = gpar(fill = NA, col = "red", lwd = 2.2))
  }
}

ht2 <- Heatmap(vaf_mat,
               name = "VAF",
               col = col_fun,
               na_col = "#f0f0f0",
               cluster_rows = FALSE,
               row_split = row_sv_class,
               cluster_row_slices = FALSE,
               row_title_rot = 0,
               cluster_columns = FALSE,
               cell_fun = cell_fun,
               column_title = sprintf(
                 "%s-linked SV breakpoints: parental vs clones\n(cell = VAF = AD/(AD+DP); red box = ecDNA-defining SV (AA_CYCLECLASS=ecDNA-like); rows: inter-chromosomal on top, intra-chromosomal below)",
                 source_barcode),
               row_names_side = "left",
               column_names_rot = 45,
               rect_gp = gpar(col = "white", lwd = 0.5),
               show_heatmap_legend = TRUE)

draw(ht2)   # <- 현재 R plot 창에 바로 표시

png(out_path("_vaf_heatmap.png"),
    width = max(800, 90 * length(sample_order)), height = max(500, 45 * nrow(vaf_mat)), res = 120)
draw(ht2)
dev.off()
message(sprintf("[done] %s", out_path("_vaf_heatmap.png")))

bp_desc <- if (reference_mode == "union") "전체 샘플 union 기준" else "parental 기준"
novel_desc <- if (reference_mode == "union") {
  sprintf("parental엔 없는 clone-only breakpoint: %d개 (포함됨, parental=NA)", n_novel)
} else {
  sprintf("매칭 안 돼서 제외된 clone-only junction: %d개", n_novel)
}
message(sprintf("Samples: %d (parental + %d clones) | Breakpoints (%s): %d | %s",
                 length(sample_order), length(clone_order), bp_desc, nrow(vaf_mat), novel_desc))
