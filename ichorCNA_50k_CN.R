


# /Users/jiwon/Desktop/research/ectres/data/ichorCNA_seg/ECTRES_ichorCNA_cna_seg_20260331.csv

library(dplyr)
library(ggplot2)
library(readr)

# =========================
# 1. 데이터 읽기
# =========================
df <- read_csv("/Users/jiwon/Desktop/research/ectres/data/ichorCNA_seg/ECTRES_ichorCNA_cna_seg_20260331.csv", show_col_types = FALSE)

# =========================
# 2. 기본 설정
# =========================
chrom_order <- c(as.character(1:22), "X", "Y")

# hg19 chromosome lengths
chr_lengths <- c(
  "1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276,
  "5"=180915260, "6"=171115067, "7"=159138663, "8"=146364022,
  "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895,
  "13"=115169878, "14"=107349540, "15"=102531392, "16"=90354753,
  "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520,
  "21"=48129895, "22"=51304566, "X"=155270560, "Y"=59373566
)

# =========================
# 3. 데이터 정리
# =========================
df2 <- df %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end),
    logR_Copy_Number = as.numeric(logR_Copy_Number),
    sample_id = as.character(sample_id),
    source_barcode = as.character(source_barcode)
  ) %>%
  filter(
    chr %in% chrom_order,
    !is.na(start),
    !is.na(end),
    !is.na(logR_Copy_Number)
  ) %>%
  mutate(
    chr = factor(chr, levels = chrom_order)
  )

# =========================
# 4. chromosome cumulative offset
# =========================
offset_df <- data.frame(
  chr = chrom_order,
  chr_len = as.numeric(chr_lengths[chrom_order])
) %>%
  mutate(offset = c(0, cumsum(chr_len)[-length(chr_len)]))

# =========================
# 5. genome-wide x축 좌표 만들기
# =========================
df2 <- df2 %>%
  left_join(offset_df, by = "chr") %>%
  mutate(
    genome_x = offset + (start + end) / 2
  )

# =========================
# 6. x축용 chr label / boundary
# =========================
chr_label_df <- offset_df %>%
  mutate(
    chr_start = offset,
    chr_end = offset + chr_len,
    chr_mid = (chr_start + chr_end) / 2
  )

chr_boundary_df <- offset_df %>%
  filter(chr != "1") %>%
  mutate(boundary = offset)

# =========================
# 7. sample 순서 정리
# =========================
sample_order <- df2 %>%
  distinct(source_barcode, sample_id) %>%
  arrange(source_barcode, sample_id) %>%
  pull(sample_id)

df2$sample_id <- factor(df2$sample_id, levels = rev(unique(sample_order)))

# =========================
# 8. heatmap 그리기
# =========================
p <- ggplot(df2, aes(x = genome_x, y = sample_id, fill = logR_Copy_Number)) +
  geom_tile(height = 0.9) +
  geom_vline(
    data = chr_boundary_df,
    aes(xintercept = boundary),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.25
  ) +
  facet_wrap(~source_barcode, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = chr_label_df$chr_mid,
    labels = chr_label_df$chr,
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 4),   # 👈 중요 (outlier clipping)
    name = "Corrected\nCopy Number"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 11),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  ) +
  labs(
    title = "Genome-wide Copy Number Heatmap by Source Barcode (hg19)",
    x = "Chromosome",
    y = "Sample ID"
  )

print(p)