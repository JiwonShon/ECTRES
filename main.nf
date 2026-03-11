nextflow.enable.dsl=2

/*
 * Minimal ichorCNA Nextflow pipeline for low-pass WGS
 * Steps:
 *   1) HMMcopy readCounter: BAM -> WIG (bin counts)
 *   2) ichorCNA runIchorCNA.R: WIG -> correctedDepth + segments + plots
 *
 * Input samplesheet CSV (header required):
 *   aliquot_barcode,patient_barcode,bam,bai
 * where 'bai' is optional; if missing the pipeline will look for <bam>.bai
 */


workflow {

    if( !params.input ) {
        error "Missing required parameter: --input <samplesheet.csv>"
    }

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.aliquot_barcode as String,
                patient_barcode: (row.patient_barcode ?: "NA") as String
            ]
            def bam = file(row.bam as String)
            def bai = row.bai ? file(row.bai as String) : file("${row.bam}.bai")
            tuple(meta, bam, bai)
        }
        .set { samples_ch }

    readcounter_out_ch = READCOUNTER(samples_ch)
    RUN_ICHORCNA(readcounter_out_ch)
}


process READCOUNTER {
    tag { meta.id }
    publishDir "${params.outdir}/wig", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    container params.hmmcopy_container

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.wig")

    script:
    """
    set -euo pipefail

    # ichorCNA wiki: BAM index should be named <bam>.bam.bai in the same directory.
    # Nextflow stages files into workdir; ensure ${bam}.bai exists.
    if [[ ! -f "${bam}.bai" ]]; then
      ln -sf "${bai}" "${bam}.bai"
    fi

    readCounter \\
      --window ${params.window} \\
      --quality ${params.mapq} \\
      --chromosome "${params.chromosomes}" \\
      "${bam}" > "${meta.id}.wig"
    """
}

process RUN_ICHORCNA {
    tag { meta.id }
    publishDir { "${params.outdir}/ichorCNA/${meta.patient_barcode}/${meta.id}" }, mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    container params.ichorcna_container

    input:
    tuple val(meta), path(wig)

    output:
    path "${meta.id}.correctedDepth.txt"
    path "${meta.id}.params.txt"
    path "${meta.id}.seg.txt", optional: true
    path "${meta.id}.cna.seg", optional: true
    path "${meta.id}", optional: true
    // The full set of outputs will be copied into outdir via publishDir.

    script:
    def normalPanelArg = params.normalPanel ? "--normalPanel \"${params.normalPanel}\"" : ""
    """
    set -euo pipefail
     
    Rscript "${params.runIchorScript}" \\
      --id "${meta.id}" \\
      --WIG "${wig}" \\
      --ploidy "${params.ploidy}" \\
      --normal "${params.normal}" \\
      --maxCN ${params.maxCN} \\
      --gcWig "${params.gcWig}" \\
      --mapWig "${params.mapWig}" \\
      ${normalPanelArg} \\
      --chrs "${params.chrs}" \\
      --chrTrain "${params.chrTrain}" \\
      --chrNormalize "${params.chrNormalize}" \\
      --genomeStyle "${params.genomeStyle}" \\
      --includeHOMD ${params.includeHOMD} \\
      --estimateNormal ${params.estimateNormal} \\
      --estimatePloidy ${params.estimatePloidy} \\
      --estimateScPrevalence ${params.estimateScPrevalence} \\
      --scStates "${params.scStates}" \\
      --txnE ${params.txnE} \\
      --txnStrength ${params.txnStrength} \\
      --rmCentromereFlankLength ${params.rmCentromereFlankLength} \\
      --fracReadsInChrYForMale ${params.fracReadsInChrYForMale} \\
      --plotFileType "${params.plotFileType}" \\
      --plotYLim "${params.plotYLim}" \\
      --outDir .

    # Make sure the key output files exist in the workdir root
    if [[ ! -f "${meta.id}.correctedDepth.txt" && -f "${meta.id}/${meta.id}.correctedDepth.txt" ]]; then
      ln -sf "${meta.id}/${meta.id}.correctedDepth.txt" "${meta.id}.correctedDepth.txt"
    fi

    if [[ ! -f "${meta.id}.params.txt" && -f "${meta.id}/${meta.id}.params.txt" ]]; then
      ln -sf "${meta.id}/${meta.id}.params.txt" "${meta.id}.params.txt"
    fi

    if [[ ! -f "${meta.id}.seg.txt" && -f "${meta.id}/${meta.id}.seg.txt" ]]; then
      ln -sf "${meta.id}/${meta.id}.seg.txt" "${meta.id}.seg.txt"
    fi

    if [[ ! -f "${meta.id}.cna.seg" && -f "${meta.id}/${meta.id}.cna.seg" ]]; then
      ln -sf "${meta.id}/${meta.id}.cna.seg" "${meta.id}.cna.seg"
    fi
    """
}
