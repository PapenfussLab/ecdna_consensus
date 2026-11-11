process REGIONS_TO_BED {
    tag "${sample_name}"
    publishDir "${params.outdir}/consensus/${sample_name}/region_beds", mode: 'copy'

    input:
    tuple val(sample_name), path(consensus_segments), path(regions)

    output:
    tuple val(sample_name), path("*.bed"), optional: true, emit: beds

    script:
    """
    Rscript ${projectDir}/bin/regions_to_bed.R \\
        --sample              ${sample_name}        \\
        --consensus_segments  ${consensus_segments} \\
        --regions             ${regions}            \\
        --min_callers         ${params.min_callers} \\
        --outdir              .
    """

    stub:
    """
    touch ${sample_name}_region1.bed
    """
}
