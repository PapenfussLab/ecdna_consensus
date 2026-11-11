process CLUSTER_CALLS {
    tag "${sample_name}"
    publishDir "${params.outdir}/consensus/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(segments)

    output:
    tuple val(sample_name), path("${sample_name}_consensus_calls.tsv"),    emit: consensus_calls
    tuple val(sample_name), path("${sample_name}_consensus_segments.tsv"), emit: consensus_segments
    tuple val(sample_name), path("${sample_name}_call_overlaps.tsv"),      emit: call_overlaps
    tuple val(sample_name), path("${sample_name}_regions.tsv"),            emit: regions

    script:
    """
    Rscript ${projectDir}/bin/cluster_calls.R \\
        --sample      ${sample_name}        \\
        --segments    ${segments}           \\
        --min_overlap ${params.min_overlap} \\
        --nested_frac ${params.nested_frac} \\
        --outdir      .
    """

    stub:
    """
    touch ${sample_name}_consensus_calls.tsv ${sample_name}_consensus_segments.tsv ${sample_name}_call_overlaps.tsv ${sample_name}_regions.tsv
    """
}
