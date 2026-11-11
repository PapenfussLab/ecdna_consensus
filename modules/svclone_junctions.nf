process SVCLONE_JUNCTIONS {
    tag "${sample_name}"
    publishDir "${params.outdir}/consensus/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(consensus_segments), path(svclone)

    output:
    tuple val(sample_name), path("${sample_name}_junction_ccf.tsv"), emit: junctions
    tuple val(sample_name), path("${sample_name}_call_ccf.tsv"),     emit: calls

    script:
    """
    Rscript ${projectDir}/bin/svclone_junctions.R \\
        --sample             ${sample_name}         \\
        --consensus_segments ${consensus_segments}  \\
        --svclone            ${svclone}             \\
        --tol                ${params.svclone_tol}  \\
        --outdir             .
    """

    stub:
    """
    touch ${sample_name}_junction_ccf.tsv ${sample_name}_call_ccf.tsv
    """
}
