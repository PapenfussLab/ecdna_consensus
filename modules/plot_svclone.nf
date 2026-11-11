process PLOT_SVCLONE {
    tag "${sample_name}"
    publishDir "${params.outdir}/plots/svclone", mode: 'copy'

    input:
    tuple val(sample_name), path(consensus_segments), path(consensus_calls), path(decoil_bed, stageAs: 'decoil_bed'), path(decoil_links, stageAs: 'decoil_links'), path(junctions)
    path gtf
    path oncogenes

    output:
    tuple val(sample_name), path("*_svclone.pdf"), optional: true, emit: plots

    script:
    """
    Rscript ${projectDir}/bin/plot_consensus.R \\
        --sample             ${sample_name}         \\
        --consensus_segments ${consensus_segments}  \\
        --consensus_calls    ${consensus_calls}     \\
        --decoil_bed         decoil_bed             \\
        --decoil_links       decoil_links           \\
        --gtf                ${gtf}                 \\
        --oncogenes          ${oncogenes}           \\
        --max_genes          ${params.max_genes}    \\
        --svclone_junctions  ${junctions}           \\
        --outdir             .
    """

    stub:
    """
    touch ${sample_name}_region1_svclone.pdf
    """
}
