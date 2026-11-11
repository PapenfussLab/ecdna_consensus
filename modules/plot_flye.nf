process PLOT_FLYE {
    tag "${sample_name}"
    publishDir "${params.outdir}/plots/flye", mode: 'copy'

    input:
    tuple val(sample_name), path(consensus_segments), path(consensus_calls), path(decoil_bed, stageAs: 'decoil_bed'), path(decoil_links, stageAs: 'decoil_links'), path(pafs)
    path gtf
    path oncogenes

    output:
    tuple val(sample_name), path("*_flye.pdf"), optional: true, emit: plots

    script:
    def paf_list = [pafs].flatten().collect { it.name }.join(',')
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
        --flye_pafs          ${paf_list}            \\
        --flye_min_overlap   ${params.flye_min_overlap} \\
        --outdir             .
    """

    stub:
    """
    touch ${sample_name}_region1_flye.pdf
    """
}
