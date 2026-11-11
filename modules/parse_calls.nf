process PARSE_CALLS {
    tag "${sample_name}"
    publishDir "${params.outdir}/consensus/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name),
          path(aa_profiles, stageAs: 'aa_profiles'),
          path(aa_beds, stageAs: 'aa_beds/*'),
          path(aa_cycles, stageAs: 'aa_cycles/*'),
          path(coral_profiles, stageAs: 'coral_profiles'),
          path(coral_beds, stageAs: 'coral_beds/*'),
          path(coral_cycles, stageAs: 'coral_cycles/*'),
          path(decoil_bed, stageAs: 'decoil_bed')

    output:
    tuple val(sample_name), path("${sample_name}_segments.tsv"), emit: segments

    script:
    """
    Rscript ${projectDir}/bin/parse_calls.R \\
        --sample                 ${sample_name}            \\
        --aa_profiles            aa_profiles               \\
        --aa_beds                aa_beds                   \\
        --aa_cycles              aa_cycles                 \\
        --coral_profiles         coral_profiles            \\
        --coral_beds             coral_beds                \\
        --coral_cycles           coral_cycles              \\
        --decoil_bed             decoil_bed                \\
        --outdir                 .
    """

    stub:
    """
    printf 'sample_id\\tchr\\tstart\\tend\\tstrand\\tsegment_order\\tcaller\\tcall_id\\trole\\n' > ${sample_name}_segments.tsv
    """
}
