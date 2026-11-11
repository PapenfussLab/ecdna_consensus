process MINIMAP2_ASSEMBLY_TO_REF {
    tag "${sample_name}-${tag_label}"
    publishDir "${params.outdir}/de_novo_assembly/flye/${sample_name}/${tag_label}", mode: 'copy'

    input:
    tuple val(sample_name), val(tag_label),
          path(circular_fasta)
    path(reference)

    output:
    tuple val(sample_name), val(tag_label),
          path("${sample_name}_${tag_label}_circular_to_ref.paf"), emit: paf

    script:
    """
    if [ ! -s ${circular_fasta} ]; then
        echo "No circular contigs to align for ${sample_name} ${tag_label}" >&2
        touch ${sample_name}_${tag_label}_circular_to_ref.paf
    else
        minimap2 \
            -x asm5 \
            --secondary=no \
            --cs \
            -t ${task.cpus} \
            ${reference} \
            ${circular_fasta} \
            > ${sample_name}_${tag_label}_circular_to_ref.paf
        echo "Alignments: \$(wc -l < ${sample_name}_${tag_label}_circular_to_ref.paf)" >&2
    fi
    """
}