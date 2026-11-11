process EXTRACT_CIRCULAR_CONTIGS {
    tag "${sample_name}-${tag_label}"
    publishDir "${params.outdir}/de_novo_assembly/flye/${sample_name}/${tag_label}", mode: 'copy'

    input:
    tuple val(sample_name), val(tag_label),
          path(assembly),
          path(assembly_info)

    output:
    tuple val(sample_name), val(tag_label),
          path("${sample_name}_${tag_label}_circular_contigs.fasta"), emit: circular_fasta

    script:
    """
    circular_ids=\$(awk -F'\t' 'NR>1 && \$4=="Y" {print \$1}' ${assembly_info})

    if [ -z "\${circular_ids}" ]; then
        echo "No circular contigs found for ${sample_name} ${tag_label}" >&2
        touch ${sample_name}_${tag_label}_circular_contigs.fasta
    else
        echo "Circular contigs: \${circular_ids}" >&2
        samtools faidx ${assembly}
        echo "\${circular_ids}" | while read id; do
            samtools faidx ${assembly} "\${id}"
        done > ${sample_name}_${tag_label}_circular_contigs.fasta
        echo "Written: ${sample_name}_${tag_label}_circular_contigs.fasta" >&2
    fi
    """
}