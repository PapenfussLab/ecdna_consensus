process EXTRACT_ECDNA {
    tag "${sample_name} - ${tag_label}"
    publishDir "${params.outdir}/de_novo_assembly/${sample_name}/${tag_label}", mode: 'copy'
    input:
    tuple val(sample_name), path(bam), path(bai), path(bed), val(tag_label)
    output:
    tuple val(sample_name), val(tag_label), path("${sample_name}_${tag_label}.ecdna.bam"),     emit: bam,   optional: true
    tuple val(sample_name), val(tag_label), path("${sample_name}_${tag_label}.ecdna.bam.bai"), emit: bai,   optional: true
    tuple val(sample_name), val(tag_label), path("${sample_name}_${tag_label}.ecdna.fastq"),   emit: fastq, optional: true
    script:
    """
    if [ ! -s ${bed} ]; then
        echo "No regions in bed file for ${tag_label}, skipping extraction"
    else
        samtools view -b -L ${bed} -@ ${task.cpus} ${bam} \\
            | samtools sort -@ ${task.cpus} -o ${sample_name}_${tag_label}.ecdna.bam
        samtools index ${sample_name}_${tag_label}.ecdna.bam
        samtools fastq -@ ${task.cpus} ${sample_name}_${tag_label}.ecdna.bam \\
            > ${sample_name}_${tag_label}.ecdna.fastq
    fi
    """
}