process FLYE {
    tag "${sample_name} - ${tag_label}"
    publishDir "${params.outdir}/de_novo_assembly/flye/${sample_name}/${tag_label}", mode: 'copy'

    input:
    tuple val(sample_name), val(tag_label), path(fastq)

    output:
    tuple val(sample_name), val(tag_label), path("assembly.fasta"),      emit: assembly
    tuple val(sample_name), val(tag_label), path("assembly_graph.gfa"),  emit: graph,         optional: true
    tuple val(sample_name), val(tag_label), path("assembly_graph.gv"),   emit: gv,            optional: true
    tuple val(sample_name), val(tag_label), path("assembly_info.txt"),   emit: info,          optional: true
    tuple val(sample_name), val(tag_label), path("flye.log"),            emit: log,           optional: true
    tuple val(sample_name), val(tag_label), path("params.json"),         emit: params_json,   optional: true

    script:
    def meta = params.flye_meta ? '--meta' : ''
    """
    flye --nano-raw ${fastq} \
        --out-dir . \
        --threads ${task.cpus} \
        --keep-haplotypes \
        ${meta}
    """
}