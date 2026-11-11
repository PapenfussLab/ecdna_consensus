params.samplesheet = 'assets/samplesheet.csv'
params.outdir      = 'results'
params.min_overlap = 0.5
params.nested_frac = 0.9
params.gtf         = null
params.oncogenes   = null
params.max_genes   = 10
params.min_callers = 2
params.run_flye    = false
params.flye_meta   = false
params.reference   = null
params.flye_min_overlap = 0.5
params.svclone_tol = 50

include { PARSE_CALLS }   from './modules/parse_calls.nf'
include { CLUSTER_CALLS } from './modules/cluster_calls.nf'
include { PLOT_CONSENSUS } from './modules/plot_consensus.nf'
include { SVCLONE_JUNCTIONS } from './modules/svclone_junctions.nf'
include { PLOT_SVCLONE } from './modules/plot_svclone.nf'
include { REGIONS_TO_BED } from './modules/regions_to_bed.nf'
include { EXTRACT_ECDNA } from './modules/extract_ecdna.nf'
include { FLYE }          from './modules/flye.nf'
include { EXTRACT_CIRCULAR_CONTIGS }   from './modules/extract_circular_contigs.nf'
include { MINIMAP2_ASSEMBLY_TO_REF }   from './modules/minimap2_assembly_to_ref.nf'
include { PLOT_FLYE }     from './modules/plot_flye.nf'

def requireFile(row, col) {
    def v = row[col]
    if (!v) error "Sample ${row.sample_name}: '${col}' is empty in ${params.samplesheet}"
    def f = file(v)
    if (!f.exists()) error "Sample ${row.sample_name}: ${col} not found: ${v}"
    f
}

def optionalFile(row, col) {
    def v = row[col]
    if (!v) return file("${projectDir}/assets/NO_FILE")
    def fs = [file(v)].flatten()
    if (fs.isEmpty() || !fs[0].exists()) error "Sample ${row.sample_name}: ${col} not found: ${v}"
    fs[0]
}

def optionalFiles(row, col) {
    def v = row[col]
    if (!v) return [file("${projectDir}/assets/NO_FILE")]
    def fs = [file(v)].flatten()
    if (fs.isEmpty() || !fs.every { it.exists() }) error "Sample ${row.sample_name}: ${col} matched no files: ${v}"
    fs
}

workflow {

    oncogenes_file = params.oncogenes ? file(params.oncogenes, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")

    sheet_lines = file(params.samplesheet).readLines()
    svclone_col = sheet_lines[0].split(',', -1).toList().indexOf('svclone_file')
    has_svclone = svclone_col >= 0 && sheet_lines.drop(1).any { l -> l.split(',', -1).size() > svclone_col && l.split(',', -1)[svclone_col].trim() }

    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def need_bam = params.run_flye
            def bam = need_bam ? requireFile(row, 'bam') : file("${projectDir}/assets/NO_FILE")
            def bai = file("${projectDir}/assets/NO_FILE")
            if (need_bam) {
                bai = file("${row.bam}.bai")
                if (!bai.exists()) error "Sample ${row.sample_name}: BAM index not found: ${bai}"
            }
            [
                sample_name : row.sample_name,
                bam         : bam,
                bai         : bai,
                aa_profiles    : optionalFile(row, 'aa_profiles'),
                aa_beds        : optionalFiles(row, 'aa_beds'),
                aa_cycles      : optionalFiles(row, 'aa_cycles'),
                coral_profiles : optionalFile(row, 'coral_profiles'),
                coral_beds     : optionalFiles(row, 'coral_beds'),
                coral_cycles   : optionalFiles(row, 'coral_cycles'),
                decoil_bed     : optionalFile(row, 'decoil_bed'),
                decoil_links   : optionalFile(row, 'decoil_links'),
                svclone     : optionalFile(row, 'svclone_file')
            ]
        }

    PARSE_CALLS(
        samples_ch.map { s -> [s.sample_name, s.aa_profiles, s.aa_beds, s.aa_cycles, s.coral_profiles, s.coral_beds, s.coral_cycles, s.decoil_bed] }
    )

    CLUSTER_CALLS(PARSE_CALLS.out.segments)

    CLUSTER_CALLS.out.consensus_calls.view()

    if (has_svclone) {
        SVCLONE_JUNCTIONS(
            CLUSTER_CALLS.out.consensus_segments
                .join(samples_ch.map { s -> [s.sample_name, s.svclone] })
                .filter { it[2].name != 'NO_FILE' }
        )
    }

    if (!params.gtf) error "Provide --gtf (gene annotation GTF) for PLOT_CONSENSUS"

    PLOT_CONSENSUS(
        CLUSTER_CALLS.out.consensus_segments
            .join(CLUSTER_CALLS.out.consensus_calls)
            .join(samples_ch.map { s -> [s.sample_name, s.decoil_bed, s.decoil_links] }),
        file(params.gtf, checkIfExists: true),
        oncogenes_file
    )

    if (has_svclone) {
        PLOT_SVCLONE(
            CLUSTER_CALLS.out.consensus_segments
                .join(CLUSTER_CALLS.out.consensus_calls)
                .join(samples_ch.map { s -> [s.sample_name, s.decoil_bed, s.decoil_links] })
                .join(SVCLONE_JUNCTIONS.out.junctions),
            file(params.gtf, checkIfExists: true),
            oncogenes_file
        )
    }

    if (params.run_flye) {

        REGIONS_TO_BED(
            CLUSTER_CALLS.out.consensus_segments.join(CLUSTER_CALLS.out.regions)
        )

        region_beds_ch = REGIONS_TO_BED.out.beds
            .flatMap { sample_name, beds ->
                [beds].flatten().collect { b ->
                    [sample_name, b.baseName.replaceFirst(/^.*_(region\d+)$/, '$1'), b]
                }
            }

        extract_in_ch = region_beds_ch
            .combine(samples_ch.map { s -> [s.sample_name, s.bam, s.bai] }, by: 0)
            .map { sample_name, tag_label, bed, bam, bai -> [sample_name, bam, bai, bed, tag_label] }

        EXTRACT_ECDNA(extract_in_ch)

        FLYE(EXTRACT_ECDNA.out.fastq)

        if (!params.reference) error "Provide --reference (GRCh38 fasta) for MINIMAP2_ASSEMBLY_TO_REF"

        EXTRACT_CIRCULAR_CONTIGS(FLYE.out.assembly.join(FLYE.out.info, by: [0, 1]))

        MINIMAP2_ASSEMBLY_TO_REF(
            EXTRACT_CIRCULAR_CONTIGS.out.circular_fasta,
            file(params.reference, checkIfExists: true)
        )

        PLOT_FLYE(
            CLUSTER_CALLS.out.consensus_segments
                .join(CLUSTER_CALLS.out.consensus_calls)
                .join(samples_ch.map { s -> [s.sample_name, s.decoil_bed, s.decoil_links] })
                .join(MINIMAP2_ASSEMBLY_TO_REF.out.paf.map { sample_name, tag_label, paf -> [sample_name, paf] }.groupTuple()),
            file(params.gtf, checkIfExists: true),
        oncogenes_file
        )
    }
}
