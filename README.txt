ecdna_consensus builds a consensus set of extrachromosomal DNA (ecDNA) regions from existing CoRAL, Decoil and AmpliconArchitect/AmpliconClassifier outputs. It optionally adds SVclone junction clonality and a local de novo assembly of long reads from each consensus region.

## Requirements

- Nextflow and micromamba (or conda).
- Build it with:

  ```bash
  bash envs/setup_env.sh /path/to/micromamba_env/ecdna_consensus
  ```

  Then set `conda_env` in `nextflow.config` to that path.

## Samplesheet

CSV with a header row. Every column except `sample_name` is optional. Leave a cell empty if that input is missing for the sample.

| Column | Required | Description |
|---|---|---|
| sample_name | yes | Sample name |
| aa_profiles | no | AmpliconClassifier amplicon classification profiles table |
| aa_beds | no | AmpliconClassifier feature interval BEDs  |
| aa_cycles | no | AmpliconClassifier annotated cycles files (AmpliconSuite) |
| coral_profiles | no | AmpliconClassifier 2.0 amplicon classification profiles table (CoRAL) |
| coral_beds | no | AmpliconClassifier 2.0 ecDNA and FAN feature interval BEDs (CoRAL)|
| coral_cycles | no | AmpliconClassifier 2.0 annotated cycles files (CoRAL)  |
| decoil_bed | no | Decoil `reconstruct.ecDNA.bed` |
| decoil_links | no | Decoil `reconstruct.links.txt`. Without it the Decoil edges are not drawn |
| svclone_file | no | SVclone output. Enables the SVclone junction step and plot for the sample |
| bam | no | Nanopore BAM. Required only with `--run_flye true`. Needs `<bam>.bai` |



```csv
sample_name,aa_profiles,aa_beds,aa_cycles,coral_profiles,coral_beds,coral_cycles,decoil_bed,decoil_links,svclone_file,bam
S1,/p/aa/S1_amplicon_classification_profiles.tsv,/p/aa/S1_amplicon*_ecDNA_*_intervals.bed,/p/aa/S1_amplicon*_annotated_cycles.txt,/p/coral/S1/*_amplicon_classification_profiles.tsv,/p/coral/S1/*_classification_bed_files/*_intervals.bed,/p/coral/S1/*_annotated_cycles_files/*_annotated_cycles.txt,/p/decoil/S1/reconstruct.ecDNA.bed,/p/decoil/S1/reconstruct.links.txt,,
```

## Usage

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --outdir results \
    --gtf /path/to/genes.gtf \
    --oncogenes /path/to/oncogenes.txt
```

With the Flye branch:

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --outdir results \
    --gtf /path/to/genes.gtf \
    --oncogenes /path/to/oncogenes.txt \
    --run_flye true \
    --reference /path/to/GRCh38.fa
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| samplesheet | assets/samplesheet.csv | Input samplesheet |
| outdir | results | Output directory |
| gtf | null | Gene annotation GTF (required) |
| oncogenes | null | Oncogene list. Used to colour genes and to filter crowded windows. One gene per line. |
| max_genes | 10 | Windows with more genes than this show only genes on the oncogene list |
| min_overlap | 0.5 | Minimum overlap threshold used for clustering calls into regions |
| nested_frac | 0.9 | Minimum fraction of a call that must lie within a larger call to be classified as nested  |
| svclone_tol | 50 | Maximum breakpoint distance (bp) to match a junction to an SVclone SV |
| run_flye | false | Run the de novo assembly branch |
| min_callers | 2 | Minimum callers supporting a region for it to be assembled |
| flye_meta | false | Run Flye in `--meta` mode |
| reference | null | GRCh38 fasta (required with run_flye) |
| flye_min_overlap | 0.5 | Minimum fraction of a contig's aligned bases inside the region for the contig to be plotted |
| conda_env | set in nextflow.config | Path to the conda environment used by all processes |

## Outputs

Under `<outdir>/`:

| File | Content |
|---|---|
| `consensus/<sample>/<sample>_segments.tsv` | Segments of every call, with caller and role (footprint, cycle or call) |
| `consensus/<sample>/<sample>_consensus_calls.tsv` | Each call with its region and nesting |
| `consensus/<sample>/<sample>_consensus_segments.tsv` | Segments with their region |
| `consensus/<sample>/<sample>_regions.tsv` | Per region: number of calls, callers, chromosomes, nested regions |
| `consensus/<sample>/<sample>_call_overlaps.tsv` | Overlapping pairs of calls: sizes, overlap, overlap score and Jaccard |
| `consensus/<sample>/<sample>_junction_ccf.tsv` | Per junction: matched SVclone event, ccf, cluster, clonal |
| `consensus/<sample>/<sample>_call_ccf.tsv` | Per call: number of junctions and median ccf |
| `consensus/<sample>/region_beds/` | Region BEDs for assembly (Flye branch) |
| `plots/consensus/*_ecdna.pdf` | One figure per region |
| `plots/svclone/*_svclone.pdf` | Figures with SVclone ccf on the edges |
| `plots/flye/*_flye.pdf` | Figures with Flye contigs |
| `de_novo_assembly/<sample>/<region>/` | Extracted region reads |
| `de_novo_assembly/flye/<sample>/<region>/` | Flye assembly, circular contigs and PAF to the reference |
