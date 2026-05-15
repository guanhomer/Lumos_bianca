# Lumos (Bianca fork)

Bianca-adapted fork of [KolmogorovLab/Lumos](https://github.com/KolmogorovLab/Lumos) for long-read somatic variant analysis on the UPPMAX Bianca cluster.

Fork comparison:

* upstream: `KolmogorovLab/Lumos`
* fork: `guanhomer/Lumos_bianca`

This fork extends the original workflow with:

* Bianca-compatible Slurm execution
* offline/container-aware execution
* aligned BAM input support through sample-local `bam/` directories
* watchdog-based automatic submission
* project-local reference and annotation organization
* UPPMAX-specific environment handling

The workflow is intended for Oxford Nanopore Technologies (ONT) long-read tumor analysis.

---

# Repository structure

The Bianca fork introduces additional directories and helper scripts compared with the upstream repository.

Example repository layout:

```text
Lumos/
├── annot/
│   ├── cosmic_genes.tsv
│   ├── hg38_cpg_cleaned.bed
│   ├── human_GRCh38_no_alt_analysis_set.trf.bed
│   └── PoN_1000G_hg38_extended.tsv.gz
│
├── bam/
│   ├── sample_part1.bam
│   ├── sample_part2.bam
│   ├── sample_part3.bam
│   └── sample_part4.bam
│
├── clair3_models_pytorch/
│   └── r1041_e82_400bps_sup_v520/
│
├── sigularity/
│   ├── google-deepsomatic-1.9.0.img
│   └── gokcekeskus-severus-v1_6.img
│
├── processes/
│   ├── merge_bams_lumos.nf
│   └── processes.nf
│
├── script/
│   ├── lumos_watchdog.sh
│   ├── download_singularity_image.sh
│   └── run_lumos.sh
│
├── reference/
│   ├── GRCh38.p14.genome.fa
│   └── GRCh38.p14.genome.fa.fai
│
├── tumorOnlyONT.nf
├── nextflow.config
├── uppmax_ryan.config
└── README.md
```

## Directory descriptions

| Path                       | Description                                   |
| -------------------------- | --------------------------------------------- |
| `annot/`                   | Annotation resources used by Lumos modules    |
| `bam/`                     | Input aligned BAM files and indexes           |
| `clair3_models_pytorch/`   | Clair3 ONT models                             |
| `reference/`               | Reference genome and indexes                  |
| `script/run_lumos.sh`      | Bianca execution wrapper                      |
| `script/lumos_watchdog.sh` | Upload monitor and automatic Slurm submission |
| `uppmax_ryan.config`       | Bianca-specific Nextflow configuration        |

---

# Workflow overview

Lumos is a Nextflow workflow for long-read somatic variant analysis.

Current functionality includes:

* structural variant calling
* copy-number analysis
* somatic SNV analysis
* methylation analysis
* phasing
* tumor-only ONT workflows

Primary tools used:

| Tool        | Purpose                    |
| ----------- | -------------------------- |
| Clair3      | Small variant calling      |
| longphase   | Read phasing               |
| Severus     | Structural variant calling |
| Wakhan      | Copy-number analysis       |
| DeepSomatic | Somatic SNV calling        |
| modkit      | DNA methylation analysis   |

The original implementation targets NIH Biowulf. This fork adapts the workflow for UPPMAX Bianca.

---

# Bianca-specific modifications

This fork modifies the upstream workflow to support execution within the Bianca secure-compute environment.

Changes include:

* Slurm scheduling defaults for UPPMAX
* Apptainer/Singularity execution
* offline Nextflow execution (`NXF_OFFLINE=true`)
* project-local cache and temporary directories
* sample-local BAM organization
* automatic monitoring/submission utilities
* compatibility with restricted network access

Depending on project allocation and storage organization, path changes may still be required.

---

# Reference data preparation

Prepare the following resources before running the workflow.

## Annotation resources

[Wakhan/scripts/cosmic.py at 0218b4ffca48c277fbdad53c04efed244b9aff21 · KolmogorovLab/Wakhan](https://github.com/KolmogorovLab/Wakhan/blob/0218b4ffca48c277fbdad53c04efed244b9aff21/scripts/cosmic.py)

Required resources:

| Resource                                   | Purpose                                  |
| ------------------------------------------ | ---------------------------------------- |
| `cosmic_genes.tsv`                         | COSMIC annotations for Wakhan            |
| `hg38_cpg_cleaned.bed`                     | CpG BED regions for methylation analysis |
| `human_GRCh38_no_alt_analysis_set.trf.bed` | VNTR/tandem repeat regions               |
| `PoN_1000G_hg38_extended.tsv.gz`           | SV panel of normals                      |

## Clair3 model

https://github.com/nanoporetech/rerio/tree/master/clair3_models

Example:

```text
clair3_models_pytorch/
└── r1041_e82_400bps_sup_v520/
```

---

# Input organization

The Bianca fork is configured for aligned BAM input.

Expected sample layout:

```text
<sample>/
└── bam/
    ├── sample_part1.bam
    ├── sample_part2.bam
    └── sample_part3.bam
```

The workflow wrapper infers the sample name from the current directory:

```bash
SAMPLE="$(basename "$PWD")"
```

The wrapper uses:

```bash
--aligned_input
--aligned_tumor_dir bam
```

instead of:

```bash
--reads_tumor
```

---

# Running the workflow

The recommended entry point is:

```text
script/run_lumos.sh
```

This wrapper configures:

* Nextflow environment variables
* container cache paths
* temporary directories
* Bianca-specific configuration
* tumor-only aligned BAM execution

## Example manual execution

```bash
module load bioinfo-tools Nextflow samtools/1.20

cd /proj/nobackup/sens2024549/Lumos

sbatch \
    -A sens2024549 \
    -n 3 \
    -t 24:00:00 \
    --wrap="bash script/run_lumos.sh" \
    -o slurm-%j.out
```

---

# Automatic upload monitoring

The repository includes:

```text
script/lumos_watchdog.sh
```

The watchdog monitors upload directories and automatically submits Lumos jobs after BAM uploads become stable.

## Expected upload structure

```text
UPLOAD_ROOT/
└── <sample>/
    └── bam/
        └── *.bam
```

Only `bam/` directories located two levels below `UPLOAD_ROOT` are monitored.

---

# Watchdog logic

The watchdog:

1. scans the upload root every 5 minutes
2. detects BAM directories
3. verifies upload completion
4. waits for file stability
5. submits `run_lumos.sh` with `sbatch`
6. records submission state

A sample is considered stable when:

* no temporary upload files exist
* BAM and index files stop changing
* file sizes remain stable for 15 minutes

Temporary upload suffixes checked:

```text
*.filepart
*.part
*.tmp
```

Monitoring aborts after 24 hours if uploads never stabilize.

---

# Watchdog state files

State directory:

```text
${UPLOAD_ROOT}/.lumos_watchdog_state/
```

## State markers

| File                  | Meaning                               |
| --------------------- | ------------------------------------- |
| `watchdog.pid`        | Prevents duplicate watchdog instances |
| `<sample>.uploading`  | Sample currently being monitored      |
| `<sample>.processing` | Workflow currently running            |
| `<sample>.done`       | Submission completed                  |

## Logs

Main watchdog log:

```text
${UPLOAD_ROOT}/watchdog.log
```

Cluster summary log:

```text
${UPLOAD_ROOT}/HPC_status.log
```

The status log records:

* `squeue`
* `projinfo`
* `projinfo -s`
* `uquota`

---

# Pipeline modes

| Mode         | Description                                  |
| ------------ | -------------------------------------------- |
| `sv_cna`     | Structural variants and copy-number analysis |
| `sv_cna_dmr` | SV, CNA, and methylation analysis            |
| `all`        | Full workflow including somatic SNV analysis |

Example:

```bash
--mode all
```

---

# Output

Default output directory:

```text
lumos_out/
```

Custom output directory:

```bash
-output-dir <directory>
```

Typical outputs:

* phased BAMs
* structural variants
* CNA results
* somatic SNVs
* methylation calls
* QC metrics
* intermediate workflow files

---

# Example runtime and storage requirements

The following example corresponds to a tumor-only ONT run using aligned BAM input.

## Example sample

| Metric               | Value        |
| -------------------- | ------------ |
| Sample               | `BT-178` |
| Estimated bases      | 176.34 Gb    |
| Input BAM size       | 130 GB       |
| Runtime              | 15 h 48 min  |
| CPU hours            | 222.5        |
| Intermediate storage | ~700 GB      |
| Final result size    | ~5.6 GB      |

## Storage considerations

Large temporary and intermediate files are generated during:

* phasing
* SV calling
* CNA analysis
* methylation pileup generation
* haplotype-specific BAM processing

The intermediate workspace may require several times the input BAM size.

For a ~130 GB BAM input, total temporary usage reached approximately:

```text
700 GB
```

Using `$SLURM_TMPDIR` or project-local scratch storage is strongly recommended.

---

# Tumor-only output structure

Example major output directories:

```text
lumos_out/
├── phased_vcf/
├── rephased_vcf/
├── severus/
├── wakhan/
├── deepsomatic/
├── methylation/
└── execution/
      ├── pipeline.report.html
      ├── pipeline.timeline.html
      └── trace.txt
```

---

# Core analysis outputs

## Phasing

| File                                                            | Description                           |
| --------------------------------------------------------------- | ------------------------------------- |
| `phased_vcf/longphase.vcf.gz`                                   | Germline variants phased by Longphase |
| `rephased_vcf/wakhan_hapcorrect/phasing_output/rephased.vcf.gz` | Wakhan-refined phasing results        |

---

## Structural variants (Severus)

| File or directory                                     | Description                       |
| ----------------------------------------------------- | --------------------------------- |
| `severus/severus_out/somatic_SVs/severus_somatic.vcf` | Somatic structural variants       |
| `severus/severus_out/all_SVs/severus_all.vcf`         | All structural variants           |
| `severus/severus_out/*/plots/`                        | Interactive breakpoint HTML plots |
| `severus/severus_out/*/breakpoint_clusters.tsv`       | Breakpoint cluster information    |

---

## Copy-number analysis (Wakhan)

Main directory:

```text
wakhan/wakhan_cna/
```

Typical contents:

| Path                  | Description                            |
| --------------------- | -------------------------------------- |
| `variation_plots/`    | Per-chromosome CNV plots               |
| `bed_output/`         | CNA segments, LOH regions, gene states |
| `vcf_output/`         | CNA calls in VCF format                |
| `coverage_data/`      | Coverage metrics and BAF data          |
| `coverage_plots/`     | Coverage visualizations                |
| `solutions_ranks.tsv` | Ranking of ploidy/purity solutions     |
| `wakhan_hapcorrect/`  | Haplotype correction outputs           |

Numeric solution directories such as:

```text
1.99_0.39_0.95
3.98_0.76_0.96
```

represent alternative ploidy/purity solutions evaluated by Wakhan.

---

## Somatic SNVs and indels (DeepSomatic)

| File                                           | Description                        |
| ---------------------------------------------- | ---------------------------------- |
| `deepsomatic/deepsomatic_out/ds.merged.vcf.gz` | Merged somatic SNV and indel calls |

---

## Methylation analysis (modkit)

| File                                   | Description                                         |
| -------------------------------------- | --------------------------------------------------- |
| `methylation/dmr.bed`                  | Differential methylation regions between haplotypes |
| `methylation/haplotagged_*.bed.gz`     | Haplotype-specific methylation pileups              |
| `methylation/haplotagged_*.stats.tsv`  | Haplotype-specific methylation statistics           |
| `methylation/pileup.bed.gz`            | Combined methylation pileup                         |
| `methylation/pileup.stats.tsv`         | Global methylation statistics                       |
| `methylation/pileup_cpg_subset.bed.gz` | CpG-filtered methylation pileup                     |

---

# Pipeline execution reports

| File                     | Description                      |
| ------------------------ | -------------------------------- |
| `execution/`             | Nextflow execution metadata      |
| `pipeline.report.html`   | Pipeline execution summary       |
| `pipeline.timeline.html` | Timeline visualization           |
| `trace.txt`              | Detailed process execution trace |

---

# Visualization entry points

## CNA visualization

```text
wakhan/wakhan_cna/*/variation_plots/CN_VARIATION_INDEX.html
```

## Structural variant visualization

```text
severus/severus_out/*/plots/severus_*.html
```

## Coverage visualization

```text
wakhan/wakhan_cna/coverage_plots/COVERAGE_INDEX.html
```

## Phase correction visualization

```text
wakhan/wakhan_cna/phasing_output/PHASE_CORRECTION_INDEX.html
```

---

# Notes

* The workflow was tested primarily with ONT data.
* Nextflow version compatibility is important.
* Bianca network restrictions may prevent automatic container downloads.
* Pre-pulling Singularity/Apptainer containers is recommended.
* Large temporary files should use `$SLURM_TMPDIR` when available.
* The repository assumes project-local storage organization on Bianca.

