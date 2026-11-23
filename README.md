# MISSion_analyses

Nextflow workflows for paired-end short-read sequencing analysis in the [MISSion project](https://gepris.dfg.de/project/425959793?lang=en). The primary workflow supports reference mapping, variant calling, isolate assembly, assembly QC, and screening of assembled sequences for antimicrobial-resistance, plasmid-associated, and virulence-associated genes. An optional Kraken2-based branch additionally performs taxon-specific read selection, abundance estimation, and metagenomic analyses.

> **Project status:** This is a project-specific **Nextflow DSL2** repository, **not an nf-core release or an nf-core-validated pipeline**. This README uses the familiar nf-core sections (Introduction, Pipeline summary, Usage, Pipeline output, Credits/Support, Citations) but documents the *actual* parameter names, source code, and publication rules.

## Introduction

The repository contains two independent entry points:

- **[`analyses_pipeline.nf`](analyses_pipeline.nf)** — the main analysis pipeline described throughout this README. It selects exactly one of two branches:
  - **`ASSEMBY_NOFILTER`**: runs when `--krakendb` and `--taxid` are **not both provided**. Reads are trimmed, mapped and assembled without taxon-specific read selection.
  - **`ASSEMBLY_FILTER`**: runs when **both** `--krakendb` and `--taxid` are provided. After trimming, Kraken2 and KrakenTools retain reads assigned to the specified taxonomy ID and its descendants; these retained reads are used for the main mapping and isolate-assembly analyses. This branch also runs extra classification/abundance and metagenomic steps.
- **[`read_filter.nf`](read_filter.nf)** — a separate host/reference-depletion-oriented workflow: fastp -> Kraken2 -> exclusion with KrakenTools -> BWA-MEM2 -> extraction of doubly unmapped read pairs. It is **not called by `analyses_pipeline.nf`** and has its own limitations (see [Standalone read-filtering workflow](#standalone-read-filtering-workflow)).


## Pipeline summary

| Analysis step | Software / implementation | Unfiltered branch | Filtered branch |
| --- | --- | :---: | :---: |
| Paired-end adapter/quality trimming and QC | fastp (`--trim_poly_g`) | Yes | Yes |
| Count trimmed reads | `zcat`, `awk` (R1 FASTQ) | — | Yes |
| Classify reads | Kraken2 | — | Yes |
| Retain reads belonging to a taxon and descendants | KrakenTools `extract_kraken_reads.py --include-children` | — | Yes |
| Reclassify retained reads for QC | Kraken2 | — | Yes |
| Abundance estimates at species and genus ranks | Bracken | — | Yes |
| Compress and count selected reads | pigz, `zcat`, `awk` | — | Yes |
| Index reference; map reads; sort/index BAM | BWA-MEM2, SAMtools | Yes | Yes |
| Summarize alignment statistics | SAMtools `flagstat`, MultiQC | Yes | Yes |
| Call per-sample variants and a joint multi-sample VCF | FreeBayes, freebayes-parallel | Yes | Yes |
| Count heterozygous/homozygous genotypes in per-sample VCFs | BCFtools | — | Yes |
| Map selected reads against an additional MLST reference | BWA-MEM2 | — | Optional (`--mlst_mapping`) |
| Assemble isolate contigs | Shovill (SPAdes backend) | Yes | Yes |
| Compress and BWA-index assemblies | pigz, BWA-MEM2 | Yes | Yes |
| Assembly assessment against FASTA/GFF reference | QUAST, MultiQC | Yes | Yes |
| Classify and select assembled contigs | Kraken2, KrakenTools | Not normally active (see limitations) | Yes |
| Screen assemblies against seven gene databases | ABRicate | Yes | Yes |
| Extract detected gene/hit sequences | Seqtk, shell utilities | Yes | Yes |
| In silico *Salmonella* typing | SISTR | Optional (`--typing`) | Optional (`--typing`) |
| Additional metagenomic assembly on **trimmed, not taxon-filtered** reads | metaSPAdes | — | Yes (graph-output issue noted below) |
| Screen metagenomic contigs for genes | ABRicate, Seqtk | — | Yes (depends on assembly outputs) |
| Profile **trimmed, not taxon-filtered** reads | MetaPhlAn | — | Optional (`--metaphlan_db`) |

## Usage

### Software and compute requirements

- A **Linux/Unix** environment with **Nextflow DSL2**, a compatible Java runtime, and common command-line tools (`bash`, `awk`, `cut`, `sed`, `gzip`, etc.). 
- **Conda** or **Mamba** to create the environments defined under [`envs/`](envs/), using `-profile mission_conda` or `-profile mission_mamba`. 
- Access to the tool-specific reference/data assets: Kraken2 and Bracken databases for taxon filtering, installed ABRicate databases, and any optional MetaPhlAn database.
- Adequate CPU, memory and disk space. Source requests include **750 GB RAM** for some Kraken2/Bracken jobs, **500 GB** for metaSPAdes, **150 GB** for MetaPhlAn, **250 GB** for BWA indexing in `read_filter.nf`, and **60 CPUs** for joint FreeBayes. These values are configurations, **not measured minimum requirements**.
- Shared visibility of input, reference, database, work and output paths on cluster worker nodes.

The configuration includes the following profiles:

| Profile | Purpose | Notes |
| --- | --- | --- |
| `mission_conda` | Enables Conda environments | Uses `--conda_cache_dir`, creation options and timeout |
| `mission_mamba` | Enables Conda with Mamba | Suitable when Mamba is installed |
| `mission_slurm` | Uses SLURM executor | `queueSize = 200`; can combine with `mission_mamba` or `mission_conda` |
| `mission_singularity` | Enables Singularity settings | Modules do **not** define container images; this profile alone does not establish a containerized setup |

### Samplesheet input

Pass the **comma-delimited CSV** samplesheet using `--config` (a **pipeline parameter**). Do **not** confuse this with Nextflow's single-dash `-c` argument for an execution configuration file. The CSV must have one header row and the following required columns:

```csv
sample,raw_fwd,raw_rev,raw_folder
SAMPLE_001,SAMPLE_001_R1.fastq.gz,SAMPLE_001_R2.fastq.gz,/data/reads
SAMPLE_002,SAMPLE_002_R1.fastq.gz,SAMPLE_002_R2.fastq.gz,/data/reads
```

| Column | Required | Interpretation in the source |
| --- | :---: | --- |
| `sample` | Yes | Sample identifier used in published paths and filenames; use a unique, filesystem-safe ID |
| `raw_fwd` | Yes | Forward (R1) FASTQ filename under `raw_folder` |
| `raw_rev` | Yes | Reverse (R2) FASTQ filename under `raw_folder` |
| `raw_folder` | Yes | Directory containing both read files; absolute paths recommended |

The code builds each input path as `raw_folder + "/" + raw_fwd` / `raw_rev`. It expects **paired-end FASTQ** files, whether plain or gzip-compressed as supported by fastp. Input filenames are not validated by a schema. A commented-out `trimmed_folder` field is not used. **Multiple lanes with repeated sample IDs are not merged**: combine them beforehand or provide unique sample IDs to avoid duplicate output names.

### References and databases

| Parameter | Required | Meaning and expected input |
| --- | --- | --- |
| `--reference_fasta` | Main workflow: **yes** | Reference genome FASTA used by BWA-MEM2, FreeBayes and QUAST; should correspond to the GFF |
| `--reference_gff` | Main workflow: **yes** | Reference annotation GFF used by QUAST |
| `--krakendb` | For filtered branch | Kraken2 database directory; must be accessible to processes |
| `--taxid` | For filtered branch | NCBI taxonomy ID to **retain**, including child taxa, e.g. `590` for *Salmonella* |
| `--krakendb_contigs` | Optional in filtered branch | Database for retained-read QC and contig classification; defaults to `--krakendb` |
| `--mlst_mapping` | Optional in filtered branch | An additional reference genome/FASTA for BWA-MEM2 mapping (not an MLST allele caller) |
| `--metaphlan_db` | Optional in filtered branch | Database setting passed to both MetaPhlAn `--bowtie2db` and `-x` (verify path vs. index-prefix syntax for your installed MetaPhlAn) |
| ABRicate databases | Required by gene screening | `argannot`, `plasmidfinder`, `megares`, `ncbi`, `card`, `vfdb`, `resfinder` installed/available to ABRicate |

**Bracken note:** Bracken uses `-r 100`, species rank (`S`) and genus rank (`G`), with threshold `0`. Prepare the appropriate distribution files for the Kraken2 database and this read length.

### Quick start: unfiltered branch

Run from the repository root:

```bash
nextflow run ./analyses_pipeline.nf \
  -profile mission_mamba \
  --config ./samplesheet.csv \
  --reference_fasta /path/to/reference.fasta \
  --reference_gff /path/to/reference.gff \
  --analysesdir ./results \
  -with-report execution_report.html
```

The pipeline selects `ASSEMBY_NOFILTER` because `--krakendb` and `--taxid` are not both supplied.

### Quick start: taxon-filtered branch

```bash
nextflow run ./analyses_pipeline.nf \
  -profile mission_slurm,mission_mamba \
  --config ./samplesheet.csv \
  --reference_fasta /path/to/reference.fasta \
  --reference_gff /path/to/reference.gff \
  --krakendb /path/to/kraken2_database \
  --taxid 590 \
  --analysesdir ./results \
  -with-report execution_report.html
```

Optional arguments for this branch include `--krakendb_contigs /path/to/other_db`, `--mlst_mapping /path/to/mlst_reference.fasta`, `--metaphlan_db /path/to/metaphlan_database`, and `--typing true` for SISTR. Omit `mission_slurm` for local execution. Add `-resume` to reuse eligible completed tasks. Standard Nextflow options such as `-with-trace`, `-with-timeline` and `-with-dag` are also available.

> **Branch-selection caveat:** Supplying *only one* of `--krakendb` and `--taxid` selects the **unfiltered** branch; it does not raise an error. The pipeline validates the presence of `--config`, `--reference_fasta` and `--reference_gff`.

### Parameters

These are the active workflow parameters or configuration parameters:

| Parameter | Default in `nextflow.config` | Effect |
| --- | --- | --- |
| `--config` | Not explicitly declared | **Required:** samplesheet CSV |
| `--reference_fasta` | `null` | **Required:** reference genome FASTA |
| `--reference_gff` | `null` | **Required:** reference annotation GFF |
| `--analysesdir` | `null` | Intended root for published results; supply explicitly |
| `--krakendb` | `null` | Kraken2 database for read classification; with `--taxid` selects filtered branch |
| `--taxid` | `null` | Taxon to retain in main filtered branch |
| `--krakendb_contigs` | `null` | Alternative Kraken2 QC/contig database |
| `--mlst_mapping` | `null` | Optional MLST-reference mapping, filtered branch only |
| `--typing` | Not declared (unset/falsey) | Runs SISTR when truthy in either branch |
| `--metaphlan_db` | `null` | Enables MetaPhlAn profiling in filtered branch |
| `--conda_cache_dir` | `conda` | Cache directory used by Conda/Mamba profiles |
| `--conda_create_options` | Empty string | Conda environment creation options |
| `--conda_create_timeout` | `20 min` | Conda environment creation timeout |
| `--singularity_cache_dir` | `singularity` | Singularity cache directory |

## Pipeline steps

### Shared analyses

1. **fastp** trims paired reads, including poly-G trimming, and writes R1/R2 gzip-compressed FASTQ plus JSON and HTML QC reports.
2. **BWA-MEM2** indexes the supplied reference and maps trimmed reads (unfiltered branch) or selected reads (filtered branch). **SAMtools** writes original and coordinate-sorted BAMs and the sorted BAM index; `flagstat` produces mapping statistics.
3. **FreeBayes** calls per-sample variants with `--gvcf`, and `freebayes-parallel` produces a joint VCF for all mapped samples. No post-calling VCF filtering is performed.
4. **Shovill** assembles paired-end reads with a fixed estimated genome size of `5M`. Contigs and the Shovill SPAdes graph are renamed, published, and contigs gzip-compressed; assembly FASTA is indexed with BWA-MEM2 for downstream tasks but the index is not published.
5. **QUAST** assesses each assembly against `--reference_fasta` and `--reference_gff`. **MultiQC** combines mapping and assembly reports across samples in separate report groups.
6. **ABRicate** screens assemblies against seven databases; **Seqtk** extracts sequences corresponding to reported hits. These databases cover more than just antimicrobial resistance (e.g. plasmid/virulence datasets).
7. **SISTR**, only when `--typing` is enabled, produces *Salmonella*-specific typing/cgMLST outputs. There is no automated check that a sample is *Salmonella*.

### Additional steps in the Kraken2-filtered branch

1. **Read counting:** `count_reads` counts records in R1 by dividing the FASTQ line count by four.
2. **Kraken2** assigns trimmed read pairs to taxa and writes classified FASTQs plus reports.
3. **KrakenTools** selects the requested taxon and descendants (`--include-children`) from the Kraken2-classified reads. Selected reads are compressed with pigz and reclassified with Kraken2 using the QC/contig database.
4. **Bracken** estimates abundance at species and genus rank from the original read-classification reports. Separate MultiQC reports are generated for primary read classification, classification QC, species/genus abundance, and contig classification.
5. **Per-sample VCF genotype counts** are calculated with BCFtools and compiled into a summary CSV.
6. **Optional additional MLST mapping:** selected reads are aligned to `--mlst_mapping`. This is a second alignment, not an MLST typing or allele-calling algorithm.
7. **Contig classification:** Kraken2 classifies Shovill contigs using `--krakendb_contigs` (or `--krakendb`); KrakenTools extracts the requested taxa from classified contigs. ABRicate analyses still use the assembled contigs, not the extracted subset.
8. **Additional metaSPAdes assembly:** uses **trimmed reads before Kraken2 selection**, not the selected reads. Its contigs are intended for an additional ABRicate/Seqtk screen. 
9. **Optional MetaPhlAn** also profiles the **trimmed, pre-selection** paired reads, when `--metaphlan_db` is provided.

## Pipeline output

### Read preprocessing, classification and abundance

| Publication destination (under `<RESULTS>`) | Expected published files | Branch |
| --- | --- | --- |
| `01_trimming/<sample>/` | `<sample>.R1.trimmed.fastq.gz`, `<sample>.R2.trimmed.fastq.gz`, `<sample>.trim_report.json`, `<sample>.trim_report.html` | Both |
| `01_trimming/<sample>/` | `<sample>_read_count.csv` for trimmed R1 | Filtered |
| `02_classification/filtered/<sample>/` | `<sample>.kraken.out.txt`, `<sample>.kraken.report.txt` (initial classification) | Filtered |
| `02_classification/qc/<sample>/` | `<sample>.kraken.out.txt`, `<sample>.kraken.report.txt` (selected-read QC) | Filtered |
| `02_classification/<sample>/` | `<sample>.S.bracken`, `<sample>.S.bracken.report`, `<sample>.G.bracken`, `<sample>.G.bracken.report` | Filtered |
| `02_read_filter/<sample>/` | `<sample>.classified.R1.fastq.gz`, `<sample>.classified.R2.fastq.gz` (retained taxon) | Filtered |
| `02_read_filter/<sample>/` | `<sample>_read_count.csv` for selected reads (**currently unreliable**, see limitations) | Filtered |

### Mapping, variant calling and variant statistics

| Publication destination | Expected published files | Branch |
| --- | --- | --- |
| `03_mapping/<mode>/<sample>/` | `<sample>.bam`, `<sample>.sort.bam`, `<sample>.sort.bam.bai`; `<mode>` = `unfiltered`, `filtered`, or optional `filtered_mlst` | Both / MLST optional |
| `03_mapping/<sample>/` | `<sample>.flagstat` (primary reference mapping only) | Both |
| `08_variantcall/<sample>/` | `<sample>.g.vcf` (per-sample FreeBayes gVCF) | Both |
| `08_variantcall/` | `all_samples.vcf.gz` (joint variant call; declared internally as `<mode>/all_samples.vcf.gz`) | Both |
| `99_reporting/vcf_statistics/<sample>/` | `<sample>.csv` (`SAMPLE,HET,HOM,HET_FRAC`) | Filtered |
| `99_reporting/` | `variant_count_compiled.csv` | Filtered |

### Assemblies, QC, taxonomy and gene screening

| Publication destination | Expected published files | Branch |
| --- | --- | --- |
| `04_assembly/<mode>/<sample>/` | `<sample>_contigs.fasta` (renamed contigs), `<sample>_contigs.gfa` (Shovill SPAdes graph), `<sample>_contigs.fasta.gz`; `<mode>` = `unfiltered` or `filtered` | Both |
| `05_qc_assembly/` | QUAST outputs declared as `<sample>/<mode>/report.html`, `<sample>/<mode>/report.tsv`, `<sample>/<mode>/icarus.html` | Both |
| `05_qc_assembly/<sample>/` | `<sample>.kraken.out.txt`, `<sample>.kraken.report.txt` for contigs | Filtered |
| `04_assembly/filtered/<sample>/` | `<sample>.classified_contigs.fasta` (Kraken2 classified output), `<sample>.classified.fasta` (KrakenTools selected contigs) | Filtered |
| `06_abricate/<sample>/` | `<sample>.<database>.<mode>.tsv` (ABRicate hits), `<sample>.<database>.<mode>.fasta` (hit sequences); modes `unfiltered`, `filtered`, and `meta` where applicable | Both; `meta` filtered only |
| `07_typing/<sample>/` | `<sample>-allele-results.json`, `<sample>-novel-alleles.fasta`, `<sample>-cgmlst-profile.csv`, `<sample>-sistr-output.tab` | Optional `--typing` |
| `04_assembly/meta/<sample>/` | metaSPAdes `contigs.fasta`, `scaffolds.fasta`; additionally intended renamed `<sample>_contigs.fasta`, `<sample>_contigs.gfa`, and compressed `<sample>_contigs.fasta.gz` | Filtered |

`<database>` is one of `argannot`, `plasmidfinder`, `megares`, `ncbi`, `card`, `vfdb`, `resfinder`. QUAST and other report content can be absent for processes that fail or are ignored.

### Metagenomics and aggregate reporting

| Publication destination | Expected published files | Branch |
| --- | --- | --- |
| `09_metagenomes/filtered/<sample>/` | `<sample>.txt` (MetaPhlAn taxonomic profile), `<sample>.bowtie.bz2` (Bowtie2 intermediate output) | Filtered, only with `--metaphlan_db` |
| `99_reporting/<analysis>_statistics/` | `<analysis>_multiqc_report.html`, `<analysis>_multiqc_report_data/` | As relevant |
| `99_reporting/` | `trimmed_read_count_compiled.csv`, `filtered_read_count_compiled.csv` (**filtered counts currently unreliable**) | Filtered |

`<analysis>` is `mapping` or `assembly` in both branches; in the filtered branch additional values are `classification`, `classification_QC`, `classification_contigs`, `abundance_species`, `abundance_genus`. Each generates a **separate** MultiQC report; no all-in-one MultiQC report is implemented.

### Representative output hierarchy

```text
<RESULTS>/
├── 01_trimming/<sample>/                  # fastp FASTQ and reports; count when filtered
├── 02_classification/                     # filtered branch: Kraken2 and Bracken
│   ├── filtered/<sample>/
│   ├── qc/<sample>/
│   └── <sample>/
├── 02_read_filter/<sample>/               # retained read pairs and read counts
├── 03_mapping/
│   ├── <mode>/<sample>/                    # BAM, sorted BAM, BAI; includes filtered_mlst
│   └── <sample>/                           # primary-mapping flagstat
├── 04_assembly/
│   ├── unfiltered/<sample>/                # unfiltered Shovill result
│   ├── filtered/<sample>/                  # filtered Shovill result + taxonomy-selected contigs
│   └── meta/<sample>/                      # metaSPAdes outputs, if produced
├── 05_qc_assembly/                        # QUAST and contig Kraken reports
├── 06_abricate/<sample>/                  # ABRicate TSV and extracted FASTA
├── 07_typing/<sample>/                    # optional SISTR
├── 08_variantcall/                        # per-sample and joint FreeBayes files
├── 09_metagenomes/filtered/<sample>/      # optional MetaPhlAn
└── 99_reporting/                          # MultiQC, read counts, variant summary
```

## Standalone read-filtering workflow

[`read_filter.nf`](read_filter.nf) is separate from the main pipeline. Its intended goal is to exclude a taxon (commonly a host) and retain read pairs not mapped to an additional reference:

1. Parse the **same comma-separated samplesheet** (`sample`, `raw_fwd`, `raw_rev`, `raw_folder`).
2. Trim the reads with fastp.
3. Build a **large-reference** BWA-MEM2 index (requests 250 GB memory).
4. Classify trimmed read pairs with Kraken2 using `--krakendb`.
5. Use KrakenTools `--exclude --include-children` to remove `--taxid` and its descendants **from the Kraken2-classified FASTQs**.
6. Compress selected reads, map them to `--reference_fasta`, and use `samtools fastq -f 12` to extract records for which **both read and mate are unmapped**.
7. Publish final compressed read pairs to `04_filtered_reads/<sample>/` as `<sample>.R1.fastq.gz` and `<sample>.R2.fastq.gz`.

Example (paths illustrative):

```bash
nextflow run ./read_filter.nf \
  -profile mission_slurm,mission_mamba \
  --config ./samplesheet.csv \
  --reference_fasta /path/to/host_reference.fasta \
  --krakendb /path/to/kraken2_database \
  --taxid 9606 \
  --analysesdir ./read_filter_results
```

## Additional repository files

- [`vcf2distance_matrix.py`](vcf2distance_matrix.py) — a standalone Python helper for sample distance matrices from VCF data; **not invoked by Nextflow**.
- [`plot_amr.Rmd`](plot_amr.Rmd), [`plot_amr_envsamples.Rmd`](plot_amr_envsamples.Rmd), [`plot_metagenomic_data.Rmd`](plot_metagenomic_data.Rmd) — standalone R Markdown analyses; **not outputs generated by this workflow**.
- [`envs/`](envs/) — Conda package specifications used by module processes.
- [`nextflow.config`](nextflow.config) — runtime parameters and site-specific profiles.

## Citations and documentation

For scientific use, record and cite the exact versions of Nextflow, tools, reference genome/annotation, taxonomy databases, and ABRicate databases used for each analysis. The repository does not supply a versioned `CITATIONS.md`. Tool documentation:

- [nf-core pipeline documentation conventions](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/README.md) and [nf-core developer documentation](https://nf-co.re/docs/developing/pipelines/template-files)
- [Nextflow documentation](https://www.nextflow.io/docs/latest/) and [`publishDir` directive](https://docs.seqera.io/nextflow/reference/process/directives/publish-dir)
- [fastp](https://github.com/OpenGene/fastp), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2), [SAMtools](http://www.htslib.org/), [FreeBayes](https://github.com/freebayes/freebayes)
- [Shovill](https://github.com/tseemann/shovill), [SPAdes and metaSPAdes](https://ablab.github.io/spades/), [QUAST](https://quast.sourceforge.net/), [MultiQC](https://multiqc.info/)
- [Kraken2](https://github.com/DerrickWood/kraken2), [KrakenTools](https://github.com/jenniferlu717/KrakenTools), [Bracken](https://github.com/jenniferlu717/Bracken), [MetaPhlAn](https://github.com/biobakery/MetaPhlAn)
- [ABRicate](https://github.com/tseemann/abricate), [Seqtk](https://github.com/lh3/seqtk), [SISTR](https://github.com/phac-nml/sistr_cmd)
