# ghga-de/QCMetrics Nextflow Pipeline

[![GitHub Actions CI Status](https://github.com/GHGA/qc/actions/workflows/nf-test.yml/badge.svg)](https://github.com/GHGA/qc/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/GHGA/qc/actions/workflows/linting.yml/badge.svg)](https://github.com/GHGA/qc/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<p align="center">
    <img title="GHGA_logo" src="docs/GHGA_short_Logo_orange.png" width=50%>
</p>

## Introduction

**ghga-de/QCMetrics** is a bioinformatics pipeline that performs basic quality control over input datasets without altering the raw data. It accepts three main input types:
1. Raw FastQ files
2. Aligned BAM/CRAM files
3. Variant called VCF/BCF files

The pipeline automatically selects the appropriate quality control tools based on your provided analysis method (e.g., WGS, RNA-Seq, Nanopore) and compiles the results into a single MultiQC report.

## Tool & Analysis Matrix

The following table details which tools are executed based on the analysis method and input data type provided in the samplesheet.

| Analysis Method | Read QC (FastQ) | Alignment QC (BAM/CRAM) | Variant QC (VCF) |
| :--- | :--- | :--- | :--- |
| **WGS / WES / TES** | FastQC, FastP, SeqFU | Mosdepth, Samtools Stats, Picard, VerifyBamID, NGS-Bits* | BCFTools Stats |
| **ATAC / ChIP-Seq** | FastQC, FastP, SeqFU | Mosdepth, Samtools Stats, Picard | BCFTools Stats |
| **RNA-Seq / smRNA** | FastQC, FastP, SeqFU | RSeQC | BCFTools Stats |
| **Nanopore** | FastQC, NanoPlot | - | BCFTools Stats |
| **PacBio** | FastPLong | - | BCFTools Stats |
| **MethylSeq** | FastQC, FastP, SeqFU | - | BCFTools Stats |

> \* *NGS-Bits SampleGender is run for WGS if `predict_sex` is enabled.*

## Usage

### 1. Prepare Samplesheet

You must create a `samplesheet.csv` containing your input data. The structure requires a `step` column to tell the pipeline which type of file you are providing:
* **step 1**: FastQ files (Read QC)
* **step 2**: BAM/CRAM files (Alignment QC)
* **step 3**: VCF files (Variant QC)

**Example `samplesheet.csv`:**

```csv
sample,fastq_1,fastq_2,bam,vcf
SAMPLE_A,read1.fq.gz,read2.fq.gz,,
SAMPLE_A,read1.fq.gz,read2.fq.gz,,
SAMPLE_B,read1.fq.gz,,,
SAMPLE_C,,,aligned.bam,
SAMPLE_D,,,,variants.vcf.gz
```

### 2. Run the Pipeline

Run the pipeline using the command below. Ensure you specify the correct method (e.g., `wgs`, `rna`, `wes`) so the pipeline loads the correct tools.

```bash
nextflow run main.nf \
   -profile docker/singularity/conda \
   --input samplesheet.csv \
   --method wgs \
   --outdir ./results
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Supported Tools

### Read QC
* [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Comprehensive quality control checks on raw sequence data.
* [**FastP**](https://github.com/OpenGene/fastp): All-in-one FASTQ preprocessor (used here for QC metrics).
* [**SeqFU**](https://telatin.github.io/seqfu2/tools/metadata.html): Sequence statstics
* [**NanoPlot**](https://github.com/wdecoster/NanoPlot): Plotting tool for long read sequencing data and alignments.
* [**FastPLong**](https://github.com/OpenGene/fastplong): Quality control for long read data (PacBio).

### Alignment QC
* [**Mosdepth**](https://github.com/brentp/mosdepth): Fast BAM/CRAM depth calculation.
* [**Samtools Stats**](http://www.htslib.org/doc/samtools.html): General statistics for alignment files.
* [**Picard CollectMultipleMetrics**](https://broadinstitute.github.io/picard/): Collects multiple classes of metrics from alignment files.
* [**RSeQC**](http://rseqc.sourceforge.net/): Quality control for RNA-seq experiments.
* [**NGS-Bits SampleGender**](https://github.com/imgag/ngs-bits): Sex determination based on coverage.
* [**VerifyBamID**](https://github.com/Griffan/VerifyBamID): A robust tool for DNA contamination estimation from sequence reads using ancestry-agnostic method.

### Variant QC
* [**BCFTools Stats**](http://samtools.github.io/bcftools/bcftools.html): Statistics for VCF/BCF files.

### Reporting
* [**MultiQC**](http://multiqc.info/): Aggregates results from all tools into a single HTML report.

## Credits

ghga-de/QCMetrics was originally written by @kubranarci.

We thank the following people for their extensive assistance in the development of this pipeline:


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).