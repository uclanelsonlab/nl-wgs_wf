# nl-wgs_wf

A comprehensive Nextflow pipeline for whole genome sequencing (WGS) analysis of germline short-read data. This pipeline is inspired by the [rare diseases pipeline from nf-core](https://nf-co.re/raredisease).

## Overview

This pipeline performs end-to-end analysis of whole genome sequencing data, including quality control, alignment, variant calling, structural variant detection, and repeat expansion analysis. It's designed for germline samples and supports both local and cloud execution.

## Features

- **Quality Control**: FASTP for adapter trimming and quality control
- **BWA-MEM2 Alignment**: Fast and accurate read alignment
- **BAM Processing**: Sorting, indexing, duplicate marking, and format conversion
- **Quality Assessment**: MultiQC report aggregation, Qualimap BAM QC, Picard metrics
- **Variant Calling**: SNV/indel detection with DeepVariant
- **Structural Variants**: Detection with Manta
- **Copy Number Variants**: Analysis with CNVpytor
- **Repeat Expansions**: Detection with ExpansionHunter and ExpansionHunterDenovo
- **Format Conversion**: BAM to CRAM conversion for storage efficiency
- **Container Support**: Docker/Singularity containers for reproducibility
- **Cloud Ready**: S3 support for input/output files

## Architecture

```mermaid
graph TD
    A[FASTQ R1/R2] --> B[FASTP]
    B --> C[BWAMEM2_MEM]
    D[Reference FASTA] --> C
    E[BWA Index Files] --> C
    
    C --> F[SAMTOOLS_SORT]
    F --> G[PICARD_MARKDUPLICATES]
    G --> H[SAMTOOLS_INDEX]
    
    H --> I[PICARD_COLLECT_MULTIPLE_METRICS]
    H --> J[PICARD_COLLECT_WGS_METRICS]
    H --> K[QUALIMAP_BAMQC]
    
    B --> L[MULTIQC]
    I --> L
    J --> L
    K --> L
    
    H --> M[EXPANSIONHUNTER]
    H --> N[MANTA_GERMLINE]
    H --> O[CNVPYTOR]
    H --> P[EXPANSIONHUNTERDENOVO_PROFILE]
    H --> Q[SAMTOOLS_BAM2CRAM]
    
    R[Variant Catalog] --> M
    S[Min Anchor MapQ] --> P
    T[Max IRR MapQ] --> P
    
    subgraph "Input Data"
        A
        D
        E
        R
        S
        T
    end
    
    subgraph "Quality Control"
        B
    end
    
    subgraph "Alignment & Processing"
        C
        F
        G
        H
    end
    
    subgraph "Quality Assessment"
        I
        J
        K
        L
    end
    
    subgraph "Variant Analysis"
        M
        N
        O
        P
    end
    
    subgraph "Output Formats"
        Q
    end
    
    style A fill:#e1f5fe
    style D fill:#e1f5fe
    style E fill:#e1f5fe
    style R fill:#e1f5fe
    style S fill:#e1f5fe
    style T fill:#e1f5fe
    style L fill:#fff3e0
    style Q fill:#c8e6c9
```

## Quick Start

### Prerequisites

- Nextflow 22.04.0 or later
- Docker or Singularity
- AWS credentials (for S3 access)

### Basic Usage

```bash
nextflow run main.nf \
    --sample_name "Sample_001" \
    --fastq_r1 "path/to/sample_R1.fastq.gz" \
    --fastq_r2 "path/to/sample_R2.fastq.gz" \
    --fasta "path/to/reference.fasta" \
    --fai "path/to/reference.fasta.fai" \
    --variant_catalog "path/to/variant_catalog.json" \
    --min_anchor_mapq 50 \
    --max_irr_mapq 40
```

## Input Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--fastq_r1` | Forward reads FASTQ file | `sample_R1.fastq.gz` |
| `--fastq_r2` | Reverse reads FASTQ file | `sample_R2.fastq.gz` |
| `--fasta` | Reference genome FASTA | `hg38.fa` |
| `--fai` | FASTA index file | `hg38.fa.fai` |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--sample_name` | Sample identifier | File basename |
| `--variant_catalog` | ExpansionHunter variant catalog | `null` |
| `--min_anchor_mapq` | Min anchor mapping quality | `null` |
| `--max_irr_mapq` | Max IRR mapping quality | `null` |
| `--outdir` | Output directory | `/mnt/workflow/pubdir` |

## Output Structure

```
outdir/
├── hg38/
│   ├── QC/
│   │   ├── *_fastp.html
│   │   ├── *_fastp.json
│   │   ├── *_multiqc_report.html
│   │   ├── *_multiqc_data/
│   │   ├── *.metrics
│   │   ├── *.wgs_metrics
│   │   └── *_bamqc/
│   ├── ALIGNMENT/
│   │   ├── *.sorted.bam
│   │   ├── *.MarkDuplicates.bam
│   │   └── *.cram
│   ├── SNV/
│   │   ├── *.deepvariant.vcf.gz
│   │   └── *.deepvariant.gvcf.gz
│   ├── SV/
│   │   ├── *.candidate_sv.vcf.gz
│   │   └── *.diploid_sv.vcf.gz
│   └── REPEATS/
│       ├── *.vcf.gz
│       ├── *.json.gz
│       └── *.locus.tsv
```

## Processes

### Quality Control
- **FASTP**: Adapter trimming, quality filtering, and QC reports

### Core Alignment
- **BWAMEM2_MEM**: BWA-MEM2 alignment with read group information
- **SAMTOOLS_SORT**: BAM file sorting and indexing
- **PICARD_MARKDUPLICATES**: Duplicate marking and removal
- **SAMTOOLS_INDEX**: BAM indexing for downstream tools

### Quality Assessment
- **PICARD_COLLECT_MULTIPLE_METRICS**: Comprehensive BAM metrics collection
- **PICARD_COLLECT_WGS_METRICS**: Whole genome sequencing metrics
- **QUALIMAP_BAMQC**: BAM quality control analysis
- **MULTIQC**: QC report aggregation from all tools

### Variant Analysis
- **DEEPVARIANT_RUNDEEPVARIANT**: SNV/indel calling (commented out)
- **MANTA_GERMLINE**: Structural variant detection
- **CNVPYTOR**: Copy number variant analysis
- **EXPANSIONHUNTER**: Repeat expansion detection
- **EXPANSIONHUNTERDENOVO_PROFILE**: De novo repeat detection

### Format Conversion
- **SAMTOOLS_BAM2CRAM**: BAM to CRAM conversion

## Configuration

### Docker Images

The pipeline uses the following Docker images (configure in `nextflow.config`):

- `fastp_docker`: FASTP quality control
- `bwa_docker`: BWA-MEM2 and Samtools
- `picard_docker`: Picard tools
- `qualimap_docker`: Qualimap BAM QC
- `multiqc_docker`: MultiQC report generation
- `deepvariant_docker`: DeepVariant
- `manta_docker`: Manta
- `expansionhunter_docker`: ExpansionHunter
- `expansionhunterdenovo_docker`: ExpansionHunterDenovo
- `cnvpytor_docker`: CNVpytor

### Resource Requirements

| Process | Memory | CPUs |
|---------|--------|------|
| FASTP | 16 GB | 8 |
| BWAMEM2_MEM | 32 GB | 16 |
| SAMTOOLS_SORT | 16 GB | 8 |
| PICARD_MARKDUPLICATES | 32 GB | 16 |
| PICARD_COLLECT_MULTIPLE_METRICS | 16 GB | 8 |
| PICARD_COLLECT_WGS_METRICS | 16 GB | 8 |
| QUALIMAP_BAMQC | 16 GB | 8 |
| MULTIQC | 4 GB | 4 |
| MANTA_GERMLINE | 16 GB | 8 |
| EXPANSIONHUNTER | 16 GB | 8 |
| CNVPYTOR | 32 GB | 16 |

## Dependencies

### BWA Index Files

The pipeline requires BWA-MEM2 index files for the reference genome:

```bash
# Generate BWA-MEM2 index
bwa-mem2 index /path/to/reference.fasta
```

This creates:
- `reference.fasta.0123` (binary sequence file)
- `reference.fasta.amb` (amb file)
- `reference.fasta.ann` (ann file)
- `reference.fasta.bwt.2bit.64` (bwt file)
- `reference.fasta.pac` (pac file)

## Troubleshooting

### Common Issues

1. **Missing BWA Index Files**
   ```
   Error: Missing required BWA index file: reference.fasta.amb
   ```
   **Solution**: Generate BWA-MEM2 index files using `bwa-mem2 index`

2. **S3 Access Issues**
   ```
   Cannot find any reads matching: s3://bucket/file.fastq.gz
   ```
   **Solution**: Ensure AWS credentials are configured and S3 permissions are set

3. **Memory Issues**
   ```
   samtools sort: couldn't allocate memory for bam_mem
   ```
   **Solution**: Increase memory allocation in `nextflow.config`

4. **MULTIQC No Files Found**
   ```
   No analysis results found. MultiQC will now exit.
   ```
   **Solution**: Check that upstream QC processes completed successfully

### Log Files

Check the `.nextflow.log` file for detailed execution logs and error messages.

## Citation

If you use this pipeline in your research, please cite:

- Nextflow: Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- BWA-MEM2: Vasimuddin, M. et al. (2019). Efficient architecture-aware acceleration of BWA-MEM for multicore systems. IEEE IPDPS.
- FASTP: Chen, S. et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.
- MultiQC: Ewels, P. et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.
- DeepVariant: Poplin, R. et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983-987.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Support

For questions and support, please open an issue on the GitHub repository.
