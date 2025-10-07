# BAM/CRAM Input Support

This pipeline now supports starting from BAM or CRAM files in addition to FASTQ files. This is useful when you already have aligned data and want to skip the alignment steps.

## Input Types

The pipeline supports three input types controlled by the `input_type` parameter:

- `fastq`: Start from paired-end FASTQ files (default behavior)
- `bam`: Start from a BAM file
- `cram`: Start from a CRAM file

## Usage Examples

### Starting from FASTQ (original behavior)
```json
{
    "input_type": "fastq",
    "fastq_r1": "path/to/sample_R1.fastq.gz",
    "fastq_r2": "path/to/sample_R2.fastq.gz",
    "sample_name": "sample_name",
    ...
}
```

### Starting from BAM
```json
{
    "input_type": "bam",
    "bam_input": "path/to/sample.bam",
    "bam_index": "path/to/sample.bam.bai",
    "sample_name": "sample_name",
    ...
}
```

### Starting from CRAM
```json
{
    "input_type": "cram",
    "cram_input": "path/to/sample.cram",
    "cram_index": "path/to/sample.cram.crai",
    "sample_name": "sample_name",
    ...
}
```

## Workflow Behavior

### FASTQ Input (`input_type: "fastq"`)
1. FASTP (quality control and trimming)
2. BWA-MEM2 (alignment)
3. SAMTOOLS_SORT (sorting)
4. PICARD_MARKDUPLICATES (duplicate marking)
5. Continue with downstream analysis...

### BAM Input (`input_type: "bam"`)
1. SAMTOOLS_SORT (re-sorting, assumes input may not be properly sorted)
2. PICARD_MARKDUPLICATES (duplicate marking)
3. Continue with downstream analysis...

### CRAM Input (`input_type: "cram"`)
1. SAMTOOLS_CRAM2BAM (convert CRAM to BAM)
2. PICARD_MARKDUPLICATES (duplicate marking)
3. Continue with downstream analysis...

## Important Notes

1. **Reference Genome**: When using BAM/CRAM input, you still need to provide the reference genome files (`fasta` and `fai`) as they are required for downstream processes.

2. **Index Files**: 
   - BAM index files (`.bai`) and CRAM index files (`.crai`) are optional but recommended for better performance
   - If not provided, the pipeline will create index files as needed
   - Providing index files can speed up processing and reduce computational overhead

3. **Duplicate Marking**: The pipeline assumes that BAM/CRAM inputs may not have duplicates marked, so it runs PICARD_MARKDUPLICATES on all input types.

4. **Sorting**: BAM inputs are re-sorted to ensure proper coordinate sorting. CRAM inputs are converted to BAM format first.

5. **MultiQC Reports**: 
   - When starting from BAM/CRAM, FASTP metrics will not be included in the MultiQC report since read trimming is skipped
   - DeepVariant visual reports are automatically included in MultiQC for all input types

6. **CRAM to CRAM**: If you start with CRAM input, the pipeline will not convert the final BAM back to CRAM to avoid redundancy.

## Example Parameter Files

See the example parameter files:
- `run_parameters_bam_example.json` - Example for BAM input
- `run_parameters_cram_example.json` - Example for CRAM input

## Required Parameters by Input Type

### For FASTQ input:
- `input_type: "fastq"`
- `fastq_r1`
- `fastq_r2`
- `sample_name`
- `fasta`
- `fai`

### For BAM input:
- `input_type: "bam"`
- `bam_input`
- `bam_index` (optional)
- `sample_name`
- `fasta`
- `fai`

### For CRAM input:
- `input_type: "cram"`
- `cram_input`
- `cram_index` (optional)
- `sample_name`
- `fasta`
- `fai`
