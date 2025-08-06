# Changelog

All notable changes to the nl-wgs_wf pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### HPC Branch Development
- HPC-optimized resource allocation and configuration
- Local ECR registry integration for faster container access
- Streamlined workflow focusing on core analysis components
- Test data integration with NA12878 sample
- Cluster-friendly configuration for SLURM/PBS schedulers
- Reduced memory and CPU requirements for HPC environments

### Added
- Initial pipeline development
- BWA-MEM2 alignment module
- Samtools sorting and indexing modules
- Picard MarkDuplicates module
- FASTP quality control and adapter trimming module
- MultiQC quality control report aggregation
- Qualimap BAM quality control
- Picard CollectMultipleMetrics and CollectWgsMetrics
- DeepVariant SNV/indel calling module
- AutoMap variant annotation module
- Docker container support for all tools
- Comprehensive error handling and validation
- Memory optimization for HPC environments
- Mermaid architecture diagram in README

### Changed
- Updated BWA-MEM2 process to use HPC-optimized memory allocation (32GB/16CPUs)
- Improved samtools sort command with memory limits
- Enhanced input/output channel handling for better data flow
- Updated process configurations for HPC resource management
- Refined error messages and validation checks
- Streamlined workflow by commenting out non-essential processes (SV, CNV, repeats)
- Updated Docker image paths to use local ECR registry
- Changed output directory to `results/` for HPC environments

### Fixed
- Syntax errors in all module files
- Input/output cardinality mismatches between processes
- Undefined variable issues in CNVPYTOR process
- Missing parameter definitions in nextflow.config
- Tab character handling in BWA-MEM2 read group strings
- Memory allocation issues in samtools sort
- File path handling for S3 and local files
- Process output file naming consistency
- FASTP module output paths and workflow integration
- Missing expansionhunterdenovo Docker parameter
- QUALIMAP_BAMQC output directory capture and MULTIQC integration

### Technical Details

#### Process Updates
- **BWAMEM2_MEM**: Added read group information, optimized memory usage
- **SAMTOOLS_SORT**: Added memory limits and temporary directory support
- **PICARD_MARKDUPLICATES**: Updated to handle new input format
- **PICARD_COLLECT_MULTIPLE_METRICS**: Added comprehensive BAM metrics collection
- **PICARD_COLLECT_WGS_METRICS**: Added whole genome sequencing metrics
- **SAMTOOLS_INDEX**: Added as separate process for BAM indexing
- **CNVPYTOR**: Fixed variable references and added error handling
- **EXPANSIONHUNTER**: Updated input/output structure
- **MANTA_GERMLINE**: Enhanced with proper reference handling
- **EXPANSIONHUNTERDENOVO_PROFILE**: Added parameter validation
- **FASTP**: Fixed output paths, added cleaned FASTQ outputs, integrated with workflow
- **QUALIMAP_BAMQC**: Added BAM quality control analysis
- **MULTIQC**: Added QC report aggregation from all tools

#### Configuration Changes
- Added S3 configuration for cloud execution
- Updated Docker image parameters
- Added memory and CPU resource allocations
- Implemented process-specific configurations
- Added parameter validation and defaults

#### File Structure
- Organized modules in separate directories
- Standardized process naming conventions
- Added version tracking for all tools
- Implemented consistent output file naming

#### HPC-Specific Changes
- **Docker Configuration**: Added `docker.runOptions = '-u $(id -u):$(id -g)'` for HPC user permissions
- **Conda Integration**: Enabled conda for additional tool support
- **Scratch Directory**: Enabled scratch for temporary file management
- **Module Binaries**: Enabled module binaries for better HPC compatibility
- **Test Data**: Pre-configured with NA12878 sample data paths
- **Resource Optimization**: Reduced memory requirements across all processes

## [1.0.0-HPC] - 2025-01-27

### Added
- HPC-optimized version of the pipeline
- Local ECR registry integration
- Test data integration with NA12878 sample
- Streamlined workflow focusing on core analysis
- Cluster-friendly resource allocation

### Changed
- Reduced memory requirements for HPC environments
- Updated Docker image paths to ECR registry
- Streamlined workflow by removing non-essential processes
- Optimized for SLURM/PBS schedulers

### Technical Details
- **Memory Optimization**: Reduced from 192GB to 32GB for DeepVariant
- **CPU Optimization**: Reduced from 48 to 16 CPUs for most processes
- **Container Registry**: Uses `224597534425.dkr.ecr.us-west-2.amazonaws.com`
- **Test Sample**: NA12878 with ERR3239334 paired-end reads
- **Output Directory**: Changed to `results/` for HPC compatibility

## [0.1.0] - 2025-07-24

### Added
- Initial pipeline structure
- Basic Nextflow workflow setup
- Core alignment processes
- Variant calling framework

### Known Issues
- Input/output cardinality mismatches between some processes
- Memory allocation issues with large datasets
- S3 file access requires proper AWS configuration
- CNVPYTOR requires reference genome resource files

### Dependencies
- Nextflow 22.04.0+
- Docker/Singularity
- BWA-MEM2 index files for reference genome
- AWS credentials for S3 access

---

## Version History

### Version 0.1.0
- **Date**: 2025-07-24
- **Status**: Development
- **Features**: Core pipeline functionality
- **Known Issues**: Several syntax and configuration issues

### Unreleased
- **Date**: Ongoing
- **Status**: Active Development
- **Features**: Comprehensive fixes and improvements
- **Target**: Production-ready pipeline

---

## Migration Guide

### From Version 0.1.0 to Current

#### Breaking Changes
- Updated input/output channel structures
- Changed process parameter names
- Modified Docker image configurations

#### Required Updates
1. Update `nextflow.config` with new parameter definitions
2. Generate BWA-MEM2 index files for reference genome
3. Configure AWS credentials for S3 access
4. Update Docker image references

#### New Features
- S3 file support
- Memory optimization
- Enhanced error handling
- Comprehensive documentation

---

## Contributing to Changelog

When adding new entries to the changelog, please follow these guidelines:

1. **Use the existing format** and structure
2. **Group changes** by type (Added, Changed, Fixed, Removed)
3. **Provide clear descriptions** of what changed
4. **Include technical details** when relevant
5. **Add migration notes** for breaking changes
6. **Update version numbers** appropriately

### Change Types

- **Added**: New features, processes, or capabilities
- **Changed**: Modifications to existing functionality
- **Fixed**: Bug fixes and error corrections
- **Removed**: Deprecated or removed features
- **Technical Details**: Implementation specifics
- **Known Issues**: Current limitations or problems 