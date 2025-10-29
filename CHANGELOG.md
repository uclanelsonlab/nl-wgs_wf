# Changelog

All notable changes to the nl-wgs_wf pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Performance optimizations for large-scale datasets
- Additional variant calling tools integration
- Enhanced reporting and visualization features
- Support for additional cloud platforms

## [1.1.0] - 2025-01-29

**Release Date**: January 29, 2025  
**Status**: Production Ready  
**Breaking Changes**: None

### Added
- **HapCUT2 Phasing Module**: Complete haplotype phasing workflow for heterozygous variants
  - `HAPCUT2_EXTRACTHAIRS`: Extract haplotype-informative reads from BAM files
  - `HAPCUT2_HAPCUT2`: Assemble haplotype blocks and phase variants
- **BCFtools Module**: VCF filtering and preprocessing for phasing
  - `BCFTOOLS_VIEW_DIPLOID`: Filter VCF for diploid genotypes (0/0, 0/1, 1/1)
- **HAPCUT2_PHASING Subworkflow**: Orchestrates the complete phasing pipeline
  - Integrates BCFtools filtering, read extraction, and haplotype phasing
  - Processes DeepVariant VCF output for phase determination
- **Subworkflow Architecture**: Modular pipeline organization
  - `FASTQ_PROCESSING`: Handle FASTQ input with FASTP, alignment, and sorting
  - `BAM_PROCESSING`: Handle BAM input with sorting
  - `CRAM_PROCESSING`: Handle CRAM input with conversion to BAM
  - `COMMON_ANALYSIS`: Shared downstream analysis for all input types
- **Docker Support**: New container parameters
  - `bcftools_docker`: BCFtools container for VCF operations
  - `hapcut2_docker`: HapCUT2 container for phasing
- **Documentation**: Comprehensive updates
  - Updated architecture mermaid diagram showing subworkflow structure
  - Added HapCUT2 and BCFtools process documentation
  - Added phasing output files to output structure documentation
  - Added resource requirements for new processes
  - Added citations for HapCUT2 and BCFtools

### Changed
- **Pipeline Architecture**: Reorganized into input-specific subworkflows
  - Main workflow now routes to FASTQ, BAM, or CRAM processing based on input type
  - All paths converge to COMMON_ANALYSIS for downstream processing
- **Output Structure**: SNV directory now includes phasing outputs
  - `*.diploid.vcf`: Filtered diploid variants
  - `*.fragment_file`: Haplotype-informative read fragments
  - `*.haplotype_output_file*`: Phased haplotype blocks

### Technical Details

#### New Processes
- **BCFTOOLS_VIEW_DIPLOID**:
  - Resources: 16GB memory, 8 CPUs, 2h runtime
  - Filters VCF using bcftools view with genotype criteria
  - Publishes to SNV directory
  
- **HAPCUT2_EXTRACTHAIRS**:
  - Resources: 16GB memory, 8 CPUs, 2h runtime
  - Extracts reads using extractHAIRS tool
  - Requires BAM, VCF, and reference genome
  
- **HAPCUT2_HAPCUT2**:
  - Resources: 16GB memory, 8 CPUs, 2h runtime
  - Performs phasing using HAPCUT2 algorithm
  - Generates haplotype block output

#### Subworkflow Organization
- **FASTQ_PROCESSING**: FASTP → BWAMEM2_MEM → SAMTOOLS_SORT
- **BAM_PROCESSING**: SAMTOOLS_SORT (re-sort for consistency)
- **CRAM_PROCESSING**: SAMTOOLS_CRAM2BAM
- **HAPCUT2_PHASING**: BCFTOOLS_VIEW_DIPLOID → HAPCUT2_EXTRACTHAIRS → HAPCUT2_HAPCUT2
- **COMMON_ANALYSIS**: All downstream QC, variant calling, and analysis

#### Configuration Changes
- Added `bcftools_docker` and `hapcut2_docker` parameters to `nextflow.config`
- Added process configurations with labels `bcftools_view_diploid` and `hapcut2`
- Added resource allocations for new processes
- Updated `parameters.json` with new Docker image parameters

#### Workflow Integration
- HapCUT2 phasing runs after DeepVariant in COMMON_ANALYSIS subworkflow
- Takes indexed BAM and DeepVariant VCF as inputs
- Integrates seamlessly with existing variant calling pipeline
- Outputs publish to SNV directory alongside DeepVariant results

## [1.0.1] - 2025-01-27

**Release Date**: January 27, 2025  
**Commit**: 7280316  
**Status**: Production Ready  
**Breaking Changes**: Output directory structure simplified (removed genome subdirectories)

### Added
- `samtools_docker` parameter for dedicated SAMtools container support
- Proper SAMtools container configuration for BAM to CRAM conversion

### Changed
- Increased Manta resources from 16GB/8CPUs to 32GB/16CPUs for better performance
- Simplified output directory structure by removing `hg38/` subdirectories
- Updated AutoMap parameter ordering for proper execution

### Fixed
- AutoMap parameter ordering issue that prevented proper execution
- SAMtools container configuration (was incorrectly using BWA container)
- Output directory structure consistency across all processes

### Technical Details

#### Process Updates
- **MANTA_GERMLINE**: Increased memory allocation to 32GB and CPU allocation to 16
- **SAMTOOLS_BAM2CRAM**: Now uses dedicated `samtools_docker` container instead of `bwa_docker`
- **AUTOMAP**: Fixed parameter order in AutoMap_v1.3.sh call

#### Configuration Changes
- Added `samtools_docker` parameter to `nextflow.config`
- Updated all `publishDir.path` configurations to remove genome-specific subdirectories
- Enhanced resource management for improved pipeline stability

## [1.0.0] - 2025-01-27

**Release Date**: January 27, 2025  
**Commit**: a19016a  
**Status**: Production Ready  
**Breaking Changes**: None from 0.1.0 (complete rewrite)

### Added
- Initial pipeline development
- BWA-MEM2 alignment module
- Samtools sorting and indexing modules
- Picard MarkDuplicates module
- Manta structural variant calling module
- CNVpytor copy number variant analysis module
- ExpansionHunter repeat expansion detection module
- ExpansionHunterDenovo de novo repeat detection module
- DeepVariant SNV/indel calling module
- AutoMap variant annotation module
- BAM to CRAM conversion module
- FASTP quality control and adapter trimming module
- MultiQC quality control report aggregation
- Qualimap BAM quality control
- Picard CollectMultipleMetrics and CollectWgsMetrics
- S3 support for input/output files
- Docker container support for all tools
- Comprehensive error handling and validation
- Memory optimization for large datasets
- Mermaid architecture diagram in README
- Complete documentation and usage examples

### Changed
- Updated BWA-MEM2 process to use optimized memory allocation
- Improved samtools sort command with memory limits
- Enhanced input/output channel handling for better data flow
- Updated process configurations for better resource management
- Refined error messages and validation checks
- Updated workflow architecture to include quality control and assessment steps
- Improved output directory structure with dedicated QC folder

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
- BWA-MEM2 FASTQ input handling and pipeline integration
- Nextflow configuration parsing errors (unquoted strings)
- publishDir configuration syntax (removed closures)
- MULTIQC input channel structure and metadata handling

### Technical Details

#### Process Updates
- **BWAMEM2_MEM**: Added read group information, optimized memory usage, fixed FASTQ input handling
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
- **DEEPVARIANT_RUNDEEPVARIANT**: Integrated SNV/indel calling
- **AUTOMAP**: Added variant annotation capabilities

#### Configuration Changes
- Added S3 configuration for cloud execution
- Updated Docker image parameters for all tools
- Added memory and CPU resource allocations
- Implemented process-specific configurations
- Added parameter validation and defaults
- Fixed configuration parsing issues
- Corrected publishDir syntax for all processes

#### File Structure
- Organized modules in separate directories
- Standardized process naming conventions
- Added version tracking for all tools
- Implemented consistent output file naming
- Updated output directory structure with QC folder

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
- Configuration parsing errors with unquoted strings
- BWA-MEM2 input handling issues
- MULTIQC integration challenges

### Dependencies
- Nextflow 22.04.0+
- Docker/Singularity
- BWA-MEM2 index files for reference genome
- AWS credentials for S3 access

---

## Version History

### Version 1.1.0
- **Date**: 2025-01-29
- **Status**: Production Ready
- **Features**: Haplotype phasing with HapCUT2, subworkflow architecture, BCFtools integration
- **Breaking Changes**: None
- **Git Tag**: v1.1.0

### Version 1.0.1
- **Date**: 2025-01-27
- **Status**: Production Ready
- **Features**: Bug fixes and resource optimizations
- **Breaking Changes**: Output directory structure simplified
- **Git Tag**: v1.0.1
- **Commit**: 7280316

### Version 1.0.0
- **Date**: 2025-01-27
- **Status**: Production Ready
- **Features**: Complete WGS analysis pipeline with all major components
- **Breaking Changes**: None from 0.1.0 (complete rewrite)
- **Git Tag**: v1.0.0
- **Commit**: a19016a

### Version 0.1.0
- **Date**: 2025-07-24
- **Status**: Development
- **Features**: Core pipeline functionality
- **Known Issues**: Several syntax and configuration issues

---

## Migration Guide

### From Version 1.0.1 to 1.1.0

#### Breaking Changes
- None. This is a feature addition with full backward compatibility.

#### Required Updates
1. Add `bcftools_docker` parameter to `nextflow.config`
2. Add `hapcut2_docker` parameter to `nextflow.config`
3. Update `parameters.json` with new Docker image parameters (if using parameter file)
4. No changes required to existing workflows - phasing runs automatically after DeepVariant

#### New Features
- Haplotype phasing with HapCUT2
- BCFtools VCF filtering
- Subworkflow architecture for better modularity
- Support for FASTQ, BAM, and CRAM inputs through dedicated subworkflows
- Enhanced documentation with updated architecture diagrams

#### Optional Configuration
- HapCUT2 processes can be disabled by setting `task.ext.when = false` in config
- Phasing outputs publish to SNV directory by default
- No additional input parameters required

### From Version 0.1.0 to 1.0.0

#### Breaking Changes
- Complete pipeline rewrite with new architecture
- Updated input/output channel structures
- Changed process parameter names
- Modified Docker image configurations
- Added new quality control and assessment processes

#### Required Updates
1. Update `nextflow.config` with new parameter definitions
2. Generate BWA-MEM2 index files for reference genome
3. Configure AWS credentials for S3 access
4. Update Docker image references
5. Ensure all string parameters are properly quoted

#### New Features
- S3 file support
- Memory optimization
- Enhanced error handling
- Comprehensive documentation
- Quality control pipeline integration
- MultiQC report generation
- BAM quality assessment tools
- Complete variant analysis suite

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