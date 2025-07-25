# Changelog

All notable changes to the nl-wgs_wf pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial pipeline development
- BWA-MEM2 alignment module
- Samtools sorting and indexing modules
- Picard MarkDuplicates module
- Manta structural variant calling module
- CNVpytor copy number variant analysis module
- ExpansionHunter repeat expansion detection module
- ExpansionHunterDenovo de novo repeat detection module
- DeepVariant SNV/indel calling module (commented out)
- BAM to CRAM conversion module
- S3 support for input/output files
- Docker container support for all tools
- Comprehensive error handling and validation
- Memory optimization for large datasets
- Mermaid architecture diagram in README

### Changed
- Updated BWA-MEM2 process to use optimized memory allocation
- Improved samtools sort command with memory limits
- Enhanced input/output channel handling for better data flow
- Updated process configurations for better resource management
- Refined error messages and validation checks

### Fixed
- Syntax errors in all module files
- Input/output cardinality mismatches between processes
- Undefined variable issues in CNVPYTOR process
- Missing parameter definitions in nextflow.config
- Tab character handling in BWA-MEM2 read group strings
- Memory allocation issues in samtools sort
- File path handling for S3 and local files
- Process output file naming consistency

### Technical Details

#### Process Updates
- **BWAMEM2_MEM**: Added read group information, optimized memory usage
- **SAMTOOLS_SORT**: Added memory limits and temporary directory support
- **PICARD_MARKDUPLICATES**: Updated to handle new input format
- **SAMTOOLS_INDEX**: Added as separate process for BAM indexing
- **CNVPYTOR**: Fixed variable references and added error handling
- **EXPANSIONHUNTER**: Updated input/output structure
- **MANTA_GERMLINE**: Enhanced with proper reference handling
- **EXPANSIONHUNTERDENOVO_PROFILE**: Added parameter validation

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