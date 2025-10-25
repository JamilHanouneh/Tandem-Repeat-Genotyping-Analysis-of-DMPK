# Tandem-Repeat-Genotyping-Analysis-of-DMPK
Nextflow pipeline for automated genotyping of tandem repeat expansions from PacBio HiFi data using TRGT, targeting 20 pathogenic loci. Includes Python/R tools for parsing, QC, and visualization of DMPK repeats and motif interruptions, aiding clinical research on repeat expansion disorders.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-276DC3.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Quick Start](#quick-start)
- [Dataset Information](#dataset-information)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [Analysis Scripts](#analysis-scripts)
- [Output Structure](#output-structure)
- [Configuration](#configuration)
- [Reproducibility](#reproducibility)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Overview

This pipeline performs comprehensive tandem repeat genotyping for 20 pathogenic repeat expansion loci associated with neurological disorders, muscular dystrophies, and other genetic diseases. It processes PacBio HiFi sequencing data from the PureTarget Repeat Expansion Panel, which targets disease loci including:

- Huntington's disease (HTT)
- Fragile X syndrome (FMR1)
- Friedreich's ataxia (FXN)
- Myotonic dystrophy (DMPK)
- Spinocerebellar ataxias (ATXN1, ATXN3, etc.)

The project includes:
- Nextflow DSL2 pipeline for TRGT genotyping
- Python scripts for VCF parsing and DMPK repeat analysis
- R scripts for quality control and statistical analysis
- Jupyter notebooks for interactive exploration
- Comprehensive visualization tools

---

## Features

- Full-length repeat sequencing: Accurately genotypes repeat expansions up to 4,000+ units
- Automated Nextflow workflow: Parallel processing with automatic retry logic
- Methylation profiling: Native detection of CpG methylation without bisulfite conversion
- Somatic mosaicism detection: Identifies cell-to-cell variation in repeat lengths
- Multi-language analysis: Python, R, and Jupyter notebook implementations
- Clinical-grade output: VCF files compatible with clinical variant databases
- Interactive visualization: Bar charts, statistical plots, and quality control dashboards
- Motif interruption detection: Identifies non-canonical repeat sequences affecting disease severity

---

## Quick Start

```
# Clone the repository
git clone https://github.com/JamilHanouneh/trgt-nextflow-pipeline.git
cd trgt-nextflow-pipeline

# Download dataset (WARNING: Large files ~50GB)
bash scripts/download_dataset.sh

# Run Nextflow pipeline
cd NextFlow
nextflow run trgt_workflow.nf -profile conda

# Analyze DMPK repeats (Python)
python analyze_dmpk_repeat_expansions.py

# Or use R analysis
Rscript parse_trgt_vcfs_for_dmpk_repeats.R

# Or explore interactively
jupyter notebook analyze_dmpk_repeat_expansions.ipynb
```

---

## Dataset Information

### PacBio PureTarget Repeat Expansion Panel

This pipeline processes data from the PacBio PureTarget Repeat Expansion Panel, which uses CRISPR-Cas9-based amplification-free enrichment.

#### IMPORTANT: Dataset NOT Included in Repository

**WARNING**: Due to file size constraints (50+ GB), the dataset is NOT included in this repository. You must download it separately.

#### Download Instructions

The dataset should be placed in the `Dataset/` directory at the project root.

**Method 1: Browser Download**

1. Navigate to: https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/
2. Download the following directories:
   - `PBMM2-BAM-Input-For-IGV-And-TRGT/` (aligned BAM files)
   - `reference/` (GRCh38 genome and repeat catalog)
   - `TRGT-VCF-files/` (pre-generated VCF files, optional)
3. Place files in the `Dataset/` directory

**Method 2: Command Line (wget)**

```
cd Dataset
wget -r -np -nH --cut-dirs=4 \
  https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/PBMM2-BAM-Input-For-IGV-And-TRGT/

wget -r -np -nH --cut-dirs=4 \
  https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/reference/
```

**Method 3: AWS CLI**

```
cd Dataset
aws s3 sync s3://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/ . --no-sign-request
```

#### Verification

```
# Check BAM file count
ls Dataset/PBMM2-BAM-Input-For-IGV-And-TRGT/*.bam | wc -l
# Expected: 20 files

# Verify reference genome
samtools faidx Dataset/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

#### Dataset Characteristics

- 20 Coriell samples from lymphoblastoid cell lines
- 20 disease loci including AR, ATN1, ATXN1, ATXN3, DMPK, FMR1, FXN, HTT, PABPN1
- Pre-aligned BAM files mapped to GRCh38 using pbmm2
- Deep coverage from 30-hour Sequel IIe sequencing runs
- Total size: Approximately 50-60 GB

---

## Installation

### Prerequisites

- Nextflow >= 20.04.0
- Conda or Mamba
- Java >= 8
- Python >= 3.8
- R >= 4.0 (for R scripts)
- Disk space: 100+ GB
- RAM: 16 GB minimum, 32 GB recommended
- CPUs: 4+ cores recommended

### Install Dependencies

**Nextflow**

```
conda install -c bioconda nextflow
```

**Python Environment**

```
conda env create -f environment.yml
conda activate trgt-pipeline
```

**Or install manually:**

```
pip install pandas matplotlib vcfpy numpy seaborn jupyter
```

**R Packages**

```
Rscript -e "install.packages(c('tidyverse', 'vcfR', 'ggplot2'))"
```

---

## Project Structure

```
trgt-nextflow-pipeline/
├── NextFlow/                              # Nextflow pipeline directory
│   ├── trgt_workflow.nf                   # Main Nextflow workflow
│   └── nextflow.config                    # Pipeline configuration
│
├── Dataset/                               # Data directory (empty - download separately)
│   ├── PBMM2-BAM-Input-For-IGV-And-TRGT/ # BAM files (not included)
│   ├── reference/                         # Reference genome (not included)
│   └── TRGT-VCF-files/                   # Pre-generated VCFs (not included)
│
├── Results/                               # Analysis results directory
│   ├── analyze_dmpk_repeat_expansions/    # Python analysis outputs
│   │   ├── plots/                         # Visualization plots
│   │   └── reports/                       # Summary reports
│   └── parse_trgt_vcfs_for_dmpk_repeats/  # R analysis outputs
│       ├── figures/                       # R-generated figures
│       └── tables/                        # Statistical tables
│
├── analyze_dmpk_repeat_expansions.py      # Python VCF parsing script
├── analyze_dmpk_repeat_expansions.ipynb   # Jupyter notebook version
├── parse_trgt_vcfs_for_dmpk_repeats.py    # Python DMPK analysis
├── parse_trgt_vcfs_for_dmpk_repeats.R     # R DMPK analysis
├── DMPK_QC_Analysis.R                     # R quality control script

```

---

## Usage

### Running the Nextflow Pipeline

```
cd NextFlow

# Basic execution
nextflow run trgt_workflow.nf -profile conda

# With custom parameters
nextflow run trgt_workflow.nf -profile conda \
  --bam_dir ../Dataset/PBMM2-BAM-Input-For-IGV-And-TRGT \
  --ref_genome ../Dataset/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --repeats_bed ../Dataset/reference/pathogenic_repeats.hg38.bed \
  --out_dir ../Results/nextflow_output \
  --trgt_cpus 8

# Resume failed runs
nextflow run trgt_workflow.nf -profile conda -resume
```

### Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bam_dir` | `../Dataset/PBMM2-BAM-Input-For-IGV-And-TRGT` | Directory containing BAM files |
| `--ref_genome` | `../Dataset/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` | Reference genome |
| `--repeats_bed` | `../Dataset/reference/pathogenic_repeats.hg38.bed` | Repeat catalog |
| `--out_dir` | `results` | Output directory |
| `--trgt_cpus` | `2` | CPU threads for TRGT |

---

## Analysis Scripts

### 1. Python VCF Parsing: analyze_dmpk_repeat_expansions.py

Parses TRGT VCF files to extract DMPK repeat lengths and classify expansions.

```
python analyze_dmpk_repeat_expansions.py
```

**Features:**
- Direct VCF parsing using vcfpy
- Extracts allele lengths from FORMAT fields
- Detects motif interruptions (CCG, CTC in CTG repeats)
- Generates bar charts with pathogenic threshold lines
- Exports results to CSV

**Output:**
- `Results/analyze_dmpk_repeat_expansions/dmpk_analysis.png`
- `Results/analyze_dmpk_repeat_expansions/dmpk_data.csv`

### 2. Python DMPK Analysis: parse_trgt_vcfs_for_dmpk_repeats.py

Alternative implementation with enhanced statistical analysis.

```
python parse_trgt_vcfs_for_dmpk_repeats.py
```

**Features:**
- Multi-sample VCF processing
- Statistical summaries (mean, median, standard deviation)
- Allele distribution analysis
- Clinical classification (Normal vs. Pathogenic)

### 3. R VCF Analysis: parse_trgt_vcfs_for_dmpk_repeats.R

R-based analysis with ggplot2 visualizations.

```
Rscript parse_trgt_vcfs_for_dmpk_repeats.R
```

**Features:**
- vcfR package for VCF parsing
- ggplot2 publication-quality plots
- Genotype-phenotype correlation analysis
- Statistical testing (t-tests, ANOVA)

**Output:**
- `Results/parse_trgt_vcfs_for_dmpk_repeats/dmpk_boxplot.pdf`
- `Results/parse_trgt_vcfs_for_dmpk_repeats/summary_statistics.csv`

### 4. R Quality Control: DMPK_QC_Analysis.R

Comprehensive quality control for DMPK genotyping.

```
Rscript DMPK_QC_Analysis.R
```

**Features:**
- Coverage depth analysis
- Read quality metrics
- Allele balance assessment
- Motif purity scoring
- Batch effect detection

**Output:**
- `Results/QC/coverage_plots.pdf`
- `Results/QC/quality_metrics.csv`

### 5. Interactive Jupyter Notebook: analyze_dmpk_repeat_expansions.ipynb

Exploratory analysis with interactive visualizations.

```
jupyter notebook analyze_dmpk_repeat_expansions.ipynb
```

**Features:**
- Step-by-step walkthrough
- Interactive plots with Plotly
- Parameter tuning
- Custom threshold testing

---

## Output Structure

### Nextflow Pipeline Output

```
Results/nextflow_output/
├── reference_index/
│   └── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
├── reindexed_bam/
│   ├── NA03697B2.pbmm2.repeats.bam
│   ├── NA03697B2.pbmm2.repeats.bam.bai
│   └── ... ( samples)
└── trgt/
    ├── NA03697B2.vcf.gz          # Genotype calls
    ├── NA03697B2.spanning.bam    # Reads spanning repeats
    ├── NA03697B2.trgt.log        # Execution log
    └── ... ( samples)
```

### Analysis Script Outputs

```
Results/
├── analyze_dmpk_repeat_expansions/
│   ├── dmpk_analysis.png
│   ├── dmpk_data.csv
│   └── motif_interruption_report.txt
├── parse_trgt_vcfs_for_dmpk_repeats/
│   ├── dmpk_boxplot.pdf
│   ├── summary_statistics.csv
│   └── allele_distribution.png
└── QC/
    ├── coverage_plots.pdf
    ├── quality_metrics.csv
    └── batch_effects.png
```

---

## Configuration

### nextflow.config

Located in `NextFlow/nextflow.config`:

```
profiles {
  conda {
    conda.enabled = true
    conda.channelPriority = 'strict'
  }
}

process {
  withName: RUN_TRGT {
    conda = 'bioconda::trgt bioconda::samtools'
    cpus = 4
    memory = '16 GB'
  }
  withName: INDEX_GENOME {
    conda = 'bioconda::samtools'
    cpus = 1
    memory = '4 GB'
  }
  withName: ENSURE_BAI {
    conda = 'bioconda::samtools'
    cpus = 1
    memory = '4 GB'
  }
}
```

### Python Configuration

Modify paths at the top of Python scripts:

```
# analyze_dmpk_repeat_expansions.py
vcf_dir = "Dataset/TRGT-VCF-files"
output_dir = "Results/analyze_dmpk_repeat_expansions"
```

### R Configuration

Modify paths in R scripts:

```
# parse_trgt_vcfs_for_dmpk_repeats.R
vcf_dir <- "Dataset/TRGT-VCF-files"
output_dir <- "Results/parse_trgt_vcfs_for_dmpk_repeats"
```

---

## Reproducibility

### Software Versions

| Tool | Version | Source |
|------|---------|--------|
| Nextflow | >= 20.04.0 | bioconda |
| TRGT | Latest | bioconda |
| samtools | >= 1.15 | bioconda |
| Python | >= 3.8 | conda-forge |
| pandas | >= 1.5.0 | conda-forge |
| vcfpy | >= 0.13.6 | bioconda |
| R | >= 4.0 | conda-forge |
| vcfR | >= 1.13.0 | CRAN |

### Workflow Execution Report

Generate comprehensive execution reports:

```
cd NextFlow
nextflow run trgt_workflow.nf -profile conda \
  -with-report execution_report.html \
  -with-timeline timeline.html \
  -with-dag flowchart.png
```

---

## Troubleshooting

### Common Issues

#### Issue: "Dataset directory not found"

**Cause**: Dataset not downloaded

**Solution**:
```
mkdir -p Dataset
cd Dataset
wget -r -np -nH --cut-dirs=4 \
  https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/PBMM2-BAM-Input-For-IGV-And-TRGT/
```

#### Issue: "No DMPK data found in VCFs"

**Cause**: VCF directory path incorrect or files missing

**Solution**:
1. Check VCF directory exists: `ls Dataset/TRGT-VCF-files/*.vcf`
2. Verify DMPK locus in VCFs: `grep DMPK Dataset/TRGT-VCF-files/*.vcf`
3. Update path in analysis scripts

#### Issue: "Module vcfpy not found"

**Cause**: Python dependencies not installed

**Solution**:
```
pip install vcfpy pandas matplotlib
```

#### Issue: "R package vcfR not available"

**Cause**: R packages not installed

**Solution**:
```
Rscript -e "install.packages('vcfR', repos='https://cran.r-project.org')"
```

#### Issue: "Nextflow pipeline fails with memory error"

**Cause**: Insufficient RAM allocation

**Solution**:
Edit `NextFlow/nextflow.config`:
```
process {
  withName: RUN_TRGT {
    memory = '32 GB'  // Increase from 16 GB
  }
}
```

---


## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Quick Contribution Guide

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Make changes and test
4. Commit: `git commit -m "Add feature X"`
5. Push: `git push origin feature/your-feature`
6. Open a Pull Request

---

## Citation

If you use this pipeline in your research, please cite:

### BibTeX

```
@software{trgt_pipeline,
  author = {Jamil Hanouneh},
  title = {TRGT Tandem Repeat Genotyping Nextflow Pipeline},
  year = {2025},
  url = {https://github.com/JamilHanouneh/trgt-nextflow-pipeline},
  version = {1.0.0}
}

@article{dolzhenko2024trgt,
  title={Characterization and visualization of tandem repeats at genome scale},
  author={Dolzhenko, Egor and others},
  journal={Nature Biotechnology},
  year={2024},
  doi={10.1038/s41587-023-02057-3}
}
```

### APA

Jamil Hanouneh. (2025). TRGT Tandem Repeat Genotyping Nextflow Pipeline (Version 1.0.0) [Computer software]. https://github.com/JamilHanouneh/trgt-nextflow-pipeline

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- PacBio for developing TRGT and providing the PureTarget dataset
- Nextflow community for DSL2 framework
- Bioconda for bioinformatics software distribution
- Coriell Institute for disease-positive cell line samples

---

## Contact

**Maintainer**: Jamil Hanouneh  
**Email**: jamil.hanouneh1997@gmail.com  
**Institution**: Friedrich-Alexander-Universität Erlangen-Nürnberg
**Issues**: https://github.com/JamilHanouneh/trgt-nextflow-pipeline/issues

---

## Roadmap

### Current Version (v1.0.0)
- Nextflow TRGT genotyping pipeline
- Python and R analysis scripts
- Jupyter notebook for interactive analysis
- DMPK-focused repeat expansion analysis

### Planned Features (v1.1.0)
- Multi-locus analysis dashboard
- Interactive HTML reports
- Automated methylation analysis
- Integration with TRVZ visualizer

### Future Versions (v2.0.0)
- Support for custom repeat catalogs
- Machine learning-based phenotype prediction
- Clinical variant database integration (ClinVar)
- Oxford Nanopore data support
```

***

## 2. .gitignore

```gitignore
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Virtual Environments
venv/
ENV/
env/
.venv

# Jupyter Notebook
.ipynb_checkpoints
*.ipynb_checkpoints

# R
.Rhistory
.RData
.Rapp.history
*.Rproj
.Rproj.user/

# Dataset directory (CRITICAL: Large files not in repo)
Dataset/
Dataset/**/*
!Dataset/.gitkeep
!Dataset/README.md

# Results and outputs
Results/
Results/**/*
!Results/.gitkeep
!Results/README.md

# Nextflow
NextFlow/work/
NextFlow/.nextflow/
NextFlow/.nextflow.log*
NextFlow/timeline.html
NextFlow/report.html
NextFlow/trace.txt
NextFlow/dag.dot

# Data files
*.bam
*.bai
*.vcf
*.vcf.gz
*.vcf.gz.tbi
*.fasta
*.fna
*.fa
*.fai
*.bed
*.zip
*.tar.gz
*.tar
*.gz

# Plot outputs
*.png
*.pdf
*.svg
*.jpg
*.jpeg

# IDE
.vscode/
.idea/
*.swp
*.swo
*~

# OS
.DS_Store
Thumbs.db
desktop.ini

# Environment files
.env
.env.local

# Testing
.pytest_cache/
.coverage
htmlcov/
.tox/
testthat/

# BUT: Keep directory structure markers
!.gitkeep
```

***

## 3. LICENSE

```
MIT License

Copyright (c) [Year] Jamil Hanouneh

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

***

## 4. CONTRIBUTING.md

```markdown
# Contributing to TRGT Nextflow Pipeline

Thank you for your interest in contributing!

## How to Contribute

### Reporting Bugs

1. Check existing issues to avoid duplicates
2. Use the bug report template
3. Include:
   - Nextflow version
   - Python/R version
   - Operating system
   - Exact error message
   - Minimal reproducible example

### Suggesting Enhancements

1. Open an issue with the enhancement label
2. Describe the use case
3. Provide examples

### Pull Request Process

#### 1. Fork and Clone

```
git clone https://github.com/JamilHanouneh/trgt-nextflow-pipeline.git
cd trgt-nextflow-pipeline
git checkout -b feature/your-feature-name
```

#### 2. Development Setup

```
conda env create -f environment.yml
conda activate trgt-pipeline
```

#### 3. Make Changes

- Follow code style guidelines
- Add tests for new features
- Update documentation


```

#### 5. Commit and Push

```
git add .
git commit -m "feat: add new feature description"
git push origin feature/your-feature-name
```

#### 6. Create Pull Request

- Use descriptive title
- Reference related issues
- Describe changes in detail

```

***

## 6. CITATION.cff

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
type: software
title: "TRGT Tandem Repeat Genotyping Nextflow Pipeline"
version: 1.0.0
date-released: 2025-10-25
authors:
  - family-names: Hanouneh
    given-names: Jamil
    orcid: ["https://orcid.org/[your-orcid]"](https://orcid.org/0009-0005-5133-4929)
    affiliation: Friedrich-Alexander-Universität Erlangen-Nürnberg
    email: jamil.hanouneh1997@gmail.com
abstract: >
  A comprehensive Nextflow DSL2 pipeline for automated genotyping of pathogenic
  tandem repeat expansions from PacBio HiFi sequencing data using TRGT. Includes
  Python and R analysis scripts for DMPK repeat expansion analysis, quality control,
  and clinical interpretation for Myotonic Dystrophy Type 1.
keywords:
  - tandem repeats
  - nextflow
  - genomics
  - pacbio
  - hifi sequencing
  - trgt
  - myotonic dystrophy
  - dmpk
  - repeat expansions
  - bioinformatics
  - vcf analysis
  - python
  - r
license: MIT
```

***


### Planned
- Multi-locus analysis dashboard
- Interactive HTML reports
- Automated methylation profiling
- TRVZ visualization integration

## [1.0.0] - 2025-10-25

### Added
- Nextflow DSL2 pipeline (trgt_workflow.nf)
- Python VCF parsing script (analyze_dmpk_repeat_expansions.py)
- Alternative Python analysis (parse_trgt_vcfs_for_dmpk_repeats.py)
- R VCF analysis script (parse_trgt_vcfs_for_dmpk_repeats.R)
- R quality control script (DMPK_QC_Analysis.R)
- Jupyter notebook for interactive analysis
- Comprehensive README with dataset download instructions
- MIT license
- Contributing guidelines
- Code of conduct
- Citation file (CFF format)

### Features
- Processes PacBio PureTarget Repeat Expansion Panel data
- Supports 20 pathogenic repeat loci
- DMPK-focused repeat expansion analysis
- Motif interruption detection
- Clinical threshold classification (50 repeats)
- Multi-language support (Python, R, Jupyter)
- Parallel sample processing with Nextflow
- Automatic BAM reindexing
- VCF, spanning BAM, and log outputs

### Documentation
- Complete README with quick start guide
- Dataset download instructions 
- Project structure documentation
- Troubleshooting section
- Analysis script documentation
```

***

## 8. Dataset/README.md

```markdown
# Dataset Directory

## Download Instructions

### Full Dataset Download

Visit: https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/

Download these directories:
1. `PBMM2-BAM-Input-For-IGV-And-TRGT/` - Aligned BAM files (required for Nextflow pipeline)
2. `reference/` - GRCh38 reference genome and repeat catalog
3. `TRGT-VCF-files/` - Pre-generated VCF files (optional for analysis scripts)

### Quick Download (Command Line)

```
cd Dataset

# Download BAM files
wget -r -np -nH --cut-dirs=4 \
  https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/PBMM2-BAM-Input-For-IGV-And-TRGT/

# Download reference
wget -r -np -nH --cut-dirs=4 \
  https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/reference/

# Download VCFs (optional)
wget -r -np -nH --cut-dirs=4 \
  https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/TRGT-VCF-files/
```

## Directory Structure

```
Dataset/
├── PBMM2-BAM-Input-For-IGV-And-TRGT/
│   ├── NA03697B2.pbmm2.repeats.bam
│   ├── NA03697B2.pbmm2.repeats.bam.bai
│   └── ... (samples)
└── TRGT-VCF-files/
    ├── NA03697B2.vcf
    └── ... (samples)
```

## Verification

```
# Check file counts
ls PBMM2-BAM-Input-For-IGV-And-TRGT/*.bam | wc -l  
ls TRGT-VCF-files/*.vcf | wc -l                     

# Verify integrity
samtools quickcheck -v PBMM2-BAM-Input-For-IGV-And-TRGT/*.bam
```
## Citation

PacBio PureTarget Repeat Expansion Panel
https://www.pacb.com/technology/puretarget/
```

***

## 9. Results/README.md

```markdown
# Results Directory
## Note
Results are generated when you run:
- `python analyze_dmpk_repeat_expansions.py`
- `Rscript parse_trgt_vcfs_for_dmpk_repeats.R`
- `Rscript DMPK_QC_Analysis.R`
```

***

## 10. Repository Metadata

**Repository Name**: `trgt-dmpk-analysis-pipeline`

**Short Description**:  
"Nextflow pipeline for TRGT tandem repeat genotyping with Python/R analysis tools for DMPK repeat expansions in Myotonic Dystrophy research"

**Topics/Tags**:
```
nextflow
pacbio
hifi-sequencing
trgt
dmpk
myotonic-dystrophy
tandem-repeats
repeat-expansions
genomics
bioinformatics
python
r-programming
vcf-analysis
clinical-genomics
jupyter-notebook
data-analysis
pipeline
workflow
dsl2
repeat-genotyping
```

***

## 11. requirements

```
pandas>=1.5.0
matplotlib>=3.5.0
vcfpy>=0.13.6
numpy>=1.21.0
seaborn>=0.11.0
jupyter>=1.0.0
plotly>=5.0.0
scipy>=1.7.0
```
