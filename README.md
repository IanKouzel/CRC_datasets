# Recent Colon Cancer Bulk RNAseq Datasets from SRA (2024-2025)

## Overview
This document provides information on the most recent bulk RNAseq datasets for colon cancer available in the NCBI Sequence Read Archive (SRA). All three datasets listed below are from 2024-2025 and contain high-quality transcriptomic data suitable for comprehensive cancer research.

---

## Dataset 1: PRJNA1304982 - A-to-I RNA Editing Across Consensus Molecular Subtypes

### Basic Information
- **BioProject ID**: PRJNA1304982
- **Registration Date**: August 12, 2025
- **Institution**: IRCCS-Istituto Tumori "Giovanni Paolo II" (Italy)
- **Funding**: "Tecnopolo per la medicina di Precisione" (Regione Puglia)
- **Data Type**: Raw sequence reads (bulk RNAseq)

### Project Description
Dissecting the biology of Consensus Molecular Subtypes (CMS) in colon cancer. The molecular subtyping of colon cancer through the four CMS groups this malignancy into biologically different layers. This project aims to dissect the intrinsic biological differences of these four groups to gain relevant information for better therapeutic management and identification of alternative targets.

### Dataset Specifications
- **Number of SRA Experiments**: 100
- **Number of BioSamples**: 100
- **Data Volume**: 
  - 981 Gbases
  - 0.33 Tbytes
- **Scope**: Multispecies
- **Relevance**: Medical

### Sequencing Details
- **Library Preparation**: TruSeq Stranded Total RNA Library Prep Gold (Illumina)
- **RNA Input**: 800 ng per sample
- **Quality Control**: 
  - RNA concentration: NanoDrop ND-1000 spectrophotometer
  - RNA quality: Agilent TapeStation 4200
- **Sequencing Provider**: Genomix4life S.R.L. (Baronissi, Salerno, Italy)

### Analysis Focus
- A-to-I (Adenosine-to-Inosine) RNA editing analysis
- Consensus Molecular Subtype classification
- Editing site identification and annotation
- Filtering using dbSNP and RepeatMasker tracks

### Access
- **SRA Link**: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=1304982
- **BioProject Link**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1304982

---

## Dataset 2: PRJNA815861 - Early-Onset vs Late-Onset Colorectal Cancer

### Basic Information
- **BioProject ID**: PRJNA815861
- **Year**: 2024
- **Institution**: University of Otago, Christchurch, New Zealand
- **Funding**: Health Research Council of New Zealand, Bowel Cancer Research Aotearoa
- **Data Type**: Raw sequence reads (bulk RNAseq)

### Project Description
Consensus molecular subtypes and gene expression in early-onset colorectal cancer. This study compares transcriptomic profiles between early-onset (≤50 years) and late-onset (>65 years) colorectal cancer patients to identify molecular differences and therapeutic targets.

### Dataset Specifications
- **Total Samples**: 300 colorectal cancer patients
  - **Early-Onset (≤50 years)**: 19 patients
  - **Late-Onset (>65 years)**: 196 patients
  - **Intermediate (51-65 years)**: 85 patients (excluded from main analysis)
- **Sample Type**: Treatment-naïve colorectal tumors
- **Collection Site**: Christchurch Hospital, Aotearoa New Zealand
- **Sample Preparation**: Surgical resections, frozen in liquid nitrogen, stored at -80°C

### Exclusion Criteria
- Chemotherapy prior to study
- Hereditary non-polyposis colorectal cancer (HNPCC)
- Familial adenomatous polyposis (FAP)

### Analysis Focus
- Consensus Molecular Subtype (CMS) classification
- Differential gene expression between age cohorts
- Identification of age-specific biomarkers
- 13 genes significantly more abundant in early-onset cohort (log2FC >1):
  - C1QTNF3, TPTEP1, STMN2, KIF25, DACH2, NIBAN3, SLC35F1, IGHD, HOXA11-AS, LINC00926, LINC02365, AP000547.3, AL360169.3

### Access
- **BioProject Link**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA815861
- **Publication**: Waddell et al., Colorectal Cancer, 2025

---

## Dataset 3: PRJNA788974 - Right vs Left-Sided Colorectal Cancer

### Basic Information
- **BioProject ID**: PRJNA788974
- **Year**: 2024
- **Data Type**: Raw sequence reads (bulk RNAseq)
- **Relevance**: Medical

### Project Description
Identifying important microbial and genomic biomarkers for differentiating right-sided versus left-sided colorectal cancer using random forest models. This project investigates molecular and microbial differences between tumors arising from different anatomical locations in the colon.

### Dataset Specifications
- **Sample Type**: Colorectal cancer tissues
- **Analysis Type**: Comparative genomics and microbiology
- **Machine Learning**: Random forest modeling of tumor laterality

### Sequencing and Analysis Details
- **Genome Mapping**: Human genome (GRCh38) using STAR (v2.73a)
- **Analysis Focus**:
  - Microbial biomarkers
  - Genomic biomarkers
  - Tumor laterality prediction
  - Right-sided vs left-sided cancer differentiation

### Related Studies
This dataset has been used in multiple publications investigating:
- Bacterial lipopolysaccharide immune response modulation
- Microbial and genomic biomarker identification
- Tumor location-specific molecular signatures

### Access
- **SRA Link**: https://www.ncbi.nlm.nih.gov/sra?term=PRJNA788974
- **BioProject Link**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA788974
- **Related Publication**: PMC10337110, PMC10447454

---

## Comparative Summary Table

| Feature | PRJNA1304982 | PRJNA815861 | PRJNA788974 |
|---------|--------------|-------------|-------------|
| **Registration** | Aug 2025 | 2024 | 2024 |
| **Samples** | 100 | 300 | Multiple |
| **Focus** | RNA editing, CMS | Age-related CMS | Tumor laterality |
| **Data Volume** | 981 Gbases | Not specified | Not specified |
| **Institution** | Italian | New Zealand | Multiple |
| **Key Analysis** | A-to-I editing | Early vs late-onset | Right vs left-sided |

---

## How to Access the Data

### Method 1: Direct SRA Access
1. Visit NCBI SRA: https://www.ncbi.nlm.nih.gov/sra
2. Search for BioProject ID (e.g., PRJNA1304982)
3. Browse available SRA experiments and runs
4. Download FASTQ files directly or use SRA Toolkit

### Method 2: Using SRA Toolkit
```bash
# Install SRA Toolkit
# Download specific runs
fastq-dump --split-files SRR_ACCESSION

# Or for faster download
fasterq-dump SRR_ACCESSION
```

### Method 3: BioProject Portal
- PRJNA1304982: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1304982
- PRJNA815861: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA815861
- PRJNA788974: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA788974

---

## Data Quality and Specifications

### Common Specifications Across Datasets
- **Organism**: Homo sapiens
- **Library Type**: Stranded total RNA
- **Sequencing Platform**: Illumina (NovaSeq/HiSeq)
- **Read Type**: Paired-end
- **Quality Control**: Q30 >90%, PF >80%

### Recommended Analysis Approaches
1. **Quality Control**: FastQC, MultiQC
2. **Alignment**: STAR, HISAT2
3. **Quantification**: featureCounts, HTSeq
4. **Differential Expression**: DESeq2, edgeR, limma
5. **Functional Analysis**: GSEA, Reactome, KEGG

---

## Citation Information

### PRJNA1304982
- **Title**: Dissecting the biology of Consensus Molecular subtypes in colon cancer
- **Institution**: IRCCS-Istituto Tumori "Giovanni Paolo II"
- **Funding**: Tecnopolo per la medicina di Precisione (Regione Puglia)

### PRJNA815861
- **Title**: Consensus molecular subtypes and gene expression in early-onset colorectal cancer
- **Authors**: Waddell et al.
- **Journal**: Colorectal Cancer, 2025
- **Funding**: Health Research Council of New Zealand, Bowel Cancer Research Aotearoa

### PRJNA788974
- **Title**: Identifying important microbial and genomic biomarkers for differentiating right- versus left-sided colorectal cancer using random forest models
- **Related Publications**: Multiple PMC articles (PMC10337110, PMC10447454)

---

## Additional Resources

### Related Databases
- **GEO (Gene Expression Omnibus)**: https://www.ncbi.nlm.nih.gov/geo
- **TCGA (The Cancer Genome Atlas)**: https://www.cancer.gov/tcga
- **COSMIC (Catalogue of Somatic Mutations in Cancer)**: https://cancer.sanger.ac.uk/cosmic

### Bioinformatics Tools
- **SRA Toolkit**: https://github.com/ncbi/sra-tools
- **Galaxy**: https://usegalaxy.org
- **Bioconductor**: https://www.bioconductor.org

### Consensus Molecular Subtypes (CMS)
- **CMS1**: Microsatellite instability (MSI) immune
- **CMS2**: Canonical (chromosomal instability)
- **CMS3**: Metabolic
- **CMS4**: Mesenchymal

---

## Notes
- All datasets are publicly available and free to access
- Data is compliant with FAIR (Findable, Accessible, Interoperable, Reusable) principles
- Researchers should cite the original publications when using these datasets
- Contact information for data providers available upon reasonable request for limited patient metadata

**Last Updated**: November 10, 2025
