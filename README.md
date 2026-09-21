# Ancestry-Specific and Multi-Ancestry Rare Variant Meta-Analysis Reveals Novel Genes for Alzheimer’s Disease  
`CARD ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/1184371261.svg)](https://doi.org/10.5281/zenodo.19072549)

**Last Updated:** September 2026  

---

## Summary
This repository accompanies the study titled **"Ancestry-Specific and Multi-Ancestry Rare Variant Meta-Analysis Reveals Novel Genes for Alzheimer’s Disease"**.  

We performed a large-scale multi-ancestry gene-based rare variant analysis of Alzheimer’s disease (AD) using whole-genome sequencing data from the Alzheimer’s Disease Sequencing Project (ADSP) and UK Biobank (UKB). The study includes 17,157 cases and 74,355 controls across diverse ancestries (EUR, AFR, AMR, AAC, AJ, EAS, and CAH).  

Analyses include ancestry-specific rare variant burden testing across seven ancestries, trans-ancestry meta-analysis, and European-only meta-analysis across multiple minor allele frequency thresholds and functional annotation masks. Meta-analyses include single-variant effect size meta-analysis (ESMA), gene-based burden meta-analysis (GENE), omnibus gene-level p-value meta-analysis (GENE_P), and p-value meta-analysis (PVMA). Significant findings were evaluated for replication in the All of Us cohort. Downstream analyses include functional annotation and phenome-wide association studies (PheWAS).

---

## Highlights
* Largest multi-ancestry rare variant gene-based analysis of AD to date using ADSP and UKB sequencing data  
* Identified X significant genes in ancestry-specific analyses (X known, Y novel)  
* Trans-ancestry meta-analysis identified X significant genes, including X novel discoveries  
* X novel AD-associated genes identified across all analyses  
* Functional annotation and PheWAS suggest brain-enriched expression and pleiotropic effects across neurological, metabolic, immune, and psychiatric traits  

---

## Citation
If you use this repository or find it helpful for your research, please cite the corresponding manuscript:

> Ancestry-Specific and Multi-Ancestry Rare Variant Meta-Analysis Reveals Novel Genes for Alzheimer’s Disease (Khani et al., 2026)  
>> GitHub DOI: 10.5281/zenodo.19072549
>> Manuscript DOI: xxx
---

## Data Statement
* Whole-genome sequencing (WGS) and array data were obtained from:
  * Alzheimer’s Disease Sequencing Project (ADSP) v5  
  * UK Biobank (UKB) v18.1  
  * All of Us (AoU) v8  

* All cohorts underwent standardized quality control and ancestry inference using the **GenoTools** pipeline  
* Access to each dataset must be requested through the respective data platforms  

---

## Repository Orientation
```
.
├── analyses
│   ├── 00_ADSP_Burden_analysis.ipynb
│   ├── 00_UKB_Burden_analysis.ipynb
│   ├── 00_ESMA_Meta_anaysis.ipynb
│   ├── 00_GENE_Meta_analysis.ipynb
│   ├── 00_GENE_P_Meta_analysis.ipynb
│   ├── 00_PVMA_Meta_analysis.ipynb
│   └── 00_All_of_Us_Burden_analysis.ipynb
└── README.md

```

---

## Analysis Overview
* Language: Python  

**File** | **Description**
--------------|---------------------------------------------------------------------------------------------------------
00_ADSP_Burden_analysis.ipynb | Gene-based rare variant burden analysis in ADSP WGS data  
00_UKB_Burden_analysis.ipynb | Gene-based rare variant burden analysis in UK Biobank WGS data  
00_All_of_Us_Burden_analysis.ipynb | Replication and burden analysis in All of Us cohort  
00_ESMA_Meta_anaysis.ipynb | Trans-ancestry and European-only single-variant effect size meta-analysis (ESMA)  
00_GENE_Meta_analysis.ipynb | Trans-ancestry and European-only gene-based burden meta-analysis (GENE)
00_GENE_P_Meta_analysis.ipynb | Trans-ancestry and European-only omnibus gene-level p-value meta-analysis (GENE_P)
00_PVMA_Meta_analysis.ipynb | Trans-ancestry and European-only p-value meta-analysis (PVMA)
