# Ancestry-Specific and Multi-Ancestry Rare Variant Meta-Analysis Reveals Novel Genes for Alzheimer’s Disease  
`CARD ❤️ Open Science 😍`

[pending DOI]

**Last Updated:** March 2026  

---

## Summary
This repository accompanies the study titled **"Ancestry-Specific and Multi-Ancestry Rare Variant Meta-Analysis Reveals Novel Genes for Alzheimer’s Disease"**.  

We performed a large-scale multi-ancestry gene-based rare variant analysis of Alzheimer’s disease (AD) using whole-genome sequencing data from the Alzheimer’s Disease Sequencing Project (ADSP) and UK Biobank (UKB). The study includes 17,157 cases and 74,355 controls across diverse ancestries (EUR, AFR, AMR, AAC, AJ, EAS, and CAH).  

Analyses include ancestry-specific burden testing, trans-ancestry meta-analysis, and European-only meta-analysis across multiple minor allele frequency thresholds and functional annotation categories. Significant findings were evaluated for replication in the All of Us cohort. Downstream analyses include functional annotation and phenome-wide association studies (PheWAS).

---

## Highlights
* Largest multi-ancestry rare variant gene-based analysis of AD to date using ADSP and UKB sequencing data  
* Identified 16 significant genes in ancestry-specific analyses (8 known, 8 novel)  
* Trans-ancestry meta-analysis identified 9 significant genes, including 4 novel discoveries  
* 10 novel AD-associated genes identified across all analyses  
* Functional annotation and PheWAS suggest brain-enriched expression and pleiotropic effects across neurological, metabolic, immune, and psychiatric traits  

---

## Citation
If you use this repository or find it helpful for your research, please cite the corresponding manuscript:

> Ancestry-Specific and Multi-Ancestry Rare Variant Meta-Analysis Reveals Novel Genes for Alzheimer’s Disease (Khani et al., 2026)  
>> GitHub DOI: xxx/zenodo.XXXXXXX  
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
│   ├── ADSP_Burden_analysis.ipynb
│   ├── All_of_Us_Burden_analysis.ipynb
│   ├── Step1_Trans_ancestry_analysis.r
│   ├── Step2_European_only_analysis.r
│   └── UKB_Burden_analysis.ipynb
└── README.md

```

---

## Analysis Overview
* Languages: Python, R  

**File** | **Description**
--------------|---------------------------------------------------------------------------------------------------------
ADSP_Burden_analysis.ipynb | Gene-based rare variant burden analysis in ADSP WGS data  
UKB_Burden_analysis.ipynb | Gene-based rare variant burden analysis in UK Biobank WGS data  
All_of_Us_Burden_analysis.ipynb | Replication and burden analysis in All of Us cohort  
Step1_Trans_ancestry_analysis.r | Trans-ancestry meta-analysis of gene-based results  
Step2_European_only_analysis.r | European-only meta-analysis of gene-based results  
