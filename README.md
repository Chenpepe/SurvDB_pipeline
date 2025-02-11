## SurvDB: Systematic Identification of Prognostic Biomarkers Across Nine Molecular Phenotypes in 33 Cancer Types

**Prognostic biomarkers** are molecular indicators used to predict disease progression, recurrence risk, or survival probability. In cancer research and clinical practice, these biomarkers assist in assessing disease severity, anticipating treatment outcomes, and guiding personalized therapies. Identifying prognostic biomarkers enables clinicians and researchers to enhance patient prognosis, optimizing treatment strategies to improve survival and quality of life.

**In [SurvDB](https://gong_lab.hzau.edu.cn/SurvDB/), we:**

1. Systematically compiled phenotype and clinical data on **9 molecular phenotypes** across **33 cancer types** from the TCGA database.
2. Applied **3 survival analysis methods**: Log-Rank test, univariate Cox regression (uni-Cox), and multivariate Cox regression (multi-Cox).
3. Identified biomarkers significantly associated with cancer prognosis for **4 clinical outcome types**: Overall survival (OS), Disease-specific survival (DSS), Progression-free interval (PFI) and Disease-free interval (DFI).

**In [SurvDB](https://gong_lab.hzau.edu.cn/SurvDB/), users can:**

1. Query survival analysis results for different clinical outcome types across various molecular phenotypes.
2. Search specific biomarkers by **cancer type**, **biomarker ID**, or **genomic position**.
3. Download comprehensive analysis results for further research.

## SurvDB pipeline

![pipeline](./pipeline.png "pipeline")

### clinical data procession

```
process_clinical_data.R
```

### phenotype data procession

```
process_pheno_data.R
```

### prognosis biomarker identification

```
do_surv.R
```

### FDR adjustment

```
do_adjust.R
```

### results filtering

```
filter_result.py
```



