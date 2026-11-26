# Multi-Omics Markers and Pathological Stage for Predicting Two-Year Survival in Lung Adenocarcinoma

This repository contains the analysis notebooks and audit files for a study of two-year survival in lung adenocarcinoma (LUAD). The study asks two questions. How well do 32 published molecular markers rank patients by mortality risk when a model is fitted within CPTAC-3? When a model is trained in TCGA and applied to other cohorts, does a ten-gene RNA score improve risk ranking beyond pathological stage and other clinical information?

The endpoint is death before day 730 among patients with known two-year status. Patients last observed alive before day 730 were excluded. CPTAC-3 also excluded known non-cancer deaths. The reported ROC AUC measures how well a model ranks patients who died ahead of those who survived to two years.

## Repository files

| File | Contents |
| --- | --- |
| [LUAD_survival_benchmark.ipynb](LUAD_survival_benchmark.ipynb) | Evaluates the published panel of 20 DNA methylation sites and 12 RNA transcripts. Includes the repeated random forest comparison within CPTAC-3 and the initial models trained on TCGA and tested on CPTAC-3. |
| [LUAD_external_transfer.ipynb](LUAD_external_transfer.ipynb) | Trains clinical and ten-gene RNA models in TCGA. Evaluates the fitted models in CPTAC-3, GSE72094, and GSE31210 and writes patient scores and metrics to CSV. |
| [LUAD_benchmark_audit.zip](LUAD_benchmark_audit.zip) | Contains the initial benchmark’s Python scripts, patient matrix, marker definitions, patient-level predictions, evaluation settings, and file checksums. Its README describes how to rebuild the data and rerun the benchmark. |

## Results

On the 85 primary CPTAC-3 patients, including 10 deaths, a class-weighted random forest using all 32 markers achieved a mean ROC AUC of **0.8301** across 20 repetitions of nested cross-validation. A forest trained after downsampling survivors achieved **0.8012** on the same splits. In the separate experiment that trained models on 267 TCGA patients and tested them on these 85 CPTAC-3 patients, pathological stage alone achieved **0.7880**.

The follow-up model uses major pathological stage to order patients across stages. Within each stage it uses substage, age, sex, and a ten-gene RNA score. The table shows the effect of adding RNA to the same clinical model. All follow-up models were fitted in TCGA and applied to each external cohort without refitting.

| External cohort | Patients | Two-year deaths | Clinical model AUC | Clinical + RNA AUC | Change in AUC |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPTAC-3 primary | 85 | 10 | 0.7427 | 0.7547 | +0.0120 |
| GSE72094 | 320 | 78 | 0.6682 | 0.6911 | +0.0229 |
| GSE31210 | 201 | 11 | 0.9493 | 0.9481 | −0.0012 |

The clearest gain from adding RNA occurred in GSE72094 (paired 95% confidence interval for the AUC change, +0.0076 to +0.0402). The initial CPTAC-3 analysis informed the follow-up model design, making the follow-up result in that cohort exploratory. Both notebooks also report a sensitivity analysis that expands the same CPTAC-3 cohort from 85 to 96 patients.

## Run the notebooks

Open either notebook in Google Colab and select **Runtime → Run all**. A CPU runtime is sufficient. Both notebooks include their analysis inputs.

- The benchmark notebook refits the external models by default and displays saved results for the within-cohort analyses. Set `REFIT_INTERNAL = True` to repeat the single-split within-cohort fits or `REFIT_STABILITY = True` to retrain both forests across all 20 repetitions.
- The external transfer notebook refits the TCGA models and writes `LUAD_transfer_metrics.csv` and `LUAD_transfer_patient_scores.csv` to the working directory.

To run the benchmark as Python scripts, extract the audit ZIP, open a terminal in the extracted directory, and run the following commands with Python 3.12.

```bash
python -m pip install -r requirements.txt
python run_eval.py
python run_stability.py
```

## Data sources

The original 32-marker panel and TCGA/CPTAC-3 analysis sets are from Jaksik et al. (2024), *Multiomics-Based Feature Extraction and Selection for the Prediction of Lung Cancer Survival*. [Article and supplementary tables](https://doi.org/10.3390/ijms25073661)

Clinical and molecular data come from the [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/) and [UCSC Xena](https://github.com/ucscXena/xena-GDC-ETL). The independent microarray cohorts are available from NCBI GEO. [GSE72094](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72094) · [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210)
