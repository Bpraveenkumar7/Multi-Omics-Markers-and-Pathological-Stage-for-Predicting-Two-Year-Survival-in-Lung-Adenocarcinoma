# Lung adenocarcinoma two-year survival benchmark

A Python feasibility analysis of Jaksik et al., International Journal of Molecular Sciences 25 (2024), 3661, DOI 10.3390/ijms25073661. All input records are from public resources. The accompanying Colab notebook includes the small dataset and saved outputs. This archive provides the CSVs and source extraction and evaluation scripts.

## Results

On the 85 paper-listed CPTAC-3 patients, the weighted random forest fitted within CPTAC-3 had an out-of-fold ROC AUC of 0.832 in the first split. Across all 20 fixed repetitions, its mean was 0.8301, range 0.7880 to 0.8640, with 16 of 20 at or above 0.816. Pooling all out-of-fold predictions across repeats gave 0.8305. The downsampled control achieved mean 0.8012 and pooled AUC 0.8029. The mean paired split difference was +0.0289. These observations are split-sensitivity checks on 10 deaths and 75 survivors. The fixed panel and all 20 seeds are retained.

In the external experiment, models were fitted entirely on the 267 TCGA cases. Frozen predictions from the molecular forest achieved AUC 0.4613, stage alone 0.7880, molecular ridge 0.4507, and stage plus molecular ridge 0.4653 on the same 85 CPTAC-3 cases. The publication's CPTAC-3 AUC (0.815) is a within-cohort reference.

## Run

Use Python 3.12.14. From the unpacked directory run the following commands. The CSV and panel definition supplied here are sufficient for model fitting.

```bash
python -m pip install -r requirements.txt
python run_eval.py
python run_stability.py
```

The evaluation script runs the locked TCGA-to-CPTAC-3 comparisons, plus five-fold inner and ten-fold outer CV within each cohort. Set BENCHMARK_WORKERS to the desired job count for that script. The stability script defaults to two workers outside the original environment and accepts the same setting. It checks the first repeat against the independently generated CPTAC-3 out-of-fold predictions.

To rebuild patient_matrix.csv from original public sources, use python fetch_sources.py in a separate copy of this folder. It downloads the publisher's three supplements, queries GDC metadata, extracts primary-tumor measurements from the Xena GDC hub, and assembles the matrix. Rebuilt results may drift with upstream revisions. The saved source metadata, sample lists, gene identifiers, and file checksums identify the snapshot analyzed here.

## Files

- patient_matrix.csv holds patient IDs, cohort, Table S1 outcomes, 12 RNA measurements, 28 methylation columns, ordinal pathological stage, original analysis-set membership, and the binary outcome. Only the designated molecular columns and stage can enter a model. Three TCGA patients lack matched 450K data, two lack stage, and one discovery-only probe is missing throughout TCGA. The analysis retained patients using training-fold median imputation.
- panels.json contains the 32 shared markers and the 40 discovery-panel markers. SHAP values in supplementary Table S3 were never used as predictors.
- external_predictions.csv contains all 96 CPTAC-3 predictions from the TCGA models. external85 and external96 files contain metrics and paired AUC intervals. External bootstrap calculations condition on fitted models and use 10,000 paired, outcome-stratified patient resamples.
- TCGA_oof.csv and CPTAC-3_oof.csv contain the patient-level predictions from the first ten-fold outer split; their corresponding settings and metrics are included.
- repeat_oof_predictions.csv contains all 1,700 held-out patient predictions per model across the 20 repetitions. repeat_auc.csv reports AUC by seed. The downsampled control uses 30% of the larger class and all deaths in each outer training fold. Inner five-fold tuning uses the reduced training set. The weighted forest trains on the full outer training fold with balanced_subsample weights. Splits are paired across both approaches. A repetition is not an independent cohort.
- protocol.txt records the original external and single-split analysis. The stability experiment was added after viewing the first within-CPTAC-3 result. The held-out external experiment was unchanged.
- source_spotchecks.json records comparisons of three RNA values and three CpG values in two patients per cohort against eight independently retrieved, checksum-verified GDC files.
- cohort_meta.json, gene_map.json, and the Xena metadata and sample manifests record the provenance of patient and measurement matching.
- SHA256SUMS.txt identifies the exact contents of each archived input and result.

## Interpretation and data sources

The published Figure 5 reports 0.839 for the full TCGA panel, 0.830 for the TCGA shared panel, and 0.815 for the CPTAC-3 shared panel. Published feature selection used TCGA outcomes before evaluation of that fixed panel. The reference paper used R/caret, inner LOOCV, downsampling, repeated nested CV and DESeq2. The available GitHub commit includes processing and feature-selection code but omits the final nested classifier/ROC driver, historical analysis-ready matrices and fold predictions. This benchmark uses Xena GDC STAR log2(TPM+1), beta values, scikit-learn, and a smaller tuning grid. The first paper-specific control is an approximation and differs in these respects. No paired confidence interval against the published scalar AUC is available. TCGA and CPTAC-3 also used different cause-of-death rules. The binary endpoint excludes patients censored before two years; it does not estimate survival for all newly diagnosed patients.

- Base paper, Jaksik et al.: https://doi.org/10.3390/ijms25073661
- Official supplements: https://www.mdpi.com/article/10.3390/ijms25073661/s1
- Author repository: https://github.com/rjaksik/LUAD_multiomic/tree/08c7ead9765309c85b0ff1d660f1ff5b41507af3
- NCI Genomic Data Commons: https://portal.gdc.cancer.gov/
- UCSC Xena GDC ETL: https://github.com/ucscXena/xena-GDC-ETL
- Xena queries: https://gdc.xenahubs.net/data/
