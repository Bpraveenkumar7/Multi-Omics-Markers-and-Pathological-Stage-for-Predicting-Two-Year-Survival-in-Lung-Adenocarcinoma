import json, os
from model_eval import *
from stability import robustness
data=pd.read_csv('patient_matrix.csv',index_col='patient_id')
panels=json.load(open('panels.json'))
df=data.loc[data.cohort.eq('CPTAC-3') & data.in_paper_analysis]
results,stats=robustness(df,panels,workers=int(os.getenv("BENCHMARK_WORKERS","2")),repetitions=20)
results.to_csv('repeat_oof_predictions.csv')
stats.to_csv('repeat_auc.csv')
orig=pd.read_csv('CPTAC-3_oof.csv',index_col='patient_id')
assert np.allclose(results[results.seed.eq(SEED)].loc[orig.index,FOREST],orig[FOREST])
print(stats.describe(),flush=True)
print('above target',sum(stats[FOREST]>=.816), 'pooled',roc_auc_score(results.death_2y,results[FOREST]),roc_auc_score(results.death_2y,results[DOWNSAMPLED]),flush=True)
