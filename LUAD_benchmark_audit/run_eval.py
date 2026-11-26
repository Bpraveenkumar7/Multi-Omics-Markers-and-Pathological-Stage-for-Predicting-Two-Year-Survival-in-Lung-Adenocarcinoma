from model_eval import *
from pathlib import Path
import json,os
root=Path(__file__).resolve().parent
data=pd.read_csv(root/'patient_matrix.csv',index_col='patient_id')
panels=json.loads((root/'panels.json').read_text())
ext,params=external(data,panels)
ext.to_csv(root/'external_predictions.csv')
(root/'external_parameters.json').write_text(json.dumps(params,indent=2))
for tag,df in [('external85',ext[ext.in_paper_analysis]),('external96',ext)]:
 for suffix,t in [('metrics',performance(df)),('auc_intervals',paired_uncertainty(df))]:
  t.to_csv(root/f'{tag}_{suffix}.csv');print(tag,suffix,t,flush=True)
for name,df in data.loc[data.in_paper_analysis].groupby('cohort',sort=False):
 methods=[FOREST,PRIMARY,STAGE,MOLECULAR]+([DISCOVERY] if name=='TCGA' else [])
 oof,p=internal(df,panels,methods,workers=int(os.getenv('BENCHMARK_WORKERS','4')))
 oof.to_csv(root/f'{name}_oof.csv')
 performance(oof).to_csv(root/f'{name}_internal_metrics.csv')
 (root/f'{name}_internal_settings.json').write_text(json.dumps(p,indent=2))
 print(name,performance(oof),flush=True)
