from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import json, time, re, hashlib
import pandas as pd
import numpy as np
root=Path(__file__).resolve().parent
def tab(part):
 return pd.read_excel(next((root/'supplement').glob(f'TableS{part} *')),header=3)
lab=tab(1).set_index('PatientID')
shap_header=tab(3)
shared=[col.split(':')[1].replace('SOX9.AS1','SOX9-AS1') for col in shap_header.columns[2:]]
included=lab.loc[shap_header.loc[shap_header['Patient ID'].isin(lab.index),'Patient ID']].copy()
imp=tab(2)
markers=imp.loc[imp.Decision=='Confirmed','Variable ID'].tolist()
def query(path,f,fields):
 u='https://api.gdc.cancer.gov/'+path+'?'+urlencode(dict(filters=json.dumps(f,separators=(',',':')),fields=fields,size=1500))
 for i in range(4):
  try:
   with urlopen(u,timeout=55) as r: data=json.load(r)['data']
   assert len(data['hits'])==data['pagination']['total'];return data['hits']
  except Exception:
   if i==3:raise
   time.sleep(1+i*2)
def match(key,v):return {'op':'in','content':{'field':key,'value':v}}
def meta(job):
 path,ids=job
 cf=match('cases.submitter_id',ids)
 if path=='cases':
  return path,query(path,cf,'submitter_id,diagnoses.ajcc_pathologic_stage,diagnoses.diagnosis_is_primary_disease,diagnoses.primary_diagnosis')
 fil={'op':'and','content':[cf,match('analysis.workflow_type',['STAR - Counts','SeSAMe Methylation Beta Estimation']),match('cases.samples.sample_type',['Primary Tumor']),match('access',['open'])]}
 fields='file_id,file_size,file_name,md5sum,platform,analysis.workflow_type,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type'
 return path,query(path,fil,fields)
def prepare():
 ids=list(lab.index)
 jobs=[(p,ids[i:i+25]) for p in ['files','cases'] for i in range(0,len(ids),25)]
 data={'files':[],'cases':[]}
 with ThreadPoolExecutor(max_workers=7) as ex:
  for p,res in ex.map(meta,jobs):data[p].extend(res)
 (root/'cohort_meta.json').write_text(json.dumps(data))
if __name__=='__main__':prepare()
