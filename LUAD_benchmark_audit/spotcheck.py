import hashlib,io,json
import numpy as np,pandas as pd
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
p=Path('.')
data=pd.read_csv('patient_matrix.csv',index_col='patient_id')
files=json.load(open('cohort_meta.json'))['files']
def work(z):
 k,id,modal,f=z
 with urlopen('https://api.gdc.cancer.gov/data/'+f['id'],timeout=120) as r:b=r.read()
 assert hashlib.md5(b).hexdigest()==f['md5sum']
 if modal=='RNA':
  d=pd.read_csv(io.BytesIO(b),sep='\t',comment='#',low_memory=False)
  v=np.log2(d.loc[d.gene_name.isin(['RAET1E','DKK1','CERS4']),'tpm_unstranded']+1)
  v.index=d.loc[v.index,'gene_name']
 else:
  d=pd.read_csv(io.BytesIO(b),sep='\t',names=['probe','value']).set_index('probe')
  v=d.loc[['cg27326750','cg22884440','cg09453760'],'value']
 diffs=(data.loc[id,list(v.index)].astype(float)-v)
 return dict(cohort=k,patient=id,assay=modal,file_id=f['id'], max_abs_difference=diffs.abs().max(),details={s:round(q,5) for s,q in diffs.items()})
jobs=[]
for k in ['TCGA','CPTAC-3']:
 n=0
 for id in data.loc[(data.cohort==k)&data.in_paper_analysis].index:
  ff=[f for f in files if any(c['submitter_id']==id for c in f['cases'])]
  a=[f for f in ff if f['file_name'].endswith('star_gene_counts.tsv')]
  b=[f for f in ff if f['file_name'].endswith('level3betas.txt') and f['platform']==('Illumina Human Methylation 450' if k=='TCGA' else 'Illumina Methylation Epic')]
  if len(a)!=1 or len(b)!=1:continue
  jobs.extend([(k,id,'RNA',a[0]),(k,id,'CpG',b[0])]);n+=1
  if n==2:break
with ThreadPoolExecutor(max_workers=5) as ex:res=list(ex.map(work,jobs))
print(json.dumps(res,indent=1))
open('source_spotchecks.json','w').write(json.dumps(res,indent=1))
