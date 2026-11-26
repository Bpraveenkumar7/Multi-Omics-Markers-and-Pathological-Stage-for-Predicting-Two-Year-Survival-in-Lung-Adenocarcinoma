from cohort import *
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict,Counter

def ask(expr):
 for attempt in range(4):
  try:
   with urlopen(Request('https://gdc.xenahubs.net/data/',data=expr.encode(),headers={'Content-Type':'text/plain'}),timeout=55) as r:return json.load(r)
  except Exception:
   if attempt==3:raise
   time.sleep(attempt+1)
def names(dataset):
 return ask('(map :value (query {:select [:value] :from [:dataset] :join [:field [:= :dataset.id :dataset_id] :code [:= :field.id :field_id]] :where [:and [:= :dataset.name '+json.dumps(dataset)+'] [:= :field.name "sampleID"]]}))')
def expression_map():
 with urlopen('https://api.gdc.cancer.gov/data/780617e8-6a28-4774-b1e6-3135e1932327',timeout=60) as r:
  e=pd.read_csv(r,comment='#',sep='\t',low_memory=False)
 df=e[e.gene_name.isin([g for g in markers if not g.startswith('cg')])]
 assert df.gene_name.is_unique and set(df.gene_name)=={g for g in markers if not g.startswith('cg')}
 return dict(zip(df.gene_name,df.gene_id))
def get_members(meta,cohort,kind):
 members=defaultdict(set)
 for f in meta['files']:
  if kind=='rna' and not f['file_name'].endswith('star_gene_counts.tsv'):continue
  if kind=='meth' and (not f['file_name'].endswith('level3betas.txt') or f['platform']!=('Illumina Human Methylation 450' if cohort=='TCGA' else 'Illumina Methylation Epic')):continue
  for c in f['cases']:
   pat=c['submitter_id']
   if lab.loc[pat,'Dataset']!=cohort:continue
   for s in c['samples']:
    if s['sample_type']=='Primary Tumor':members[pat].add(s['submitter_id'])
 return members
def getmatrix(job):
 c,kind,ds,m,ids=job
 available=set(names(ds))
 valid={p:sorted(v & available) for p,v in m.items() if v & available}
 missing=list(set(lab[lab.Dataset==c].index)-set(valid))
 xs=sorted({s for samps in valid.values() for s in samps})
 print(c,kind,'matched cases',len(valid),'samples',len(xs),'missing',missing[:8],flush=True)
 meta=ask('(query {:select [:name :text] :from [:dataset] :where [:= :name '+json.dumps(ds)+']})')
 print(ds,str(meta)[:1800],flush=True)
 pr=list(ids.values())
 matrix=ask('(fetch [{:table '+json.dumps(ds)+' :columns '+json.dumps(pr)+' :samples '+json.dumps(xs)+'}])')
 a=pd.DataFrame(matrix,index=list(ids),columns=xs).apply(pd.to_numeric,errors='coerce').T
 print(ds,'observed by marker',a.count().to_dict(),flush=True)
 return ds,pd.DataFrame({p:a.loc[v].median() for p,v in valid.items()}).T,valid,meta
def run():
 info=json.loads((root/'cohort_meta.json').read_text())
 genes=expression_map()
 probes={g:g for g in markers if g.startswith('cg')}
 jobs=[]
 for c,p,mt in [('TCGA','TCGA-LUAD','methylation450'),('CPTAC-3','CPTAC-3','methylation_epic')]:
  for kind,ds,ids in [('rna',p+'.star_tpm.tsv',genes),('meth',p+'.'+mt+'.tsv',probes)]:
   jobs.append((c,kind,ds,get_members(info,c,kind),ids))
 with ThreadPoolExecutor(max_workers=4) as ex:
  for ds,t,v,meta in ex.map(getmatrix,jobs):
   t.to_csv(root/(ds+'.csv'));(root/(ds+'.samples.json')).write_text(json.dumps(v));(root/(ds+'.meta.json')).write_text(json.dumps(meta))
 (root/'gene_map.json').write_text(json.dumps(genes))
if __name__=='__main__':run()
