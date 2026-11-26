from cohort import *
def build():
 sets=[]
 for c,p,m in [('TCGA','TCGA-LUAD','methylation450'),('CPTAC-3','CPTAC-3','methylation_epic')]:
  ids=lab[lab.Dataset==c].index
  omics=[pd.read_csv(root/(p+'.'+ds+'.tsv.csv'),index_col=0) for ds in ('star_tpm',m)]
  sets.append(pd.concat(omics,axis=1).reindex(ids))
 X=pd.concat(sets)
 stages={}
 for case in json.loads((root/'cohort_meta.json').read_text())['cases']:
  eligible=[d for d in case['diagnoses'] if d.get('diagnosis_is_primary_disease')]
  levels=set()
  for d in eligible:
   hit=re.fullmatch(r'Stage (IV|III|II|I)[A-C]?[1-3]?',d.get('ajcc_pathologic_stage') or '')
   if hit:levels.add({'I':1,'II':2,'III':3,'IV':4}[hit[1]])
  assert len(levels)<=1
  stages[case['submitter_id']]=next(iter(levels),np.nan)
 frame=lab.rename(columns={'Dataset':'cohort','Survival (2-year)':'outcome'}).join(X)
 frame['stage']=frame.index.to_series().map(stages)
 frame['in_paper_analysis']=frame.index.isin(included.index)
 frame['death_2y']=(frame.outcome=='Died').astype(int)
 assert set(frame.outcome)=={'Died','Survived'} and frame.index.is_unique
 for c,d in frame.groupby('cohort'):
  print(c,len(d),len(d.loc[d.in_paper_analysis]),'stage missing',d.stage.isna().sum(),'shared molecular completely missing',d[shared].isna().all(axis=1).sum(),'CpG block missing',d[[x for x in shared if x.startswith('cg')]].isna().all(axis=1).sum(),'RNA block missing',d[[x for x in shared if not x.startswith('cg')]].isna().all(axis=1).sum())
 frame.to_csv(root/'patient_matrix.csv',index_label='patient_id')
 (root/'panels.json').write_text(json.dumps({'shared':shared,'discovery':markers}))
if __name__=='__main__':build()
