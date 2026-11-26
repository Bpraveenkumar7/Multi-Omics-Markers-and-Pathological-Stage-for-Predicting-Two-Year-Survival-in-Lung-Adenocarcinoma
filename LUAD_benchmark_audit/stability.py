from model_eval import (SEED, FOREST, DOWNSAMPLED, get_scores, predictions)
from sklearn.model_selection import StratifiedKFold
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

def compare_fold(k, seed, tr, te, cohort, panels):
    train, holdout = cohort.iloc[tr], cohort.iloc[te]
    full, p_full = predictions(FOREST, train, holdout, panels)
    # The paper kept all minority cases and a random 30% of majority cases.
    # Our control applies this rule to each outer training fold before inner tuning.
    events = train.loc[train.death_2y.eq(1)]
    survivors = train.loc[train.death_2y.eq(0)]
    n_keep = round(0.3 * len(survivors))
    sampled = pd.concat([events, survivors.sample(n=n_keep,random_state=seed + k)]).sort_index()
    small, p_small = predictions(DOWNSAMPLED, sampled, holdout, panels)
    out = pd.DataFrame({'seed':seed, 'fold':k, 'death_2y':holdout.death_2y,
                        FOREST:full, DOWNSAMPLED:small,
                        'weighted_mtry':str(p_full['model__max_features']),
                        'downsampled_mtry':str(p_small['model__max_features'])}, index=holdout.index)
    return out

def robustness(cohort, panels, workers=2, repetitions=20):
    jobs = []
    for seed in range(SEED,SEED + repetitions):
        folds = StratifiedKFold(10, shuffle=True, random_state=seed).split(cohort,cohort.death_2y)
        for k,(tr,te) in enumerate(folds,1):
            jobs.append(delayed(compare_fold)(k,seed,tr,te,cohort,panels))
    records = pd.concat(Parallel(n_jobs=workers)(jobs))
    score_rows = []
    for seed, frame in records.groupby('seed'):
        assert len(frame)==len(cohort) and frame.index.is_unique
        metrics = {name:get_scores(frame.death_2y,frame[name])['ROC AUC']
                   for name in [FOREST,DOWNSAMPLED]}
        score_rows.append({'seed':seed, **metrics, 'paired_difference':metrics[FOREST]-metrics[DOWNSAMPLED]})
    return records,pd.DataFrame(score_rows).set_index('seed')
