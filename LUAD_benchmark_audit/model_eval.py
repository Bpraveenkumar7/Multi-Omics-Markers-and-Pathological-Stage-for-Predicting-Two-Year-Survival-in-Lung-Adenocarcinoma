"""Two-year LUAD model comparison on the paper-listed patients."""
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20261125
PRIMARY = 'Stage + molecular ridge'
FOREST = 'Molecular RF'
MOLECULAR = 'Molecular ridge'
STAGE = 'Stage alone'
DISCOVERY = 'Discovery-panel RF'
DOWNSAMPLED = 'Downsampled molecular RF'

def columns_for(name, panels):
    if name == STAGE:
        return ['stage']
    if name == DISCOVERY:
        return list(panels['discovery'])
    return list(panels['shared']) + (['stage'] if name == PRIMARY else [])

def fit_model(name, train, panels):
    columns = columns_for(name, panels)
    # A feature with no measured value anywhere in this training set cannot enter the model.
    columns = [c for c in columns if train[c].notna().any()]
    steps = [('impute', SimpleImputer(strategy='median'))]
    if name in [FOREST, DISCOVERY, DOWNSAMPLED]:
        weight = None if name == DOWNSAMPLED else 'balanced_subsample'
        model = RandomForestClassifier(n_estimators=300, class_weight=weight,
                                       n_jobs=1, random_state=SEED)
        grid = {'model__max_features': ['sqrt', 0.5, 1.0]}
    else:
        steps.append(('standardize', StandardScaler()))
        model = LogisticRegression(solver='lbfgs', max_iter=2000)
        grid = {'model__C': [0.01, 0.1, 1.0, 10.0]}
    steps.append(('model', model))
    inner = StratifiedKFold(5, shuffle=True, random_state=SEED)
    search = GridSearchCV(Pipeline(steps), grid, scoring='roc_auc', cv=inner, n_jobs=1,
                          error_score='raise')
    search.fit(train[columns], train.death_2y)
    return search, columns

def predictions(name, train, test, panels):
    model, columns = fit_model(name, train, panels)
    return model.predict_proba(test[columns])[:, 1], model.best_params_

def get_scores(y, risk):
    return {'ROC AUC':roc_auc_score(y, risk), 'Average precision':average_precision_score(y, risk),
            'Brier score':brier_score_loss(y, risk)}

def external(data, panels):
    train = data.loc[data.cohort.eq('TCGA')]
    test = data.loc[data.cohort.eq('CPTAC-3')]
    scores = test[['death_2y', 'in_paper_analysis']].copy()
    settings = {}
    for name in [FOREST, STAGE, MOLECULAR, PRIMARY]:
        scores[name], settings[name] = predictions(name, train, test, panels)
    return scores, settings

def performance(scores):
    rows = {name:get_scores(scores.death_2y, scores[name]) for name in scores
            if name not in ['death_2y', 'in_paper_analysis', 'fold']}
    return pd.DataFrame.from_dict(rows, orient='index')

def paired_uncertainty(scores, repeats=10000):
    # For fixed predictions, ROC AUC is the average of all event/non-event pair comparisons.
    # Pair matrices preserve ties. Patient resamples and their multiplicities are shared by models.
    models = [c for c in scores if c not in ['death_2y', 'in_paper_analysis', 'fold']]
    values = scores[models].to_numpy()
    events, survivors = values[scores.death_2y.eq(1)], values[scores.death_2y.eq(0)]
    gap = events[:,None,:] - survivors[None,:,:]
    wins = (gap > 0).astype(float) + 0.5 * (gap == 0)
    rng = np.random.default_rng(SEED)
    w_event = rng.multinomial(len(events), np.full(len(events), 1/len(events)), size=repeats)
    w_surv = rng.multinomial(len(survivors), np.full(len(survivors), 1/len(survivors)), size=repeats)
    auc_samples = np.einsum('ri,ijm,rj->rm',w_event,wins,w_surv,optimize=True)/(len(events)*len(survivors))
    points = wins.mean(axis=(0,1))
    assert np.allclose(points,[roc_auc_score(scores.death_2y,scores[m]) for m in models])
    intervals = {m: [points[i], *np.quantile(auc_samples[:,i],[.025,.975])]
                 for i,m in enumerate(models)}
    for other in [FOREST, STAGE, MOLECULAR]:
        a, b = models.index(PRIMARY), models.index(other)
        intervals[f'{PRIMARY} minus {other}'] = [points[a]-points[b], *np.quantile(auc_samples[:,a]-auc_samples[:,b],[.025,.975])]
    return pd.DataFrame.from_dict(intervals, orient='index', columns=['estimate','lower_95','upper_95'])

def fold_prediction(name, fold, train_rows, test_rows, cohort, panels):
    pred, params = predictions(name, cohort.iloc[train_rows], cohort.iloc[test_rows], panels)
    return name, fold, cohort.index[test_rows], pred, params

def internal(cohort, panels, names, workers=2):
    folds = list(StratifiedKFold(10, shuffle=True, random_state=SEED).split(cohort,cohort.death_2y))
    jobs = (delayed(fold_prediction)(name,k,tr,te,cohort,panels)
            for name in names for k,(tr,te) in enumerate(folds,1))
    fitted = Parallel(n_jobs=workers)(jobs)
    result = cohort[['death_2y']].copy()
    settings = []
    for name,k,ids,risk,params in fitted:
        result.loc[ids,'fold'] = k
        result.loc[ids,name] = risk
        settings.append({'model':name,'fold':k,**params})
    assert result[names].notna().all().all()
    return result, settings
