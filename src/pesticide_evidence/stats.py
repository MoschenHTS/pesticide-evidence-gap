"""Count-model utilities: Poisson pseudo-maximum-likelihood (PPML) with cluster-robust errors."""

import numpy as np
import pandas as pd

AUTHORIZED = "authorized"
# the three groups entered as indicators in the model; authorized is the reference
GROUPS = ["withdrawn: hazard/treaty", "lapsed: other reason", "absent from Northern registers"]
ABSENT = GROUPS[2]
GROUP_KEYS = ["withdrawn", "lapsed", "absent"]


def round5(x):
    """Float rounded to five decimals, as stored in the results files."""
    return round(float(x), 5)


def ppml(y, X, cluster, iters=300, init=None):
    """Poisson pseudo-ML by iteratively reweighted least squares; returns coefficients and cluster-robust standard errors."""
    beta, cov = ppml_cov(y, X, cluster, iters=iters, init=init)
    return beta, np.sqrt(np.diag(cov))


def ppml_cov(y, X, cluster, iters=300, init=None):
    """As `ppml`, returning the full cluster-robust covariance matrix."""
    k = X.shape[1]
    beta = np.zeros(k) if init is None else np.array(init, dtype=float)
    if init is None:
        beta[0] = np.log(max(y.mean(), 0.1))
    ridge = 1e-8 * np.eye(k)
    for _ in range(iters):
        mu = np.exp(np.clip(X @ beta, -30, 30))
        z = X @ beta + (y - mu) / mu
        XtW = X.T * mu
        new = np.linalg.solve(XtW @ X + ridge, XtW @ z)
        if np.abs(new - beta).max() < 1e-10:
            beta = new
            break
        beta = new
    mu = np.exp(np.clip(X @ beta, -30, 30))
    bread = np.linalg.inv((X.T * mu) @ X + ridge)
    scores = pd.DataFrame(X * (y - mu)[:, None]).groupby(np.asarray(cluster)).sum().values
    return beta, bread @ (scores.T @ scores) @ bread


def pct_change(x):
    """Percentage change from a log coefficient."""
    return (np.exp(x) - 1) * 100


def pct(b, s):
    """Percentage change and 95% interval from a log coefficient and its standard error."""
    return (np.exp(b) - 1) * 100, (np.exp(b - 1.96 * s) - 1) * 100, (np.exp(b + 1.96 * s) - 1) * 100


def analysis_sample(table, min_south_t=50, outcome="en_zh_tox"):
    """Ingredients with enough Southern use, a known target class and complete literature counts."""
    keep = (table.South >= min_south_t) & table.known & table.target.notna() & table[outcome].notna() & table.oa_zh_tox.notna()
    return table[keep].copy()


def design(data, terms=None, extra=(), drop=()):
    """Primary model design with named columns.

    Covariates: ln Southern use, ln Northern use + 1, ln market age, target class, crop-group shares.
    `terms` are column names of `data` entered as the effects of interest (default: the three
    non-authorized groups); `extra` are further covariates entered before them; `drop` names covariates left out.
    """
    for key, group in zip(GROUP_KEYS, GROUPS):
        data = data.assign(**{key: (data.grp4 == group).astype(float)})
    if terms is None:
        terms = GROUP_KEYS
    target = pd.get_dummies(data.target, drop_first=True).astype(float)
    crops = [c for c in data.columns if c.startswith("crop_")][:-1]
    cols = {"const": np.ones(len(data)), "ln_south": np.log(data.South.values), "ln_north": np.log(data.North.values + 1), "ln_age": np.log(data.age.values)}
    for c in extra:
        cols[c] = data[c].values.astype(float)
    for c in terms:
        cols[c] = data[c].values.astype(float)
    for c in target:
        cols["target_" + str(c)] = target[c].values
    for c in crops:
        cols[c] = data[c].fillna(0).values
    for c in drop:
        cols.pop(c, None)
    return np.column_stack(list(cols.values())), list(cols)


def fit(data, outcome="en_zh_tox", terms=None, **kwargs):
    """PPML fit; returns {term: (ratio, lower, upper)} for the effects of interest and the full coefficient table."""
    X, names = design(data, terms=terms, **kwargs)
    b, s = ppml(data[outcome].values.astype(float), X, np.arange(len(data)))
    keys = GROUP_KEYS if terms is None else terms
    out = {
        k: (float(np.exp(b[names.index(k)])), float(np.exp(b[names.index(k)] - 1.96 * s[names.index(k)])), float(np.exp(b[names.index(k)] + 1.96 * s[names.index(k)])))
        for k in keys
    }
    return out, dict(zip(names, b)), dict(zip(names, s))
