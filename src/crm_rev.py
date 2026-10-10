"""
crm_rev.py - extra analyses asked for in the review of the paper.

Every function reads the main table `master` built in section 4 of the notebook
(columns: group, s_n5_class, composition, elements, formula, has_* flags, metric columns).
Nothing here changes the earlier results; all outputs are new files whose names start with R.

R1  cut-off sweep for the common support (and the support without 4-element metal-only compounds)
R2  cluster bootstrap for the main gaps and for the Mantel-Haenszel odds ratios
R3  application shortlists matched by compound type on the common support
R4  stricter lower-risk alternatives
R5  robustness checks run on the common support
R6  element ratios with compound-type weights, 1000 bootstrap samples
R7  gap decomposition in both weight orders, 1000 bootstrap samples
R8  stability margins controlled for the number of MP entries in each chemical system
"""
from __future__ import annotations

import itertools
import math
import numpy as np
import pandas as pd

import crm_audit as ca

PGM = ["Ru", "Rh", "Pd", "Os", "Ir", "Pt"]
REE = set(ca.REE)
KEY_METRICS = ["crm_eu_atom", "crm_us_atom", "has_ree", "has_f_in_core", "has_pgm"]
COUNTRY_METRICS = ["hhi_new", "china_frac", "conc_frac"]


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def metrics_present(master, extra=True):
    m = [c for c in KEY_METRICS if c in master.columns]
    if extra:
        m += [c for c in COUNTRY_METRICS if c in master.columns and master[c].notna().any()]
    return m


def support_types(master, ref="MP_known", min_ref=100, strata="s_n5_class", drop=()):
    """Compound types with at least `min_ref` compounds in the reference group."""
    ct = master.loc[master["group"] == ref, strata].value_counts()
    ok = [s for s in ct.index[ct >= min_ref] if s not in set(drop)]
    return sorted(ok)


def matched_gap(T, R, metric, strata="s_n5_class"):
    """Gap between target and re-weighted reference with the analytic 95% CI used in the paper
    (eqn 7). T and R must already be limited to the compared types."""
    T = T.dropna(subset=[metric]); R = R.dropna(subset=[metric])
    common = sorted(set(T[strata]) & set(R[strata]))
    T = T[T[strata].isin(common)]
    if not len(T):
        return None
    tv = T[metric].astype(float).values
    pi = T[strata].value_counts(normalize=True)
    g = R[R[strata].isin(common)].assign(_v=lambda x: x[metric].astype(float)) \
        .groupby(strata)["_v"].agg(["mean", "var", "size"]).loc[pi.index]
    ref_m = float((pi * g["mean"]).sum())
    var_ref = float((pi ** 2 * g["var"].fillna(0) / g["size"]).sum())
    var_t = tv.var(ddof=1) / len(tv) if len(tv) > 1 else 0.0
    d, se = tv.mean() - ref_m, math.sqrt(var_t + var_ref)
    return dict(metric=metric, target=tv.mean(), ref_matched=ref_m, diff=d,
                ci_low=d - 1.96 * se, ci_high=d + 1.96 * se, n_target=len(tv), n_ref=int(g["size"].sum()))


def ln_key(comp_str):
    """Near-copy template used in the paper: every rare earth replaced by one symbol."""
    c = ca.parse_composition(comp_str)
    return ";".join(sorted(f"{'Ln' if e in REE else e}:{v:g}" for e, v in c.items()))


def family_key(elements_str):
    """Coarser cluster: the set of elements, with rare earths -> Ln and Pt-group metals -> Pg,
    stoichiometry ignored. Puts whole substitution families into one cluster."""
    out = []
    for e in elements_str.split(","):
        out.append("Ln" if e in REE else ("Pg" if e in PGM else e))
    return "-".join(sorted(out))


# --------------------------------------------------------------------------------------
# R1  cut-off sweep
# --------------------------------------------------------------------------------------
def cutoff_sweep(master, ref="MP_known", target="GNoME", cutoffs=(100, 300, 500, 1000),
                 strata="s_n5_class"):
    metrics = metrics_present(master)
    n_all = int((master["group"] == target).sum())
    variants = [(f"known >= {c} per type", c, ()) for c in cutoffs]
    variants.append(("known >= 100, without 4-element metal-only", 100, ("4_intermetallic",)))
    rows = []
    for label, c, drop in variants:
        ok = support_types(master, ref, c, strata, drop)
        T = master[(master["group"] == target) & master[strata].isin(ok)]
        R = master[(master["group"] == ref) & master[strata].isin(ok)]
        for m in metrics:
            r = matched_gap(T, R, m, strata)
            if r:
                rows.append(dict(variant=label, types_kept=" ".join(ok), coverage=len(T) / n_all, **r))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# R2  cluster bootstrap
# --------------------------------------------------------------------------------------
def add_cluster_keys(df):
    """Adds the two cluster keys once (slow step, about a minute on the full table)."""
    df = df.copy()
    comp = pd.Series(df["composition"].unique())
    df["_ln_key"] = df["composition"].map(dict(zip(comp, comp.map(ln_key))))
    els = pd.Series(df["elements"].unique())
    df["_fam_key"] = df["elements"].map(dict(zip(els, els.map(family_key))))
    return df


def _cluster_codes(df, how):
    col = "_ln_key" if how == "ln_template" else "_fam_key"
    if col in df.columns:
        keys = df[col]
    else:
        keys = df["composition"].map(ln_key) if how == "ln_template" else df["elements"].map(family_key)
    return pd.factorize(keys)[0]


def cluster_boot_gap(master, metric, cluster="ln_template", ref="MP_known", target="GNoME",
                     min_ref=100, strata="s_n5_class", n_boot=1000, seed=0):
    """Matched gap on the common support with a cluster (Poisson) bootstrap: each cluster of
    near-copies gets one random weight, so near-copies move together."""
    rng = np.random.default_rng(seed)
    ok = support_types(master, ref, min_ref, strata)
    d = master[master["group"].isin([target, ref]) & master[strata].isin(ok)].dropna(subset=[metric])
    d = d.assign(_cl=_cluster_codes(d, cluster))
    out = {}
    stypes = sorted(set(d[strata]))
    pi = d[d["group"] == target][strata].value_counts(normalize=True).reindex(stypes).fillna(0).values
    for grp in (target, ref):
        g = d[d["group"] == grp]
        s_code = pd.Categorical(g[strata], categories=stypes).codes
        # one row per (cluster, type): sums and counts
        key = g["_cl"].values * len(stypes) + s_code
        u, inv = np.unique(key, return_inverse=True)
        sums = np.bincount(inv, weights=g[metric].astype(float).values)
        cnts = np.bincount(inv).astype(float)
        cl_of_u = u // len(stypes); st_of_u = u % len(stypes)
        cl_u, cl_inv = np.unique(cl_of_u, return_inverse=True)
        out[grp] = (sums, cnts, st_of_u, cl_inv, len(cl_u))
    def stat(wT, wR):
        sT, nT, stT, ciT, _ = out[target]; sR, nR, stR, ciR, _ = out[ref]
        a, b = wT[ciT], wR[ciR]
        t_mean = (sT * a).sum() / (nT * a).sum()
        rs = np.bincount(stR, weights=sR * b, minlength=len(stypes))
        rn = np.bincount(stR, weights=nR * b, minlength=len(stypes))
        r_mean = float(np.nansum(pi * np.where(rn > 0, rs / np.where(rn > 0, rn, 1), np.nan)))
        return t_mean - r_mean
    point = stat(np.ones(out[target][4]), np.ones(out[ref][4]))
    bs = np.array([stat(rng.poisson(1.0, out[target][4]), rng.poisson(1.0, out[ref][4]))
                   for _ in range(n_boot)])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return dict(metric=metric, cluster=cluster, diff=point, ci_low=lo, ci_high=hi,
                clusters_target=out[target][4], clusters_ref=out[ref][4],
                compounds_target=int((d["group"] == target).sum()), compounds_ref=int((d["group"] == ref).sum()))


def _mh(a, b, c, dd):
    n = a + b + c + dd
    ok = (n > 0) & ((a + b) > 0) & ((c + dd) > 0)
    R = (a * dd / np.where(ok, n, 1))[ok].sum(); S = (b * c / np.where(ok, n, 1))[ok].sum()
    return R / S if S > 0 else np.nan


def cluster_boot_mh(G, flag="has_f_in_core", margin_col="e_decomp_mp", thr_meV=25.0,
                    cluster="ln_template", strata="s_n_class", n_boot=1000, seed=0):
    """Mantel-Haenszel odds ratio of a small margin, with clusters of near-copies resampled together."""
    rng = np.random.default_rng(seed)
    d = G.dropna(subset=[margin_col]).copy()
    d["near"] = (-1000 * d[margin_col].astype(float)) < thr_meV
    d["_cl"] = _cluster_codes(d, cluster)
    stypes = sorted(set(d[strata])); s_code = pd.Categorical(d[strata], categories=stypes).codes
    f, y = d[flag].values.astype(bool), d["near"].values
    cell = np.select([f & y, f & ~y, ~f & y], [0, 1, 2], 3)
    key = (d["_cl"].values * len(stypes) + s_code) * 4 + cell
    u, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    cell_u = u % 4; st_u = (u // 4) % len(stypes); cl_u = (u // 4) // len(stypes)
    clu, cl_inv = np.unique(cl_u, return_inverse=True)
    def stat(w):
        ww = cnt * w[cl_inv]
        tab = np.bincount(st_u * 4 + cell_u, weights=ww, minlength=len(stypes) * 4).reshape(len(stypes), 4)
        return _mh(tab[:, 0], tab[:, 1], tab[:, 2], tab[:, 3])
    point = stat(np.ones(len(clu)))
    bs = np.array([stat(rng.poisson(1.0, len(clu))) for _ in range(n_boot)])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    share_f, share_nf = y[f].mean(), y[~f].mean()
    return dict(flag=flag, margin_col=margin_col, threshold_meV=thr_meV, cluster=cluster,
                or_mh=point, ci_low=lo, ci_high=hi, clusters=len(clu), compounds=len(d),
                share_with=share_f, share_without=share_nf, risk_difference=share_f - share_nf)


# --------------------------------------------------------------------------------------
# R3  shortlists matched by compound type
# --------------------------------------------------------------------------------------
def shortlists_type_matched(eg, master, ref="MP_known", min_ref=100, conf_min=0.9,
                            strata="s_n5_class", conf_col="AI-experts confidence (-)"):
    g = master[master["group"] == "GNoME"]
    ok = support_types(master, ref, min_ref, strata)
    R = master[(master["group"] == ref) & master[strata].isin(ok)]
    rows = []
    for app in sorted(eg["application"].unique()):
        e = eg[eg["application"] == app]
        sets = [("all", e)]
        if conf_col in e.columns:
            sets.append((f"confidence>={conf_min}", e[pd.to_numeric(e[conf_col], errors="coerce") >= conf_min]))
        for lab, es in sets:
            sub = g[g["gnome_id"].isin(set(es["Material Id"].astype(str)))]
            T = sub[sub[strata].isin(ok)]
            for m in ["crm_eu_atom", "has_ree", "has_f_in_core", "has_pgm"]:
                r = matched_gap(T, R, m, strata)
                if r:
                    rows.append(dict(application=app, subset=lab, n_shortlist=len(sub),
                                     share_on_support=len(T) / max(len(sub), 1), **r))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# R4  stricter alternatives
# --------------------------------------------------------------------------------------
def alternatives_strict(egj, app_perf, perf_fn, conf_col, risk="crm_eu_atom", high_risk=0.5,
                        rel_tol=0.10, abs_tol=0.10, min_risk_drop=0.2, conf_min=0.5, trim_top=0.01,
                        same_type=False):
    """Same rule as the paper (eqn 14). With same_type=True the alternative must also have the
    same number of elements and the same chemistry (anion or metal-only)."""
    rows = []
    for (app, sub), d in egj.groupby(["application", "subset"]):
        if app not in app_perf:
            continue
        d = d.copy()
        d["perf"] = perf_fn(d, app)
        if conf_col in d.columns:
            d = d[pd.to_numeric(d[conf_col], errors="coerce") >= conf_min]
        d = d.dropna(subset=["perf", risk]).drop_duplicates("gnome_id")
        if len(d) < 30:
            continue
        d = d[d["perf"] <= d["perf"].quantile(1 - trim_top)]
        anion = d["elements"].str.split(",").apply(lambda e: bool(set(e) & {"O", "S", "Se", "Te", "N", "F", "Cl", "Br", "I", "H"}))
        typ = (d["nelements"].astype(int).clip(upper=5).astype(str) + np.where(anion, "_anion", "_metal")).values
        mode = app_perf[app][1]
        p, r = d["perf"].values.astype(float), d[risk].values.astype(float)
        tol = rel_tol * np.abs(p) if mode == "higher" else np.full(len(p), abs_tol)
        hr = r >= high_risk
        n_alt = np.zeros(len(p), int)
        for i in np.where(hr)[0]:
            m = (np.abs(p - p[i]) <= tol[i]) & (r <= r[i] - min_risk_drop)
            if same_type:
                m &= typ == typ[i]
            n_alt[i] = int(m.sum())
        rows.append(dict(application=app, subset=sub, n=len(d), high_risk_n=int(hr.sum()),
                         share_with_1plus_alt=(n_alt[hr] >= 1).mean() if hr.any() else np.nan,
                         median_alts=float(np.median(n_alt[hr])) if hr.any() else np.nan,
                         conf_min=conf_min, same_type=same_type, tolerance=rel_tol))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# R5  robustness on the common support
# --------------------------------------------------------------------------------------
def robustness_common_support(master, dedupe_fn, ref="MP_known", min_ref=100):
    ok = support_types(master, ref, min_ref, "s_n5_class")
    base = master[master["s_n5_class"].isin(ok)]
    metrics = ["crm_eu_atom", "has_ree", "has_f_in_core", "has_pgm"]
    out = []
    def run(label, d, strata="s_n5_class"):
        T = d[d["group"] == "GNoME"]; R = d[d["group"] == ref]
        for m in metrics:
            r = matched_gap(T, R, m, strata)
            if r:
                out.append(dict(check=label, matched_on=strata, coverage=len(T) / (master.group == "GNoME").sum(), **r))
    run("common support (main)", base)
    run("duplicates removed", dedupe_fn(base, "e_decomp_mp" if "e_decomp_mp" in base else None))
    run("radioactive removed", base[~base["has_radioactive"].astype(bool)])
    if "crystal_system" in base.columns:
        x = base[base["crystal_system"].astype(str) != "unknown"].copy()
        x["s_n5_xtal"] = x["s_n5_class"] + "_" + x["crystal_system"].astype(str)
        okx = x.loc[x.group == ref, "s_n5_xtal"].value_counts()
        x = x[x["s_n5_xtal"].isin(okx.index[okx >= min_ref])]
        run("crystal system added (known >= 100 per cell)", x, "s_n5_xtal")
    return pd.DataFrame(out)


# --------------------------------------------------------------------------------------
# R6  element ratios with compound-type weights
# --------------------------------------------------------------------------------------
def element_ratios_type_matched(master, ref="MP_known", min_ref=100, strata="s_n5_class",
                                n_boot=1000, seed=0, min_count=20):
    rng = np.random.default_rng(seed)
    ok = support_types(master, ref, min_ref, strata)
    d = master[master[strata].isin(ok)]
    T, R = d[d["group"] == "GNoME"], d[d["group"] == ref]
    pi = T[strata].value_counts(normalize=True); pr = R[strata].value_counts(normalize=True)
    w = R[strata].map(lambda s: pi.get(s, 0.0) / pr[s]).values.astype(float)  # 0 for types GNoME lacks
    MT, els = ca.element_indicator(T); MR, _ = ca.element_indicator(R, els)
    MT, MR = MT.astype(np.float32), MR.astype(np.float32)
    eps = 0.5
    def ratio(wt, wr):
        pt = (wt @ MT + eps) / (wt.sum() + 2 * eps)
        prr = ((wr * w) @ MR + eps * w.mean()) / ((wr * w).sum() + 2 * eps * w.mean())
        return np.log2(pt / prr), pt, prr
    lr, pt, prr = ratio(np.ones(len(MT), np.float32), np.ones(len(MR), np.float32))
    boots = np.empty((n_boot, len(els)))
    for b in range(n_boot):
        boots[b] = ratio(rng.poisson(1.0, len(MT)).astype(np.float32), rng.poisson(1.0, len(MR)).astype(np.float32))[0]
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    res = pd.DataFrame({"element": els, "share_gnome": pt, "share_known_matched": prr,
                        "ratio": 2 ** lr, "log2_ratio": lr, "ratio_ci_low": 2 ** lo, "ratio_ci_high": 2 ** hi,
                        "count_gnome": MT.sum(0), "count_known": MR.sum(0)})
    res = res[(res.count_gnome + res.count_known) >= min_count]
    return res.sort_values("log2_ratio", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------------------
# R7  decomposition in both weight orders
# --------------------------------------------------------------------------------------
def decompose_both(master, metric, ref="MP_known", min_ref=100, strata="s_n5_class", n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    ok = support_types(master, ref, min_ref, strata)
    d = master[master["group"].isin(["GNoME", ref]) & master[strata].isin(ok)].dropna(subset=[metric])
    codes, uniq = pd.factorize(d[strata]); K = len(uniq)
    isT = (d["group"] == "GNoME").values; v = d[metric].astype(float).values
    def parts(wt_obs, wr_obs):
        nT = np.bincount(codes[isT], weights=wt_obs, minlength=K); nR = np.bincount(codes[~isT], weights=wr_obs, minlength=K)
        sT = np.bincount(codes[isT], weights=v[isT] * wt_obs, minlength=K); sR = np.bincount(codes[~isT], weights=v[~isT] * wr_obs, minlength=K)
        okk = (nT > 0) & (nR > 0)
        wT, wR = nT[okk] / nT[okk].sum(), nR[okk] / nR[okk].sum()
        mT, mR = sT[okk] / nT[okk], sR[okk] / nR[okk]
        gap = float((wT * mT).sum() - (wR * mR).sum())
        mix_a = float(((wT - wR) * mR).sum())        # order used in the paper
        mix_b = float(((wT - wR) * mT).sum())        # other order
        return gap, mix_a / gap if gap else np.nan, mix_b / gap if gap else np.nan
    nT_, nR_ = int(isT.sum()), int((~isT).sum())
    g, a, b = parts(np.ones(nT_), np.ones(nR_))
    bs = np.array([parts(rng.poisson(1.0, nT_), rng.poisson(1.0, nR_)) for _ in range(n_boot)])
    ci = lambda x: (np.nanpercentile(x, 2.5), np.nanpercentile(x, 97.5))
    avg = (bs[:, 1] + bs[:, 2]) / 2
    return dict(metric=metric, gap=g, mix_share_paper_order=a, ci_paper_order=ci(bs[:, 1]),
                mix_share_other_order=b, ci_other_order=ci(bs[:, 2]),
                mix_share_average=(a + b) / 2, ci_average=ci(avg), n_boot=n_boot)


# --------------------------------------------------------------------------------------
# R8  margins controlled for MP data density
# --------------------------------------------------------------------------------------
def mp_entries_per_system(master, G, mp_groups=("MP_known", "MP_theoretical")):
    """For each GNoME compound: number of MP entries whose elements are all inside its chemical
    system (its possible competing phases)."""
    mp = master[master["group"].isin(mp_groups)]
    cnt = mp["elements"].map(lambda s: "-".join(sorted(s.split(",")))).value_counts().to_dict()
    def count(s):
        el = sorted(s.split(","))
        return sum(cnt.get("-".join(c), 0) for k in range(1, len(el) + 1) for c in itertools.combinations(el, k))
    u = pd.Series(G["elements"].unique())
    return G["elements"].map(dict(zip(u, u.map(count))))


def margins_density_control(G, n_mp, flags=("has_f_in_core", "has_ree"), margin_col="e_decomp_mp",
                            thresholds=(10.0, 25.0), n_bins=4):
    d = G.dropna(subset=[margin_col]).copy()
    d["n_mp"] = n_mp.loc[d.index].values
    d["density_bin"] = pd.qcut(d["n_mp"].rank(method="first"), n_bins, labels=False).astype(str)
    d["s_dens"] = d["s_n_class"] + "_d" + d["density_bin"]
    d["margin_meV"] = -1000 * d[margin_col].astype(float)
    rows = []
    for f in flags:
        rows.append(dict(flag=f, measure="median MP entries in chemical system (with / without)",
                         value_with=float(d.loc[d[f], "n_mp"].median()), value_without=float(d.loc[~d[f], "n_mp"].median())))
        for thr in thresholds:
            d["lost"] = d["margin_meV"] < thr; d["stable_before"] = True
            for strata, lab in (("s_n_class", "type only (as in paper)"), ("s_dens", f"type x MP-density quartile")):
                o = ca.mantel_haenszel_or(d.assign(group="GNoME"), f, "GNoME", strata)
                rows.append(dict(flag=f, measure=f"odds ratio, margin < {thr:g} meV/atom, strata: {lab}",
                                 value_with=o.get("or_mh"), value_without=np.nan,
                                 ci_low=o.get("ci_low"), ci_high=o.get("ci_high")))
    return pd.DataFrame(rows)
