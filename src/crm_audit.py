"""
crm_audit.py — Does AI materials discovery over-produce critical-element compounds,
and are those "discoveries" artifacts of known errors in the training data?

Three tests:
  A. Criticality: are AI-proposed (GNoME) compounds more supply-risky than known
     compounds, after matching on number of elements?
  B. Link: are the elements that AI over-produces the same elements whose energies
     change most between GGA(+U) and the more accurate r2SCAN method?
  C. Survival: do AI-proposed compounds with those elements lose stability under
     r2SCAN more often than comparable known compounds?
Plus an optional "artifact risk" classifier, evaluated with grouped cross-validation.
"""
from __future__ import annotations

import os
import math
import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------------------------
# Element groups (edit here; every analysis reads from these lists)
# --------------------------------------------------------------------------------------
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
               "Ho", "Er", "Tm", "Yb", "Lu"]
REE = LANTHANIDES + ["Sc", "Y"]

# Lanthanides whose 4f electrons are frozen in the core ("_3"/"_2" PAW potentials) in
# Materials Project calculations. Verified in pymatgen's MP input sets: MPRelaxSet.yaml
# (GGA) and MPSCANRelaxSet.yaml (r2SCAN) both use Pr_3, Nd_3, Sm_3, Tb_3, Dy_3, Ho_3,
# Er_3, Tm_3, Lu_3 and Yb_2/Yb_3; La, Ce, Eu and Gd use 4f-in-valence potentials.
# Pm has no MP potential choice listed; it is kept with its neighbours.
F_IN_CORE = ["Pr", "Nd", "Pm", "Sm", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]

# Elements that get a Hubbard U correction in MP oxides/fluorides only (mixing scheme).
U_ELEMENTS = ["V", "Cr", "Mn", "Fe", "Co", "Ni", "Mo", "W"]

# Radioactive elements with no stable isotope (not usable as real materials at scale).
RADIOACTIVE = ["Tc", "Pm", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
               "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]

# EU 2023 Critical Raw Materials list, mapped to elements.
# Clear one-to-one cases only. Raw materials that do not map cleanly onto an element
# in a compound (bauxite->Al, fluorspar->F, natural graphite->C, silicon metal->Si,
# coking coal, feldspar, phosphate rock, helium) are kept apart in CRM_AMBIGUOUS so
# you can test both versions. VERIFY against the official EU list before publishing.
CRM_EU2023 = (["Sb", "As", "Ba", "Be", "Bi", "B", "Co", "Cu", "Ga", "Ge", "Hf",
               "Li", "Mg", "Mn", "Ni", "Nb", "P", "Sr", "Ta", "Ti", "W", "V",
               "Ru", "Rh", "Pd", "Os", "Ir", "Pt"] + REE)
CRM_AMBIGUOUS = ["Al", "F", "C", "Si", "He"]

THERMO_GGA = "GGA_GGA+U"
THERMO_MIX = "GGA_GGA+U_R2SCAN"
THERMO_R2 = "R2SCAN"


# --------------------------------------------------------------------------------------
# 1. Data loading (Materials Project via mp-api). Everything is cached to disk.
# --------------------------------------------------------------------------------------
SUMMARY_FIELDS = ["material_id", "formula_pretty", "composition_reduced", "elements",
                  "nelements", "nsites", "theoretical", "energy_above_hull",
                  "is_stable", "formation_energy_per_atom"]
THERMO_FIELDS = ["material_id", "thermo_type", "energy_above_hull", "is_stable",
                 "formation_energy_per_atom"]


def _cached(path, builder):
    """Load a parquet file if it exists, otherwise build it and save it."""
    if path and os.path.exists(path):
        return pd.read_parquet(path)
    df = builder()
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        df.to_parquet(path, index=False)
    return df


def _is_el(sym):
    from pymatgen.core.periodic_table import Element
    return Element.is_valid_symbol(sym)


def _summary_docs_to_df(docs):
    rows = []
    for d in docs:
        d = d if isinstance(d, dict) else d.model_dump()
        comp = d.get("composition_reduced")
        if hasattr(comp, "as_dict"):
            comp = comp.as_dict()
        comp = {str(k): float(v) for k, v in (comp or {}).items() if _is_el(str(k))}
        rows.append({
            "material_id": str(d["material_id"]),
            "formula": d.get("formula_pretty"),
            "composition": comp,
            "elements": sorted(str(getattr(e, "symbol", e)) for e in (d.get("elements") or [])),
            "nelements": d.get("nelements"),
            "nsites": d.get("nsites"),
            "theoretical": d.get("theoretical"),
            "e_hull_default": d.get("energy_above_hull"),
            "is_stable_default": d.get("is_stable"),
            "ef_default": d.get("formation_energy_per_atom"),
        })
    df = pd.DataFrame(rows)
    # parquet cannot store dicts with varying keys well -> store composition as JSON-like str
    df["composition"] = df["composition"].apply(lambda c: ";".join(f"{k}:{v}" for k, v in c.items()))
    df["elements"] = df["elements"].apply(lambda e: ",".join(e))
    return df


def parse_elements(s):
    return s.split(",") if isinstance(s, str) and s else []


def parse_composition(s):
    out = {}
    if isinstance(s, str) and s:
        for part in s.split(";"):
            k, v = part.split(":")
            out[k] = float(v)
    return out


def fetch_summary(api_key, include_gnome, cache_path=None, chunk_size=1000):
    def build():
        from mp_api.client import MPRester
        with MPRester(api_key, use_document_model=False, mute_progress_bars=False) as mpr:
            docs = mpr.materials.summary.search(
                fields=SUMMARY_FIELDS, deprecated=False,
                include_gnome=include_gnome, chunk_size=chunk_size)
        return _summary_docs_to_df(docs)
    return _cached(cache_path, build)


def fetch_thermo(api_key, thermo_type, cache_path=None, chunk_size=1000):
    def build():
        from mp_api.client import MPRester
        with MPRester(api_key, use_document_model=False, mute_progress_bars=False) as mpr:
            docs = mpr.materials.thermo.search(
                thermo_types=[thermo_type], fields=THERMO_FIELDS, chunk_size=chunk_size)
        rows = []
        for d in docs:
            d = d if isinstance(d, dict) else d.model_dump()
            rows.append({
                "material_id": str(d["material_id"]),
                "thermo_type": str(getattr(d.get("thermo_type"), "value", d.get("thermo_type"))),
                "e_hull": d.get("energy_above_hull"),
                "is_stable": d.get("is_stable"),
                "ef": d.get("formation_energy_per_atom"),
            })
        return pd.DataFrame(rows)
    return _cached(cache_path, build)


def label_groups(df_with_gnome, df_without_gnome):
    """GNoME = present only when GNoME is included. Known = MP non-theoretical
    (has an experimental ICSD record). Theoretical = other computed proposals."""
    ids_without = set(df_without_gnome["material_id"])
    df = df_with_gnome.copy()
    df["group"] = np.where(~df["material_id"].isin(ids_without), "GNoME",
                   np.where(df["theoretical"].astype(bool), "MP_theoretical", "MP_known"))
    return df


GNOME_CSV_URL = ("https://storage.googleapis.com/gdm_materials_discovery/"
                 "gnome_data/stable_materials_summary.csv")


def load_gnome_csv(url=GNOME_CSV_URL, cache_path=None):
    """DeepMind's own GNoME summary table (PBE energies from GNoME's hull).
    Column names are detected loosely because they may change between releases."""
    def build():
        return pd.read_csv(url)
    raw = _cached(cache_path, build)
    cols = {c.lower().replace(" ", "_"): c for c in raw.columns}
    def pick(*cands):
        for c in cands:
            if c in cols:
                return cols[c]
        raise KeyError(f"None of {cands} in columns: {list(raw.columns)}")
    out = pd.DataFrame({
        "gnome_id": raw[pick("materialid", "material_id")],
        "reduced_formula": raw[pick("reduced_formula", "reduced_formula_")],
        "nsites": raw[pick("nsites", "n_sites")],
        "gnome_e_decomp": raw[pick("decomposition_energy_per_atom")],
    })
    return out


def match_gnome_csv(df_meta, gnome_csv):
    """Attach GNoME's own decomposition energy to MP GNoME entries by matching
    reduced formula + number of sites. Ambiguous matches (several rows) are dropped."""
    from pymatgen.core import Composition
    g = df_meta[df_meta["group"] == "GNoME"].copy()
    g["reduced_formula"] = g["formula"].apply(lambda f: Composition(f).reduced_formula)
    c = gnome_csv.copy()
    c["reduced_formula"] = c["reduced_formula"].apply(lambda f: Composition(f).reduced_formula)
    c = c.drop_duplicates(["reduced_formula", "nsites"], keep=False)
    m = g.merge(c, on=["reduced_formula", "nsites"], how="inner")
    return m[["material_id", "gnome_e_decomp"]]


def coverage_table(df_meta, thermo_by_type):
    """How many compounds in each group have each thermo type. Read this first:
    it tells you which comparisons are actually possible."""
    rows = []
    for g, d in df_meta.groupby("group"):
        row = {"group": g, "n": len(d)}
        ids = set(d["material_id"])
        for t, th in thermo_by_type.items():
            row[t] = int(th["material_id"].isin(ids).sum())
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# 2. Compound-level criticality metrics
# --------------------------------------------------------------------------------------
def element_hhi_table():
    """Element HHI (production, reserves) from pymatgen (Gaultois et al., 2013).
    NOTE: these values reflect ~2011 USGS data; consider updating with current USGS."""
    from pymatgen.analysis.hhi import HHIModel
    h = HHIModel()
    return {el: (float(p), float(r)) for el, (p, r) in h.symbol_hhip_hhir.items()}


def add_compound_metrics(df, hhi=None, crm=None):
    """Adds: hhi_p / hhi_r (atom-fraction weighted), hhi_p_max, crm_frac, flags."""
    hhi = hhi or element_hhi_table()
    crm = set(crm or CRM_EU2023)
    f_core, ree, uel, rad = set(F_IN_CORE), set(REE), set(U_ELEMENTS), set(RADIOACTIVE)
    out = {k: [] for k in ["hhi_p", "hhi_r", "hhi_p_max", "crm_frac", "hhi_missing"]}
    for s in df["composition"].values:
        comp = parse_composition(s)
        tot = sum(comp.values()) or 1.0
        hp = hr = known = crm_amt = 0.0
        hmax = 0.0
        missing = False
        for el, amt in comp.items():
            frac = amt / tot
            if el in hhi:
                hp += frac * hhi[el][0]
                hr += frac * hhi[el][1]
                known += frac
                hmax = max(hmax, hhi[el][0])
            else:
                missing = True
            if el in crm:
                crm_amt += frac
        out["hhi_p"].append(hp / known if known else np.nan)
        out["hhi_r"].append(hr / known if known else np.nan)
        out["hhi_p_max"].append(hmax if known else np.nan)
        out["crm_frac"].append(crm_amt)
        out["hhi_missing"].append(missing)
    df = df.copy()
    for k, v in out.items():
        df[k] = v
    els = df["elements"].apply(lambda s: set(parse_elements(s)))
    df["has_ree"] = els.apply(lambda e: bool(e & ree))
    df["has_f_in_core"] = els.apply(lambda e: bool(e & f_core))
    df["has_u_element"] = els.apply(lambda e: bool(e & uel))
    df["has_radioactive"] = els.apply(lambda e: bool(e & rad))
    df["chemsys"] = df["elements"].apply(lambda s: "-".join(sorted(parse_elements(s))))
    return df


# --------------------------------------------------------------------------------------
# 3. Test A — criticality and element enrichment, matched on number of elements
# --------------------------------------------------------------------------------------
def _strata_weights(target, ref, strata_col):
    """Weights that make `ref` have the same distribution of `strata_col` as `target`
    (direct standardisation). Strata absent from ref are dropped from target."""
    t = target[strata_col].value_counts(normalize=True)
    r = ref[strata_col].value_counts(normalize=True)
    common = t.index.intersection(r.index)
    t = t[common] / t[common].sum()
    w = ref[strata_col].map(lambda s: t.get(s, 0.0) / r[s] if s in common else 0.0)
    return w.values, set(common)


def compare_metric(df, metric, target="GNoME", ref="MP_known", strata_col="nelements",
                   n_boot=500, seed=0):
    """Difference in mean `metric` between target and ref, with ref reweighted to
    the target's nelements mix. Bootstrap CI (resampling compounds in both groups)."""
    rng = np.random.default_rng(seed)
    T = df[df["group"] == target].dropna(subset=[metric])
    R = df[df["group"] == ref].dropna(subset=[metric])
    w, common = _strata_weights(T, R, strata_col)
    T = T[T[strata_col].isin(common)]
    tv, rv = T[metric].astype(float).values, R[metric].astype(float).values
    point_t, point_r = tv.mean(), np.average(rv, weights=w)
    boots = []
    for _ in range(n_boot):
        ti = rng.integers(0, len(tv), len(tv))
        ri = rng.integers(0, len(rv), len(rv))
        wr = w[ri]
        if wr.sum() == 0:
            continue
        boots.append(tv[ti].mean() - np.average(rv[ri], weights=wr))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"metric": metric, "target": target, "ref": ref,
            "mean_target": point_t, "mean_ref_matched": point_r,
            "mean_ref_raw": rv.mean(), "diff": point_t - point_r,
            "ci_low": lo, "ci_high": hi, "n_target": len(tv), "n_ref": len(rv)}


def element_indicator(df, elements=None):
    lists = df["elements"].apply(parse_elements)
    elements = elements or sorted({e for l in lists for e in l})
    idx = {e: i for i, e in enumerate(elements)}
    M = np.zeros((len(df), len(elements)), dtype=bool)
    for r, l in enumerate(lists):
        for e in l:
            if e in idx:
                M[r, idx[e]] = True
    return M, elements


def element_enrichment(df, target="GNoME", ref="MP_known", strata_col="nelements",
                       n_boot=200, seed=0, min_count=20):
    """For each element: share of target compounds containing it vs the
    nelements-matched share in ref. Returns log2 ratio with bootstrap CI."""
    rng = np.random.default_rng(seed)
    T = df[df["group"] == target]
    R = df[df["group"] == ref]
    w, common = _strata_weights(T, R, strata_col)
    T = T[T[strata_col].isin(common)]
    MT, els = element_indicator(T)
    MR, _ = element_indicator(R, els)
    eps = 0.5  # pseudo-count to keep log finite
    def ratio(mt, mr, wr):
        pt = (mt.sum(0) + eps) / (len(mt) + 2 * eps)
        pr = ((mr * wr[:, None]).sum(0) + eps * wr.mean()) / (wr.sum() + 2 * eps * wr.mean())
        return np.log2(pt / pr), pt, pr
    lr, pt, pr = ratio(MT, MR, w)
    boots = np.empty((n_boot, len(els)))
    for b in range(n_boot):
        ti = rng.integers(0, len(MT), len(MT))
        ri = rng.integers(0, len(MR), len(MR))
        boots[b] = ratio(MT[ti], MR[ri], w[ri])[0]
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    res = pd.DataFrame({"element": els, "share_target": pt, "share_ref_matched": pr,
                        "log2_enrichment": lr, "ci_low": lo, "ci_high": hi,
                        "count_target": MT.sum(0), "count_ref": MR.sum(0)})
    res = res[(res["count_target"] + res["count_ref"]) >= min_count]
    return res.sort_values("log2_enrichment", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------------------
# 4. Test B — link enrichment to a per-element "data reliability" signal
# --------------------------------------------------------------------------------------
def element_functional_sensitivity(df_meta, thermo_a, thermo_b, groups=("MP_known",),
                                   value="ef", min_count=10):
    """Per element: median |value_b - value_a| over compounds (in `groups`) that have
    both thermo types. Large = energies for this element change a lot between
    GGA(+U) and r2SCAN = less trustworthy training labels."""
    a = thermo_a[["material_id", value]].rename(columns={value: "a"})
    b = thermo_b[["material_id", value]].rename(columns={value: "b"})
    m = df_meta[df_meta["group"].isin(groups)][["material_id", "elements"]]
    j = m.merge(a, on="material_id").merge(b, on="material_id").dropna(subset=["a", "b"])
    j["absdiff"] = (j["b"] - j["a"]).abs()
    rows = {}
    for els, d in zip(j["elements"], j["absdiff"]):
        for e in parse_elements(els):
            rows.setdefault(e, []).append(d)
    out = pd.DataFrame([{"element": e, "sensitivity": float(np.median(v)), "n_pairs": len(v)}
                        for e, v in rows.items()])
    return out[out["n_pairs"] >= min_count].reset_index(drop=True)


def _rank_residual(y, x):
    ry, rx = stats.rankdata(y), stats.rankdata(x)
    beta = np.polyfit(rx, ry, 1)
    return ry - np.polyval(beta, rx)


def link_test(enrich, sens, data_counts=None, n_perm=10000, seed=0):
    """Spearman correlation between element enrichment and element sensitivity.
    If data_counts (element -> number of known compounds) is given, also a partial
    correlation that removes the effect of how much data each element has.
    Permutation p-values (two-sided)."""
    rng = np.random.default_rng(seed)
    j = enrich.merge(sens, on="element")
    if data_counts is not None:
        j["log_count"] = np.log10(j["element"].map(data_counts).fillna(1).clip(lower=1))
    x, y = j["log2_enrichment"].values, j["sensitivity"].values
    rho = stats.spearmanr(x, y).statistic
    perm = np.array([stats.spearmanr(x, rng.permutation(y)).statistic for _ in range(n_perm)])
    out = {"n_elements": len(j), "spearman_rho": rho,
           "p_perm": (np.sum(np.abs(perm) >= abs(rho)) + 1) / (n_perm + 1)}
    if data_counts is not None:
        rx, ry = _rank_residual(x, j["log_count"]), _rank_residual(y, j["log_count"])
        prho = stats.pearsonr(rx, ry).statistic
        pperm = np.array([stats.pearsonr(rx, rng.permutation(ry)).statistic for _ in range(n_perm)])
        out.update({"partial_rho_given_datacount": prho,
                    "p_perm_partial": (np.sum(np.abs(pperm) >= abs(prho)) + 1) / (n_perm + 1)})
    return out, j


def leave_one_out_link(enrich, sens, data_counts=None):
    """Robustness: drop each element once and recompute rho. A result driven by
    one element is fragile."""
    rows = []
    for e in enrich.merge(sens, on="element")["element"]:
        r, _ = link_test(enrich[enrich.element != e], sens[sens.element != e],
                         data_counts, n_perm=1)
        rows.append({"dropped": e, "rho": r["spearman_rho"],
                     "partial_rho": r.get("partial_rho_given_datacount")})
    return pd.DataFrame(rows)


def family_test(enrich, family=F_IN_CORE):
    """Mann-Whitney: are family elements more enriched than all other elements?"""
    inf = enrich[enrich.element.isin(family)]["log2_enrichment"]
    oth = enrich[~enrich.element.isin(family)]["log2_enrichment"]
    if len(inf) < 2 or len(oth) < 2:
        return {"n_family": len(inf), "note": "too few elements"}
    u = stats.mannwhitneyu(inf, oth, alternative="greater")
    return {"n_family": len(inf), "n_other": len(oth),
            "median_family": float(inf.median()), "median_other": float(oth.median()),
            "U": float(u.statistic), "p_one_sided": float(u.pvalue)}


# --------------------------------------------------------------------------------------
# 5. Test C — do proposals survive r2SCAN?
# --------------------------------------------------------------------------------------
def build_survival_table(df_meta, thermo_before, thermo_after, thr=0.0):
    """Joins 'before' (GGA-based hull) and 'after' (r2SCAN hull) energies above hull.
    stable_before: e_hull_before <= thr. lost: stable_before and e_hull_after > thr."""
    b = thermo_before[["material_id", "e_hull"]].rename(columns={"e_hull": "e_before"})
    a = thermo_after[["material_id", "e_hull"]].rename(columns={"e_hull": "e_after"})
    j = df_meta.merge(b, on="material_id").merge(a, on="material_id")
    j = j.dropna(subset=["e_before", "e_after"])
    j["stable_before"] = j["e_before"] <= thr + 1e-9
    j["lost"] = j["stable_before"] & (j["e_after"] > thr + 1e-9)
    return j


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, c - h, c + h)


def loss_rates(surv, flag_col="has_f_in_core", groups=("GNoME", "MP_theoretical", "MP_known")):
    """Share of 'stable before' compounds that lose stability, by group x flag."""
    s = surv[surv["stable_before"]]
    rows = []
    for g in groups:
        for f in (True, False):
            d = s[(s["group"] == g) & (s[flag_col] == f)]
            p, lo, hi = wilson(int(d["lost"].sum()), len(d))
            rows.append({"group": g, flag_col: f, "n": len(d), "lost": int(d["lost"].sum()),
                         "rate": p, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def mantel_haenszel_or(surv, flag_col="has_f_in_core", group="GNoME", strata_col="nelements"):
    """Odds ratio of losing stability for flagged vs not, within `group`, pooled over
    nelements strata (Mantel-Haenszel) with Robins-Breslow-Greenland 95% CI."""
    s = surv[(surv["stable_before"]) & (surv["group"] == group)]
    R = S = P_R = PS_QR = Q_S = 0.0
    for _, d in s.groupby(strata_col):
        a = ((d[flag_col]) & d["lost"]).sum(); b = ((d[flag_col]) & ~d["lost"]).sum()
        c = ((~d[flag_col]) & d["lost"]).sum(); dd = ((~d[flag_col]) & ~d["lost"]).sum()
        n = a + b + c + dd
        if n == 0 or (a + b) == 0 or (c + dd) == 0:
            continue
        r, s_ = a * dd / n, b * c / n
        p, q = (a + dd) / n, (b + c) / n
        R += r; S += s_
        P_R += p * r; PS_QR += p * s_ + q * r; Q_S += q * s_
    if R == 0 or S == 0:
        return {"or_mh": np.nan, "note": "not enough events"}
    orr = R / S
    se = math.sqrt(P_R / (2 * R * R) + PS_QR / (2 * R * S) + Q_S / (2 * S * S))
    return {"group": group, "flag": flag_col, "or_mh": orr,
            "ci_low": math.exp(math.log(orr) - 1.96 * se),
            "ci_high": math.exp(math.log(orr) + 1.96 * se)}


# --------------------------------------------------------------------------------------
# 6. Optional: artifact-risk classifier with honest evaluation
# --------------------------------------------------------------------------------------
def composition_features(df, elements):
    X = np.zeros((len(df), len(elements)))
    idx = {e: i for i, e in enumerate(elements)}
    for r, s in enumerate(df["composition"].values):
        comp = parse_composition(s)
        tot = sum(comp.values()) or 1
        for e, v in comp.items():
            if e in idx:
                X[r, idx[e]] = v / tot
    extra = df[["nelements", "nsites"]].astype(float).values
    return np.hstack([X, extra])


def evaluate_classifier(surv, group="GNoME", seed=0, n_splits=5, holdout_family=REE):
    """Predict 'lost' among stable_before compounds. Reports three evaluations:
      random K-fold (optimistic), GroupKFold by chemical system (honest),
      and train-without / test-on a held-out element family (hardest).
    Baselines: always-majority and the one-line rule 'contains f-in-core lanthanide'."""
    from sklearn.model_selection import KFold, GroupKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score, average_precision_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    s = surv[(surv["stable_before"]) & (surv["group"] == group)].reset_index(drop=True)
    y = s["lost"].astype(int).values
    if y.sum() < 10 or (len(y) - y.sum()) < 10:
        return pd.DataFrame([{"note": "too few positive/negative examples"}])
    els = sorted({e for l in s["elements"].apply(parse_elements) for e in l})
    X = composition_features(s, els)
    models = {
        "logistic": lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
        "grad_boost": lambda: HistGradientBoostingClassifier(random_state=seed, max_iter=300,
                                                             learning_rate=0.05),
    }
    rule = s["has_f_in_core"].astype(float).values
    rows = []

    def score(name, split, yt, p):
        if len(np.unique(yt)) < 2:
            return
        rows.append({"model": name, "split": split, "roc_auc": roc_auc_score(yt, p),
                     "pr_auc": average_precision_score(yt, p), "base_rate": yt.mean(),
                     "n_test": len(yt)})

    splits = {
        "random_kfold": list(KFold(n_splits, shuffle=True, random_state=seed).split(X)),
        "group_by_chemsys": list(GroupKFold(n_splits).split(X, y, s["chemsys"])),
    }
    fam = s["elements"].apply(lambda l: bool(set(parse_elements(l)) & set(holdout_family))).values
    splits["holdout_family"] = [(np.where(~fam)[0], np.where(fam)[0])]

    for split_name, folds in splits.items():
        for mname, mk in models.items():
            yt_all, p_all = [], []
            for tr, te in folds:
                if len(np.unique(y[tr])) < 2:
                    continue
                m = mk().fit(X[tr], y[tr])
                yt_all.append(y[te]); p_all.append(m.predict_proba(X[te])[:, 1])
            if yt_all:
                score(mname, split_name, np.concatenate(yt_all), np.concatenate(p_all))
        te_all = np.concatenate([te for _, te in folds])
        score("rule_f_in_core", split_name, y[te_all], rule[te_all])
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# 7. Stand-in data for testing the pipeline without an API key
# --------------------------------------------------------------------------------------
def make_synthetic(n_known=6000, n_theo=4000, n_gnome=8000, effect=True, seed=0):
    """Fake data with the same columns as the real pipeline. If effect=True, GNoME
    over-samples f-in-core lanthanides AND those compounds lose stability under
    r2SCAN more often. If effect=False, there is no planted signal: tests should
    come back null. ONLY for checking the code; never report these numbers."""
    rng = np.random.default_rng(seed)
    hhi = element_hhi_table()
    pool = [e for e in hhi if e not in RADIOACTIVE and e not in ("He",)]
    base_w = np.ones(len(pool))
    common = {"O": 8, "Li": 3, "Na": 3, "Mg": 3, "Al": 3, "Si": 3, "P": 3, "S": 3,
              "Fe": 3, "Cu": 2, "Zn": 2, "Ca": 3, "K": 2, "Ti": 2, "Mn": 2, "Co": 2, "Ni": 2}
    for i, e in enumerate(pool):
        base_w[i] *= common.get(e, 1.0)
    gn_w = base_w.copy()
    if effect:
        for i, e in enumerate(pool):
            if e in F_IN_CORE:
                gn_w[i] *= 3.0
    sens_el = {e: max(0.005, (0.12 if (effect and e in F_IN_CORE) else 0.04) + rng.normal(0, 0.01))
               for e in pool}

    def gen(n, group, w, nel_p):
        rows = []
        p = w / w.sum()
        for k in range(n):
            ne = rng.choice([1, 2, 3, 4, 5], p=nel_p)
            els = list(rng.choice(pool, size=ne, replace=False, p=p))
            comp = {e: float(rng.integers(1, 5)) for e in els}
            rows.append({"material_id": f"{group}-{k}", "formula": "".join(els),
                         "composition": ";".join(f"{e}:{v}" for e, v in comp.items()),
                         "elements": ",".join(sorted(els)), "nelements": ne,
                         "nsites": int(sum(comp.values()) * rng.integers(1, 4)),
                         "theoretical": group != "known", "group_raw": group})
        return rows
    rows = (gen(n_known, "known", base_w, [0.05, 0.35, 0.4, 0.15, 0.05]) +
            gen(n_theo, "theo", base_w, [0.02, 0.3, 0.45, 0.18, 0.05]) +
            gen(n_gnome, "gnome", gn_w, [0.01, 0.1, 0.35, 0.39, 0.15]))
    df = pd.DataFrame(rows)
    df["group"] = df["group_raw"].map({"known": "MP_known", "theo": "MP_theoretical",
                                       "gnome": "GNoME"})
    # energies
    ef_gga = rng.normal(-1.5, 0.6, len(df))
    s = df["elements"].apply(lambda l: np.mean([sens_el[e] for e in parse_elements(l)])).values
    shift = rng.normal(0, s)
    fcore = df["elements"].apply(lambda l: bool(set(parse_elements(l)) & set(F_IN_CORE))).values
    is_gn = (df["group"] == "GNoME").values
    e_before = np.where(is_gn, np.abs(rng.normal(0, 0.005, len(df))) * (rng.random(len(df)) < 0.3),
                        np.abs(rng.normal(0, 0.03, len(df))) * (rng.random(len(df)) < 0.5))
    extra = 0.02 + (0.08 if effect else 0.0) * (fcore & is_gn)
    e_after = np.clip(e_before + rng.normal(extra, 0.03, len(df)), 0, None) * (rng.random(len(df)) < 0.85)
    tg = pd.DataFrame({"material_id": df.material_id, "thermo_type": THERMO_GGA,
                       "e_hull": e_before, "ef": ef_gga})
    tr = pd.DataFrame({"material_id": df.material_id, "thermo_type": THERMO_R2,
                       "e_hull": e_after, "ef": ef_gga + shift})
    return df.drop(columns="group_raw"), tg, tr