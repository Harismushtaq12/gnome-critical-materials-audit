"""
crm_v2.py — extended audit (Levels 1-4). Builds on crm_audit.py (v1).

Level 1  make the main result solid: two critical lists (EU 2023, US 2022), atom and
         mass fractions, duplicates removed, chemistry/crystal matching, a "uniform
         element choice" expectation, null calibration.
Level 2  explain the cause: gap decomposition (chemistry mix vs within-type choice),
         inheritance from computed data, near-copy templates, stability margins,
         r2SCAN re-check from DeepMind's own file.
Level 3  real-world relevance: Energy-GNoME application shortlists; non-GNoME
         high-throughput baselines from DeepMind's external hull file (OQMD, WBM, ...).
Level 4  own calculations: re-check a sample with a universal ML potential.
"""
from __future__ import annotations

import io
import os
import json
import math
import zipfile
from math import comb

import numpy as np
import pandas as pd
from scipy import stats

import crm_audit as ca

# --------------------------------------------------------------------------------------
# Lists and constants
# --------------------------------------------------------------------------------------
CRM_EU = list(ca.CRM_EU2023)                 # EU 2023 (element-mappable entries)
CRM_EU_AMBIG = list(ca.CRM_AMBIGUOUS)         # Al, F, C, Si, He
# US 2022 Final List of Critical Minerals (USGS, Federal Register 87 FR 10381).
CRM_US = ["Al", "Sb", "As", "Ba", "Be", "Bi", "Ce", "Cs", "Cr", "Co", "Dy", "Er", "Eu",
          "F", "Gd", "Ga", "Ge", "C", "Hf", "Ho", "In", "Ir", "La", "Li", "Lu", "Mg",
          "Mn", "Nd", "Ni", "Nb", "Pd", "Pt", "Pr", "Rh", "Rb", "Ru", "Sm", "Sc", "Ta",
          "Te", "Tb", "Tm", "Sn", "Ti", "W", "V", "Yb", "Y", "Zn", "Zr"]
US_AMBIG = ["Al", "F", "C"]                   # bauxite/aluminum, fluorspar, graphite
CRM_US_STRICT = [e for e in CRM_US if e not in US_AMBIG]
ANIONS = {"O", "S", "Se", "Te", "N", "F", "Cl", "Br", "I", "H"}
NOBLE = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}

GNOME_BASE = "https://storage.googleapis.com/gdm_materials_discovery/"
GNOME_SUMMARY_URL = GNOME_BASE + "gnome_data/stable_materials_summary.csv"
GNOME_R2SCAN_URL = GNOME_BASE + "gnome_data/stable_materials_r2scan.csv"
GNOME_BYID_ZIP_URL = GNOME_BASE + "gnome_data/by_id.zip"
GNOME_EXTERNAL_URL = GNOME_BASE + "external_data/external_materials_summary.csv"

EG_BASE = "https://raw.githubusercontent.com/paolodeangelis/Energy-GNoME/main/data/final/"
EG_FILES = ([("cathode", ion, f"cathodes/insertion/{ion}")
             for ion in ["Al", "Ca", "Cs", "K", "Li", "Mg", "Na", "Rb", "Y"]] +
            [("perovskite", "pure", "perovskites/pure_models"),
             ("perovskite", "mixed", "perovskites/mixed_models")] +
            [("thermoelectric", t, f"thermoelectrics/{t}")
             for t in ["300K", "430K", "560K", "690K", "820K", "950K"]])


# --------------------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------------------
def cached_parquet(path, builder):
    if path and os.path.exists(path):
        return pd.read_parquet(path)
    df = builder()
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        df.to_parquet(path, index=False)
    return df


def find_col(df, *keys, required=True):
    """Find a column by loose name match (case, spaces, underscores ignored)."""
    norm = {c: c.lower().replace(" ", "").replace("_", "").replace("(", "").replace(")", "")
            for c in df.columns}
    for k in keys:
        k = k.lower().replace(" ", "").replace("_", "")
        for c, n in norm.items():
            if n == k:
                return c
    for k in keys:
        k = k.lower().replace(" ", "").replace("_", "")
        for c, n in norm.items():
            if k in n:
                return c
    if required:
        raise KeyError(f"none of {keys} in {list(df.columns)}")
    return None


def comps_to_strings(comps):
    """pymatgen Compositions -> ('El:amt;...', 'A,B,C', nelements, reduced formula)."""
    cs, es, ns, rf = [], [], [], []
    for c in comps:
        d = {e.symbol: float(a) for e, a in c.items()}
        cs.append(";".join(f"{k}:{v:g}" for k, v in d.items()))
        es.append(",".join(sorted(d)))
        ns.append(len(d))
        rf.append(c.reduced_formula)
    return cs, es, ns, rf


def parse_formulas(formulas):
    from pymatgen.core import Composition
    out = []
    for f in formulas:
        try:
            out.append(Composition(str(f)))
        except Exception:
            out.append(None)
    return out


# --------------------------------------------------------------------------------------
# 1. Data loading
# --------------------------------------------------------------------------------------
def load_mp_compact(api_key, cache_dir, local_lake="/root/mp_lake"):
    """All Materials Project materials in one fast download (data lakehouse), reduced
    to the columns we need and saved as one small parquet file in cache_dir."""
    out = os.path.join(cache_dir, "mp_compact.parquet")

    def build():
        import glob
        import pyarrow as pa
        import pyarrow.dataset as pads
        from mp_api.client import MPRester
        files = [f for f in glob.glob(f"{local_lake}/**/*.parquet", recursive=True)
                 if "summary" in f and "_delta_log" not in f]
        if not files:
            with MPRester(api_key, local_dataset_cache=local_lake) as mpr:
                _ = mpr.materials.summary.search(fields=["material_id"], include_gnome=True)
            files = [f for f in glob.glob(f"{local_lake}/**/*.parquet", recursive=True)
                     if "summary" in f and "_delta_log" not in f]
        dset = pads.dataset(files, format="parquet")
        names = dset.schema.names
        want = [c for c in ["material_id", "deprecated", "theoretical", "nsites",
                            "formula_pretty", "energy_above_hull",
                            "formation_energy_per_atom", "symmetry"] if c in names]
        tab = dset.to_table(columns=want)
        df = tab.drop(["symmetry"]).to_pandas() if "symmetry" in want else tab.to_pandas()
        if "symmetry" in want:
            sym = tab.column("symmetry").combine_chunks()
            if pa.types.is_struct(sym.type):
                cs = sym.field("crystal_system").to_pandas()
            else:
                cs = pd.Series([(x or {}).get("crystal_system") for x in sym.to_pylist()])
            df["crystal_system"] = cs.astype(str).str.lower().values
        df["material_id"] = df["material_id"].astype(str)
        return df.drop_duplicates("material_id").reset_index(drop=True)
    return cached_parquet(out, build)


def mp_to_master(mp):
    """Materials Project compact table -> common schema with groups."""
    mp = mp[mp["deprecated"].astype(str) != "True"].copy()
    comps = parse_formulas(mp["formula_pretty"])
    ok = [c is not None for c in comps]
    mp, comps = mp[ok].copy(), [c for c in comps if c is not None]
    cs, es, ns, rf = comps_to_strings(comps)
    out = pd.DataFrame({
        "material_id": mp["material_id"].values, "formula": rf, "composition": cs,
        "elements": es, "nelements": ns, "nsites": mp["nsites"].values,
        "theoretical": mp["theoretical"].astype(bool).values,
        "crystal_system": mp.get("crystal_system", pd.Series(["unknown"] * len(mp))).values,
        "e_hull_mp": mp["energy_above_hull"].values})
    out["group"] = np.where(out["theoretical"], "MP_theoretical", "MP_known")
    return out


def load_gnome_raw(cache_dir, url=GNOME_SUMMARY_URL):
    return cached_parquet(os.path.join(cache_dir, "gnome_summary_raw.parquet"),
                          lambda: pd.read_csv(url))


def gnome_to_master(graw):
    comps = parse_formulas(graw[find_col(graw, "reducedformula")])
    keep = np.array([c is not None for c in comps])
    g = graw[keep].reset_index(drop=True)
    comps = [c for c in comps if c is not None]
    cs, es, ns, rf = comps_to_strings(comps)

    def opt(*k):
        c = find_col(g, *k, required=False)
        return g[c].values if c else np.full(len(g), np.nan)
    out = pd.DataFrame({
        "material_id": "gnome-" + g[find_col(g, "materialid")].astype(str).values,
        "gnome_id": g[find_col(g, "materialid")].astype(str).values,
        "formula": rf, "composition": cs, "elements": es, "nelements": ns,
        "nsites": g[find_col(g, "nsites")].values, "theoretical": True,
        "crystal_system": pd.Series(opt("crystalsystem")).astype(str).str.lower().values,
        "is_train": opt("istrain"),
        "e_decomp": opt("decompositionenergyperatom"),
        "e_decomp_all": opt("decompositionenergyperatomall"),
        "e_decomp_rel": opt("decompositionenergyperatomrelative"),
        "e_decomp_mp": opt("decompositionenergyperatommp"),
        "e_decomp_mp_oqmd": opt("decompositionenergyperatommpoqmd")})
    out["group"] = "GNoME"
    return out


def load_external(cache_dir, url=GNOME_EXTERNAL_URL):
    return cached_parquet(os.path.join(cache_dir, "gnome_external_raw.parquet"),
                          lambda: pd.read_csv(url))


def external_to_master(xraw):
    """DeepMind's external hull file (MP / OQMD / WBM / ... entries). Source is read
    from the id prefix when possible."""
    fcol = find_col(xraw, "reducedformula", "composition", "formula")
    comps = parse_formulas(xraw[fcol])
    keep = np.array([c is not None for c in comps])
    x = xraw[keep].reset_index(drop=True)
    comps = [c for c in comps if c is not None]
    cs, es, ns, rf = comps_to_strings(comps)
    idc = find_col(x, "materialid", "entryid", "id", required=False)
    ids = x[idc].astype(str) if idc else pd.Series([f"ext{i}" for i in range(len(x))])
    src = ids.str.extract(r"^([A-Za-z]+)", expand=False).fillna("other").str.lower()
    nsc = find_col(x, "nsites", required=False)
    out = pd.DataFrame({
        "material_id": "ext-" + ids.values, "formula": rf, "composition": cs,
        "elements": es, "nelements": ns,
        "nsites": x[nsc].values if nsc else np.nan, "theoretical": True,
        "crystal_system": "unknown", "source": src.values})
    out["group"] = "EXT_" + out["source"].str.upper()
    return out


def load_gnome_r2scan(cache_dir, url=GNOME_R2SCAN_URL):
    return cached_parquet(os.path.join(cache_dir, "gnome_r2scan_raw.parquet"),
                          lambda: pd.read_csv(url))


def load_energy_gnome(cache_dir):
    def build():
        import urllib.request
        rows = []
        for app, sub, path in EG_FILES:
            with urllib.request.urlopen(EG_BASE + path + "/candidates.json") as r:
                obj = json.loads(r.read().decode())
            d = pd.DataFrame(obj)
            d["application"], d["subset"] = app, sub
            rows.append(d)
        d = pd.concat(rows, ignore_index=True)
        for c in d.columns:            # keep every column; make text columns parquet-safe
            if d[c].dtype == object and c != "Material Id":
                num = pd.to_numeric(d[c], errors="coerce")
                d[c] = num if num.notna().mean() > 0.9 else d[c].astype(str)
        return d.astype({"Material Id": str})
    return cached_parquet(os.path.join(cache_dir, "energy_gnome_full.parquet"), build)


# --------------------------------------------------------------------------------------
# 2. Metrics
# --------------------------------------------------------------------------------------
def element_masses():
    from pymatgen.core import Element
    out = {}
    for el in Element:
        try:
            out[el.symbol] = float(el.atomic_mass)
        except Exception:
            pass
    return out


def load_hhi(override_csv=None):
    """Element HHI (production, reserves). Default: pymatgen (Gaultois et al. 2013,
    ~2011 USGS data). Override with a CSV: element,hhi_p[,hhi_r] (0-10000 scale)."""
    hhi = ca.element_hhi_table()
    if override_csv and os.path.exists(override_csv):
        o = pd.read_csv(override_csv)
        for _, r in o.iterrows():
            p = float(r["hhi_p"])
            rr = float(r["hhi_r"]) if "hhi_r" in o.columns and pd.notna(r["hhi_r"]) else hhi.get(r["element"], (p, p))[1]
            hhi[str(r["element"])] = (p, rr)
    return hhi


def add_metrics(df, hhi=None):
    """All per-compound measures, in one pass."""
    hhi = hhi or load_hhi()
    mass = element_masses()
    sets = {"eu": set(CRM_EU), "eu_amb": set(CRM_EU) | set(CRM_EU_AMBIG),
            "us": set(CRM_US_STRICT), "us_amb": set(CRM_US)}
    cols = {f"crm_{k}_{w}": [] for k in sets for w in ("atom", "mass")}
    cols.update({"hhi_p": [], "hhi_r": [], "hhi_p_mass": [], "hhi_p_max": []})
    for s in df["composition"].values:
        comp = ca.parse_composition(s)
        tot_a = sum(comp.values()) or 1.0
        tot_m = sum(v * mass.get(e, 0.0) for e, v in comp.items()) or 1.0
        for k, st in sets.items():
            a = sum(v for e, v in comp.items() if e in st)
            m = sum(v * mass.get(e, 0.0) for e, v in comp.items() if e in st)
            cols[f"crm_{k}_atom"].append(a / tot_a)
            cols[f"crm_{k}_mass"].append(m / tot_m)
        known = [(e, v) for e, v in comp.items() if e in hhi]
        ka = sum(v for _, v in known)
        km = sum(v * mass.get(e, 0.0) for e, v in known)
        cols["hhi_p"].append(sum(v * hhi[e][0] for e, v in known) / ka if ka else np.nan)
        cols["hhi_r"].append(sum(v * hhi[e][1] for e, v in known) / ka if ka else np.nan)
        cols["hhi_p_mass"].append(sum(v * mass.get(e, 0) * hhi[e][0] for e, v in known) / km if km else np.nan)
        cols["hhi_p_max"].append(max((hhi[e][0] for e, _ in known), default=np.nan))
    df = df.copy()
    for k, v in cols.items():
        df[k] = v
    df["crm_frac"] = df["crm_eu_atom"]           # name used by v1 functions
    els = df["elements"].str.split(",").apply(set)
    df["has_ree"] = els.apply(lambda e: bool(e & set(ca.REE)))
    df["has_f_in_core"] = els.apply(lambda e: bool(e & set(ca.F_IN_CORE)))
    df["has_u_element"] = els.apply(lambda e: bool(e & set(ca.U_ELEMENTS)))
    df["has_radioactive"] = els.apply(lambda e: bool(e & set(ca.RADIOACTIVE)))
    df["has_pgm"] = els.apply(lambda e: bool(e & {"Ru", "Rh", "Pd", "Os", "Ir", "Pt"}))
    df["chem_class"] = np.where(els.apply(lambda e: bool(e & ANIONS)), "anion", "intermetallic")
    df["chemsys"] = df["elements"].str.replace(",", "-")
    df["s_n"] = df["nelements"].astype(str)
    df["s_n_class"] = df["s_n"] + "_" + df["chem_class"]
    df["s_n_class_xtal"] = df["s_n_class"] + "_" + df["crystal_system"].astype(str)
    return df


def dedupe_by_formula(df, prefer_col=None):
    """One row per (group, reduced formula). Keeps the lowest prefer_col if given."""
    d = df
    if prefer_col and prefer_col in d.columns:
        d = d.sort_values(prefer_col, na_position="last")
    return d.drop_duplicates(["group", "formula"]).reset_index(drop=True)


# --------------------------------------------------------------------------------------
# 3. Level 1 analyses
# --------------------------------------------------------------------------------------
HEADLINE_METRICS = ["crm_eu_atom", "crm_eu_mass", "crm_us_atom", "crm_us_mass",
                    "has_ree", "has_f_in_core", "has_pgm", "hhi_p", "hhi_p_mass"]


def compare_grid(df, metrics, refs, strata_list, target="GNoME", n_boot=200, seed=0,
                 label=""):
    rows = []
    for ref in refs:
        if (df["group"] == ref).sum() == 0:
            continue
        for st in strata_list:
            for m in metrics:
                d = df[df["group"].isin([target, ref])]
                d = d[~d[st].astype(str).str.contains("unknown|nan|None", regex=True)] \
                    if st == "s_n_class_xtal" else d
                try:
                    r = ca.compare_metric(d.assign(**{m: d[m].astype(float)}), m, target, ref,
                                          strata_col=st, n_boot=n_boot, seed=seed)
                except (ZeroDivisionError, ValueError, IndexError):
                    continue          # no overlapping strata (e.g. no crystal data)
                rows.append(dict(check=label, ref=ref, matched_on=st, metric=m,
                                 target=r["mean_target"], ref_matched=r["mean_ref_matched"],
                                 ref_raw=r["mean_ref_raw"], diff=r["diff"],
                                 ci_low=r["ci_low"], ci_high=r["ci_high"],
                                 n_target=r["n_target"], n_ref=r["n_ref"]))
    return pd.DataFrame(rows)


def uniform_expectation(df, groups, family, label):
    """If each compound picked its k elements uniformly at random from the element
    pool, how often would it contain at least one member of `family`?
    Compared with what each group actually shows."""
    pool = sorted({e for s in df["elements"] for e in s.split(",")} - NOBLE)
    E, F = len(pool), len(set(family) & set(pool))
    p_none = {k: comb(E - F, k) / comb(E, k) for k in range(1, 12)}
    rows = []
    fam = set(family)
    for g in groups:
        d = df[df["group"] == g]
        if not len(d):
            continue
        k = d["nelements"].clip(upper=11).astype(int)
        expected = (1 - k.map(p_none)).mean()
        observed = d["elements"].str.split(",").apply(lambda e: bool(set(e) & fam)).mean()
        rows.append(dict(family=label, group=g, observed=observed, uniform_expected=expected,
                         ratio=observed / expected, pool_size=E, family_size=F))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# 4. Level 2 analyses
# --------------------------------------------------------------------------------------
def decompose_gap(df, metric, target, ref, strata_col, n_boot=200, seed=0):
    """Splits the gap in `metric` between target and ref into:
       mix    = from having a different share of each stratum (e.g. more metal-only)
       within = from differences inside the same stratum.
    gap = sum_s wT(s)mT(s) - sum_s wR(s)mR(s)
        = sum_s wT(s)[mT(s)-mR(s)]  (within)  +  sum_s [wT(s)-wR(s)] mR(s)  (mix).
    Only strata present in both groups are used."""
    rng = np.random.default_rng(seed)
    d = df[df["group"].isin([target, ref])].dropna(subset=[metric])
    codes, uniq = pd.factorize(d[strata_col])
    isT = (d["group"] == target).values
    v = d[metric].astype(float).values
    K = len(uniq)

    def parts(idxT, idxR):
        nT = np.bincount(codes[idxT], minlength=K).astype(float)
        nR = np.bincount(codes[idxR], minlength=K).astype(float)
        sT = np.bincount(codes[idxT], weights=v[idxT], minlength=K)
        sR = np.bincount(codes[idxR], weights=v[idxR], minlength=K)
        ok = (nT > 0) & (nR > 0)
        wT, wR = nT[ok] / nT[ok].sum(), nR[ok] / nR[ok].sum()
        mT, mR = sT[ok] / nT[ok], sR[ok] / nR[ok]
        within = float((wT * (mT - mR)).sum())
        mix = float(((wT - wR) * mR).sum())
        return within + mix, mix, within
    iT, iR = np.where(isT)[0], np.where(~isT)[0]
    tot, mix, within = parts(iT, iR)
    bs = np.array([parts(rng.choice(iT, len(iT)), rng.choice(iR, len(iR)))
                   for _ in range(n_boot)])
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    share = bs[:, 1] / np.where(bs[:, 0] == 0, np.nan, bs[:, 0])
    return dict(metric=metric, strata=strata_col, gap=tot, mix=mix, within=within,
                mix_share=mix / tot if tot else np.nan,
                mix_share_ci=(np.nanpercentile(share, 2.5), np.nanpercentile(share, 97.5)),
                gap_ci=(lo[0], hi[0]), mix_ci=(lo[1], hi[1]), within_ci=(lo[2], hi[2]))


def amplification(df, n_boot=50, seed=0):
    """Per element: enrichment of computed MP proposals vs known (an existing
    computational bias) against enrichment of GNoME vs known. Slope > 1 = GNoME
    amplifies a bias that already existed in computed data."""
    e1 = ca.element_enrichment(df, "MP_theoretical", "MP_known", n_boot=n_boot, seed=seed)
    e2 = ca.element_enrichment(df, "GNoME", "MP_known", n_boot=n_boot, seed=seed)
    j = e1[["element", "log2_enrichment"]].merge(
        e2[["element", "log2_enrichment"]], on="element", suffixes=("_mp_theo", "_gnome"))
    rho = stats.spearmanr(j["log2_enrichment_mp_theo"], j["log2_enrichment_gnome"]).statistic
    ts = stats.theilslopes(j["log2_enrichment_gnome"], j["log2_enrichment_mp_theo"])
    return dict(n_elements=len(j), spearman_rho=rho, theil_slope=ts.slope,
                slope_ci=(ts.low_slope, ts.high_slope)), j


def ln_templates(df, groups=("GNoME", "MP_known", "MP_theoretical")):
    LN = set(ca.REE)

    def key(s):
        c = ca.parse_composition(s)
        return ";".join(sorted(f"{'Ln' if e in LN else e}:{v:g}" for e, v in c.items()))
    rows = []
    for g in groups:
        d = df[(df["group"] == g) & df["has_ree"]]
        if not len(d):
            continue
        sizes = d["composition"].map(key).value_counts()
        rows.append(dict(group=g, ree_compounds=len(d), templates=len(sizes),
                         copies_per_template=len(d) / len(sizes),
                         share_templates_5plus=(sizes >= 5).mean(),
                         max_copies=int(sizes.max())))
    return pd.DataFrame(rows)


def margin_table(g, flags=("has_f_in_core", "has_ree", "has_u_element"),
                 margin_col="e_decomp_rel", thresholds=(0.010, 0.025)):
    """GNoME stability margins from DeepMind's own columns. margin = how far below the
    hull of all OTHER entries (meV/atom). 'Fragile' = margin below threshold.
    Compared between compounds with and without each flag (MH odds ratio over
    number of elements and chemistry type)."""
    d = g.dropna(subset=[margin_col]).copy()
    d["margin_meV"] = -1000 * d[margin_col].astype(float)
    rows = []
    for f in flags:
        for thr in thresholds:
            d["stable_before"] = True
            d["lost"] = d["margin_meV"] < thr * 1000
            for val in (True, False):
                s = d[d[f] == val]
                p, lo, hi = ca.wilson(int(s["lost"].sum()), len(s))
                rows.append(dict(flag=f, threshold_meV=thr * 1000, flag_value=val, n=len(s),
                                 median_margin_meV=s["margin_meV"].median(),
                                 fragile_share=p, ci_low=lo, ci_high=hi))
            o = ca.mantel_haenszel_or(d.assign(group="GNoME"), f, "GNoME", "s_n_class")
            rows.append(dict(flag=f, threshold_meV=thr * 1000, flag_value="odds_ratio",
                             n=len(d), fragile_share=o.get("or_mh"),
                             ci_low=o.get("ci_low"), ci_high=o.get("ci_high")))
    return pd.DataFrame(rows), d


def r2scan_check(r2raw, gmaster, thresholds=(0.001, 0.025)):
    """Uses DeepMind's r2SCAN file. Finds the id and a decomposition/hull column by
    name; if none is found it returns the column list so the analysis can be adapted."""
    idc = find_col(r2raw, "materialid", "material id", "id", required=False)
    cand = [c for c in r2raw.columns if any(k in c.lower() for k in
            ("decomposition", "hull", "e_above", "ehull"))]
    info = dict(columns=list(r2raw.columns), id_col=idc, stability_cols=cand, n=len(r2raw))
    if idc is None or not cand:
        return info, None, None
    col = cand[0]
    j = gmaster.merge(r2raw[[idc, col]].rename(columns={idc: "gnome_id", col: "e_r2"})
                      .astype({"gnome_id": str}), on="gnome_id", how="inner")
    info.update(used_col=col, matched=len(j))
    rows = []
    for thr in thresholds:
        j["stable_before"] = True
        j["lost"] = j["e_r2"].astype(float) > thr
        for f in ("has_f_in_core", "has_ree", "has_u_element"):
            for val in (True, False):
                s = j[j[f] == val]
                p, lo, hi = ca.wilson(int(s["lost"].sum()), len(s))
                rows.append(dict(flag=f, threshold_eV=thr, flag_value=val, n=len(s),
                                 unstable_share=p, ci_low=lo, ci_high=hi))
            o = ca.mantel_haenszel_or(j.assign(group="GNoME"), f, "GNoME", "s_n_class")
            rows.append(dict(flag=f, threshold_eV=thr, flag_value="odds_ratio", n=len(j),
                             unstable_share=o.get("or_mh"), ci_low=o.get("ci_low"),
                             ci_high=o.get("ci_high")))
    return info, pd.DataFrame(rows), j


# --------------------------------------------------------------------------------------
# 5. Level 3 analyses
# --------------------------------------------------------------------------------------
def energy_gnome_table(eg, master, metrics=("crm_eu_atom", "crm_us_atom", "has_ree",
                       "has_f_in_core", "has_pgm", "hhi_p"), conf_min=0.9, n_boot=200, seed=0):
    """Criticality of GNoME candidates shortlisted for real applications, compared with
    all GNoME and with known materials (matched on number of elements)."""
    g = master[master["group"] == "GNoME"]
    conf = "AI-experts confidence (-)"
    rows, joined = [], []
    for app in sorted(eg["application"].unique()):
        e = eg[eg["application"] == app]
        sets = [("all", e)]
        if conf in e.columns:
            sets.append((f"confidence>={conf_min}", e[e[conf].astype(float) >= conf_min]))
        for lab, es in sets:
            ids = set(es["Material Id"].astype(str))
            sub = g[g["gnome_id"].isin(ids)].assign(group=f"EG_{app}")
            if not len(sub):
                rows.append(dict(application=app, subset=lab, n_listed=len(ids), n_matched=0))
                continue
            joined.append(sub.assign(subset=lab))
            for m in metrics:
                both = pd.concat([sub, master[master["group"] == "MP_known"]])
                r = ca.compare_metric(both.assign(**{m: both[m].astype(float)}), m,
                                      f"EG_{app}", "MP_known", "s_n", n_boot=n_boot, seed=seed)
                rows.append(dict(application=app, subset=lab, n_listed=len(ids),
                                 n_matched=len(sub), metric=m, shortlist=r["mean_target"],
                                 all_gnome=float(g[m].astype(float).mean()),
                                 known_matched=r["mean_ref_matched"], diff_vs_known=r["diff"],
                                 ci_low=r["ci_low"], ci_high=r["ci_high"]))
    return pd.DataFrame(rows), (pd.concat(joined) if joined else pd.DataFrame())


# --------------------------------------------------------------------------------------
# 6. Level 4: re-check a sample with a universal ML potential
# --------------------------------------------------------------------------------------
def pick_level4_sample(master, n_per_group=40, seed=0, max_elements=3):
    """Ternary-or-smaller GNoME compounds; half with f-in-core lanthanides, half with
    no rare earth, same chemistry-type mix. Systems where GGA and GGA+U would be mixed
    (U-element together with O or F) are excluded to keep energies comparable."""
    rng = np.random.default_rng(seed)
    g = master[(master["group"] == "GNoME") & (master["nelements"] <= max_elements)].copy()
    els = g["elements"].str.split(",").apply(set)
    mixed_u = els.apply(lambda e: bool(e & set(ca.U_ELEMENTS)) and bool(e & {"O", "F"}))
    g = g[~mixed_u & ~g["has_radioactive"]]
    a = g[g["has_f_in_core"]]
    b = g[~g["has_ree"]]
    picks = []
    for cls in ("anion", "intermetallic"):
        aa, bb = a[a["chem_class"] == cls], b[b["chem_class"] == cls]
        k = min(n_per_group // 2, len(aa), len(bb))
        if k:
            picks.append(aa.sample(k, random_state=int(rng.integers(1e9))))
            picks.append(bb.sample(k, random_state=int(rng.integers(1e9))))
    return pd.concat(picks).reset_index(drop=True) if picks else g.head(0)


def subsystems(chemsys):
    els = chemsys.split("-")
    out = []
    for r in range(1, len(els) + 1):
        for c in __import__("itertools").combinations(sorted(els), r):
            out.append("-".join(c))
    return out


def gnome_structure_reader(url=GNOME_BYID_ZIP_URL, local_zip=None):
    """Returns a function gnome_id -> pymatgen Structure, reading single CIFs from the
    remote zip with HTTP range requests (no full download). Falls back to a local zip."""
    import fsspec
    from pymatgen.core import Structure
    fobj = open(local_zip, "rb") if local_zip else fsspec.open(url, "rb", block_size=2**20).open()
    z = zipfile.ZipFile(fobj)
    index = {}
    for n in z.namelist():
        base = os.path.basename(n)
        stem = base.rsplit(".", 1)[0]
        if base.lower().endswith(".cif"):
            index[stem] = n
    if not index:
        raise RuntimeError(f"No CIF files in zip. First names: {z.namelist()[:10]}")

    def read(gid):
        name = index.get(str(gid))
        if name is None:
            return None
        return Structure.from_str(z.read(name).decode(), fmt="cif")
    read.n_files = len(index)
    return read


def mp_structure_fetcher(api_key, e_hull_max=0.05):
    """Returns a function chemsys -> list of (material_id, Structure) for all MP
    entries in the chemical system and its subsystems near the hull."""
    from mp_api.client import MPRester
    from pymatgen.core import Structure

    def fetch(chemsys):
        with MPRester(api_key, use_document_model=False, mute_progress_bars=True) as mpr:
            docs = mpr.materials.summary.search(
                chemsys=subsystems(chemsys), energy_above_hull=(0, e_hull_max),
                deprecated=False, fields=["material_id", "structure"])
        out = []
        for d in docs:
            s = d["structure"]
            out.append((str(d["material_id"]), s if hasattr(s, "lattice") else Structure.from_dict(s)))
        return out
    return fetch


def make_calculator(name="mace", device=None):
    if name == "emt":                       # only for testing the code
        from ase.calculators.emt import EMT
        return EMT()
    if device is None:
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        except Exception:
            device = "cpu"
    if name == "mace":
        from mace.calculators import mace_mp
        return mace_mp(model="medium", device=device, default_dtype="float32")
    if name == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator
        return CHGNetCalculator(use_device=device)
    raise ValueError(name)


def relax_energy(structure, calc, fmax=0.05, steps=300):
    from pymatgen.io.ase import AseAtomsAdaptor
    try:
        from ase.filters import FrechetCellFilter as CellFilter
    except Exception:
        from ase.constraints import ExpCellFilter as CellFilter
    from ase.optimize import FIRE
    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.calc = calc
    opt = FIRE(CellFilter(atoms), logfile=None)
    opt.run(fmax=fmax, steps=steps)
    return float(atoms.get_potential_energy()), len(atoms), int(opt.nsteps)


def run_level4(sample, gnome_reader, mp_fetcher, calc, cache_path, max_mp_per_system=60,
               log_every=10):
    """Relaxes each sampled GNoME compound and the MP phases of its chemical system
    with the same ML potential, then measures how far the GNoME compound lies below
    (negative) or above (positive) the ML-potential hull of the MP phases.
    Energies are checkpointed to cache_path so an interrupted run can resume."""
    from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
    done = pd.read_parquet(cache_path) if os.path.exists(cache_path) else \
        pd.DataFrame(columns=["key", "formula", "energy", "natoms", "steps", "ok"])
    E = {r.key: r for r in done.itertuples()}
    new_rows = []

    def energy_of(key, struct):
        if key in E:
            return E[key]
        try:
            e, n, st = relax_energy(struct, calc)
            row = dict(key=key, formula=struct.composition.formula, energy=e, natoms=n,
                       steps=st, ok=bool(np.isfinite(e)))
        except Exception as ex:
            row = dict(key=key, formula=struct.composition.formula, energy=np.nan,
                       natoms=len(struct), steps=-1, ok=False)
        new_rows.append(row)
        E[key] = pd.Series(row)
        return E[key]

    def checkpoint():
        if new_rows:
            cur = pd.DataFrame([dict(key=k, formula=v.formula, energy=v.energy, natoms=v.natoms,
                                     steps=v.steps, ok=v.ok) for k, v in E.items()])
            os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
            cur.to_parquet(cache_path, index=False)
            new_rows.clear()

    results = []
    for i, r in enumerate(sample.itertuples()):
        s = gnome_reader(r.gnome_id)
        if s is None:
            results.append(dict(gnome_id=r.gnome_id, status="no_structure"))
            continue
        g = energy_of(f"gnome:{r.gnome_id}", s)
        mp = mp_fetcher(r.chemsys)[:max_mp_per_system]
        entries = []
        for mid, ms in mp:
            m = energy_of(f"mp:{mid}", ms)
            if m.ok:
                entries.append(PDEntry(ms.composition, float(m.energy)))
        elems = set(r.elements.split(","))
        have = {e.symbol for en in entries for e in en.composition.elements
                if len(en.composition.elements) == 1}
        if not g.ok or have != elems:
            results.append(dict(gnome_id=r.gnome_id, status="incomplete_references"))
        else:
            pdg = PhaseDiagram(entries)
            cand = PDEntry(s.composition, float(g.energy))
            _, e_ml = pdg.get_decomp_and_e_above_hull(cand, allow_negative=True)
            results.append(dict(gnome_id=r.gnome_id, status="ok", e_hull_ml=float(e_ml),
                                n_mp_phases=len(entries)))
        if (i + 1) % log_every == 0:
            checkpoint()
            ok = sum(1 for x in results if x.get("status") == "ok")
            print(f"  {i + 1}/{len(sample)} compounds done ({ok} ok)")
    checkpoint()
    res = pd.DataFrame(results).merge(sample, on="gnome_id", how="left")
    return res


def level4_summary(res, thresholds=(0.0, 0.025, 0.05)):
    d = res[res["status"] == "ok"].copy()
    rows = []
    for thr in thresholds:
        d["stable_before"] = True
        d["lost"] = d["e_hull_ml"] > thr
        for val in (True, False):
            s = d[d["has_f_in_core"] == val]
            p, lo, hi = ca.wilson(int(s["lost"].sum()), len(s))
            rows.append(dict(threshold_eV=thr, has_f_in_core=val, n=len(s),
                             not_stable_share=p, ci_low=lo, ci_high=hi,
                             median_e_hull_ml=s["e_hull_ml"].median()))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# 7. Demo data in the REAL file layouts (to test the whole notebook without a key)
# --------------------------------------------------------------------------------------
def make_demo_sources(seed=0, n_mp=6000, n_gnome=8000):
    """Fake inputs with the same columns as the real files, with a planted effect:
    GNoME favours f-in-core lanthanides and metal-only compounds, and its
    f-in-core compounds have thinner stability margins. ONLY for testing code."""
    rng = np.random.default_rng(seed)
    hhi = ca.element_hhi_table()
    pool = [e for e in hhi if e not in ca.RADIOACTIVE and e not in NOBLE]
    w = np.array([8 if e in ANIONS else 1.0 for e in pool])
    wg = np.array([2 if e in ANIONS else (3.0 if e in ca.F_IN_CORE else 1.0) for e in pool])
    xtal = ["cubic", "tetragonal", "orthorhombic", "hexagonal", "trigonal", "monoclinic", "triclinic"]

    def formulas(n, weights, nel_p):
        out = []
        p = weights / weights.sum()
        for _ in range(n):
            k = rng.choice([1, 2, 3, 4, 5], p=nel_p)
            els = rng.choice(pool, size=k, replace=False, p=p)
            out.append("".join(f"{e}{int(rng.integers(1, 5))}" for e in els))
        return out
    f_mp = formulas(n_mp, w, [0.05, 0.35, 0.4, 0.15, 0.05])
    theo = rng.random(n_mp) < 0.65
    mp = pd.DataFrame({"material_id": [f"mp-{i}" for i in range(n_mp)], "deprecated": False,
                       "theoretical": theo, "nsites": rng.integers(2, 60, n_mp),
                       "formula_pretty": f_mp, "energy_above_hull": np.abs(rng.normal(0, .05, n_mp)),
                       "formation_energy_per_atom": rng.normal(-1.5, .6, n_mp),
                       "crystal_system": rng.choice(xtal, n_mp)})
    f_g = formulas(n_gnome, wg, [0.0, 0.1, 0.4, 0.38, 0.12])
    fc = np.array([any(e in f for e in ca.F_IN_CORE) for f in f_g])  # crude flag
    margin = np.abs(rng.normal(0.06, 0.04, n_gnome)) * np.where(fc, 0.5, 1.0)
    graw = pd.DataFrame({
        "Composition": f_g, "MaterialId": [f"{i:010x}" for i in range(n_gnome)],
        "Reduced Formula": f_g, "Elements": "", "NSites": rng.integers(4, 80, n_gnome),
        "Crystal System": rng.choice(xtal, n_gnome), "Is Train": rng.random(n_gnome) < 0.3,
        "Decomposition Energy Per Atom": -margin, "Decomposition Energy Per Atom All": -margin,
        "Decomposition Energy Per Atom Relative": -margin,
        "Decomposition Energy Per Atom MP": -margin - 0.02,
        "Decomposition Energy Per Atom MP OQMD": -margin - 0.01})
    r2 = pd.DataFrame({"MaterialId": graw["MaterialId"].sample(3000, random_state=seed).values})
    fc_r2 = graw.set_index("MaterialId").loc[r2["MaterialId"], "Reduced Formula"].apply(
        lambda f: any(e in f for e in ca.F_IN_CORE)).values
    r2["Decomposition Energy Per Atom"] = rng.normal(0.0, 0.03, len(r2)) + 0.03 * fc_r2
    ext_f = formulas(4000, w * np.array([1.5 if e in ca.REE else 1 for e in pool]),
                     [0.05, 0.3, 0.45, 0.15, 0.05])
    xraw = pd.DataFrame({"MaterialId": [(f"mp-{i}" if i % 2 else f"{i:012x}") for i in range(4000)],
                         "Reduced Formula": ext_f, "NSites": rng.integers(2, 40, 4000)})
    ids = graw["MaterialId"].values
    eg = pd.DataFrame({"Material Id": rng.choice(ids, 1500, replace=False),
                       "AI-experts confidence (-)": rng.random(1500),
                       "application": rng.choice(["cathode", "perovskite", "thermoelectric"], 1500),
                       "subset": "demo"})
    thermo_rows = []
    for t, shift in ((ca.THERMO_GGA, 0.0), (ca.THERMO_R2, 0.05)):
        thermo_rows.append(pd.DataFrame({"material_id": mp["material_id"], "thermo_type": t,
                                         "e_hull": np.clip(mp["energy_above_hull"] + rng.normal(0, shift, n_mp), 0, None),
                                         "is_stable": False,
                                         "ef": mp["formation_energy_per_atom"] + rng.normal(0, shift, n_mp)}))
    thermo = {ca.THERMO_GGA: thermo_rows[0], ca.THERMO_R2: thermo_rows[1]}
    return mp, graw, r2, xraw, eg, thermo