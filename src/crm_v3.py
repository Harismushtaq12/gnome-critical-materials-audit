"""
crm_v3.py — two extensions on top of crm_v2.py.

A. Country exposure. Uses USGS Mineral Commodity Summaries 2025 world production data
   (production by country, latest year) to measure, per compound, how much of its
   supply depends on single countries. Gives an up-to-date HHI as well.

B. Lower-risk alternatives. For GNoME candidates shortlisted for real applications
   (Energy-GNoME), asks: for a high-risk candidate, is there another candidate with
   similar predicted performance but lower supply risk? And how much performance is
   lost if we only keep low-risk candidates?
"""
from __future__ import annotations

import io
import os
import re
import json
import zipfile
import urllib.request

import numpy as np
import pandas as pd

import crm_audit as ca
import crm_v2 as v2

# ======================================================================================
# A. Country exposure
# ======================================================================================
USGS_WORLD_ITEM = "6798fd34d34ea8c18376e8ee"   # MCS 2025 World Production, Capacity, Reserves
CHINA = "China"

# element -> ordered list of (commodity keyword, preferred TYPE keyword or None).
# The first commodity found in the file is used. "proxy" marks elements that borrow
# another commodity's country pattern (stated in the paper's methods).
ELEMENT_SOURCES = {
    "Al": [("BAUXITE", "mine"), ("ALUMINUM", None)], "Sb": [("ANTIMONY", "mine")],
    "As": [("ARSENIC", None)], "Ba": [("BARITE", "mine")], "Be": [("BERYLLIUM", "mine")],
    "Bi": [("BISMUTH", None)], "B": [("BORON", None)], "Br": [("BROMINE", None)],
    "Cd": [("CADMIUM", None)], "Cr": [("CHROMIUM", None)], "Co": [("COBALT", "mine")],
    "Cu": [("COPPER", "mine")], "F": [("FLUORSPAR", None)], "Ga": [("GALLIUM", None)],
    "Ge": [("GERMANIUM", None)], "Au": [("GOLD", "mine")], "C": [("GRAPHITE", None)],
    "In": [("INDIUM", None)], "I": [("IODINE", None)], "Fe": [("IRON ORE", None)],
    "Pb": [("LEAD", "mine")], "Li": [("LITHIUM", "mine")],
    "Mg": [("MAGNESIUM METAL", None), ("MAGNESIUM", None)], "Mn": [("MANGANESE", "mine")],
    "Hg": [("MERCURY", "mine")], "Mo": [("MOLYBDENUM", "mine")], "Ni": [("NICKEL", "mine")],
    "Nb": [("NIOBIUM", "mine")], "N": [("NITROGEN", None)], "P": [("PHOSPHATE", "mine")],
    "Pt": [("PLATINUM", "platinum")], "Pd": [("PLATINUM", "palladium")],
    "Ru": [("PLATINUM", "platinum")], "Rh": [("PLATINUM", "platinum")],      # proxy
    "Ir": [("PLATINUM", "platinum")], "Os": [("PLATINUM", "platinum")],      # proxy
    "K": [("POTASH", None)], "Re": [("RHENIUM", None)], "Se": [("SELENIUM", None)],
    "Si": [("SILICON", None)], "Ag": [("SILVER", "mine")], "Na": [("SALT", None)],
    "Cl": [("SALT", None)], "Sr": [("STRONTIUM", None)], "S": [("SULFUR", None)],
    "Ta": [("TANTALUM", "mine")], "Te": [("TELLURIUM", None)], "Sn": [("TIN", "mine")],
    "Ti": [("TITANIUM MINERAL", "ilmenite"), ("TITANIUM MINERAL", None), ("TITANIUM", None)],
    "W": [("TUNGSTEN", "mine")], "V": [("VANADIUM", "mine")], "Zn": [("ZINC", "mine")],
    "Zr": [("ZIRCONIUM", None)], "Hf": [("ZIRCONIUM", None)],                 # proxy
    "Y": [("YTTRIUM", None), ("RARE EARTH", "mine")],
}
for _e in ca.LANTHANIDES + ["Sc"]:
    ELEMENT_SOURCES.setdefault(_e, [("RARE EARTH", "mine")])               # proxy per element
PROXY_ELEMENTS = {"Ru", "Rh", "Ir", "Os", "Hf", "Sc"} | set(ca.LANTHANIDES)


def fetch_usgs_world(cache_dir, item_id=USGS_WORLD_ITEM, local_file=None):
    """Download the MCS 2025 world data table (CSV, possibly inside a zip) from
    ScienceBase, or read a file you uploaded yourself (local_file)."""
    out = os.path.join(cache_dir, "usgs_mcs2025_world_raw.parquet")
    if os.path.exists(out):
        return pd.read_parquet(out)
    if local_file:
        raw = _read_any_table(open(local_file, "rb").read(), local_file)
    else:
        meta = json.loads(urllib.request.urlopen(
            f"https://www.sciencebase.gov/catalog/item/{item_id}?format=json", timeout=60).read())
        files = list(meta.get("files", []))
        for fc in meta.get("facets", []) or []:
            files += fc.get("files", [])
        print("ScienceBase files:", [f.get("name") for f in files])
        pick = [f for f in files if re.search(r"world", f.get("name", ""), re.I)
                and f.get("name", "").lower().endswith((".csv", ".zip", ".xlsx"))]
        pick = pick or [f for f in files if f.get("name", "").lower().endswith((".csv", ".zip"))]
        if not pick:
            raise RuntimeError("No CSV/ZIP found; download the world data file by hand and set USGS_LOCAL_FILE.")
        f = pick[0]
        url = f.get("downloadUri") or f.get("url")
        print("downloading", f.get("name"), "from", url)
        raw = _read_any_table(urllib.request.urlopen(url, timeout=300).read(), f.get("name", ""))
    raw = raw.astype(str)
    raw.to_parquet(out, index=False)
    return raw


def _read_any_table(blob, name):
    name = name.lower()
    if name.endswith(".zip"):
        z = zipfile.ZipFile(io.BytesIO(blob))
        cands = [n for n in z.namelist() if n.lower().endswith(".csv")]
        world = [n for n in cands if "world" in n.lower()] or cands
        print("zip members used:", world[:3])
        return pd.concat([_read_csv_bytes(z.read(n)) for n in world], ignore_index=True)
    if name.endswith(".xlsx"):
        return pd.read_excel(io.BytesIO(blob))
    return _read_csv_bytes(blob)


def _read_csv_bytes(b):
    for enc in ("utf-8", "latin-1"):
        try:
            return pd.read_csv(io.BytesIO(b), encoding=enc, dtype=str)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(io.BytesIO(b), encoding="latin-1", dtype=str, on_bad_lines="skip")


def _num(x):
    s = re.sub(r"[,\s]", "", str(x))
    s = re.sub(r"[^\d.\-eE]", "", s) if re.search(r"\d", s) else ""
    try:
        return float(s)
    except ValueError:
        return np.nan


def parse_usgs_world(raw):
    """-> long table: commodity, type, country, value (latest production year)."""
    cols = {c: re.sub(r"[^a-z0-9]", "", c.lower()) for c in raw.columns}
    def find(*keys):
        for k in keys:
            for c, n in cols.items():
                if k in n:
                    return c
        return None
    c_comm, c_country, c_type = find("commodity"), find("country"), find("type")
    prod = [(int(m.group()), c) for c, n in cols.items()
            if ("prod" in n) and (m := re.search(r"20\d\d", n))]
    if not (c_comm and c_country and prod):
        raise KeyError(f"Unexpected USGS columns: {list(raw.columns)}")
    year, c_val = max(prod)
    long = pd.DataFrame({"commodity": raw[c_comm].astype(str).str.upper().str.strip(),
                         "type": raw[c_type].astype(str).str.strip() if c_type else "",
                         "country": raw[c_country].astype(str).str.strip(),
                         "value": raw[c_val].map(_num)})
    long.attrs["year_column"] = c_val
    return long


def element_country_shares(long):
    """Per element: country shares of world production, plus summary
    (new HHI 0-10000, top country, top share, China share)."""
    shares, summary = [], []
    is_world = lambda s: s.str.contains("world", case=False)
    is_other = lambda s: s.str.contains("other", case=False)
    for el, sources in ELEMENT_SOURCES.items():
        used = None
        for kw, tkw in sources:
            d = long[long["commodity"].str.contains(kw, regex=False)]
            if tkw and len(d) and d["type"].str.contains(tkw, case=False).any():
                d = d[d["type"].str.contains(tkw, case=False)]
            elif len(d) and d["type"].nunique() > 1:
                mine = d["type"].str.contains("mine", case=False)
                d = d[mine] if mine.any() else d[d["type"] == d["type"].iloc[0]]
            d = d.dropna(subset=["value"])
            if len(d[~is_world(d["country"])]):
                used = (kw, tkw, d)
                break
        if used is None:
            summary.append(dict(element=el, source=None))
            continue
        kw, tkw, d = used
        world = d.loc[is_world(d["country"]), "value"]
        cty = d[~is_world(d["country"])].groupby("country", as_index=False)["value"].sum()
        total = max(world.max() if len(world) else 0, cty["value"].sum())
        if total <= 0:
            summary.append(dict(element=el, source=kw))
            continue
        cty["share"] = cty["value"] / total
        named = cty[~is_other(cty["country"])]
        hhi = float((100 * named["share"]) .pow(2).sum())
        top = named.sort_values("share", ascending=False).iloc[0]
        china = named.loc[named["country"].str.contains(CHINA, case=False), "share"].sum()
        summary.append(dict(element=el, source=kw, type_hint=tkw, proxy=el in PROXY_ELEMENTS,
                            hhi_new=hhi, top_country=top["country"], top_share=float(top["share"]),
                            china_share=float(china), n_countries=len(named)))
        shares.append(cty.assign(element=el)[["element", "country", "share"]])
    return pd.DataFrame(summary), (pd.concat(shares, ignore_index=True) if shares else pd.DataFrame())


def add_country_metrics(df, summary, conc_threshold=0.5):
    """Per compound (atom-weighted over elements with data):
    hhi_new, china_frac, top_share_max (weakest link), conc_frac = atom share from
    elements where one country produces >= conc_threshold, coverage."""
    s = summary.dropna(subset=["hhi_new"]).set_index("element")
    H, CH, TOP = s["hhi_new"].to_dict(), s["china_share"].to_dict(), s["top_share"].to_dict()
    out = {k: [] for k in ("hhi_new", "china_frac", "top_share_max", "conc_frac", "country_coverage")}
    for comp in df["composition"].values:
        c = ca.parse_composition(comp)
        tot = sum(c.values()) or 1.0
        kn = {e: v for e, v in c.items() if e in H}
        ka = sum(kn.values())
        out["country_coverage"].append(ka / tot)
        if ka == 0:
            for k in ("hhi_new", "china_frac", "top_share_max", "conc_frac"):
                out[k].append(np.nan)
            continue
        out["hhi_new"].append(sum(v * H[e] for e, v in kn.items()) / ka)
        out["china_frac"].append(sum(v * CH[e] for e, v in kn.items()) / ka)
        out["top_share_max"].append(max(TOP[e] for e in kn))
        out["conc_frac"].append(sum(v for e, v in kn.items() if TOP[e] >= conc_threshold) / ka)
    df = df.copy()
    for k, v in out.items():
        df[k] = v
    return df


def country_profile(df, shares, groups, ref="MP_known", target="GNoME", strata_col="s_n_class", top=12):
    """Average dependence on each country (share of a group's supply, atom-weighted over
    elements with data). The reference group is re-weighted to the target's strata mix."""
    S = shares.pivot_table(index="element", columns="country", values="share", aggfunc="sum").fillna(0)
    els = list(S.index)
    idx = {e: i for i, e in enumerate(els)}
    t_mix = df.loc[df["group"] == target, strata_col].value_counts(normalize=True)
    res = {}
    for g in groups:
        d = df[df["group"] == g]
        if not len(d):
            continue
        if g != target:
            r = d[strata_col].value_counts(normalize=True)
            w = d[strata_col].map(lambda s: t_mix.get(s, 0) / r[s]).values
        else:
            w = np.ones(len(d))
        acc = np.zeros(len(els)); wsum = 0.0
        for comp, wi in zip(d["composition"].values, w):
            if wi == 0:
                continue
            c = ca.parse_composition(comp)
            kn = {e: v for e, v in c.items() if e in idx}
            ka = sum(kn.values())
            if ka == 0:
                continue
            for e, v in kn.items():
                acc[idx[e]] += wi * v / ka
            wsum += wi
        res[g] = (acc / wsum) @ S.values if wsum else np.zeros(S.shape[1])
    prof = pd.DataFrame(res, index=S.columns)
    prof = prof.loc[prof.max(axis=1).sort_values(ascending=False).index[:top]]
    return prof


# ======================================================================================
# B. Lower-risk alternatives among application shortlists
# ======================================================================================
# application -> (performance column, mode, parameter). mode "higher": bigger is better,
# compared with a relative tolerance; "target": closeness to a target value (eV),
# compared with an absolute tolerance.
APP_PERF = {
    "cathode": ("Gravimetric energy (Wh/kg)", "higher", None),
    "thermoelectric": ("Average zT (-)", "higher", None),
    "perovskite": ("Average Band Gap (eV)", "target", 1.34),   # single-junction solar optimum
}
CONF_COL = "AI-experts confidence (-)"


def eg_join(eg, master, risk_cols=("crm_eu_atom", "crm_us_atom", "has_ree", "hhi_p",
                                   "hhi_new", "china_frac", "top_share_max", "conc_frac")):
    g = master[master["group"] == "GNoME"]
    keep = ["gnome_id", "formula", "elements", "nelements"] + [c for c in risk_cols if c in g.columns]
    j = eg.merge(g[keep], left_on="Material Id", right_on="gnome_id", how="inner")
    return j


def _performance(d, app):
    col, mode, par = APP_PERF[app]
    if col not in d.columns:
        return None
    x = pd.to_numeric(d[col], errors="coerce")
    return x if mode == "higher" else -(x - par).abs()


def alternatives(egj, risk="crm_eu_atom", high_risk=0.5, rel_tol=0.10, abs_tol=0.10,
                 min_risk_drop=0.20, conf_min=0.5):
    """For each application subset:
      - share of high-risk candidates (risk >= high_risk) that have at least one
        alternative with performance within tolerance AND risk lower by >= min_risk_drop;
      - share of Pareto-optimal candidates (best performance-risk trade-off) that are high-risk;
      - cost of choosing low-risk only: best and top-10 mean performance, all vs low-risk.
    Returns (summary table, per-candidate table with Pareto flags)."""
    rows, per = [], []
    for (app, sub), d in egj.groupby(["application", "subset"]):
        if app not in APP_PERF:
            continue
        d = d.copy()
        d["perf"] = _performance(d, app)
        if d["perf"] is None or d["perf"].isna().all():
            continue
        if CONF_COL in d.columns:
            d = d[pd.to_numeric(d[CONF_COL], errors="coerce") >= conf_min]
        d = d.dropna(subset=["perf", risk]).drop_duplicates("gnome_id")
        if len(d) < 5:
            continue
        mode = APP_PERF[app][1]
        p, r = d["perf"].values.astype(float), d[risk].values.astype(float)
        order = np.argsort(-p)
        ps, rs = p[order], r[order]
        prefmin = np.minimum.accumulate(rs)                       # lowest risk among >= perf
        hr = r >= high_risk
        thr = p - (rel_tol * np.abs(p) if mode == "higher" else abs_tol)
        k = np.searchsorted(-ps, -thr, side="right") - 1         # last index with ps >= thr
        has_alt = (k >= 0) & (prefmin[np.clip(k, 0, None)] <= r - min_risk_drop)
        # Pareto front: no other candidate with >= perf and < risk
        front_sorted = rs <= np.minimum.accumulate(np.r_[np.inf, rs[:-1]]) - 1e-12
        front = np.zeros(len(d), bool); front[order] = front_sorted
        low = ~hr
        best_all, best_low = ps[0], (p[low].max() if low.any() else np.nan)
        top_all = np.sort(p)[::-1][:10].mean()
        top_low = np.sort(p[low])[::-1][:10].mean() if low.sum() >= 1 else np.nan
        rows.append(dict(application=app, subset=sub, n=len(d), high_risk_share=hr.mean(),
                         high_risk_with_alternative=has_alt[hr].mean() if hr.any() else np.nan,
                         pareto_n=int(front.sum()), pareto_high_risk_share=hr[front].mean(),
                         best_perf_all=best_all, best_perf_low_risk=best_low,
                         top10_mean_all=top_all, top10_mean_low_risk=top_low,
                         top10_loss=((top_all - top_low) / abs(top_all) if mode == "higher"
                                     else top_all - top_low),   # % for higher-better; eV for target
                         loss_unit="fraction" if mode == "higher" else "eV from 1.34",
                         perf_metric=APP_PERF[app][0] + (" (closeness to 1.34 eV)" if mode == "target" else "")))
        per.append(d.assign(pareto=front, high_risk=hr, has_alternative=has_alt))
    return pd.DataFrame(rows), (pd.concat(per, ignore_index=True) if per else pd.DataFrame())


# ======================================================================================
# Demo inputs in the assumed real layouts (testing only)
# ======================================================================================
def make_demo_usgs(seed=0):
    rng = np.random.default_rng(seed)
    countries = ["China", "Australia", "Chile", "Congo (Kinshasa)", "Russia", "South Africa",
                 "Brazil", "United States", "Indonesia", "Canada", "Other countries"]
    comms = sorted({kw for v in ELEMENT_SOURCES.values() for kw, _ in v})
    rows = []
    for cm in comms:
        tp = "Mine production" if cm not in ("GALLIUM", "GERMANIUM", "SILICON") else "Refinery production"
        types = [("Mine production, platinum", 1), ("Mine production, palladium", 1)] if cm == "PLATINUM" else [(tp, 1)]
        for t, _ in types:
            china_bias = 8 if cm in ("RARE EARTH", "GALLIUM", "GERMANIUM", "GRAPHITE", "TUNGSTEN") else 1
            vals = rng.gamma(1.0, 100, len(countries)); vals[0] *= china_bias
            for c, v in zip(countries, vals):
                rows.append({"Source": "MCS2025", "COMMODITY": cm.title() if cm != "RARE EARTH" else "Rare earths",
                             "COUNTRY": c, "TYPE": t, "UNIT_MEAS": "metric tons",
                             "PROD_2023": f"{v * 0.95:,.0f}", "PROD_EST_ 2024": f"{v:,.0f}"})
            rows.append({"Source": "MCS2025", "COMMODITY": cm.title(), "COUNTRY": "World total (rounded)", "TYPE": t,
                         "UNIT_MEAS": "metric tons", "PROD_2023": "", "PROD_EST_ 2024": f"{vals.sum():,.0f}"})
    return pd.DataFrame(rows)


def add_demo_properties(eg, master, seed=0):
    """Gives the demo shortlist realistic property columns (weakly linked to risk)."""
    rng = np.random.default_rng(seed)
    g = master[master.group == "GNoME"].set_index("gnome_id")["crm_eu_atom"]
    risk = eg["Material Id"].map(g).fillna(0.5).values
    eg = eg.copy()
    eg["Gravimetric energy (Wh/kg)"] = np.where(eg.application == "cathode", rng.gamma(4, 120, len(eg)) * (1 + 0.3 * risk), np.nan)
    eg["Average zT (-)"] = np.where(eg.application == "thermoelectric", rng.gamma(3, 0.12, len(eg)) * (1 + 0.3 * risk), np.nan)
    eg["Average Band Gap (eV)"] = np.where(eg.application == "perovskite", rng.uniform(0.5, 3, len(eg)), np.nan)
    eg["subset"] = np.where(eg.application == "cathode", rng.choice(["Li", "Na"], len(eg)), "demo")
    return eg