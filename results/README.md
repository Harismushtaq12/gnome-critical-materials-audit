# Results

Tables (CSV) and figures (PNG) from the paper run (real data, 3 October 2026).
"Known" means experimentally known materials in the Materials Project (`MP_known`).
"Matched" means the reference group is re-weighted to GNoME's mix of compound types.

## Groups

| Name | What it is |
|---|---|
| `GNoME` | DeepMind's GNoME stable materials (554,054) |
| `MP_known` | Materials Project entries marked as experimentally known (49,235) |
| `MP_theoretical` | Materials Project entries that are computed only (105,142) |
| `EXT_CODED` | Entries from other computed databases in DeepMind's external file, with coded IDs (45,063) |
| `EXT_MP` | Entries in DeepMind's external file with Materials Project IDs (9,421) |

## Measures

| Column | Meaning |
|---|---|
| `crm_eu_atom`, `crm_us_atom` | Share of atoms that are on the EU 2023 / US 2022 critical list |
| `crm_eu_mass`, `crm_us_mass` | The same, by mass |
| `crm_eu_amb_atom`, `crm_us_amb_atom` | The same, with Al, F, C, Si, He also counted as critical |
| `has_ree` | Contains a rare earth (Sc, Y, La–Lu) |
| `has_f_in_core` | Contains Pr, Nd, Pm, Sm, Tb, Dy, Ho, Er, Tm, Yb or Lu, whose 4f electrons are frozen in standard Materials Project calculations |
| `has_pgm` | Contains a platinum-group metal (Ru, Rh, Pd, Os, Ir, Pt) |
| `hhi_p` | Supply concentration (Herfindahl–Hirschman index, 2011 production data) |
| `hhi_new`, `china_frac`, `conc_frac`, `top_share_max` | Supply concentration, China share, share of highly concentrated elements and largest single-country share, from USGS 2025 data |
| `s_n`, `s_n_class`, `s_n5_class` | Matching groups: number of elements; plus chemistry type (anion compound or metal-only); with 5+ elements grouped |

## Files

| File | Section | What it shows |
|---|---|---|
| `HEADLINE.csv` | 22 | The main numbers quoted in the paper |
| `L1_testA_grid.csv` | 6 | GNoME vs each reference group, every measure and matching level |
| `L1_testA_capped.csv` | 6b | The same, with 5+-element compounds grouped |
| `L1_testA_common_support.csv` | 6b | **Main comparison**: only where at least 100 reference compounds exist (72.6% of GNoME) |
| `L1_element_enrichment.csv` | 7 | Which elements GNoME uses more or less than known materials |
| `L1_uniform_expectation.csv` | 8 | Observed use vs random element choice |
| `L2_gap_decomposition.csv`, `L2_gap_decomposition_common_support.csv` | 9 | How much of each gap comes from the type of chemistry |
| `L2_amplification.csv`, `L2_is_train.csv` | 10 | Is the bias inherited from training data? |
| `L2_templates.csv` | 11 | Near-copy compounds (rare earth swapped for another) |
| `L2_stability_margins.csv` | 12 | How close to the stability limit compounds are |
| `L2_r2scan_gnome.csv`, `L2_r2scan_mp.csv` | 13 | Do compounds stay stable with the r²SCAN method? |
| `L2_testB.csv` | 14 | Link between over-use and energy sensitivity (H2) |
| `L3_energy_gnome.csv` | 15 | Critical share in application shortlists |
| `country_element_table.csv`, `L3_country_grid.csv`, `L3_country_profile.csv`, `L3_country_common_support.csv` | 4b, 16 | Country exposure (USGS 2025) |
| `L3_alternatives.csv` | 17 | First design of the alternatives analysis (kept for the record, not used) |
| `L3_alternatives_v2.csv` | 17 | **Final design**: lower-risk alternatives with similar predicted performance |
| `X_list_check.csv`, `X_big_region.csv`, `X_r2scan_subset.csv`, `X_margins_mp_oqmd.csv` | 21 | Extra checks for the limits |
| `null_calibration.csv` | 20 | False-alarm rate of each test on fake data with no effect |
| `run_log.csv` | 22 | Settings of each run |
| `fig1_headline.png` … `fig8_alternatives.png` | — | Figures 1–8 |

Note: in the paper run, the extra checks (section 21) and the final alternatives design were run after the first
full run, in the same session. The notebook in this repository runs everything in one pass.
