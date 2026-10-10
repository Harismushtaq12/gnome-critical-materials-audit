# Rare-earth and platinum-group metals in AI-predicted crystals: an audit of the GNoME dataset

Code and results for an audit of Google DeepMind's **GNoME** dataset of 554,054 AI-predicted stable crystals.
The study asks two things:

1. Do GNoME's materials contain critical, supply-risky elements (for example rare earths and platinum-group metals)
   more often than known materials?
2. If so, where does the excess come from, and does it reach real application shortlists?

Everything uses public data and runs in Google Colab on a normal laptop. No lab work and no GPU are needed.

## Main findings

Fair comparison means GNoME is compared with known materials that have the same number of elements and the same
chemistry type, using only compound types with at least 100 known examples. It covers 72.6% of GNoME.

| Measure | GNoME | Known materials |
|---|---|---|
| Share of atoms on the EU 2023 critical list | 0.59 | 0.52 |
| Share of atoms on the US 2022 critical list | 0.61 | 0.44 |
| Contains a rare earth | 72% | 45% |
| Contains a lanthanide with frozen 4f electrons in standard DFT | 60% | 27% |
| Contains a platinum-group metal | 41% | 19% |
| Supply concentration (HHI, USGS 2025, 0 to 10,000) | 4,174 | 3,723 |

Other results (all in `results/`):

- The rare-earth excess remains against other computed databases; the platinum-group excess is shared with them.
- Known materials use rare earths less often than a random choice of elements would (0.59×), and GNoME more often
  (1.24×), so the gap reflects both experimental practice and an excess in the dataset.
- Materials with frozen-4f lanthanides are closer to the stability limit (odds ratio 1.83; 1.62 for all rare earths).
- Application shortlists (battery cathodes, solar perovskites, thermoelectrics) also contain more critical elements,
  but with compound-type matching the gaps are much smaller. In most shortlists, candidates with fewer critical
  elements reach similar predicted performance; lithium cathodes are the exception (median 442 vs 186 Wh/kg).

### Checks added in revision (10 October 2026)

Level 5 of the notebook (files `R1` to `R8` in `results/`) tests the main result in eight ways. It holds without the
thin 4-element metal-only type (rare-earth gap +0.301), with near-copies resampled together (95% CI 0.215 to 0.328),
and after removing duplicates or radioactive elements. Some results are smaller with stricter matching: the
largest element ratios fall sharply, shortlist gaps shrink, and for cathodes close alternatives of the same compound type are
often missing.

The limits of the study are listed in the paper and in [`docs/LIMITS.md`](docs/LIMITS.md).

## Interactive results page

Explore the results in your browser (periodic table, applications, countries):
**https://harismushtaq12.github.io/gnome-critical-materials-audit/**

The page is the single file `docs/index.html`. It needs no installation.

## How to reproduce

1. Open `critical_element_audit.ipynb` in Google Colab (File → Upload notebook, or the GitHub tab).
2. **Test run:** set `DEMO_MODE = True` in section 2 and choose Runtime → Run all.
   It uses fake data and takes a few minutes. It only checks that the code works.
3. **Real run:**
   - Get a free API key from [materialsproject.org](https://materialsproject.org) (Dashboard → API key).
   - In Colab, open the key icon (Secrets), add `MP_API_KEY`, and allow notebook access.
   - Set `DEMO_MODE = False` and choose Run all.
   - The first run downloads about 2 GB and takes a few hours. Downloads and the main table are saved to
     Google Drive (`crm_audit_cache`), so later runs are much faster.
4. All tables (CSV) and figures (PNG) are written to `results/` and downloaded as a zip at the end.
5. **Revision checks only:** keep `USE_CHECKPOINT = True`, run sections 1 to 5 (they reload the saved table from
   Google Drive), then run the Level 5 cells (R1 to R8). The last cell prints the values used in the paper and saves
   `results/R0_values_for_paper.csv`.

The notebook writes its own code files (`crm_audit.py`, `crm_v2.py`, `crm_v3.py`, and `crm_rev.py` for Level 5) as it runs, so it runs on its own.
The same files are in `src/` for easy reading.

Results can change slightly if the online databases are updated. The paper run used the versions available on
3 October 2026; the revision run of 10 October 2026 reused the same saved data, and all its main result files are
identical to those of the paper run.

## What is in this repository

| Path | What it is |
|---|---|
| `critical_element_audit.ipynb` | The full analysis, from download to final tables, including the Level 5 revision checks |
| `src/` | The analysis code as plain Python files (the same code the notebook writes) |
| `results/` | All result tables and figures from the paper run. See [`results/README.md`](results/README.md) |
| `docs/LIMITS.md` | Known limits of the study and the checks done on them |
| `docs/index.html` | The interactive results page (served by GitHub Pages) |
| `requirements.txt` | Python packages used |

## Data sources

The raw data is **not** included. The notebook downloads it from the original sources.

| Data | Source | Licence |
|---|---|---|
| GNoME stable materials, r²SCAN re-check, external hull entries | Merchant et al., *Nature* 624, 80–85 (2023), [doi:10.1038/s41586-023-06735-9](https://doi.org/10.1038/s41586-023-06735-9); [github.com/google-deepmind/materials_discovery](https://github.com/google-deepmind/materials_discovery) | CC BY-NC 4.0 |
| Materials Project | Jain et al., *APL Materials* 1, 011002 (2013), [doi:10.1063/1.4812323](https://doi.org/10.1063/1.4812323) | CC BY 4.0 |
| Energy-GNoME application shortlists | De Angelis et al., *Energy and AI* 22, 100605 (2025), [doi:10.1016/j.egyai.2025.100605](https://doi.org/10.1016/j.egyai.2025.100605) | CC BY 4.0 |
| World production by country | USGS Mineral Commodity Summaries 2025 Data Release, [doi:10.5066/P13XCP3R](https://doi.org/10.5066/P13XCP3R) | Public domain (US Government) |
| Critical raw material lists | European Commission (2023), *Study on the Critical Raw Materials for the EU 2023*; US Final List of Critical Minerals 2022 (87 FR 10381) | Public |
| Supply concentration (HHI, 2011) | Gaultois et al., *Chem. Mater.* 25, 2911–2920 (2013), as packaged in pymatgen | — |

## How to cite

If you use this code or these results, please cite the archived version (DOI on Zenodo; see the "Cite this
repository" button on GitHub) and the paper, once published.

## Licence

- Code: MIT licence (see `LICENSE`).
- Result tables and figures in `results/`: CC BY 4.0. They are summary statistics. They contain no rows of the
  original GNoME data, which remains under its own CC BY-NC 4.0 licence.

## Author

Mohammad Haris Mushtaq
