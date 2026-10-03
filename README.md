# Does AI materials discovery lean on critical raw materials?

Code and results for an audit of Google DeepMind's **GNoME** dataset, about 554,000 computed stable materials.
The study asks two things:

1. Are GNoME's materials made of critical, supply-risky elements (for example rare earths and platinum-group metals)
   more often than known materials?
2. If so, why, and does it matter for real applications?

Everything uses public data and runs in Google Colab on a normal laptop. No lab work and no GPU are needed.

## Main findings (paper run, 3 October 2026)

Fair comparison means GNoME is compared with known materials that have the same number of elements and the same
chemistry type. It covers 72.6% of GNoME.

| Measure | GNoME | Known materials |
|---|---|---|
| Share of atoms on the EU 2023 critical list | 0.59 | 0.52 |
| Share of atoms on the US 2022 critical list | 0.61 | 0.44 |
| Contains a rare earth | 72% | 45% |
| Contains a lanthanide with frozen 4f electrons in standard DFT | 60% | 27% |
| Contains a platinum-group metal | 41% | 19% |
| Supply concentration (HHI, 2011 data) | 4,808 | 3,171 |

Other results (all in `results/`):

- GNoME uses rare earths more often than random element choice would (1.24×). Known materials use them less often (0.59×).
- Rare-earth-containing GNoME materials are closer to the stability limit (odds ratio 1.83).
- In application shortlists (battery cathodes, solar perovskites, thermoelectrics), critical-element shares are
  about 1.5 to 3 times those of known materials. For 93–100% of high-risk candidates, a lower-risk candidate with
  similar predicted performance exists. Lithium cathodes are the exception: there, high-risk candidates are predicted to perform better
  (median 442 vs 186 Wh/kg).

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

The notebook writes its own code files (`crm_audit.py`, `crm_v2.py`, `crm_v3.py`) at the start, so it runs on its own.
The same files are in `src/` for easy reading.

Results can change slightly if the online databases are updated. The paper run used the versions available on
3 October 2026.

## What is in this repository

| Path | What it is |
|---|---|
| `critical_element_audit.ipynb` | The full analysis, from download to final tables, in 22 sections |
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
