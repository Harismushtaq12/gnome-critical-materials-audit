# Limits of the study

What this study can and cannot show, and the checks done on each limit. Checks marked R1–R8 were added in revision
(Level 5 of the notebook, run on 10 October 2026 with the same saved data as the paper run).

1. **No lab work.** The study judges which elements materials contain, not whether they can be made.
   Whether GNoME's rare-earth compounds are real needs 4f-explicit calculations or experiments.
2. **Known materials are not a neutral reference.** They reflect what chemists have chosen to make, which may be
   shaped by cost, availability and research interest. Known materials contain rare earths less often than a random
   choice of elements would give (0.59), and GNoME more often (1.24), so the gap has two sides.
3. **The fair comparison covers 72.6% of GNoME**, and 58% of its weight is the 4-element metal-only type, which rests
   on 276 known compounds. *Check (R1):* without that type the gaps are similar (rare earths +0.301, frozen-4f
   +0.315, Pt-group +0.145, EU share +0.100; 30.6% of GNoME). Cut-offs of 300, 500 and 1000 known compounds per type
   give the same result. 5+-element metal-only compounds can only be compared with other computed databases
   (rare earths 90% vs 93%; only 275 `EXT_CODED` compounds).
4. **Near-copies.** Many compounds differ only by swapping one rare earth for another, so the simple intervals are too
   narrow. *Check (R2):* with near-copies resampled together, the rare-earth gap is 0.268 (95% CI 0.215 to 0.328;
   0.207 to 0.342 with broader groups). Many GNoME entries may also be ordered versions of disordered phases; the
   study does not correct for this.
5. **Other databases are not separated.** DeepMind's external file does not say which entries come from which database.
6. **r²SCAN subset is not random.** DeepMind re-checked 117,043 GNoME materials. *Check:* this subset has fewer
   rare-earth compounds (66% vs 80%) and smaller compounds (3.8 vs 4.3 elements) than the rest. The r²SCAN result
   applies to this subset only. r²SCAN also uses frozen 4f electrons, so it cannot test the 4f question.
7. **Stability margins.** The main margins are measured against Materials Project phases. *Checks:* against
   Materials Project plus OQMD, margins are smaller, but frozen-4f compounds are still more often near the hull
   (odds ratio 1.71 instead of 1.83). With the number of MP entries per chemical system in the strata (R8), the
   odds ratio is 1.91.
8. **Matching method changes the size of some results.** *Checks:* with compound-type matching, the largest element
   ratios fall sharply (R6; Ir and Rh fall to about 2, osmium rises), shortlist gaps are much smaller (R3; no clear rare-earth
   excess for thermoelectrics), and with the other weight order the type mix explains 61% of the Pt-group gap (R7).
9. **Country data and lists.** USGS reports rare earths as one group and Ru, Rh, Os, Ir together with platinum, so
   these elements share one pattern. Data are for mining, not refining. USGS gives no world germanium figure
   (germanium is in 3.4% of GNoME and 5.4% of known compounds). *Check:* counting Al, F, C, Si and He as critical
   makes the gap smaller but it stays clearly above zero (EU +0.043, US +0.132).
10. **Performance values are predictions.** The alternatives analysis uses Energy-GNoME's machine-learning
    predictions. *Check (R4):* with stricter rules (model confidence of at least 0.9 and an alternative of the same
    compound type), 87–94% of high-risk thermoelectric candidates still have an alternative, but only 0–46% of
    cathode candidates.
11. **Exploratory analyses.** Three hypotheses (H1–H3) guided the analysis (H1 supported; H2 and H3 not).
    Everything else is exploratory. The alternatives analysis was redesigned once after the first run: the first
    design accepted any candidate at least as good as the high-risk one and kept extreme predictions, and its results
    (alternatives for 96–100% of high-risk candidates, no loss in 10 of 17 shortlists) were judged too optimistic.
    Both versions are saved.
12. **Statistical checks.** On fake data with no effect, the tests gave false alarms 4–7% of the time (100 runs),
    close to the expected 5%.
13. **One dataset.** All results refer to the August 2024 release of the GNoME data, and the audit has so far been
    applied only to this dataset.
