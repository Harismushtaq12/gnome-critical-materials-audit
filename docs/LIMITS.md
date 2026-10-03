# Limits of the study

What this study can and cannot show, and the checks done on each limit.

1. **No lab work.** The study judges which elements materials contain, not whether they can be made.
   Whether GNoME's rare-earth compounds are real needs 4f-explicit calculations or experiments.
2. **The fair comparison covers 72.6% of GNoME.** Large metal-only compounds (4 and 5+ elements) have almost no known
   counterparts. *Check:* compared with other computed databases (`EXT_CODED`), GNoME is still higher with 4 elements
   (rare earths 80% vs 71%), but about the same with 5+ elements (90% vs 93%; only 275 `EXT_CODED` compounds).
3. **Other databases are not separated.** DeepMind's external file does not say which entries come from which database.
4. **r²SCAN subset is not random.** DeepMind re-checked 117,043 GNoME materials. *Check:* this subset has fewer
   rare-earth compounds (67% vs 80%) and smaller compounds (3.8 vs 4.3 elements) than the rest. The r²SCAN result
   applies to this subset only. r²SCAN also uses frozen 4f electrons, so it cannot test the 4f question.
5. **Stability margins.** The main margins are measured against Materials Project phases. *Check:* against
   Materials Project plus OQMD, margins are smaller, but rare-earth compounds are still more fragile
   (odds ratio 1.71 instead of 1.83).
6. **Country data and lists.** USGS reports rare earths as one group and Ru, Rh, Os, Ir together with platinum, so
   these elements share one pattern. USGS gives no world germanium figure (germanium is in 3.4% of GNoME and 5.4% of
   known compounds). *Check:* counting Al, F, C, Si and He as critical makes the gap smaller but it stays clearly
   above zero (EU +0.043, US +0.132).
7. **Performance values are predictions.** The alternatives analysis uses Energy-GNoME's machine-learning predictions.
8. **Exploratory analyses.** Hypotheses H1–H3 were fixed before the real data (H1 supported; H2 and H3 not).
   Everything else is exploratory. The alternatives analysis was redesigned once after the first run, because the
   first design counted near-identical candidates. Both versions are saved.
9. **Statistical checks.** On fake data with no effect, the tests gave false alarms 4–7% of the time (100 runs),
   close to the expected 5%.
