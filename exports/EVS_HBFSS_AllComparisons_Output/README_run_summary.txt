SEQUENCE manuscript analysis run summary

Analysis name: EVS_HBFSS_AllComparisons_Output
Runtime: ~12.3 min
Output directory: /root/REAPER98632/exports/EVS_HBFSS_AllComparisons_Output

Core rules
  DESeq2 adjusted p value threshold: 0.2
  Strong effect boundary |LFC|: 1
  Primary EVS cutoff method: fourier_shared_cutoff
  Backup EVS cutoff method: fixed_top_n
  Crossing stability window: 3

Comparison cutoffs
  RT0_ZT6 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 4870 | quantile = 0.8457 | cutoff_value = 0.00120729
  RT2_ZT8 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 3372 | quantile = 0.8931 | cutoff_value = 0.00116006
  RT4_ZT10 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 4599 | quantile = 0.8547 | cutoff_value = 0.00098214
  RT8_ZT14 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 4334 | quantile = 0.8628 | cutoff_value = 0.00119603

Export families
  figures/
    comparison-level regime crossing panels
    volcano plots
    dispersion plots
    PCA panels
  tables/
    manuscript gene tables
    display class count tables
    dataset summary tables
    cutoff manifest and crossing summaries
    EVS loading rank tables
