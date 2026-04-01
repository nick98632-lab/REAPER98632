SEQUENCE manuscript analysis run summary

Analysis name: EVS_HBFSS_AllComparisons_Output
Runtime: ~11.9 min
Output directory: /root/REAPER98632/exports/EVS_HBFSS_AllComparisons_Output

Core rules
  DESeq2 adjusted p value threshold: 0.2
  Strong effect boundary |LFC|: 1
  Primary EVS cutoff method: fourier_shared_cutoff
  Backup EVS cutoff method: fixed_top_n
  Crossing stability window: 3

Comparison cutoffs
  RT0_ZT6 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 2130 | quantile = 0.9325 | cutoff_value = 8.98183
  RT2_ZT8 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 1840 | quantile = 0.9417 | cutoff_value = 9.2649
  RT4_ZT10 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 1840 | quantile = 0.9419 | cutoff_value = 9.54883
  RT8_ZT14 | method = first_stable_crossing | selected_reason = first_stable_crossing | rank = 1548 | quantile = 0.951 | cutoff_value = 9.98193

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
