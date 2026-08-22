#!/usr/bin/env Rscript

# =============================================================================
# ZERO-ARM NB-EXCESS DIAGNOSTIC
# =============================================================================
#
# PURPOSE
# -------
# This is a SENSITIVITY CHECK. It does NOT silently change the primary method.
#
# The current primary analysis uses:
#
#   E_ig = max(V_pool,i - mu_ig, 0)
#
# where V_pool,i is the pooled within-group variance across all 8 arms.
#
# Therefore, if a feature has raw counts:
#
#   (0,0,0,0,0)
#
# in arm g, then mu_ig = 0, but V_pool,i may still be >0 because that same
# feature varies in other arms. The current formula can consequently assign
# positive NB-excess mass to an arm in which the feature was never observed.
#
# This diagnostic asks:
#
#   1. How often does that happen?
#   2. How much of each arm's NB-excess mass comes from all-zero features?
#   3. Does forcing E_ig = 0 ONLY when the raw count sum in arm g is zero
#      materially change the cumulative PC1-NB divergence geometry?
#   4. Does that sensitivity correction change shared c1/c2?
#   5. Does it change the final weighted-Pareto top-k cutoff?
#
# IMPORTANT
# ---------
# The PC1 rankings and raw-count variance geometry are NOT changed.
# Only the NB-excess mass assigned to all-zero arm/feature combinations is
# zero-masked in the sensitivity analysis.
#
# Run from the project directory:
#
#   cd /root/REAPER98632
#   Rscript ZERO_ARM_NB_DIAGNOSTIC.R
#
# This script sources the current EMPERICALCUTOFF.R first, so the diagnostic
# uses the exact same data loading, normalization, PC1 ranking, divergence
# model, and weighted-Pareto functions as the primary analysis.
#
# =============================================================================


# =============================================================================
# LOAD THE CURRENT PRIMARY ANALYSIS
# =============================================================================

PRIMARY_SCRIPT <- "EMPERICALCUTOFF.R"

if (!file.exists(PRIMARY_SCRIPT)) {
  stop(
    "Could not find ",
    PRIMARY_SCRIPT,
    " in the current working directory."
  )
}

message("Running current primary analysis first...")
source(PRIMARY_SCRIPT, local = FALSE)

required_objects <- c(
  "count_mat",
  "group_labels",
  "normalized_counts",
  "pooled",
  "group_results",
  "COMPARISONS",
  "C1",
  "C2",
  "GLOBAL_K",
  "global_scan",
  "global_opt",
  "pair_optima",
  "fit_shared_knots",
  "smooth_divergence_for_display",
  "scan_pair_cutoffs",
  "select_weighted_pareto_optimum",
  "aggregate_global_cutoff_scan",
  "rank_cutoff_from_k",
  "OUT_ROOT"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Primary script did not create required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}


# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

DIAG_ROOT <- file.path(
  OUT_ROOT,
  "ZeroArm_NB_Diagnostic"
)

dir.create(
  DIAG_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# HELPERS
# =============================================================================

safe_pct <- function(num, den) {
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }

  100 * num / den
}


get_group_raw_sum <- function(group_name) {
  idx <- which(
    as.character(group_labels) ==
      group_name
  )

  if (length(idx) < 1L) {
    stop(
      "No samples found for group ",
      group_name
    )
  }

  x <- count_mat[
    ,
    idx,
    drop = FALSE
  ]

  out <- rowSums(
    x,
    na.rm = TRUE
  )

  names(out) <- rownames(
    count_mat
  )

  out
}


get_group_nonzero_n <- function(group_name) {
  idx <- which(
    as.character(group_labels) ==
      group_name
  )

  x <- count_mat[
    ,
    idx,
    drop = FALSE
  ]

  out <- rowSums(
    x > 0,
    na.rm = TRUE
  )

  names(out) <- rownames(
    count_mat
  )

  out
}


region_from_rank <- function(
    rank,
    c1,
    c2) {

  ifelse(
    rank < c1,
    "Remainder",
    ifelse(
      rank <= c2,
      "Divergence",
      "Leading edge"
    )
  )
}


make_zero_masked_group_result <- function(
    group_name,
    group_result) {

  df <- group_result$data

  raw_sum_map <- get_group_raw_sum(
    group_name
  )

  raw_nonzero_n_map <- get_group_nonzero_n(
    group_name
  )

  raw_sum_ranked <- unname(
    raw_sum_map[
      df$feature_id
    ]
  )

  raw_nonzero_n_ranked <- unname(
    raw_nonzero_n_map[
      df$feature_id
    ]
  )

  if (
    any(!is.finite(raw_sum_ranked)) ||
    any(!is.finite(raw_nonzero_n_ranked))
  ) {
    stop(
      "Failed to align raw-count diagnostics for group ",
      group_name
    )
  }

  all_zero <- (
    raw_sum_ranked == 0
  )

  E_original <- df$nb_excess_variance

  # Sensitivity correction ONLY:
  # an all-zero arm/feature combination cannot contribute positive
  # arm-specific NB-excess mass.
  E_masked <- E_original
  E_masked[all_zero] <- 0

  E_total <- sum(
    E_masked,
    na.rm = TRUE
  )

  if (
    !is.finite(E_total) ||
    E_total <= 0
  ) {
    stop(
      "Zero-masked NB excess mass is undefined for group ",
      group_name
    )
  }

  q_masked <- E_masked / E_total

  F_E_masked <- cumsum(
    q_masked
  )

  D_masked <- (
    F_E_masked -
    df$cumulative_pc1_mass
  )

  display_D_masked <-
    smooth_divergence_for_display(
      rank = df$rank,
      D = D_masked
    )

  df$raw_group_sum <- raw_sum_ranked
  df$raw_nonzero_sample_n <- raw_nonzero_n_ranked
  df$all_zero_in_arm <- all_zero

  df$nb_excess_variance_original <-
    E_original

  df$nb_excess_variance_masked <-
    E_masked

  df$nb_excess_removed_by_zero_mask <-
    E_original - E_masked

  # Replace the fields consumed by fit_shared_knots().
  df$nb_excess_variance <-
    E_masked

  df$nb_excess_variance_mass <-
    q_masked

  df$cumulative_nb_mass <-
    F_E_masked

  df$cumulative_divergence <-
    D_masked

  df$display_D <-
    display_D_masked

  out <- group_result
  out$data <- df

  out
}


# =============================================================================
# ZERO-COUNT DIAGNOSTIC TABLES
# =============================================================================

message(
  "Quantifying all-zero arm/feature combinations..."
)

arm_summary_rows <- list()
feature_rows <- list()

for (g in names(group_results)) {
  current_df <- group_results[[g]]$data

  raw_sum_map <- get_group_raw_sum(
    g
  )

  raw_nonzero_n_map <- get_group_nonzero_n(
    g
  )

  raw_sum <- unname(
    raw_sum_map[
      current_df$feature_id
    ]
  )

  raw_nonzero_n <- unname(
    raw_nonzero_n_map[
      current_df$feature_id
    ]
  )

  all_zero <- (
    raw_sum == 0
  )

  E_current <-
    current_df$nb_excess_variance

  E_total <- sum(
    E_current,
    na.rm = TRUE
  )

  E_zero <- sum(
    E_current[
      all_zero
    ],
    na.rm = TRUE
  )

  positive_E_zero <- (
    all_zero &
    E_current > 0
  )

  arm_summary_rows[[g]] <-
    data.frame(
      group = g,

      n_features =
        nrow(current_df),

      n_samples =
        sum(
          as.character(group_labels) ==
            g
        ),

      all_zero_feature_n =
        sum(all_zero),

      all_zero_with_positive_pooled_variance_n =
        sum(
          all_zero &
          current_df$pooled_normalized_variance >
            0
        ),

      all_zero_with_positive_current_E_n =
        sum(
          positive_E_zero
        ),

      current_total_nb_excess =
        E_total,

      current_nb_excess_from_all_zero =
        E_zero,

      pct_current_nb_excess_from_all_zero =
        safe_pct(
          E_zero,
          E_total
        ),

      current_nb_mass_from_all_zero =
        sum(
          current_df$nb_excess_variance_mass[
            all_zero
          ],
          na.rm = TRUE
        ),

      max_rank_of_all_zero_feature =
        if (
          any(all_zero)
        ) {
          max(
            current_df$rank[
              all_zero
            ]
          )
        } else {
          NA_integer_
        },

      all_zero_features_beyond_current_c2 =
        sum(
          all_zero &
          current_df$rank >
            C2
        ),

      stringsAsFactors = FALSE
    )

  if (any(all_zero)) {
    feature_rows[[g]] <-
      data.frame(
        group = g,
        feature_id =
          current_df$feature_id[
            all_zero
          ],
        rank =
          current_df$rank[
            all_zero
          ],
        current_region =
          region_from_rank(
            current_df$rank[
              all_zero
            ],
            C1,
            C2
          ),
        raw_group_sum =
          raw_sum[
            all_zero
          ],
        raw_nonzero_sample_n =
          raw_nonzero_n[
            all_zero
          ],
        normalized_group_mean =
          current_df$normalized_group_mean[
            all_zero
          ],
        pooled_normalized_variance =
          current_df$pooled_normalized_variance[
            all_zero
          ],
        current_nb_excess_variance =
          E_current[
            all_zero
          ],
        current_nb_excess_mass =
          current_df$nb_excess_variance_mass[
            all_zero
          ],
        abs_pc1_loading =
          current_df$abs_pc1_loading[
            all_zero
          ],
        pc1_variance_contribution =
          current_df$pc1_variance_contribution[
            all_zero
          ],
        stringsAsFactors = FALSE
      )
  }
}

arm_summary <- dplyr::bind_rows(
  arm_summary_rows
)

zero_feature_table <- dplyr::bind_rows(
  feature_rows
) %>%
  dplyr::arrange(
    group,
    dplyr::desc(
      current_nb_excess_variance
    )
  )

write.csv(
  arm_summary,
  file.path(
    DIAG_ROOT,
    "Table_ZeroArm_ByArm.csv"
  ),
  row.names = FALSE
)

write.csv(
  zero_feature_table,
  file.path(
    DIAG_ROOT,
    "Table_ZeroArm_Features.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# ZERO PATTERNS WITHIN EACH CONTROL/TREATMENT COMPARISON
# =============================================================================

comparison_rows <- list()

for (
  comparison_name in
  names(COMPARISONS)
) {
  mapping <-
    COMPARISONS[[
      comparison_name
    ]]

  control_group <-
    unname(
      mapping[[
        "control"
      ]]
    )

  treatment_group <-
    unname(
      mapping[[
        "treatment"
      ]]
    )

  c_zero <- (
    get_group_raw_sum(
      control_group
    ) == 0
  )

  t_zero <- (
    get_group_raw_sum(
      treatment_group
    ) == 0
  )

  ids <- names(
    c_zero
  )

  t_zero <- t_zero[
    ids
  ]

  # Current global top-k union.
  c_df <-
    group_results[[
      control_group
    ]]$data

  t_df <-
    group_results[[
      treatment_group
    ]]$data

  c_top <- tail(
    c_df$feature_id,
    GLOBAL_K
  )

  t_top <- tail(
    t_df$feature_id,
    GLOBAL_K
  )

  selected_union <-
    union(
      c_top,
      t_top
    )

  comparison_rows[[
    comparison_name
  ]] <-
    data.frame(
      comparison =
        comparison_name,

      control_group =
        control_group,

      treatment_group =
        treatment_group,

      both_arms_all_zero_n =
        sum(
          c_zero &
          t_zero
        ),

      control_only_all_zero_n =
        sum(
          c_zero &
          !t_zero
        ),

      treatment_only_all_zero_n =
        sum(
          !c_zero &
          t_zero
        ),

      neither_all_zero_n =
        sum(
          !c_zero &
          !t_zero
        ),

      current_global_k =
        GLOBAL_K,

      selected_union_n =
        length(
          selected_union
        ),

      selected_union_both_arms_all_zero_n =
        sum(
          selected_union %in%
            ids[
              c_zero &
              t_zero
            ]
        ),

      selected_union_control_all_zero_n =
        sum(
          selected_union %in%
            ids[
              c_zero
            ]
        ),

      selected_union_treatment_all_zero_n =
        sum(
          selected_union %in%
            ids[
              t_zero
            ]
        ),

      stringsAsFactors = FALSE
    )
}

comparison_summary <-
  dplyr::bind_rows(
    comparison_rows
  )

write.csv(
  comparison_summary,
  file.path(
    DIAG_ROOT,
    "Table_ZeroArm_ByComparison.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# SENSITIVITY: FORCE E=0 ONLY FOR ALL-ZERO ARM/FEATURE COMBINATIONS
# =============================================================================

message(
  "Recomputing cumulative divergence after zero-masking NB excess..."
)

masked_results <- vector(
  "list",
  length(
    group_results
  )
)

names(
  masked_results
) <- names(
  group_results
)

for (
  g in names(
    group_results
  )
) {
  masked_results[[g]] <-
    make_zero_masked_group_result(
      group_name = g,
      group_result =
        group_results[[g]]
    )
}

message(
  "Refitting shared c1/c2 under zero-mask sensitivity analysis..."
)

masked_knot_fit <-
  fit_shared_knots(
    masked_results
  )

MASKED_C1 <-
  masked_knot_fit$c1

MASKED_C2 <-
  masked_knot_fit$c2

message(
  "Current c1/c2 = ",
  C1,
  " / ",
  C2
)

message(
  "Zero-masked c1/c2 = ",
  MASKED_C1,
  " / ",
  MASKED_C2
)


# =============================================================================
# SENSITIVITY: RE-RUN WEIGHTED PARETO CUTOFF WITH MASKED c1/c2
# =============================================================================

message(
  "Re-running weighted Pareto cutoff using zero-masked c1/c2..."
)

masked_pair_scans <- list()
masked_pair_optima <- list()

for (
  comparison_name in
  names(COMPARISONS)
) {
  mapping <-
    COMPARISONS[[
      comparison_name
    ]]

  control_group <-
    unname(
      mapping[[
        "control"
      ]]
    )

  treatment_group <-
    unname(
      mapping[[
        "treatment"
      ]]
    )

  raw_scan <-
    scan_pair_cutoffs(
      control_df =
        masked_results[[
          control_group
        ]]$data,

      treatment_df =
        masked_results[[
          treatment_group
        ]]$data,

      c1 =
        MASKED_C1,

      c2 =
        MASKED_C2,

      comparison_name =
        comparison_name,

      control_group =
        control_group,

      treatment_group =
        treatment_group
    )

  opt <-
    select_weighted_pareto_optimum(
      raw_scan,
      good_col =
        "good_n",
      cost_col =
        "remainder_cross_n"
    )

  masked_pair_scans[[
    comparison_name
  ]] <- opt$scan

  masked_pair_optima[[
    comparison_name
  ]] <- opt
}

masked_global_scan_raw <-
  aggregate_global_cutoff_scan(
    masked_pair_scans
  )

masked_global_opt <-
  select_weighted_pareto_optimum(
    masked_global_scan_raw,
    good_col =
      "good_n",
    cost_col =
      "remainder_cross_n"
  )

masked_global_scan <-
  masked_global_opt$scan

MASKED_GLOBAL_K <-
  masked_global_opt$selected_k


# =============================================================================
# SENSITIVITY TABLES
# =============================================================================

pair_sensitivity_rows <- list()

for (
  comparison_name in
  names(COMPARISONS)
) {
  current_k <-
    pair_optima[[
      comparison_name
    ]]$selected_k

  masked_k <-
    masked_pair_optima[[
      comparison_name
    ]]$selected_k

  current_selected <-
    pair_optima[[
      comparison_name
    ]]$scan %>%
    dplyr::filter(
      k ==
        current_k
    ) %>%
    dplyr::slice(1L)

  masked_selected <-
    masked_pair_optima[[
      comparison_name
    ]]$scan %>%
    dplyr::filter(
      k ==
        masked_k
    ) %>%
    dplyr::slice(1L)

  pair_sensitivity_rows[[
    comparison_name
  ]] <-
    data.frame(
      comparison =
        comparison_name,

      current_c1 =
        C1,

      zero_masked_c1 =
        MASKED_C1,

      delta_c1 =
        MASKED_C1 -
        C1,

      current_c2 =
        C2,

      zero_masked_c2 =
        MASKED_C2,

      delta_c2 =
        MASKED_C2 -
        C2,

      current_pair_k =
        current_k,

      zero_masked_pair_k =
        masked_k,

      delta_pair_k =
        masked_k -
        current_k,

      current_good_n =
        current_selected$good_n[
          1L
        ],

      zero_masked_good_n =
        masked_selected$good_n[
          1L
        ],

      current_remainder_cross_n =
        current_selected$remainder_cross_n[
          1L
        ],

      zero_masked_remainder_cross_n =
        masked_selected$remainder_cross_n[
          1L
        ],

      stringsAsFactors = FALSE
    )
}

pair_sensitivity <-
  dplyr::bind_rows(
    pair_sensitivity_rows
  )

write.csv(
  pair_sensitivity,
  file.path(
    DIAG_ROOT,
    "Table_ZeroArm_Cutoff_Sensitivity.csv"
  ),
  row.names = FALSE
)

current_global_selected <-
  global_scan %>%
  dplyr::filter(
    k ==
      GLOBAL_K
  ) %>%
  dplyr::slice(1L)

masked_global_selected <-
  masked_global_scan %>%
  dplyr::filter(
    k ==
      MASKED_GLOBAL_K
  ) %>%
  dplyr::slice(1L)

overall_sensitivity <-
  data.frame(
    metric = c(
      "c1",
      "c2",
      "leading_edge_size",
      "global_weighted_k",
      "global_cutoff_rank",
      "global_good_n",
      "global_remainder_cross_n"
    ),

    current = c(
      C1,
      C2,
      nrow(
        group_results[[1L]]$data
      ) - C2,
      GLOBAL_K,
      rank_cutoff_from_k(
        nrow(
          group_results[[1L]]$data
        ),
        GLOBAL_K
      ),
      current_global_selected$good_n[
        1L
      ],
      current_global_selected$remainder_cross_n[
        1L
      ]
    ),

    zero_masked = c(
      MASKED_C1,
      MASKED_C2,
      nrow(
        masked_results[[1L]]$data
      ) - MASKED_C2,
      MASKED_GLOBAL_K,
      rank_cutoff_from_k(
        nrow(
          masked_results[[1L]]$data
        ),
        MASKED_GLOBAL_K
      ),
      masked_global_selected$good_n[
        1L
      ],
      masked_global_selected$remainder_cross_n[
        1L
      ]
    ),

    stringsAsFactors = FALSE
  )

overall_sensitivity$delta <-
  overall_sensitivity$zero_masked -
  overall_sensitivity$current

write.csv(
  overall_sensitivity,
  file.path(
    DIAG_ROOT,
    "Table_ZeroArm_Overall_Sensitivity.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# CLEAN DIAGNOSTIC FIGURE
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

diag_theme <- function(
    base_size = 12) {

  theme_minimal(
    base_size =
      base_size
  ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          size =
            base_size + 1.4
        ),

      plot.subtitle =
        element_text(
          size =
            base_size - 0.4,
          color =
            "#444444"
        ),

      axis.title =
        element_text(
          face =
            "bold"
        ),

      panel.grid.minor =
        element_blank(),

      panel.grid.major.x =
        element_blank(),

      legend.position =
        "bottom",

      legend.title =
        element_blank(),

      plot.margin =
        margin(
          8,
          10,
          8,
          10
        )
    )
}


GROUP_COLORS <- c(
  RT0  = "#377EB8",
  ZT6  = "#00A087",
  RT2  = "#4DAF4A",
  ZT8  = "#984EA3",
  RT4  = "#FF7F00",
  ZT10 = "#E64B35",
  RT8  = "#A65628",
  ZT14 = "#F39B7F"
)

METHOD_COLORS <- c(
  "Current pooled E" =
    "#2C7BB6",

  "Zero-masked E" =
    "#D01C8B"
)


# -------------------------------------------------------------------------
# Panel A: how much NB-excess mass comes from all-zero features?
# -------------------------------------------------------------------------

pA <-
  ggplot(
    arm_summary,
    aes(
      x = group,
      y =
        pct_current_nb_excess_from_all_zero,
      fill = group
    )
  ) +
  geom_col(
    width = 0.72
  ) +
  geom_text(
    aes(
      label =
        paste0(
          all_zero_with_positive_current_E_n,
          " sites"
        )
    ),
    vjust = -0.35,
    size = 3.0
  ) +
  scale_fill_manual(
    values =
      GROUP_COLORS
  ) +
  labs(
    title =
      "A. All-zero features receiving NB-excess mass",

    subtitle =
      "Current pooled-variance formulation",

    x =
      NULL,

    y =
      "% of arm NB-excess mass"
  ) +
  diag_theme() +
  theme(
    legend.position =
      "none"
  )


# -------------------------------------------------------------------------
# Panel B: current vs zero-masked cumulative divergence.
# -------------------------------------------------------------------------

median_current <-
  data.frame(
    rank =
      group_results[[1L]]$data$rank,

    D =
      apply(
        do.call(
          cbind,
          lapply(
            group_results,
            function(z) {
              z$data$cumulative_divergence
            }
          )
        ),
        1L,
        median,
        na.rm = TRUE
      ),

    method =
      "Current pooled E",

    stringsAsFactors =
      FALSE
  )

median_masked <-
  data.frame(
    rank =
      masked_results[[1L]]$data$rank,

    D =
      apply(
        do.call(
          cbind,
          lapply(
            masked_results,
            function(z) {
              z$data$cumulative_divergence
            }
          )
        ),
        1L,
        median,
        na.rm = TRUE
      ),

    method =
      "Zero-masked E",

    stringsAsFactors =
      FALSE
  )

median_D <-
  bind_rows(
    median_current,
    median_masked
  )

boundary_df <-
  data.frame(
    rank = c(
      C1,
      C2,
      MASKED_C1,
      MASKED_C2
    ),

    method = c(
      "Current pooled E",
      "Current pooled E",
      "Zero-masked E",
      "Zero-masked E"
    ),

    boundary = c(
      "c1",
      "c2",
      "c1",
      "c2"
    ),

    stringsAsFactors =
      FALSE
  )

pB <-
  ggplot(
    median_D,
    aes(
      x = rank,
      y = D,
      color = method
    )
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dotted",
    color = "#777777"
  ) +
  geom_line(
    linewidth = 1.15
  ) +
  geom_vline(
    data = boundary_df,
    aes(
      xintercept = rank,
      color = method,
      linetype = boundary
    ),
    linewidth = 0.8,
    show.legend = TRUE
  ) +
  scale_color_manual(
    values =
      METHOD_COLORS
  ) +
  scale_linetype_manual(
    values = c(
      c1 = "dashed",
      c2 = "longdash"
    )
  ) +
  labs(
    title =
      "B. Divergence sensitivity to all-zero masking",

    subtitle =
      paste0(
        "Current c1/c2 = ",
        C1,
        "/",
        C2,
        "   |   Masked = ",
        MASKED_C1,
        "/",
        MASKED_C2
      ),

    x =
      "PC1 rank",

    y =
      "Median cumulative divergence"
  ) +
  diag_theme()


# -------------------------------------------------------------------------
# Panel C: weighted Pareto cutoff sensitivity.
# -------------------------------------------------------------------------

current_frontier <-
  global_scan %>%
  filter(
    is_pareto
  ) %>%
  arrange(
    remainder_cross_n,
    good_n
  ) %>%
  mutate(
    method =
      "Current pooled E"
  )

masked_frontier <-
  masked_global_scan %>%
  filter(
    is_pareto
  ) %>%
  arrange(
    remainder_cross_n,
    good_n
  ) %>%
  mutate(
    method =
      "Zero-masked E"
  )

frontiers <-
  bind_rows(
    current_frontier,
    masked_frontier
  )

selected_pts <-
  bind_rows(
    current_global_selected %>%
      mutate(
        method =
          "Current pooled E",
        selected_k =
          GLOBAL_K
      ),

    masked_global_selected %>%
      mutate(
        method =
          "Zero-masked E",
        selected_k =
          MASKED_GLOBAL_K
      )
  )

pC <-
  ggplot(
    frontiers,
    aes(
      x =
        remainder_cross_n,
      y =
        good_n,
      color =
        method
    )
  ) +
  geom_line(
    linewidth = 1.25
  ) +
  geom_point(
    size = 1.4,
    alpha = 0.55
  ) +
  geom_point(
    data =
      selected_pts,
    aes(
      x =
        remainder_cross_n,
      y =
        good_n,
      color =
        method
    ),
    shape = 21,
    fill = "white",
    stroke = 1.4,
    size = 4.7
  ) +
  geom_text(
    data =
      selected_pts,
    aes(
      x =
        remainder_cross_n,
      y =
        good_n,
      label =
        paste0(
          "k*=",
          selected_k
        ),
      color =
        method
    ),
    nudge_y = 140,
    size = 3.4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values =
      METHOD_COLORS
  ) +
  labs(
    title =
      "C. Final cutoff sensitivity",

    subtitle =
      "Weighted Pareto optimum using remainder crossings as contamination",

    x =
      "Opposite-arm remainder crossings",

    y =
      "Joint + permissible disjoint sites"
  ) +
  diag_theme()


# -------------------------------------------------------------------------
# Save three-panel diagnostic.
# -------------------------------------------------------------------------

FIG_PATH <-
  file.path(
    DIAG_ROOT,
    "Figure_ZeroArm_NB_Sensitivity.png"
  )

grDevices::png(
  filename =
    FIG_PATH,
  width = 15,
  height = 12,
  units = "in",
  res = 360,
  bg = "white"
)

grid.newpage()

pushViewport(
  viewport(
    layout =
      grid.layout(
        nrow = 3,
        ncol = 1
      )
  )
)

print(
  pA,
  vp =
    viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
)

print(
  pB,
  vp =
    viewport(
      layout.pos.row = 2,
      layout.pos.col = 1
    )
)

print(
  pC,
  vp =
    viewport(
      layout.pos.row = 3,
      layout.pos.col = 1
    )
)

dev.off()


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

message("============================================================")
message("ZERO-ARM NB-EXCESS DIAGNOSTIC COMPLETE")
message("")

for (
  i in seq_len(
    nrow(
      arm_summary
    )
  )
) {
  message(
    arm_summary$group[i],
    ": all-zero=",
    arm_summary$all_zero_feature_n[i],
    "; all-zero with current E>0=",
    arm_summary$all_zero_with_positive_current_E_n[i],
    "; % NB-excess mass from all-zero=",
    signif(
      arm_summary$pct_current_nb_excess_from_all_zero[i],
      4
    )
  )
}

message("")
message(
  "CURRENT shared c1/c2: ",
  C1,
  " / ",
  C2
)

message(
  "ZERO-MASKED shared c1/c2: ",
  MASKED_C1,
  " / ",
  MASKED_C2
)

message(
  "Delta c1/c2: ",
  MASKED_C1 - C1,
  " / ",
  MASKED_C2 - C2
)

message("")
message(
  "CURRENT global weighted k*: ",
  GLOBAL_K
)

message(
  "ZERO-MASKED global weighted k*: ",
  MASKED_GLOBAL_K
)

message(
  "Delta global k*: ",
  MASKED_GLOBAL_K -
    GLOBAL_K
)

message("")
message(
  "Diagnostic outputs: ",
  DIAG_ROOT
)

message(
  "Figure: ",
  FIG_PATH
)

message("============================================================")
