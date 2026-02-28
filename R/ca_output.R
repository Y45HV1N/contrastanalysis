# =========================================================================
# INTERNAL: APA-formatted output
# =========================================================================
#' @keywords internal
.ca_apa_output <- function(results, alpha, alpha_tost) {
  .fmt_p <- function(p) {
    if (p < .001) return("< .001")
    formatted <- sub("^0", "", sprintf("%.3f", p))
    paste("=", formatted)
  }
  .fmt_num <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)

  lines <- character(0)
  r <- results$interest

  lines <- c(lines, paste0(
    "Contrast of interest: b = ", .fmt_num(r$b), ", SE = ", .fmt_num(r$se),
    ", t(", r$df, ") = ", .fmt_num(r$t), ", p ", .fmt_p(r$p),
    ", d = ", .fmt_num(r$d),
    ", ", round((1 - alpha) * 100), "% CI [",
    .fmt_num(r$ci_d[1]), ", ", .fmt_num(r$ci_d[2]), "]."
  ))

  for (j in seq_along(results$residuals)) {
    rr <- results$residuals[[j]]
    equiv_label <- if (rr$p_tost < alpha_tost) "equivalent to zero" else "not equivalent to zero"
    ci_pct <- round(rr$ci_level_tost * 100)

    if (!is.null(results$delta_type) && results$delta_type == "share_signal") {
      bounds_desc <- paste0(
        "equivalence bounds = +/-", .fmt_num(rr$delta, 3),
        " (from ", .fmt_num(results$delta_share_pct[j], 1),
        "% share, d_int_min = ", .fmt_num(results$d_interest_min, 2), ")")
    } else {
      bounds_desc <- paste0("equivalence bounds = +/-", .fmt_num(rr$delta, 2))
    }

    lines <- c(lines, paste0(
      "Residual contrast ", j, " (TOST): b = ", .fmt_num(rr$b),
      ", SE = ", .fmt_num(rr$se),
      ", d = ", .fmt_num(rr$d),
      ", ", ci_pct, "% CI [",
      .fmt_num(rr$ci_d_tost[1]), ", ", .fmt_num(rr$ci_d_tost[2]), "]",
      ", ", bounds_desc,
      ", t1(", rr$df, ") = ", .fmt_num(rr$t1),
      ", p1 ", .fmt_p(rr$p1),
      ", t2(", rr$df, ") = ", .fmt_num(rr$t2),
      ", p2 ", .fmt_p(rr$p2),
      ", p_TOST ", .fmt_p(rr$p_tost),
      " (", equiv_label, ")."
    ))
  }

  paste(lines, collapse = "\n")
}


# =========================================================================
# INTERNAL: Variance decomposition
# =========================================================================
#' @keywords internal
.ca_variance_decomposition <- function(results) {
  if (!is.null(results$model)) {
    at <- results$anova
    ss_all <- at[["Sum Sq"]]
    n_pred <- nrow(at) - 1
    ss_interest  <- ss_all[1]
    ss_residuals <- if (n_pred > 1) ss_all[2:n_pred] else numeric(0)
    ss_error     <- ss_all[nrow(at)]
    ss_model     <- sum(ss_all[1:n_pred])

    pct_interest  <- (ss_interest / ss_model) * 100
    pct_residuals <- if (length(ss_residuals) > 0) (ss_residuals / ss_model) * 100 else numeric(0)
  } else {
    ss_interest  <- results$ss_interest
    ss_residuals <- results$ss_residuals
    ss_model     <- ss_interest + sum(ss_residuals)

    pct_interest  <- (ss_interest / ss_model) * 100
    pct_residuals <- if (length(ss_residuals) > 0) (ss_residuals / ss_model) * 100 else numeric(0)
  }

  list(
    ss_interest   = ss_interest,
    ss_residuals  = ss_residuals,
    ss_model      = ss_model,
    pct_interest  = pct_interest,
    pct_residuals = pct_residuals
  )
}


# =========================================================================
# INTERNAL: Print summary
# =========================================================================
#' @keywords internal
.ca_print_summary <- function(results) {
  cat("\n")
  cat("=================================================================\n")
  cat("      CONTRAST ANALYSIS WITH EQUIVALENCE TESTING\n")
  cat("=================================================================\n\n")

  cat("Design:  ", results$design, "-subjects\n", sep = "")
  cat("Groups:  ", paste(results$hyp_levels, collapse = ", "), " (k = ", results$k, ")\n", sep = "")
  cat("N:       ", results$N, "\n", sep = "")
  cat("Alpha:   ", results$alpha, "\n", sep = "")
  if (results$bonferroni_applied) {
    cat("NOTE: Bonferroni correction applied to TOST alpha: ", results$alpha_tost, "\n", sep = "")
  }
  cat("Orthogonal contrasts: ", ifelse(results$is_orthogonal, "Yes", "No"), "\n", sep = "")

  if (results$design == "between") {
    cat("Effect size metric: d = b / sigma (b = contrast coefficient, sigma = sqrt(MSE))\n")
  } else {
    cat("Effect size metric: d_z = M / SD (M and SD of within-subject contrast scores)\n")
  }
  cat("\n")

  # --- Hypothesized pattern ---
  h <- results$hypothesis
  ordered_idx   <- order(h)
  ordered_names <- names(h)[ordered_idx]
  ordered_vals  <- h[ordered_idx]
  val_groups    <- split(ordered_names, ordered_vals)
  pattern_parts <- vapply(val_groups, function(g) {
    if (length(g) == 1) g else paste0("[", paste(g, collapse = " = "), "]")
  }, character(1))
  cat("Hypothesized pattern: ", paste(pattern_parts, collapse = " < "), "\n\n", sep = "")

  # --- Contrast weights table ---
  cat("----- Contrast Weights -----\n")
  cm <- results$contrast_matrix
  # Header row
  col_labels <- c("Interest", paste0("Residual ", seq_len(results$n_residuals)))
  header <- sprintf("  %-12s", "")
  for (cl in col_labels) header <- paste0(header, sprintf("%10s", cl))
  cat(header, "\n")
  # Data rows
  for (i in seq_along(results$hyp_levels)) {
    row_str <- sprintf("  %-12s", results$hyp_levels[i])
    for (j in seq_len(ncol(cm))) {
      row_str <- paste0(row_str, sprintf("%10.4f", cm[i, j]))
    }
    cat(row_str, "\n")
  }
  if (results$n_residuals > 0) {
    cat("\n  Note: Residual contrasts are one possible orthonormal basis for the\n")
    cat("  subspace orthogonal to the predicted pattern. Individual residual\n")
    cat("  vectors are not uniquely defined, but their joint test is invariant:\n")
    cat("  if any systematic deviation from the predicted pattern exists, at\n")
    cat("  least one residual will detect it.\n")
  }
  cat("\n")

  # --- Contrast of interest ---
  cat("----- Contrast of Interest (NHST) -----\n")
  r <- results$interest
  sig_label <- ifelse(r$p < results$alpha, "SIGNIFICANT", "not significant")
  cat(sprintf("  b = %.3f, SE = %.3f, t(%d) = %.3f, p %s [%s]\n",
              r$b, r$se, r$df, r$t,
              ifelse(r$p < .001, "< .001", sprintf("= %.3f", r$p)),
              sig_label))
  cat(sprintf("  d = %.3f, %d%% CI [%.3f, %.3f]\n\n",
              r$d, round((1 - results$alpha) * 100), r$ci_d[1], r$ci_d[2]))

  # --- Residual contrasts ---
  if (results$n_residuals > 0) {
    cat("----- Residual Contrasts (TOST Equivalence Tests) -----\n")

    if (!is.null(results$delta_type) && results$delta_type == "share_signal") {
      cat("  Delta type: share_signal (% of model variance)\n")
      cat("  Reference d_interest_min: ", results$d_interest_min, "\n", sep = "")
      cat("  Equivalence bounds (converted to d-scale):\n")
      for (j in seq_len(results$n_residuals)) {
        cat(sprintf("    Residual %d: max share = %.1f%% => d bound = +/-%.4f\n",
                    j, results$delta_share_pct[j], results$delta_vec[j]))
      }
    } else {
      cat("  Equivalence bounds (standardized): ",
          paste(sprintf("+/-%.3f", results$delta_vec), collapse = ", "), "\n", sep = "")
    }
    cat("\n")

    for (j in seq_along(results$residuals)) {
      rr <- results$residuals[[j]]
      equiv_label <- ifelse(rr$p_tost < results$alpha_tost,
                            "EQUIVALENT (within bounds)",
                            "NOT equivalent")
      ci_pct <- round(rr$ci_level_tost * 100)

      cat(sprintf("  Residual %d:\n", j))
      # Show weights inline
      wts <- results$residual_basis[, j]
      wt_str <- paste(sprintf("%s: %+.3f", results$hyp_levels, wts), collapse = ", ")
      cat("    Weights: ", wt_str, "\n", sep = "")
      cat(sprintf("    NHST: b = %.3f, SE = %.3f, t(%d) = %.3f, p %s\n",
                  rr$b, rr$se, rr$df, rr$t,
                  ifelse(rr$p_nhst < .001, "< .001", sprintf("= %.3f", rr$p_nhst))))
      cat(sprintf("    TOST: d = %.3f, %d%% CI [%.3f, %.3f]\n",
                  rr$d, ci_pct, rr$ci_d_tost[1], rr$ci_d_tost[2]))
      cat(sprintf("          p1 %s, p2 %s, p_TOST %s [%s]\n\n",
                  ifelse(rr$p1 < .001, "< .001", sprintf("= %.3f", rr$p1)),
                  ifelse(rr$p2 < .001, "< .001", sprintf("= %.3f", rr$p2)),
                  ifelse(rr$p_tost < .001, "< .001", sprintf("= %.3f", rr$p_tost)),
                  equiv_label))
    }
  }

  # --- Variance decomposition ---
  cat("----- Variance Decomposition -----\n")
  v <- results$variance
  cat(sprintf("  Contrast of Interest: %.1f%%\n", v$pct_interest))
  for (j in seq_along(v$pct_residuals)) {
    cat(sprintf("  Residual %d:           %.1f%%\n", j, v$pct_residuals[j]))
  }
  cat("\n")

  # --- Share_signal interpretation note ---
  if (!is.null(results$delta_type) && results$delta_type == "share_signal") {
    cat("----- Share-Based Equivalence Interpretation -----\n")
    cat("  Equivalence bounds were derived from a maximum tolerable share of\n")
    cat("  model variance, assuming the contrast of interest is at least\n")
    cat("  d_interest_min = ", results$d_interest_min, " in standardized units.\n", sep = "")
    cat("  Formula: delta_d = d_interest_min * sqrt(p / (1 - p)),\n")
    cat("  where p = max_share / 100.\n")
    cat("  This interpretation holds when the observed interest effect is at\n")
    cat("  least d_interest_min; the conclusion logic already enforces this\n")
    cat("  implicitly (interest must be significant for SUPPORTED).\n\n")
  }

  # --- Conclusion ---
  cat("=================================================================\n")
  cat("  CONCLUSION: ", results$conclusion, "\n", sep = "")
  cat("=================================================================\n")
  cat(strwrap(results$conclusion_text, width = 65, prefix = "  "), sep = "\n")
  cat("\n")

  # --- APA output ---
  cat("----- APA-Formatted Results -----\n")
  cat(results$apa, "\n\n")
}
