#' Contrast Analysis with Equivalence Testing
#'
#' Decomposes between-group variance into a theory-driven contrast of
#' interest and orthogonal residual contrasts, then tests whether the
#' predicted pattern is significant (NHST) and whether deviations from it
#' are negligibly small (TOST equivalence tests). Integrates Campbell (2024)
#' integrated with Campbell (2024) and the Intersection-Union Test (IUT;
#' Berger, 1982).
#'
#' @param data A data frame.
#' @param dv Name of the dependent variable column (string). Set `NULL`
#'   for within-subjects designs (wide format).
#' @param conditions Name of the grouping variable column (string). Set
#'   `NULL` for within-subjects designs.
#' @param hypothesis A named numeric vector specifying the predicted
#'   pattern. Names must match the levels of `conditions` (between) or
#'   column names (within). Values encode the predicted relative ordering
#'   (e.g., `c(Control = 0, Low = 1, High = 2)`). Only the relative
#'   spacing matters; the vector is centered and normalized internally.
#' @param design Either `"between"` (default) or `"within"`.
#' @param id Name of the participant ID column (string). Required for
#'   within-subjects designs.
#' @param delta Equivalence bound(s) for residual contrasts. A scalar
#'   (applied to all residuals) or a vector of length *k* - 2 (one per
#'   residual). Units depend on `delta_type`:
#'   \itemize{
#'     \item `"dz"`: standardized effect size (*d* units).
#'     \item `"share_signal"`: maximum tolerable percentage of model
#'       variance (between 0 and 100, exclusive).
#'   }
#' @param delta_type Either `"dz"` (default) or `"share_signal"`. See
#'   **Details**.
#' @param d_interest_min A single positive number: the smallest meaningful
#'   standardized effect (SESOI) for the contrast of interest. Required
#'   when `delta_type = "share_signal"`. Ignored when `delta_type = "dz"`.
#' @param alpha Significance level (default 0.05).
#' @param confirm If `TRUE`, display the setup (hypothesis, weights,
#'   bounds) and pause for confirmation before running the analysis.
#' @param bonferroni If `TRUE` (default) and contrasts are non-orthogonal,
#'   apply Bonferroni correction to the TOST alpha. Auto-generated
#'   contrasts are always orthogonal, so this rarely activates.
#' @param graph1 If `TRUE` (default), display the forest plot.
#' @param graph1.only If `TRUE`, return *only* the forest plot as a ggplot
#'   object, without printing text output.
#' @param graph2 If `TRUE`, display the variance decomposition plot.
#' @param graph2.only If `TRUE`, return *only* the variance plot as a
#'   ggplot object, without printing text output.
#'
#' @details
#' ## Effect Sizes
#'
#' - **Between-subjects:** `d = b / sigma`, where `b` is the contrast
#'   regression coefficient and `sigma = sqrt(MSE)`.
#' - **Within-subjects:** `d_z = M / SD`, where `M` and `SD` are the mean
#'   and standard deviation of within-subject contrast scores.
#'
#' ## Equivalence Bound Types
#'
#' When `delta_type = "dz"` (the default), `delta` is specified directly
#' in standardized *d* units.
#'
#' When `delta_type = "share_signal"`, `delta` is specified as a maximum
#' tolerable percentage of model variance. The function converts this to a
#' *d*-scale bound using:
#'
#' \deqn{\delta_d = d_{\text{interest, min}} \times \sqrt{\frac{p}{1 - p}}}
#'
#' where `p = delta / 100`. This conversion requires `d_interest_min` to
#' be specified. The bound is pre-specified (not data-dependent), preserving
#' valid Type I error control.
#'
#' ## Conclusions
#'
#' - **SUPPORTED:** Contrast of interest is significant AND all residuals
#'   are equivalent to zero.
#' - **PARTIALLY SUPPORTED:** Contrast of interest is significant, BUT at
#'   least one residual is not equivalent to zero.
#' - **NOT SUPPORTED:** Contrast of interest is not significant.
#'
#' ## Alpha Correction
#'
#' Per the Intersection-Union Test (IUT; Berger, 1982), no alpha
#' correction is needed across orthogonal contrasts. Within each TOST,
#' each one-sided test is conducted at alpha (not alpha/2), following
#' Lakens (2017).
#'
#' @return If `graph1.only = TRUE`, a ggplot object (forest plot). If
#'   `graph2.only = TRUE`, a ggplot object (variance plot). Otherwise, the
#'   function prints a formatted summary and invisibly returns a list with:
#'   \describe{
#'     \item{interest}{List: b, se, t, df, p, d, se_d, ci_b, ci_d.}
#'     \item{residuals}{List of lists (one per residual), each with NHST
#'       and TOST results.}
#'     \item{model}{Fitted `lm` object (between) or `NULL` (within).}
#'     \item{design}{`"between"` or `"within"`.}
#'     \item{hypothesis}{The hypothesis vector supplied.}
#'     \item{hyp_levels}{Condition names in hypothesis order.}
#'     \item{w_interest}{Normalized contrast weights.}
#'     \item{residual_basis}{Matrix of residual contrast weights.}
#'     \item{n_residuals}{Number of residual contrasts (k - 2).}
#'     \item{delta_vec}{Equivalence bounds in d units (after conversion).}
#'     \item{delta_type}{`"dz"` or `"share_signal"`.}
#'     \item{delta_share_pct}{Share percentages specified (`NA` if dz).}
#'     \item{d_interest_min}{Reference d for interest (`NA` if dz).}
#'     \item{alpha, alpha_tost}{Alpha levels used.}
#'     \item{is_orthogonal}{Logical: are contrasts orthogonal?}
#'     \item{bonferroni_applied}{Logical: was Bonferroni applied?}
#'     \item{contrast_matrix}{Full contrast matrix (interest + residuals).}
#'     \item{k}{Number of groups/conditions.}
#'     \item{conclusion}{`"SUPPORTED"`, `"PARTIALLY SUPPORTED"`, or
#'       `"NOT SUPPORTED"`.}
#'     \item{conclusion_text}{Plain-English summary.}
#'     \item{residuals_equiv}{Logical vector: equivalence verdict per
#'       residual.}
#'     \item{apa}{APA-formatted results string.}
#'     \item{variance}{List: ss_interest, ss_residuals, ss_model,
#'       pct_interest, pct_residuals.}
#'     \item{plot_forest}{ggplot object or NULL.}
#'     \item{plot_variance}{ggplot object or NULL.}
#'   }
#'
#' @references
#' Berger, R. L. (1982). Multiparameter hypothesis testing and acceptance
#'   sampling. *Technometrics*, *24*(4), 295--300.
#'
#' Campbell, H. (2024). Equivalence testing for linear regression.
#'   *Psychological Methods*, *29*(1), 88--98.
#'   https://doi.org/10.1037/met0000596
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for *t* tests,
#'   correlations, and meta-analyses. *Social Psychological and Personality
#'   Science*, *8*(4), 355--362.
#'
#' @examples
#' # Between-subjects: linear dose-response
#' set.seed(42)
#' dat <- data.frame(
#'   group = rep(c("A","B","C","D"), each = 50),
#'   score = c(rnorm(50,2,3), rnorm(50,4,3), rnorm(50,6,3), rnorm(50,8,3))
#' )
#' result <- contrast_analysis(
#'   dat, "score", "group",
#'   hypothesis = c(A=0, B=1, C=2, D=3),
#'   delta = 0.3
#' )
#'
#' # Within-subjects: Control < [DrugA = DrugB]
#' set.seed(77)
#' Sigma <- matrix(c(9,4.5,4.5, 4.5,9,4.5, 4.5,4.5,9), 3)
#' scores <- MASS::mvrnorm(50, mu = c(5,8,8), Sigma = Sigma)
#' dat_w <- data.frame(pid=1:50, Ctrl=scores[,1],
#'                     DrugA=scores[,2], DrugB=scores[,3])
#' result_w <- contrast_analysis(
#'   dat_w, dv=NULL, conditions=NULL,
#'   hypothesis = c(Ctrl=0, DrugA=1, DrugB=1),
#'   design = "within", id = "pid", delta = 0.3
#' )
#'
#' # Using share_signal: residuals < 5% of model variance
#' result_s <- contrast_analysis(
#'   dat, "score", "group",
#'   hypothesis = c(A=0, B=1, C=2, D=3),
#'   delta = 5, delta_type = "share_signal",
#'   d_interest_min = 0.8
#' )
#'
#' @export
#' @importFrom MASS Null
#' @importFrom stats lm anova as.formula pt qt sd rnorm
#' @importFrom ggplot2 ggplot aes geom_vline geom_segment geom_point
#'   annotate labs theme_minimal theme element_blank element_text margin
#'   geom_col geom_text scale_y_continuous expansion
contrast_analysis <- function(
    data,
    dv,
    conditions,
    hypothesis,
    design      = "between",
    id          = NULL,
    delta       = 0.3,
    delta_type  = "dz",
    d_interest_min = NULL,
    alpha       = 0.05,
    confirm     = FALSE,
    bonferroni  = TRUE,
    graph1      = TRUE,
    graph1.only = FALSE,
    graph2      = FALSE,
    graph2.only = FALSE
) {

  # =========================================================================
  # 1. INPUT VALIDATION
  # =========================================================================
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  design <- match.arg(design, c("between", "within"))

  # --- Validate hypothesis ---
  if (!is.numeric(hypothesis) || is.null(names(hypothesis))) {
    stop("`hypothesis` must be a named numeric vector.\n",
         "  Names = condition/group labels.\n",
         "  Values = predicted relative pattern (e.g., c(Control=0, Low=1, High=2)).")
  }
  k <- length(hypothesis)
  if (k < 3) stop("Need at least 3 groups/conditions for contrast analysis.")

  # --- Design-specific validation ---
  if (design == "between") {
    if (!dv %in% names(data)) stop(paste0("DV column '", dv, "' not found in data."))
    if (!conditions %in% names(data)) stop(paste0("Conditions column '", conditions, "' not found in data."))
    if (!is.numeric(data[[dv]])) stop("DV must be numeric (continuous).")
    data[[conditions]] <- as.factor(data[[conditions]])
    obs_levels <- levels(data[[conditions]])
    hyp_levels <- names(hypothesis)
    if (!all(hyp_levels %in% obs_levels)) {
      stop("Hypothesis names must match factor levels of `conditions`.\n",
           "  Hypothesis levels: ", paste(hyp_levels, collapse = ", "), "\n",
           "  Data levels:       ", paste(obs_levels, collapse = ", "))
    }
    if (!all(obs_levels %in% hyp_levels)) {
      warning("Some factor levels in data are not in hypothesis and will be dropped: ",
              paste(setdiff(obs_levels, hyp_levels), collapse = ", "))
      data <- data[data[[conditions]] %in% hyp_levels, , drop = FALSE]
      data[[conditions]] <- droplevels(data[[conditions]])
    }
  } else {
    if (is.null(id)) stop("For within-subjects designs, `id` (participant ID column) is required.")
    if (!id %in% names(data)) stop(paste0("ID column '", id, "' not found in data."))
    hyp_levels <- names(hypothesis)
    missing_cols <- hyp_levels[!hyp_levels %in% names(data)]
    if (length(missing_cols) > 0) {
      stop("Condition columns not found in data: ", paste(missing_cols, collapse = ", "))
    }
    for (col in hyp_levels) {
      if (!is.numeric(data[[col]])) stop(paste0("Column '", col, "' must be numeric."))
    }
  }

  # --- Validate delta_type ---
  delta_type <- match.arg(delta_type, c("dz", "share_signal"))

  # --- Validate delta ---
  n_residuals <- k - 2
  if (length(delta) == 1) {
    delta_vec <- rep(delta, n_residuals)
  } else if (length(delta) == n_residuals) {
    delta_vec <- delta
  } else {
    stop("`delta` must be length 1 (same for all residuals) or length ", n_residuals,
         " (one per residual contrast).")
  }

  # --- Process delta based on delta_type ---
  delta_share_pct <- NA_real_
  d_interest_min_used <- NA_real_

  if (delta_type == "share_signal") {
    if (is.null(d_interest_min) || !is.numeric(d_interest_min) || length(d_interest_min) != 1) {
      stop("`d_interest_min` must be a single positive number when delta_type = 'share_signal'.\n",
           "  This is the smallest meaningful standardized effect (SESOI) for the\n",
           "  contrast of interest.")
    }
    if (d_interest_min <= 0) stop("`d_interest_min` must be positive.")
    if (any(delta_vec <= 0) || any(delta_vec >= 100)) {
      stop("When delta_type = 'share_signal', delta values must be between 0 and 100 (exclusive),\n",
           "  representing the maximum tolerable % of model variance for each residual.\n",
           "  Example: delta = 5 means each residual may account for at most 5% of model variance.")
    }
    delta_share_pct <- delta_vec
    d_interest_min_used <- d_interest_min
    p_share <- delta_vec / 100
    delta_vec <- d_interest_min * sqrt(p_share / (1 - p_share))
  } else {
    if (any(delta_vec <= 0)) stop("All delta values must be positive.")
  }

  # =========================================================================
  # 2. CONSTRUCT CONTRAST OF INTEREST
  # =========================================================================
  h_raw <- hypothesis[hyp_levels]
  w_interest <- h_raw - mean(h_raw)
  w_interest_norm <- w_interest / sqrt(sum(w_interest^2))
  names(w_interest_norm) <- hyp_levels

  # =========================================================================
  # 3. CONSTRUCT RESIDUAL CONTRASTS (orthogonal complement)
  # =========================================================================
  ones <- rep(1, k)
  constraint_mat <- cbind(ones, as.numeric(w_interest_norm))
  residual_basis <- MASS::Null(constraint_mat)

  if (n_residuals > 0) {
    colnames(residual_basis) <- paste0("Residual_", seq_len(n_residuals))
    rownames(residual_basis) <- hyp_levels
  }

  contrast_matrix <- cbind(Interest = w_interest_norm, residual_basis)

  # =========================================================================
  # 4. CHECK ORTHOGONALITY
  # =========================================================================
  cross_products <- crossprod(contrast_matrix)
  diag(cross_products) <- 0
  is_orthogonal <- all(abs(cross_products) < 1e-10)

  alpha_tost <- alpha
  bonferroni_applied <- FALSE
  if (!is_orthogonal && bonferroni && n_residuals > 1) {
    alpha_tost <- alpha / n_residuals
    bonferroni_applied <- TRUE
  }

  # =========================================================================
  # 5. CONFIRM WITH USER (if requested)
  # =========================================================================
  if (confirm) {
    cat("\n===== CONTRAST ANALYSIS: CONFIRMATION =====\n\n")
    cat("HYPOTHESIZED PATTERN:\n")
    ordered_idx <- order(h_raw)
    ordered_names <- hyp_levels[ordered_idx]
    ordered_vals  <- h_raw[ordered_idx]
    val_groups <- split(ordered_names, ordered_vals[order(ordered_vals)])
    pattern_parts <- vapply(val_groups, function(g) {
      if (length(g) == 1) g else paste0("[", paste(g, collapse = " = "), "]")
    }, character(1))
    cat("  ", paste(pattern_parts, collapse = " < "), "\n")
    cat("  (Groups in brackets are predicted to be equal.)\n\n")

    cat("CONTRAST OF INTEREST (weights, normalized):\n")
    for (i in seq_along(hyp_levels)) {
      cat(sprintf("  %-20s: %+.4f\n", hyp_levels[i], w_interest_norm[i]))
    }

    cat("\nRESIDUAL CONTRASTS (", n_residuals, " total):\n", sep = "")
    if (n_residuals > 0) {
      for (j in seq_len(n_residuals)) {
        cat(sprintf("  Residual %d:\n", j))
        for (i in seq_along(hyp_levels)) {
          cat(sprintf("    %-20s: %+.4f\n", hyp_levels[i], residual_basis[i, j]))
        }
      }
    }

    cat("\nORTHOGONALITY: ",
        ifelse(is_orthogonal, "All contrasts are orthogonal (no correction needed).",
               "WARNING: Contrasts are NOT orthogonal."),
        "\n")
    if (bonferroni_applied) {
      cat("  Bonferroni correction applied to residual TOST alpha: ",
          round(alpha_tost, 4), "\n")
    }

    cat("\nEQUIVALENCE BOUNDS:\n")
    if (delta_type == "share_signal") {
      cat("  Delta type: share_signal (% of model variance)\n")
      cat("  Reference d_interest_min: ", d_interest_min_used, "\n")
      for (j in seq_len(n_residuals)) {
        cat(sprintf("  Residual %d: max share = %.1f%% => converted d bound = %.4f\n",
                    j, delta_share_pct[j], delta_vec[j]))
      }
    } else {
      cat("  Delta type: dz (standardized effect size)\n")
      for (j in seq_len(n_residuals)) {
        cat(sprintf("  Residual %d: delta = %.3f\n", j, delta_vec[j]))
      }
    }

    cat("\nDesign: ", design, "-subjects\n", sep = "")
    cat("Alpha: ", alpha, "\n")
    cat("\nProceed? [Enter to continue, Esc to abort] ")
    readline()
    cat("\n")
  }

  # =========================================================================
  # 6. FIT MODEL & EXTRACT RESULTS
  # =========================================================================
  if (design == "between") {
    results <- .ca_between(data, dv, conditions, hyp_levels, w_interest_norm,
                           residual_basis, n_residuals, delta_vec,
                           alpha, alpha_tost, k)
  } else {
    results <- .ca_within(data, id, hyp_levels, w_interest_norm,
                          residual_basis, n_residuals, delta_vec,
                          alpha, alpha_tost, k)
  }

  # Attach metadata
  results$design           <- design
  results$hypothesis       <- hypothesis
  results$hyp_levels       <- hyp_levels
  results$w_interest       <- w_interest_norm
  results$residual_basis   <- residual_basis
  results$n_residuals      <- n_residuals
  results$delta_vec        <- delta_vec
  results$delta_type       <- delta_type
  results$delta_share_pct  <- delta_share_pct
  results$d_interest_min   <- d_interest_min_used
  results$alpha            <- alpha
  results$alpha_tost       <- alpha_tost
  results$is_orthogonal    <- is_orthogonal
  results$bonferroni_applied <- bonferroni_applied
  results$contrast_matrix  <- contrast_matrix
  results$k                <- k

  # =========================================================================
  # 7. QUALITATIVE CONCLUSION
  # =========================================================================
  interest_sig <- results$interest$p < alpha

  if (n_residuals > 0) {
    residuals_equiv <- vapply(results$residuals, function(r) r$p_tost < alpha_tost, logical(1))
    all_residuals_equiv <- all(residuals_equiv)
  } else {
    residuals_equiv <- logical(0)
    all_residuals_equiv <- TRUE
  }

  if (interest_sig && all_residuals_equiv) {
    conclusion <- "SUPPORTED"
    conclusion_text <- paste0(
      "The data fully support the hypothesized pattern. The contrast of interest ",
      "is statistically significant, and all residual contrasts are equivalent to ",
      "zero (within the specified equivalence bounds).")
  } else if (interest_sig && !all_residuals_equiv) {
    conclusion <- "PARTIALLY SUPPORTED"
    conclusion_text <- paste0(
      "The data partially support the hypothesized pattern. The contrast of interest ",
      "is statistically significant, but not all residual contrasts could be declared ",
      "equivalent to zero. The hypothesized pattern does not fully account for the ",
      "between-condition differences.")
  } else {
    conclusion <- "NOT SUPPORTED"
    conclusion_text <- paste0(
      "The data do not support the hypothesized pattern. The contrast of interest ",
      "is not statistically significant.")
  }

  results$conclusion      <- conclusion
  results$conclusion_text <- conclusion_text
  results$residuals_equiv <- residuals_equiv

  # =========================================================================
  # 8. APA OUTPUT, VARIANCE, GRAPHS
  # =========================================================================
  results$apa      <- .ca_apa_output(results, alpha, alpha_tost)
  results$variance <- .ca_variance_decomposition(results)

  p1 <- NULL
  if (graph1 || graph1.only) {
    p1 <- .ca_forest_plot(results)
    results$plot_forest <- p1
  }

  p2 <- NULL
  if (graph2 || graph2.only) {
    p2 <- .ca_variance_plot(results)
    results$plot_variance <- p2
  }

  # =========================================================================
  # 9. OUTPUT LOGIC
  # =========================================================================
  if (graph1.only) return(p1)
  if (graph2.only) return(p2)

  .ca_print_summary(results)
  if (!is.null(p1)) print(p1)
  if (!is.null(p2)) print(p2)
  invisible(results)
}
