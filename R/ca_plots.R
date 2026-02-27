# =========================================================================
# INTERNAL: Forest plot (graph1)
# =========================================================================
#' @keywords internal
.ca_forest_plot <- function(results) {

  alpha      <- results$alpha
  alpha_tost <- results$alpha_tost
  n_res      <- results$n_residuals
  delta_vec  <- results$delta_vec
  delta_max  <- max(delta_vec)

  # --- Build data frame ---
  r <- results$interest
  plot_df <- data.frame(
    label    = r$label,
    d        = r$d,
    ci_low   = r$ci_d[1],
    ci_high  = r$ci_d[2],
    type     = "interest",
    res_idx  = NA_integer_,
    res_delta = NA_real_,
    stringsAsFactors = FALSE
  )

  if (n_res > 0) {
    for (j in seq_len(n_res)) {
      rr <- results$residuals[[j]]
      plot_df <- rbind(plot_df, data.frame(
        label    = rr$label,
        d        = rr$d,
        ci_low   = rr$ci_d_tost[1],
        ci_high  = rr$ci_d_tost[2],
        type     = "residual",
        res_idx  = j,
        res_delta = delta_vec[j],
        stringsAsFactors = FALSE
      ))
    }
  }

  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

  # --- Determine colors ---
  colors <- character(nrow(plot_df))
  for (i in seq_len(nrow(plot_df))) {
    if (plot_df$type[i] == "interest") {
      excludes_zero <- plot_df$ci_low[i] > 0 | plot_df$ci_high[i] < 0
      colors[i] <- ifelse(excludes_zero, "#2E8B57", "#CD5C5C")
    } else {
      dj <- plot_df$res_delta[i]
      within_bounds <- plot_df$ci_low[i] >= -dj & plot_df$ci_high[i] <= dj
      colors[i] <- ifelse(within_bounds, "#2E8B57", "#CD5C5C")
    }
  }
  plot_df$color <- colors

  # --- Caption ---
  ci_pct_interest <- round((1 - alpha) * 100)
  ci_pct_tost <- if (n_res > 0) round(results$residuals[[1]]$ci_level_tost * 100) else NA

  ci_note <- paste0(
    "Contrast of interest: ", ci_pct_interest, "% CI; ",
    "Residual contrasts: ", ci_pct_tost, "% CI (TOST-consistent)"
  )

  if (n_res > 0 && length(unique(delta_vec)) > 1) {
    ci_note <- paste0(ci_note,
                      "\nShaded region = largest delta (+/-", round(delta_max, 3),
                      "); coloring uses per-residual delta")
  }

  if (!is.null(results$delta_type) && results$delta_type == "share_signal") {
    ci_note <- paste0(ci_note,
                      "\nEquivalence bounds derived from share_signal: max share = ",
                      paste(paste0(round(results$delta_share_pct, 1), "%"), collapse = ", "),
                      "; d_interest_min = ", round(results$d_interest_min, 2))
  }

  # --- Build plot ---
  p <- ggplot(plot_df, aes(x = d, y = label)) +
    annotate("rect",
             xmin = -delta_max, xmax = delta_max,
             ymin = -Inf, ymax = Inf,
             fill = "#B0C4DE", alpha = 0.25) +
    geom_vline(xintercept = c(-delta_max, delta_max),
               linetype = "solid", color = "#4682B4", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    {
      if (n_res > 0 && length(unique(delta_vec)) > 1) {
        lapply(unique(delta_vec), function(dj) {
          list(
            geom_vline(xintercept = c(-dj, dj),
                       linetype = "dotted", color = "#4682B4", linewidth = 0.3)
          )
        })
      }
    } +
    geom_segment(aes(x = ci_low, xend = ci_high, y = label, yend = label),
                 color = plot_df$color, linewidth = 1.2) +
    geom_point(aes(x = d, y = label), color = plot_df$color,
               size = 3.5, shape = 18) +
    labs(
      x = "Standardized Effect (d)",
      y = NULL,
      caption = ci_note
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y  = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 8, hjust = 0, color = "grey40"),
      plot.margin  = margin(10, 15, 10, 10)
    )

  return(p)
}


# =========================================================================
# INTERNAL: Variance decomposition plot (graph2)
# =========================================================================
#' @keywords internal
.ca_variance_plot <- function(results) {

  v     <- results$variance
  n_res <- results$n_residuals

  labels <- c("Contrast\nof Interest",
              if (n_res > 0) paste0("Residual\n", seq_len(n_res)) else NULL)
  pcts   <- c(v$pct_interest, v$pct_residuals)
  types  <- c("interest", rep("residual", n_res))

  bar_df <- data.frame(
    label = factor(labels, levels = labels),
    pct   = pcts,
    type  = types,
    stringsAsFactors = FALSE
  )

  p <- ggplot(bar_df, aes(x = label, y = pct)) +
    geom_col(data = bar_df[bar_df$type == "interest", , drop = FALSE],
             fill = "grey30", color = "grey30", width = 0.6) +
    {
      if (n_res > 0) {
        geom_col(data = bar_df[bar_df$type == "residual", , drop = FALSE],
                 fill = NA, color = "grey30", linetype = "dashed", width = 0.6)
      }
    } +
    geom_text(aes(label = paste0(sprintf("%.1f", pct), "%")),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.08)),
                       labels = function(x) paste0(x, "%")) +
    labs(
      x     = NULL,
      y     = "Percentage of Model Variance",
      title = "Variance Decomposition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x  = element_text(face = "bold", size = 10),
      plot.title   = element_text(face = "bold", size = 13, hjust = 0.5)
    )

  return(p)
}
