# =========================================================================
# INTERNAL: Between-subjects analysis
# =========================================================================
#' @keywords internal
.ca_between <- function(data, dv, conditions, hyp_levels, w_interest,
                        residual_basis, n_residuals, delta_vec,
                        alpha, alpha_tost, k) {

  data$.C_Interest <- w_interest[as.character(data[[conditions]])]

  if (n_residuals > 0) {
    for (j in seq_len(n_residuals)) {
      cname <- paste0(".C_Res", j)
      data[[cname]] <- residual_basis[as.character(data[[conditions]]), j]
    }
  }

  predictors <- ".C_Interest"
  if (n_residuals > 0) {
    predictors <- c(predictors, paste0(".C_Res", seq_len(n_residuals)))
  }
  fml <- as.formula(paste(dv, "~", paste(predictors, collapse = " + ")))
  model <- lm(fml, data = data)
  smry  <- summary(model)
  coefs <- smry$coefficients
  mse   <- smry$sigma^2
  sigma <- smry$sigma
  df_resid <- model$df.residual
  N <- nrow(data)

  # --- Contrast of interest: NHST ---
  b_int  <- coefs[".C_Interest", "Estimate"]
  se_int <- coefs[".C_Interest", "Std. Error"]
  t_int  <- coefs[".C_Interest", "t value"]
  p_int  <- coefs[".C_Interest", "Pr(>|t|)"]
  d_int  <- b_int / sigma

  ci_b_int <- b_int + c(-1, 1) * qt(1 - alpha/2, df_resid) * se_int
  se_d_int <- se_int / sigma
  ci_d_int <- d_int + c(-1, 1) * qt(1 - alpha/2, df_resid) * se_d_int

  interest <- list(
    b = b_int, se = se_int, t = t_int, df = df_resid, p = p_int,
    d = d_int, se_d = se_d_int,
    ci_b = ci_b_int, ci_d = ci_d_int,
    label = "Contrast of Interest"
  )

  # --- Residual contrasts: NHST + TOST ---
  residuals_list <- list()
  if (n_residuals > 0) {
    for (j in seq_len(n_residuals)) {
      cname <- paste0(".C_Res", j)
      b_res  <- coefs[cname, "Estimate"]
      se_res <- coefs[cname, "Std. Error"]
      t_res  <- coefs[cname, "t value"]
      p_res  <- coefs[cname, "Pr(>|t|)"]
      d_res  <- b_res / sigma
      se_d_res <- se_res / sigma

      ci_d_res <- d_res + c(-1, 1) * qt(1 - alpha/2, df_resid) * se_d_res

      delta_raw <- delta_vec[j] * sigma

      t1 <- (b_res + delta_raw) / se_res
      t2 <- (b_res - delta_raw) / se_res
      p1 <- pt(t1, df_resid, lower.tail = FALSE)
      p2 <- pt(t2, df_resid, lower.tail = TRUE)
      p_tost <- max(p1, p2)

      ci_level_tost <- 1 - 2 * alpha_tost
      ci_d_tost <- d_res + c(-1, 1) * qt(1 - alpha_tost, df_resid) * se_d_res

      residuals_list[[j]] <- list(
        b = b_res, se = se_res, t = t_res, df = df_resid,
        p_nhst = p_res, d = d_res, se_d = se_d_res,
        ci_d_nhst = ci_d_res,
        t1 = t1, t2 = t2, p1 = p1, p2 = p2, p_tost = p_tost,
        delta = delta_vec[j], delta_raw = delta_raw,
        ci_d_tost = ci_d_tost, ci_level_tost = ci_level_tost,
        label = paste0("Residual ", j)
      )
    }
  }

  anova_table <- anova(model)

  list(interest = interest, residuals = residuals_list, model = model,
       anova = anova_table, mse = mse, sigma = sigma, N = N,
       df_resid = df_resid)
}


# =========================================================================
# INTERNAL: Within-subjects analysis
# =========================================================================
#' @keywords internal
.ca_within <- function(data, id, hyp_levels, w_interest,
                       residual_basis, n_residuals, delta_vec,
                       alpha, alpha_tost, k) {

  N <- nrow(data)
  cond_mat <- as.matrix(data[, hyp_levels, drop = FALSE])

  cs_interest <- as.numeric(cond_mat %*% w_interest)

  cs_residuals <- list()
  if (n_residuals > 0) {
    for (j in seq_len(n_residuals)) {
      cs_residuals[[j]] <- as.numeric(cond_mat %*% residual_basis[, j])
    }
  }

  # --- Contrast of interest: one-sample t-test ---
  m_int  <- mean(cs_interest)
  sd_int <- sd(cs_interest)
  se_int <- sd_int / sqrt(N)
  t_int  <- m_int / se_int
  df_int <- N - 1
  p_int  <- 2 * pt(abs(t_int), df_int, lower.tail = FALSE)

  d_z_int   <- m_int / sd_int
  se_dz_int <- sqrt(1/N + d_z_int^2 / (2*N))
  ci_dz_int <- d_z_int + c(-1, 1) * qt(1 - alpha/2, df_int) * se_dz_int

  interest <- list(
    b = m_int, se = se_int, t = t_int, df = df_int, p = p_int,
    d = d_z_int, se_d = se_dz_int,
    ci_b = m_int + c(-1, 1) * qt(1 - alpha/2, df_int) * se_int,
    ci_d = ci_dz_int,
    label = "Contrast of Interest"
  )

  # --- Residual contrasts: one-sample TOST ---
  residuals_list <- list()
  if (n_residuals > 0) {
    for (j in seq_len(n_residuals)) {
      cs_j <- cs_residuals[[j]]
      m_j  <- mean(cs_j)
      sd_j <- sd(cs_j)
      se_j <- sd_j / sqrt(N)
      t_j  <- m_j / se_j
      df_j <- N - 1
      p_j  <- 2 * pt(abs(t_j), df_j, lower.tail = FALSE)

      d_z_j   <- m_j / sd_j
      se_dz_j <- sqrt(1/N + d_z_j^2 / (2*N))
      ci_dz_nhst <- d_z_j + c(-1, 1) * qt(1 - alpha/2, df_j) * se_dz_j

      delta_raw <- delta_vec[j] * sd_j

      t1 <- (m_j + delta_raw) / se_j
      t2 <- (m_j - delta_raw) / se_j
      p1 <- pt(t1, df_j, lower.tail = FALSE)
      p2 <- pt(t2, df_j, lower.tail = TRUE)
      p_tost <- max(p1, p2)

      ci_level_tost <- 1 - 2 * alpha_tost
      ci_dz_tost <- d_z_j + c(-1, 1) * qt(1 - alpha_tost, df_j) * se_dz_j

      residuals_list[[j]] <- list(
        b = m_j, se = se_j, t = t_j, df = df_j,
        p_nhst = p_j, d = d_z_j, se_d = se_dz_j,
        ci_d_nhst = ci_dz_nhst,
        t1 = t1, t2 = t2, p1 = p1, p2 = p2, p_tost = p_tost,
        delta = delta_vec[j], delta_raw = delta_raw,
        ci_d_tost = ci_dz_tost, ci_level_tost = ci_level_tost,
        label = paste0("Residual ", j)
      )
    }
  }

  ss_interest <- N * m_int^2
  ss_residuals <- if (n_residuals > 0) {
    vapply(cs_residuals, function(cs) N * mean(cs)^2, numeric(1))
  } else {
    numeric(0)
  }

  list(interest = interest, residuals = residuals_list, model = NULL,
       ss_interest = ss_interest, ss_residuals = ss_residuals,
       N = N, df_resid = N - 1)
}
