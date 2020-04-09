

PPP_model <- function(FQ = NULL, t = NULL, alpha = NULL, truncation = NULL, truncation_type = "lp", dispersion = 1, Status = 0, Comment = "OK") {
  obj <- list(FQ = NULL, t = NULL, alpha = NULL, truncation = NULL, truncation_type = "lp", dispersion = 1, Status = 2, Comment = "")
  class(obj) <- "PPP_model"
  if (is.null(FQ) && is.null(t) && is.null(alpha)) {
    return(obj)
  }
  if (!is.numeric(FQ)) {
    warning("FQ must be numeric")
    obj$Comment = "FQ must be numeric"
    return(obj)
  }
  if (length(FQ) != 1) {
    warning("FQ must have length 1")
    obj$Comment = "FQ must have length 1"
    return(obj)
  }
  if (FQ < 0) {
    warning("FQ must not be negative")
    obj$Comment = "FQ must not be negative"
    return(obj)
  }
  if (!is.numeric(t)) {
    warning("t must be numeric.")
    obj$Comment <- "t must be numeric."
    obj$Status <- 2
    return(obj)
  }
  k <-length(t)
  if (!is.numeric(alpha)) {
    warning("alpha must be numeric.")
    obj$Comment <- "alpha must be numeric."
    obj$Status <- 2
    return(obj)
  }
  if (k<1) {
    warning("t must have lenght >= 1.")
    obj$Comment <- "t must have lenght >= 1."
    obj$Status <- 2
    return(obj)
  }
  if (length(alpha) != k) {
    warning("t and alpha must have the same lenght.")
    obj$Comment <- "t and alpha must have the same lenght."
    obj$Status <- 2
    return(obj)
  }
  if (sum(t > 0, na.rm = T) < k) {
    warning("Elements of t must be positive.")
    obj$Comment <- "Elements of t must be positive."
    obj$Status <- 2
    return(obj)
  }
  if (sum(alpha >= 0, na.rm = T) < k) {
    warning("Entries of alpha must be non-negative.")
    obj$Comment <- "Entries of alpha must be non-negative."
    obj$Status <- 2
    return(obj)
  }
  if (alpha[k] <= 0) {
    warning("Last entry alpha must be positive.")
    obj$Comment <- "Last entry alpha must be positive."
    obj$Status <- 2
    return(obj)
  }
  if (k > 1 && min(diff(t)) <= 0) {
    warning("t must be increasing.")
    obj$Comment <- "t must be increasing."
    obj$Status <- 2
    return(obj)
  }


    obj$FQ <- FQ
    obj$t <- t
    obj$alpha <- alpha
    obj$truncation <- truncation
    obj$truncation_type <- truncation_type
    obj$dispersion <- dispersion
    obj$Status <- Status
    obj$Comment <- Comment

    return(obj)
}


print.PPP_model <- function(x, ...) {
  if (x$dispersion == 1) {
    fq_dist <- "Poisson"
  } else if (x$dispersion > 1) {
    fq_dist <- "Negative Binomial"
  } else {
    fq_dist <- "Binomial"
  }
  cat("\nPanjer & Piecewise Pareto model\n\n")
  cat("Collective model with a ", fq_dist, " distribution for the claim count and a Piecewise Pareto distributed severity.", sep = "")
  cat("\n\n", fq_dist, " Distribution:", sep = "")
  cat("\nExpected Frequency:   ", x$FQ, sep = "")
  cat("\nDispersion:           ", x$dispersion, sep = "")
  cat("\n\nPiecewise Pareto Distribution:")
  cat("\nThresholds:           ", x$t, sep = "")
  cat("\nAlphas:               ", x$alpha, sep = "")
  if (!is.null(x$truncation)) {
    cat("\nTruncation:           ", x$truncation, sep = "")
    cat("\nTruncation Type:      '", x$truncation_type,"'", sep = "")
  }
}


PPP_Model_Exp_Layer_Loss_s <- function(Cover, Attachment_Point, PPP_Model) {
  return(PPP_Model$FQ * PiecewisePareto_Layer_Mean(Cover, Attachment_Point, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type))
}
PPP_Model_Exp_Layer_Loss <- Vectorize(PPP_Model_Exp_Layer_Loss_s, c("Cover", "Attachment_Point"))


PPP_Model_Layer_Var_s <- function(Cover, Attachment_Point, PPP_Model) {
  E_N <- PPP_Model$FQ
  Var_N <- E_N * PPP_Model$dispersion
  E_X <- PiecewisePareto_Layer_Mean(Cover, Attachment_Point, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type)
  Var_X <- PiecewisePareto_Layer_Var(Cover, Attachment_Point, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type)
  return(E_N * Var_X + Var_N * E_X^2)
}
PPP_Model_Layer_Var <- Vectorize(PPP_Model_Layer_Var_s, c("Cover", "Attachment_Point"))


PPP_Model_Layer_Sd_s <- function(Cover, Attachment_Point, PPP_Model) {
  return(sqrt(PPP_Model_Layer_Var(Cover, Attachment_Point, PPP_Model)))
}
PPP_Model_Layer_Sd <- Vectorize(PPP_Model_Layer_Sd_s, c("Cover", "Attachment_Point"))


PPP_Model_Excess_Frequency <- function(x, PPP_Model) {
  return(PPP_Model$FQ * (1 - pPiecewisePareto(x, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type)))
}


PPP_Model_Simulate <- function(n, PP_Model) {
  #return(PP_Model$FQ * PiecewisePareto_Layer_Mean(Cover, Attachment_Point, PP_Model$t, PP_Model$alpha, truncation = PP_Model$truncation, truncation_type = PP_Model$truncation_type))
}


