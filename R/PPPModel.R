#' PPP_Model (Collective Panjer & Piecewise Pareto Model) Object
#'
#' @description Constructor function for the PPP_Model object
#'
#' @param FQ Numerical. Expected claim count of the collective model.
#' @param t Numeric vector. Vector containing the thresholds of the Piecewise Pareto distribution.
#' @param alpha Numeric vector. Vector containing the alphas of the Piecewise Pareto distribution.
#' @param truncation Numerical. If \code{truncation} is not \code{NULL} and \code{truncation > max(t)}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#' @param dispersion Numerical. Dispersion of the Panjer distribution (i.e. variance to mean ratio).
#' @param Status Numerical indicator if a function returns a PPP_Model object: 0 = success, 1 = some information has been ignored, 2 = no solution found
#' @param Comment Charakter. An optional comment.

#' @examples
#' PPPM <- PPP_Model(2, c(1000,2000), c(1,2), dispersion = 2)
#' PPPM
#'
#' @export

PPP_Model <- function(FQ = NULL, t = NULL, alpha = NULL, truncation = NULL, truncation_type = "lp", dispersion = 1, Status = 0, Comment = "OK") {
  #
  obj <- list(FQ = FQ, t = t, alpha = alpha, truncation = truncation, truncation_type = truncation_type, dispersion = dispersion, Status = Status, Comment = Comment)
  class(obj) <- "PPP_Model"

  if (!is.valid.PPP_Model(obj)) {
    obj <- list(FQ = NULL, t = NULL, alpha = NULL, truncation = NULL, truncation_type = "lp", dispersion = 1, Status = 2, Comment = is.valid.PPP_Model(obj, comment = TRUE))
    class(obj) <- "PPP_Model"
  }
  return(obj)

}

#' Print a PPP_Model Object(Collective Panjer & Piecewise Pareto Model) Object
#'
#' @description Print method for PPP_Model objects
#'
#' @param x PPP_Model object.
#' @param ... Other arguments, all currently ignored.
#'
#' @export

print.PPP_Model <- function(x, ...) {
  if (!is.positive.finite.number(x$dispersion)) {
    fq_dist <- "Panjer"
  } else if (x$dispersion == 1) {
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
  if (is.positive.finite.number(x$dispersion)  && x$dispersion != 1) {
    cat("\nDispersion:           ", x$dispersion, sep = "")
    if (is.positive.finite.number(x$FQ) && x$dispersion > 1) {
      cat(" (i.e. contagion = ", (x$dispersion - 1)/x$FQ, ")", sep = "")
    }
  }
  cat("\n\nPiecewise Pareto Distribution:")
  cat("\nThresholds:      ", x$t, sep = "   ")
  cat("\nAlphas:           ", x$alpha, sep = "   ")
  if (!is.null(x$truncation)) {
    cat("\nTruncation:           ", x$truncation, sep = "")
    cat("\nTruncation Type:      '", x$truncation_type,"'", sep = "")
  } else {
    cat("\nThe distribution is not truncated.")
  }
  cat("\n\nStatus:              ", x$Status)
  cat("\nComments:            ", x$Comment)
  if (!is.valid.PPP_Model(x)) {
    cat("\n\nThe model is not valid.\n")
    cat(is.valid.PPP_Model(x, comment = TRUE))
  }
  cat("\n\n")

}

#' Check if an object is a PPP_Model
#'
#' @description Checks if the class of an object is 'PPP_Model'
#'
#' @param x Object to be checked.

#' @examples
#' PPPM <- PPP_Model(2, c(1000,2000), c(1,2), dispersion = 2)
#' PPPM
#' is.valid.PPP_Model(PPPM)
#'
#' PPPM$alpha <- 2
#' is.valid.PPP_Model(PPPM)
#' is.PPP_Model(PPPM)
#'
#' @export

is.PPP_Model <- function(x) {
  if (class(x) == "PPP_Model") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' Check if an object is a valid PPP_Model
#'
#' @description Checks if an object is a PPP_Model object and whether it is valid for the use in functions like \code{Layer_Mean}
#'
#' @param x Object to be checked.
#' @param comment If FALSE then the function returns a boolean indicating whether \code{x} is a valid PPP_Model. If TRUE then the function returns a comment instead.

#' @examples
#' PPPM <- PPP_Model(2, c(1000,2000), c(1,2), dispersion = 2)
#' PPPM
#' is.valid.PPP_Model(PPPM)
#' is.valid.PPP_Model(PPPM, comment = TRUE)
#'
#' PPPM$alpha <- 2
#' is.valid.PPP_Model(PPPM)
#' is.valid.PPP_Model(PPPM, comment = TRUE)
#'
#' @export

is.valid.PPP_Model <- function(x, comment = FALSE) {
  if (class(x) != "PPP_Model" || typeof(x) != "list") {
    if (!comment) {
      return(FALSE)
    } else {
      return("Object does not have class PPP_Model.")
    }
  }
  required_elements <- c("FQ", "t", "alpha", "truncation", "truncation_type", "dispersion", "Status", "Comment")
  available <- required_elements %in% names(x)
  if (sum(!available) > 0) {
    if (!comment) {
      return(FALSE)
    } else {
      return(paste("Not all required list elements available. Missing elements:", paste(required_elements[!available], collapse = ", ")))
    }
  }


  if (!is.nonnegative.finite.number(x$FQ)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("FQ must be a nonnegative number.")
    }
  }
  if (!valid.parameters.PiecewisePareto(x$t, x$alpha, x$truncation, x$truncation_type)) {
    if (!comment) {
      return(FALSE)
    } else {
      return(valid.parameters.PiecewisePareto(x$t, x$alpha, x$truncation, x$truncation_type, comment = TRUE))
    }
  }
  if (!is.positive.finite.number(x$dispersion)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("dispersion must be a positive number.")
    }
  }


  if (!comment) {
    return(TRUE)
  } else {
    return("OK")
  }
}



PPP_Model_Exp_Layer_Loss_s <- function(Cover, AttachmentPoint, PPP_Model) {
  if (!is.valid.PPP_Model(PPP_Model)) {
    warning(is.valid.PPP_Model(PPP_Model, comment = TRUE))
    return(NaN)
  } else {
    return(PPP_Model$FQ * PiecewisePareto_Layer_Mean(Cover, AttachmentPoint, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type))
  }
}

#' Expected Loss of a Reinsurance Layer
#'
#' @description  Calculates the expected loss of a reinsurance layer for a PPP_Model. This function is deprecated. Use \code{Layer_Mean} instead.
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param PPP_Model PPP_Model object.
#'
#' @return The expected loss of the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{PPP_Model}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Example1_AP
#' Example1_Cov
#' Example1_EL
#' Layer_Mean(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

PPP_Model_Exp_Layer_Loss <- function(Cover, AttachmentPoint, PPP_Model) {
  .Deprecated("Layer_Mean")
  PPP_Model_Exp_Layer_Loss_v(Cover, AttachmentPoint, PPP_Model)
}




PPP_Model_Exp_Layer_Loss_v <- Vectorize(PPP_Model_Exp_Layer_Loss_s, c("Cover", "AttachmentPoint"))



PPP_Model_Layer_Var_s <- function(Cover, AttachmentPoint, PPP_Model) {
  if (!is.valid.PPP_Model(PPP_Model)) {
    warning(is.valid.PPP_Model(PPP_Model, comment = TRUE))
    return(NaN)
  } else {
    E_N <- PPP_Model$FQ
    Var_N <- E_N * PPP_Model$dispersion
    E_X <- PiecewisePareto_Layer_Mean(Cover, AttachmentPoint, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type)
    Var_X <- PiecewisePareto_Layer_Var(Cover, AttachmentPoint, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type)
    return(E_N * Var_X + Var_N * E_X^2)
  }
}

#' Variance of a Reinsurance Layer
#'
#' @description  Calculates the variance of the loss in a reinsurance layer for a PPP_Model. This function is deprecated. Use \code{Layer_Var} instead.
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param PPP_Model PPP_Model object.
#'
#' @return The variance of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{PPP_Model}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Var(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

PPP_Model_Layer_Var <- function(Cover, AttachmentPoint, PPP_Model) {
  .Deprecated("Layer_Var")
  PPP_Model_Layer_Var_v(Cover, AttachmentPoint, PPP_Model)
}

PPP_Model_Layer_Var_v <- Vectorize(PPP_Model_Layer_Var_s, c("Cover", "AttachmentPoint"))





PPP_Model_Layer_Sd_s <- function(Cover, AttachmentPoint, PPP_Model) {
  if (!is.valid.PPP_Model(PPP_Model)) {
    warning(is.valid.PPP_Model(PPP_Model, comment = TRUE))
    return(NaN)
  } else {
    return(sqrt(PPP_Model_Layer_Var_v(Cover, AttachmentPoint, PPP_Model)))
  }
}


#' Standard Deviation of a Reinsurance Layer
#'
#' @description  Calculates the standard deviation of the loss in a reinsurance layer for a PPP_Model
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param PPP_Model PPP_Model object.
#'
#' @return The standard deviation of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{PPP_Model}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Sd(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

PPP_Model_Layer_Sd <- function(Cover, AttachmentPoint, PPP_Model) {
  .Deprecated("Layer_Sd")
  PPP_Model_Layer_Sd_v(Cover, AttachmentPoint, PPP_Model)
}

PPP_Model_Layer_Sd_v <- Vectorize(PPP_Model_Layer_Sd_s, c("Cover", "AttachmentPoint"))






PPP_Model_Excess_Frequency_s <- function(x, PPP_Model) {
  if (!is.valid.PPP_Model(PPP_Model)) {
    warning(is.valid.PPP_Model(PPP_Model, comment = TRUE))
    return(NaN)
  } else if (!is.atomic(x) || !is.numeric(x) || length(x) != 1 || is.na(x)) {
    warning("x must be a number.")
    return(NaN)
  } else {
    return(PPP_Model$FQ * (1 - pPiecewisePareto(x, PPP_Model$t, PPP_Model$alpha, truncation = PPP_Model$truncation, truncation_type = PPP_Model$truncation_type)))
  }
}

#' Expected Frequency in Excess of a Threshold
#'
#' @description  Calculates the expected frequency in excess of a threshold for a PPP_Model
#'
#' @param x Numeric. Threshold.
#' @param PPP_Model PPP_Model object.
#'
#' @return The expected frequency in excess of \code{x} for the given \code{PPP_Model}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Excess_Frequency(PPPM, c(-Inf, 0, 1000, 2000, 3000, Inf))
#'
#' @export

PPP_Model_Excess_Frequency <- function(x, PPP_Model) {
  .Deprecated("Excess_Frequency")
  PPP_Model_Excess_Frequency_v(x, PPP_Model)
}


PPP_Model_Excess_Frequency_v <- Vectorize(PPP_Model_Excess_Frequency_s, c("x"))





#' Simulate Losses with a PPP_Model
#'
#' @description  Simulates losses of a PPP_Model
#'
#' @param n Integer. Number of Simulations.
#' @param PPP_Model PPP_Model object.
#'
#' @return A matrix where row k contains the simulated losses of the kth simulation.
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(c(1000, 2000, 3000), c(2000, 1000, 500),
#'                                            truncation = 10000, truncation_type = "wd")
#' PPPM
#' Simulate_Losses(PPPM, 100)
#'
#' @export

PPP_Model_Simulate <- function(n, PPP_Model) {
  .Deprecated("Simulate_Losses")
  if (!is.valid.PPP_Model(PPP_Model)) {
    warning(is.valid.PPP_Model(PPP_Model, comment = TRUE))
    return(NaN)
  }
  if (!is.positive.finite.number(n)) {
    warning("n must be a positive number.")
    return(NaN)
  } else {
    n <- ceiling(n)
  }
  claim_count <- rPanjer(n, PPP_Model$FQ, PPP_Model$dispersion)
  claims <- rPiecewisePareto(sum(claim_count), PPP_Model$t, PPP_Model$alpha, PPP_Model$truncation, PPP_Model$truncation_type)
  result <- matrix(NaN, nrow = n, ncol = max(claim_count))
  result[col(result) <= claim_count] <- claims
  return(result)
}


