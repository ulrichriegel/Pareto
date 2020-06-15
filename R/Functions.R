
#' Layer Mean of the Pareto Distribution

#' @description  Calculates the expected loss of a Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Numeric. Pareto alpha.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The expected loss of the (truncated) Pareto distribution with parameters \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' Pareto_Layer_Mean(4000, 1000, 2)
#' Pareto_Layer_Mean(4000, 1000, alpha = 2, t = 1000)
#' Pareto_Layer_Mean(4000, 1000, alpha = 2, t = 5000)
#' Pareto_Layer_Mean(4000, 1000, alpha = 2, t = 1000, truncation = 5000)
#' Pareto_Layer_Mean(9000, 1000, alpha = 2, t = 1000, truncation = 5000)
#'
#' @export

Pareto_Layer_Mean <- function(Cover, AttachmentPoint, alpha, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point in zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  if (is.infinite(Cover)) {
    if (alpha <= 1) {
      return(Inf)
    } else if (t <= AttachmentPoint) {
      EP <- -(t / AttachmentPoint)^alpha / (1 - alpha) * AttachmentPoint
    } else {
      EP <- t - AttachmentPoint
      EP <- EP - t / (1 - alpha)
    }
    return(EP)
  } else {
    # Calculation ignoring truncation
    if (t <= AttachmentPoint) {
      if (alpha == 0) {
        EP <- Cover
      } else if (alpha == 1) {
        EP <- t * (log(Cover + AttachmentPoint) - log(AttachmentPoint))
      } else {
        EP <- t / (1 - alpha) * (((Cover + AttachmentPoint) / t)^(1 - alpha) - (AttachmentPoint / t)^(1 - alpha))
      }
    } else if (t >= AttachmentPoint + Cover) {
      EP <- Cover
    } else {
      EP <- t - AttachmentPoint
      if (alpha == 0) {
        EP <- Cover
      } else if (alpha == 1) {
        EP <- EP + t * (log(Cover + AttachmentPoint) - log(t))
      } else {
        EP <- EP + t / (1 - alpha) * (((Cover + AttachmentPoint) / t)^(1 - alpha) - 1)
      }
    }

    if (is.positive.finite.number(truncation)) {
      # then Cover + AttachmentPoint <= truncation
      if (alpha < 1e-6) {
        IndefitineIntegral <- function(x) {
          return(x - 1 / log(truncation / t) * (x * log(x / t) - x))
        }
        if (t <= AttachmentPoint) {
           EP <- IndefitineIntegral(Cover + AttachmentPoint) - IndefitineIntegral(AttachmentPoint)
        } else if (t >= Cover + AttachmentPoint) {
          EP <- Cover
        } else {
          EP <- t - AttachmentPoint + IndefitineIntegral(Cover + AttachmentPoint) - IndefitineIntegral(t)
        }
      } else {
        # Adjustment for truncation
        FQ_at_truncation <- (t / truncation)^alpha
        EP <- (EP - FQ_at_truncation * Cover) / (1 - FQ_at_truncation)
      }
    }
    return(EP)
  }

}





# #' Calculates the second moment of the xs loss of a Pareto distribution in a reinsurance layer
# #'
# #' Not visible to the user. Used in Pareto_Layer_Var.
# #' @param Cover Numeric. Cover of the reinsurance layer. Use Inf for unlimited layers.
# #' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
# #' @param alpha Numeric. Pareto alpha.


Pareto_Layer_Second_Moment_simple <- function(Cover, AttachmentPoint, alpha) {
  if (!is.positive.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a positive number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (Cover == 0) {
    return(0)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }

  if (is.infinite(Cover)) {
    if (alpha <= 2) {
      # warning("alpha must be > 2 for unlimited covers!")
      return(Inf)
    } else {
      SM <- 2 * AttachmentPoint^2 * (1/(alpha-2) - 1/(alpha-1))
    }
    return(SM)
  } else {
    if (alpha == 0) {
      SM <- Cover^2
    } else if (alpha == 1) {
      SM <- 2 * AttachmentPoint^2 * (Cover/AttachmentPoint - log(1+Cover/AttachmentPoint))
    } else if (alpha == 2) {
      SM <- 2 * AttachmentPoint^2 * (-Cover/(Cover+AttachmentPoint) + log(1+Cover/AttachmentPoint))
    } else {
      SM <- 2 * AttachmentPoint^2 * (((1+Cover/AttachmentPoint)^(2-alpha)-1) / (2-alpha) - ((1+Cover/AttachmentPoint)^(1-alpha)-1) / (1-alpha))
    }
    return(SM)
  }

}

#' Layer Variance of the Pareto Distribution
#'
#' @description Calculates the variance of a Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Numeric. Pareto alpha.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The variance of the(truncated) Pareto distribution with parameters \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' Pareto_Layer_Var(4000, 1000, 2)
#' Pareto_Layer_Var(4000, 1000, alpha = 2, t = 1000)
#' Pareto_Layer_Var(4000, 1000, alpha = 2, t = 5000)
#' Pareto_Layer_Var(4000, 1000, alpha = 2, t = 1000, truncation = 5000)
#' Pareto_Layer_Var(9000, 1000, alpha = 2, t = 1000, truncation = 5000)
#'
#' @export

Pareto_Layer_Var <- function(Cover, AttachmentPoint, alpha, t = NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point is zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  ExitPoint <- Cover + AttachmentPoint
  if (t >= ExitPoint) {
    return(0)
  }
  AttachmentPoint <- max(AttachmentPoint, t)
  Cover <- ExitPoint - AttachmentPoint

  if (is.infinite(Cover) && alpha <= 2) {
    return(Inf)
  }

  # Second moment of layer loss if t = AttachmentPoint
  SM <- Pareto_Layer_Second_Moment_simple(Cover, AttachmentPoint, alpha)
  if (is.positive.finite.number(truncation) && alpha > 0) {
    # probability of truncation if t = AttachmentPoint
    p <- 1- pPareto(truncation, AttachmentPoint, alpha)
    # consider truncation in second moment
    SM <- (SM - p * Cover^2) / (1 - p)
  }
  # consider thresholds t < AttachmentPoint
  p <- 1 - pPareto(AttachmentPoint, t, alpha, truncation = truncation)
  SM <- p * SM

  if (is.positive.finite.number(truncation) && alpha == 0) {
    IndefiniteIntegral <- function(x) {
      if (x <= t) {
        return(x^2 + 0.5 * t^2 / log(truncation / t))
      } else {
        return(x^2 - x^2 / log(truncation / t) * log(x / t) + 0.5 * x^2 / log(truncation / t))
      }
    }
    SM <- IndefiniteIntegral(Cover + AttachmentPoint) - IndefiniteIntegral(AttachmentPoint) - 2 * AttachmentPoint * Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, t, truncation = truncation)
  }


  # calculate Variance
  Result <- SM - Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, t, truncation = truncation)^2
  return(Result)
}


#' Second Layer Moment of the Pareto Distribution
#'
#' @description Calculates the second moment of a Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Numeric. Pareto alpha.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The second moment of the (truncated) Pareto distribution with parameters \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' Pareto_Layer_SM(4000, 1000, 2)
#' Pareto_Layer_SM(4000, 1000, alpha = 2, t = 1000)
#' Pareto_Layer_SM(4000, 1000, alpha = 2, t = 5000)
#' Pareto_Layer_SM(4000, 1000, alpha = 2, t = 1000, truncation = 5000)
#' Pareto_Layer_SM(9000, 1000, alpha = 2, t = 1000, truncation = 5000)
#'
#' @export


Pareto_Layer_SM <- function(Cover, AttachmentPoint, alpha, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point is zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  Var <- Pareto_Layer_Var(Cover, AttachmentPoint, alpha, t, truncation = truncation)
  if (is.infinite(Var)) {
    return(Inf)
  }
  Result <- Var + Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, t, truncation = truncation)^2
  return(Result)
}




#' Simulation of the Pareto Distribution
#'
#' @description Generates random deviates of a Pareto distribution
#'
#' @param n Numeric. Number of observations.
#' @param t Numeric vector. Thresholds of the Pareto distributions
#' @param alpha Numeric vector. Pareto alphas of the Pareto distributions.
#' @param truncation NULL or Numeric vector. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation} (resampled Pareto)
#'
#' @return A vector of \code{n} samples from the (truncated) Pareto distribution with parameters \code{t} and \code{alpha}
#'
#' @examples
#' rPareto(100, 1000, 2)
#' rPareto(100, 1000, 2, truncation = 2000)
#' rPareto(100, t = c(1, 10, 100, 1000, 10000), alpha = c(1, 2, 4, 8, 16))
#'
#' @export

rPareto <- function(n, t, alpha, truncation = NULL) {
  if (!is.positive.finite.number(n)) {
    warning("n must be a positive number.")
    return(NaN)
  }
  n <- ceiling(n)
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (n %% length(t) != 0) {
    warning("n is not a multiple of length(t).")
    t <- rep(t, ceiling(n / length(t)))[1:n]
  }
  if (!is.nonnegative.finite.vector(alpha)) {
    warning("alpha must be a non-negative vector.")
    return(NaN)
  }
  if (n %% length(alpha) != 0) {
    warning("n is not a multiple of length(alpha).")
    alpha <- rep(alpha, ceiling(n / length(alpha)))[1:n]
  }

  if (!is.null(truncation)) {
    if (!is.positive.vector(truncation)) {
      warning("truncation must be NULL or a positive vector (elements may be 'Inf').")
      return(NaN)
    }
    if (n %% length(truncation) != 0) {
      warning("n is not a multiple of length(truncation).")
      truncation <- rep(truncation, ceiling(n / length(truncation)))[1:n]
    }
    if (sum(truncation <= t) > 0) {
      warning("truncation must be > t.")
      return(NaN)
    }
  }


  FinvPareto <- function(x,t,alpha) {
    return(t/(1-x)^(1/alpha))
  }
  FinvTruncParetoAlphaZero <- function(x, t, truncation) {
    return(t * (truncation / t)^x)
  }


  u <- 0
  o <- 1

  if (min(alpha > 0)) {
    if (!is.null(truncation)) {
      o <- 1 - (t / truncation)^alpha
    }

    return(FinvPareto(stats::runif(n, u, o),t,alpha))
  } else {
    if (length(alpha) < n) {
      alpha <- rep(alpha, ceiling(n / length(alpha)))[1:n]
    }
    if (length(t) < n) {
      t <- rep(t, ceiling(n / length(t)))[1:n]
    }
    if (!is.null(truncation)) {
      o <- ifelse(alpha == 0, 1, 1 - (t / truncation)^alpha)
      if (length(truncation) < n) {
        truncation <- rep(truncation, ceiling(n / length(truncation)))[1:n]
      }
    } else {
      o <- rep(1, n)
      truncation <- rep(Inf, n)
    }
    u <- rep(0, n)
    index_alpha_pos <- alpha > 0
    index_alpha_zero_untruncated <- alpha == 0 & is.infinite(truncation)
    index_alpha_zero_truncated <- alpha == 0 & !is.infinite(truncation)

    sim_losses <- numeric(n)
    sim_losses[index_alpha_pos] <- FinvPareto(stats::runif(sum(index_alpha_pos), u[index_alpha_pos], o[index_alpha_pos]), t[index_alpha_pos], alpha[index_alpha_pos])
    sim_losses[index_alpha_zero_untruncated] <- Inf
    sim_losses[index_alpha_zero_truncated] <- FinvTruncParetoAlphaZero(stats::runif(sum(index_alpha_zero_truncated)), t[index_alpha_zero_truncated], truncation[index_alpha_zero_truncated])

    return(sim_losses)
  }
}




#' Pareto Extrapolation
#'
#' @description Uses a Pareto distribution to derive the expected loss of a layer from the expected loss of another layer
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param Cover_1 Numeric. Cover of the layer from which we extrapolate. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint_1 Numeric. Attachment point of the layer from which we extrapolate.
#' @param Cover_2 Numeric. Cover of the layer to which we extrapolate. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint_2 Numeric. Attachment point of the layer to which we extrapolate.
#' @param alpha Numeric. Pareto alpha used for the extrapolation.
#' @param ExpLoss_1 Numeric. Expected loss of the layer from which we extrapolate. If \code{NULL} (default) then the function provides only the ratio between the expected losses of the layers.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > AttachmentPoint_1}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The expected loss of the layer \code{Cover_2} xs \code{AttachmentPoint_2} given that  \code{Cover_1} xs \code{AttachmentPoint_1} has expected
#' loss \code{ExpLoss_1} and assuming a (truncated) Pareto distribution with parameters \code{t} and \code{alpha}. If missing then \code{ExpLoss_1 == 1} is assumed.
#'
#' @examples
#' Pareto_Extrapolation(1000, 1000, 2000, 2000, 2, ExpLoss_1 = 100)
#' Pareto_Extrapolation(1000, 1000, 2000, 2000, 2) * 100
#' Pareto_Extrapolation(1000, 1000, 2000, 2000, 2, truncation = 5000, ExpLoss_1 = 100)
#' Pareto_Extrapolation(1000, 1000, 2000, 2000, 2, truncation = 5000) * 100
#'
#' @export

Pareto_Extrapolation <- function(Cover_1, AttachmentPoint_1, Cover_2, AttachmentPoint_2, alpha, ExpLoss_1 = NULL, truncation = NULL) {
  if (!is.positive.number(Cover_1)) {
    warning("Cover_1 must be a positive number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(AttachmentPoint_1)) {
    warning("AttachmentPoint_1 must be a positive number.")
    return(NaN)
  }
  if (!is.nonnegative.number(Cover_2)) {
    warning("Cover_2 must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(AttachmentPoint_2)) {
    warning("AttachmentPoint_2 must be a positive number.")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(ExpLoss_1)) {
    ExpLoss_1 <- 1
  } else {
    if (!is.nonnegative.finite.number(ExpLoss_1)) {
      warning("ExpLoss_1 must be NULL or a non-negative number.")
      return(NaN)
    }
  }

  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be positive ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= AttachmentPoint_1) {
      warning("truncation must be greater than AttachmentPoint_1")
      return(NaN)
    }
  }

  Smaller_AP <- min(AttachmentPoint_1, AttachmentPoint_2)

  if (ExpLoss_1 == 0) {
    # frequency is zero
    return(0)
  }
  if (is.infinite(Pareto_Layer_Mean(Cover_1, AttachmentPoint_1, alpha, t = Smaller_AP, truncation = truncation))) {
    warning("Pareto_Layer_Mean of layer 1 must be finite.")
    return(NaN)
  }

  Result <- ExpLoss_1 * Pareto_Layer_Mean(Cover_2, AttachmentPoint_2, alpha, t = Smaller_AP, truncation = truncation) / Pareto_Layer_Mean(Cover_1, AttachmentPoint_1, alpha, t = Smaller_AP, truncation = truncation)
  return(Result)

}



#' Pareto Alpha Between Two Layers
#'
#' @description Finds the Pareto alpha between two layers
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param Cover_1 Numeric. Cover of the first layer.
#' @param AttachmentPoint_1 Numeric. Attachment point of the first layer.
#' @param ExpLoss_1 Numeric. Expected loss of the first layer.
#' @param Cover_2 Numeric. Cover of the second layer.
#' @param AttachmentPoint_2 Numeric. Attachment point of the second layer.
#' @param ExpLoss_2 Numeric. Expected loss of the second layer.
#' @param max_alpha Numeric. Upper limit for the alpha that is returned.
#' @param tolerance Numeric. Accuracy of the result.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The Pareto alpha between the layer \code{Cover_1} xs \code{AttachmentPoint_1} with expected loss \code{ExpLoss_1}
#' and the layer  \code{Cover_2} xs \code{AttachmentPoint_2} with expected loss \code{ExpLoss_2}
#'
#' @examples
#' Pareto_Find_Alpha_btw_Layers(100, 100, 100, 200, 200, 50)
#' Pareto_Find_Alpha_btw_Layers(100, 100, 100, 200, 200, 50, truncation = 500)
#'
#' @export

Pareto_Find_Alpha_btw_Layers <- function(Cover_1, AttachmentPoint_1, ExpLoss_1, Cover_2, AttachmentPoint_2, ExpLoss_2, max_alpha = 100, tolerance = 1e-10, truncation = NULL) {
  if (!is.positive.number(Cover_1) || !is.positive.number(Cover_2)) {
    warning("Covers must be positive ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(AttachmentPoint_1) || !is.positive.finite.number(AttachmentPoint_2)) {
    warning("Attachment points must be positive.")
    return(NaN)
  }
  if (!is.positive.finite.number(ExpLoss_1)) {
    warning("ExpLoss_1 must be positive.")
    return(NaN)
  }
  if (!is.positive.finite.number(ExpLoss_2)) {
    warning("ExpLoss_2 must be positive.")
    return(NaN)
  }
  min_alpha <- 0
  if (!is.positive.finite.number(max_alpha)) {
    warning("max_alpha must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(tolerance)) {
    warning("tolerance must be a positive number.")
    return(NaN)
  }

  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("tuncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= max(AttachmentPoint_1, AttachmentPoint_2)) {
      warning("tuncation must be NULL or greater than both attachment points!")
      return(NaN)
    } else {
      Cover_1 <- min(truncation - AttachmentPoint_1, Cover_1)
      Cover_2 <- min(truncation - AttachmentPoint_2, Cover_2)
      min_alpha <- tolerance
    }
  }

  f <- function(alpha) {
    Pareto_Extrapolation(Cover_1, AttachmentPoint_1, Cover_2, AttachmentPoint_2, alpha, ExpLoss_1, truncation = truncation) - ExpLoss_2
  }

  Result <- NA
  if (!is.infinite(Cover_1) && !is.infinite(Cover_2)) {
    try(Result <- stats::uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = TRUE)
    if (is.na(Result)) {
      if (AttachmentPoint_1 > AttachmentPoint_2 && AttachmentPoint_1 + Cover_1 >= AttachmentPoint_2 + Cover_2 && f(max_alpha) < 0) {
        Result <- max_alpha
      } else if (AttachmentPoint_1 < AttachmentPoint_2 && AttachmentPoint_1 + Cover_1 <= AttachmentPoint_2 + Cover_2 && f(max_alpha) > 0) {
        Result <- max_alpha
      }
    }
  } else {
    try(Result <- stats::uniroot(f, c(1 + tolerance, max_alpha), tol = tolerance)$root, silent = TRUE)
    if (is.na(Result)) {
      if (Cover_1 == Inf && Cover_2 == Inf) {
        if (AttachmentPoint_1 < AttachmentPoint_2 && f(max_alpha) > 0) {
          Result <- max_alpha
        }
        if (AttachmentPoint_1 > AttachmentPoint_2 && f(max_alpha) < 0) {
          Result <- max_alpha
        }
      } else if (Cover_1 == Inf) {
        if (AttachmentPoint_1 > AttachmentPoint_2 && f(max_alpha) < 0) {
          Result <- max_alpha
        }
      } else if (Cover_2 == Inf) {
        if (AttachmentPoint_1 < AttachmentPoint_2 && f(max_alpha) > 0) {
          Result <- max_alpha
        }
      }
    }
  }
  if (is.na(Result)) {
    warning("Did not find a solution!")
    return(NaN)
  }
  return(Result)

}


#' Pareto Alpha Between a Frequency and a Layer
#'
#' @description Finds the Pareto alpha between an excess frequency and the expected loss of a layer
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param Threshold Numeric. Threshold
#' @param Frequency Numeric. Expected frequency in excess of \code{Thershold}
#' @param Cover Numeric. Cover of the second layer.
#' @param AttachmentPoint Numeric. Attachment point of the layer.
#' @param ExpLoss Numeric. Expected loss of the layer.
#' @param max_alpha Numeric. Upper limit for the alpha that is returned.
#' @param tolerance Numeric. Accuracy of the result.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The Pareto alpha between the expected number of claims \code{Frequency} excess \code{Threshold}
#' and the layer  \code{Cover} xs \code{AttachmentPoint} with expected loss \code{ExpLoss}
#'
#' @examples
#' Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 1000, 1000, 500)
#' Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 1000, 1000, 500, truncation = 5000)
#'
#' @export

Pareto_Find_Alpha_btw_FQ_Layer <- function(Threshold, Frequency, Cover, AttachmentPoint, ExpLoss, max_alpha = 100, tolerance = 1e-10, truncation = NULL) {
  if (!is.positive.finite.number(Threshold)) {
    warning("Threshold must be positive.")
    return(NaN)
  }
  if (!is.positive.finite.number(Frequency)) {
    warning("Frequency must be a positive number.")
    return(NaN)
  }
  if (!is.positive.number(Cover)) {
    warning("Cover must be positive ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be positive.")
    return(NaN)
  }
  if (!is.positive.finite.number(ExpLoss)) {
    warning("ExpLoss must be positive.")
    return(NaN)
  }
  if (Threshold > AttachmentPoint && Threshold < AttachmentPoint + Cover) {
    warning("Threshold must be <= AttachmentPoint or >= Cover + AttachmentPoint")
    return(NaN)
  }
  if (!is.positive.finite.number(max_alpha)) {
    warning("max_alpha must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(tolerance)) {
    warning("tolerance must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= AttachmentPoint || truncation <= Threshold) {
      warning("Threshold and AttachmentPoint must be < truncation")
      return(NaN)
    }
    Cover <- min(Cover, truncation - AttachmentPoint)
  }

  f <- function(alpha) {
    if (AttachmentPoint < Threshold) {
      FQ_Factor <- 1 / (1 - pPareto(Threshold, AttachmentPoint, alpha, truncation = truncation))
    } else {
      FQ_Factor <- 1 - pPareto(AttachmentPoint, Threshold, alpha, truncation = truncation)
    }

    Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, truncation = truncation) * FQ_Factor * Frequency - ExpLoss
  }

  Result <- NA
  if (is.infinite(Cover)) {
    min_alpha <- 1 + tolerance
  } else if (!is.null(truncation)) {
    min_alpha <- tolerance
  } else {
    min_alpha <- 0
  }
  while (is.infinite(f(max_alpha))) {
    max_alpha <- max_alpha / 2
  }
  try(Result <- stats::uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = TRUE)
  if (is.na(Result)) {
    if (AttachmentPoint >= Threshold && f(max_alpha) > 0) {
      Result <- max_alpha
    } else if (AttachmentPoint + Cover <= Threshold && f(max_alpha) < 0) {
      Result <- max_alpha
    }
  }

  if (is.na(Result)) {
    warning("Did not find a solution!")
    return(NaN)
  }
  return(Result)

}


#' Pareto Alpha Between Two Frequencies
#'
#' @description Finds the Pareto alpha between two excess frequencies
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param Threshold_1 Numeric. Threshold 1
#' @param Frequency_1 Numeric. Expected frequency in excess of \code{Threshold_1}
#' @param Threshold_2 Numeric. Threshold 2
#' @param Frequency_2 Numeric. Expected frequency in excess of \code{Threshold_2}
#' @param max_alpha Numeric. Upper limit for the alpha that is returned.
#' @param tolerance Numeric. Accuracy of the result.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The Pareto alpha between the expected number of claims \code{Frequency_1} excess \code{Threshold_1}
#' and the expected number of claims \code{Frequency_2} excess \code{Threshold_2}
#'
#' @examples
#' Pareto_Find_Alpha_btw_FQs(1000, 1, 2000, 0.5)
#' Pareto_Find_Alpha_btw_FQs(1000, 1, 2000, 0.5, truncation = 5000)
#'
#' @export

Pareto_Find_Alpha_btw_FQs <- function(Threshold_1, Frequency_1, Threshold_2, Frequency_2, max_alpha = 100, tolerance = 1e-10, truncation = NULL) {
  if (!is.positive.finite.number(Threshold_1)) {
    warning("Threshold_1 must be positive!")
    return(NaN)
  }
  if (!is.positive.finite.number(Frequency_1)) {
    warning("Frequency_1 must be positive!")
    return(NaN)
  }
  if (!is.positive.finite.number(Threshold_2)) {
    warning("Threshold_2 must be positive!")
    return(NaN)
  }
  if (!is.positive.finite.number(Frequency_2)) {
    warning("Frequency_2 must be positive!")
    return(NaN)
  }
  if (!is.positive.finite.number(max_alpha)) {
    warning("max_alpha must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(tolerance)) {
    warning("tolerance must be a positive number.")
    return(NaN)
  }

  if (Threshold_1 > Threshold_2) {
    temp1 <- Threshold_1
    temp2 <- Frequency_1
    Threshold_1 <- Threshold_2
    Frequency_1 <- Frequency_2
    Threshold_2 <- temp1
    Frequency_2 <- temp2
  }
  if (Threshold_2 == Threshold_1) {
    warning("Thresholds must not be equal")
    return(NaN)
  }
  if (Frequency_2 > Frequency_1) {
    warning("Frequency of larger threshold must be less than or equal to frequency at lower threshold")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= Threshold_2) {
      warning("Thresholds must be less than truncation")
      return(NaN)
    }
  }

  if (is.null(truncation) || is.infinite(truncation)) {
    alpha <- log(Frequency_1 / Frequency_2) / log(Threshold_2 / Threshold_1)
  } else {
    f <- function(alpha) {
      1 - pPareto(Threshold_2, Threshold_1, alpha, truncation = truncation) - Frequency_2 / Frequency_1
    }

    Result <- NA
    min_alpha <- tolerance
    if (f(max_alpha) > 0) {
      Result = max_alpha
    } else {
      try(Result <- stats::uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = TRUE)
    }
    if (is.na(Result)) {
      warning("Did not find a solution!")
      return(NaN)
    }
    alpha <- Result
  }
  return(alpha)

}



#' Layer Mean of the Piecewise Pareto Distribution
#'
#' @description Calculates the expected loss of a piecewise Pareto distribution in a reinsurance layer
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param Cover Numeric. Cover of the reinsurance layer.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return The expected loss of the (truncated) piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' PiecewisePareto_Layer_Mean(4000, 1000, t, alpha)
#' PiecewisePareto_Layer_Mean(4000, 1000, t, alpha, truncation = 5000)
#' PiecewisePareto_Layer_Mean(4000, 1000, t, alpha, truncation = 5000, truncation_type = "lp")
#' PiecewisePareto_Layer_Mean(4000, 1000, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

PiecewisePareto_Layer_Mean <- function(Cover, AttachmentPoint, t, alpha, truncation = NULL, truncation_type = "lp") {
  if(!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }
  n <- length(t)
  if (!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }

  if (n == 1) {
    return(Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, t, truncation))
  }
  if (!is.null(truncation)) {
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    Cover <- min(Cover, truncation - AttachmentPoint)
    if (is.infinite(truncation)) {
      truncation <- NULL
    }
  }
  if (Cover == 0) {
    return(0)
  }
  if (is.infinite(Cover) && alpha[n] <= 1) {
    return(Inf)
  }

  factors_t <- t[2:n] / t[1:(n-1)]
  excess_prob <- numeric(n)
  excess_prob[1] <- 1
  excess_prob[2:n] <- cumprod((1/factors_t)^alpha[1:(n-1)])

  k1 <- sum(t <= AttachmentPoint)
  if (Cover == Inf) {
    k2 <- n
  } else {
    k2 <- sum(t < AttachmentPoint + Cover)
  }

  Cover_orig <- Cover
  if (k1 == 0 && k2 == 0) {
    return(Cover)
  } else if (k1 == 0) {
    Result <- t[1] - AttachmentPoint
    AttachmentPoint <- t[1]
    if (is.numeric(Cover)) {
      Cover <- Cover - Result
    }
    k1 <- 1
  } else {
    Result <- 0
  }

  Att <- numeric(n)
  Exit <- numeric(n)

  Att[k1:k2] <- pmax(t[k1:k2], AttachmentPoint)
  if (k2 < n) {
    Exit[k1:k2] <- pmin(t[(k1+1):(k2+1)], AttachmentPoint + Cover)
  } else if (k1 < k2) {
    if (is.numeric(Cover)) {
      Exit[k1:k2] <- c(t[(k1+1):k2], AttachmentPoint + Cover)
    } else {
      Exit[k1:k2] <- c(t[(k1+1):k2], 0)
    }
  } else {
    if (is.numeric(Cover)) {
      Exit[k1] <- AttachmentPoint + Cover
    } else {
      Exit[k1] <- 0
    }
  }

  if (!is.infinite(Cover)) {
    for (i in k1:k2) {
      if (!is.null(truncation) && truncation_type == "lp" && i == n) {
        Result <- Result + Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = truncation) * excess_prob[i]
      } else {
        Result <- Result + Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
      }
    }
  } else if (is.infinite(Cover)) {
    if (k1 < k2) {
      for (i in k1:(k2-1)) {
        Result <- Result + Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
      }
    }
    Result <- Result + Pareto_Layer_Mean(Inf, Att[n], alpha[n], t[n]) * excess_prob[n]
  }
  if (!is.null(truncation) && truncation_type == "wd") {
    p <- (1 - pPareto(truncation, t[n], alpha[n])) * excess_prob[n]
    Result <- (Result - p * Cover_orig) / (1 - p)
  }

  return(Result)
}



#' Second Layer Moment of the Piecewise Pareto Distribution
#'
#' @description Calculates the second moment of a piecewise Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return The second moment of the (truncated) piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' PiecewisePareto_Layer_SM(4000, 1000, t, alpha)
#' PiecewisePareto_Layer_SM(4000, 1000, t, alpha, truncation = 5000)
#' PiecewisePareto_Layer_SM(4000, 1000, t, alpha, truncation = 5000, truncation_type = "lp")
#' PiecewisePareto_Layer_SM(4000, 1000, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

PiecewisePareto_Layer_SM <- function(Cover, AttachmentPoint, t, alpha, truncation = NULL, truncation_type = "lp") {
  if(!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }
  n <- length(t)

  if (!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }


  if (n == 1) {
    return(Pareto_Layer_SM(Cover, AttachmentPoint, alpha, t, truncation))
  }

  if (!is.null(truncation)) {
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    Cover <- min(Cover, truncation - AttachmentPoint)
    if (is.infinite(truncation)) {
      truncation <- NULL
    }
  }
  if (Cover == 0) {
    return(0)
  }
  if (is.infinite(Cover) && alpha[n] <= 2) {
    return(Inf)
  }

  factors_t <- t[2:n] / t[1:(n-1)]
  excess_prob <- numeric(n)
  excess_prob[1] <- 1
  excess_prob[2:n] <- cumprod((1/factors_t)^alpha[1:(n-1)])

  prob <- c(- diff(excess_prob), excess_prob[n])

  k1 <- sum(t <= AttachmentPoint)
  if (Cover == Inf) {
    k2 <- n
  } else {
    k2 <- sum(t < AttachmentPoint + Cover)
  }

  AttachmentPoint_orig <- AttachmentPoint
  Cover_orig <- Cover

  if (k1 == 0 && k2 == 0) {
    return(Cover^2)
  } else if (k1 == 0) {
    AttachmentPoint <- t[1]
    if (is.numeric(Cover)) {
      Cover <- Cover - (AttachmentPoint - AttachmentPoint_orig)
    }
    k1 <- 1
  }

  Result <- 0


  Att <- numeric(n)
  Exit <- numeric(n)

  Att[k1:k2] <- pmax(t[k1:k2], AttachmentPoint)
  if (k2 < n) {
    Exit[k1:k2] <- pmin(t[(k1+1):(k2+1)], AttachmentPoint + Cover)
  } else if (k1 < k2) {
    Exit[k1:k2] <- c(t[(k1+1):k2], AttachmentPoint + Cover)
  } else {
    Exit[k1] <- AttachmentPoint + Cover
  }

  for (i in k1:k2) {
    if (!is.null(truncation) && truncation_type == "lp" && i == n) {
      Result <- Result + Pareto_Layer_SM(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = truncation) * prob[i]
      Result <- Result + (Att[i] - AttachmentPoint_orig)^2 * prob[i] * (1 - pPareto(Att[n], t[n], alpha[n]))
      Result <- Result + 2 * (Att[i] - AttachmentPoint_orig) * Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = truncation) * prob[i]
    } else if (i == k2) {
      Result <- Result + Pareto_Layer_SM(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
      Result <- Result + (Att[i] - AttachmentPoint_orig)^2 * excess_prob[i] * (t[i] / Att[i])^alpha[i]
      Result <- Result + 2 * (Att[i] - AttachmentPoint_orig) * Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
    } else {
      Result <- Result + Pareto_Layer_SM(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = Exit[i]) * prob[i]
      Result <- Result + (Att[i] - AttachmentPoint_orig)^2 * prob[i] * (t[i] / Att[i])^alpha[i]
      Result <- Result + 2 * (Att[i] - AttachmentPoint_orig) * Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = Exit[i]) * prob[i]
    }
  }
  if (!is.null(truncation) && truncation_type == "wd") {
    p <- (1 - pPareto(truncation, t[n], alpha[n])) * excess_prob[n]
    Result <- (Result - p * Cover_orig^2) / (1 - p)
  }


  return(Result)
}

#' Layer Variance of the Piecewise Pareto Distribution
#'
#' @description Calculate the variance of a piecewise Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return The variance of the (truncated) piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' PiecewisePareto_Layer_Var(4000, 1000, t, alpha)
#' PiecewisePareto_Layer_SM(4000, 1000, t, alpha) - PiecewisePareto_Layer_Mean(4000, 1000, t, alpha)^2
#' PiecewisePareto_Layer_Var(4000, 1000, t, alpha, truncation = 5000)
#' PiecewisePareto_Layer_Var(4000, 1000, t, alpha, truncation = 5000, truncation_type = "lp")
#' PiecewisePareto_Layer_Var(4000, 1000, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

PiecewisePareto_Layer_Var <- function(Cover, AttachmentPoint, t, alpha, truncation = NULL, truncation_type = "lp") {

  if(!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }
  n <- length(t)

  if (!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }

  SM <- PiecewisePareto_Layer_SM(Cover, AttachmentPoint, t, alpha, truncation, truncation_type)
  if (is.infinite(SM)) {
    return(Inf)
  }
  Result <- SM - PiecewisePareto_Layer_Mean(Cover, AttachmentPoint, t, alpha, truncation, truncation_type)^2
  return(Result)
}



#' Simulation of the Piecewise Pareto Distribution
#'
#' @description Generates random deviates of a piecewise Pareto distribution
#'
#' @param n Numeric. Number of simulations
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#' @param scale_pieces Numeric vector. If not \code{NULL} then the density of the i-th Pareto piece (on the interval (\code{t[i], t[i+1])}) is scaled with the factor \code{const * scale_pieces[i]} (where \code{const} is a normalization constant)
#'
#' @return A vector of \code{n} samples from the (truncated) piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' rPiecewisePareto(100, t, alpha)
#' rPiecewisePareto(100, t, alpha, truncation = 5000)
#' rPiecewisePareto(100, t, alpha, truncation = 5000, truncation_type = "lp")
#' rPiecewisePareto(100, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export


rPiecewisePareto <- function(n, t, alpha, truncation = NULL, truncation_type = "lp", scale_pieces = NULL) {

  if(!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }
  if (!is.null(scale_pieces) && !is.TRUEorFALSE(scale_pieces)) {
    warning("scale_pieces must be NULL or TRUE or FALSE.")
    return(NaN)
  }

  k <- length(t)

  if (k == 1) {
    return(rPareto(n, t, alpha, truncation))
  }
  if (!is.null(truncation) && !is.null(scale_pieces)) {
    warning("either truncation or scale_pieces must be NULL")
    return(NaN)
  }
  if (!is.null(scale_pieces)) {
    if (!is.nonnegative.finite.vector(scale_pieces)) {
      warning("scale_pieces must be NULL or a positive number.")
      return(NaN)
    }
    if (length(scale_pieces) != length(t)) {
      warning("t and scale_pieces must have the same length.")
      return(NaN)
    }
    if (sum(scale_pieces) <= 0) {
      warning("scale_pieces must have a positive entry")
      return(NaN)
    }
  }

  factors_t <- t[2:k] / t[1:(k-1)]
  excess_prob <- numeric(k)
  excess_prob[1] <- 1
  excess_prob[2:k] <- cumprod((1/factors_t)^alpha[1:(k-1)])
  prob_for_pieces <- c(-diff(excess_prob), excess_prob[k])

  if (!is.null(truncation) && truncation_type == "wd") {
    prob_for_pieces[k] <- excess_prob[k] - excess_prob[k] * (t[k] / truncation)^alpha[k]
    prob_for_pieces <- prob_for_pieces / sum(prob_for_pieces)
  }
  if (!is.null(scale_pieces)) {
    prob_for_pieces <- prob_for_pieces * scale_pieces
    prob_for_pieces <- prob_for_pieces / sum(prob_for_pieces)
  }

  Simulated_Pieces <- c(sample(1:k, n, replace = T, prob = prob_for_pieces))
  NumberOfSimulations_for_Pieces <- numeric(k)
  for (i in 1:k) {
    NumberOfSimulations_for_Pieces[i] <- sum(Simulated_Pieces == i)
  }

  FinvPareto <- function(x,t,alpha) {
    return(t/(1-x)^(1/alpha))
  }

  Result <- numeric(n)

  for (i in 1:k) {
    if (i == k) {
      if (!is.null(truncation)) {
        CDF <- 1 - (t[k] / truncation)^alpha[i]
      } else {
        CDF <- 1
      }
    } else {
      CDF <- 1 - (t[i] / t[i+1])^alpha[i]
    }
    Result[Simulated_Pieces == i] <- FinvPareto(stats::runif(NumberOfSimulations_for_Pieces[i], min = 0, max = CDF), t[i], alpha[i])
  }

  return(Result)
}




#' Match a Tower of Expected Layers Losses
#'
#' @description Matches the expected losses of a tower of reinsurance layers using a piecewise Pareto severity
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param Attachment_Points Numeric vector. Vector containing the attachment points of consecutive layers in increasing order
#' @param Expected_Layer_Losses Numeric vector. Vector containing the expected losses of layers xs the attachment points.
#' @param Unlimited_Layers Logical. If \code{TRUE}, then \code{Expected_Layer_Losses[i]} contains the expected loss of \code{Inf} xs \code{Attachment_Points[i]}. If \code{FALSE} then \code{Expected_Layer_Losses[i]} contains the expected loss of the layer \code{Attachment_Points[i+1]} xs \code{Attachment_Points[i]}
#' @param Frequencies Numeric vector. Expected frequencies excess the attachment points. The vector may contain NAs. If \code{NULL} then the function calculates frequencies.
#' @param FQ_at_lowest_AttPt Numerical. Expected frequency excess \code{Attachment_Points[1]}. Overrules first entry in Frequencies.
#' @param FQ_at_highest_AttPt Numerical. Expected frequency excess \code{Attachment_Points[k]}. Overrules last entry in Frequencies.
#' @param TotalLoss_Frequencies Numeric vector. \code{TotalLoss_Frequencies[i]} is the frequency of total losses to layer \code{i} (i.e. \code{Attachment_Points[i+1] - Attachment_Points[i]} xs \code{Attachment_Points[i]}).    \code{TotalLoss_Frequencies[i]} is the frequency for losses larger than or equal to \code{Attachment_Points[i+1]}, whereas \code{Frequencies[i]} is the frequency of losses larger than \code{Attachment_Points[i]}.    \code{TotalLoss_Frequencies[i] > Frequencies[i+1]} means that there is a point mass of the severity at \code{Attachment_Points[i+1]}.
#' @param minimize_ratios Logical. If \code{TRUE} then ratios between alphas are minimized.
#' @param Use_unlimited_Layer_for_FQ Logical. Only relevant if no frequency is provided for the highest attachment point by the user. If \code{TRUE} then the frequency is calculated using the Pareto alpha between the last two layers.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#' @param dispersion Numerical. Dispersion of the claim count distribution in the resulting PPP_Model.
#' @param tolerance Numeric. Numerical tolerance.
#' @param alpha_max Numerical. Maximum alpha to be used for the matching.
#' @param merge_tolerance Numerical. Consecutive Pareto pieces are merged if the alphas deviate by less than merge_tolerance.
#' @param RoL_tolerance Numerical. Consecutive layers are merged if RoL decreases less than factor \code{1 - RoL_tolerance}.

#' @return A PPP_Model object that contains the information about a collective model with a Panjer distributed claim count and a Piecewise Pareto distributed severity. The object contains the following elements: \itemize{
#' \item \code{FQ} Numerical. Frequency in excess of the lowest threshold of the piecewise Pareto distribution
#' \item \code{t} Numeric vector. Vector containing the thresholds for the piecewise Pareto distribution
#' \item \code{alpha} Numeric vector. Vector containing the Pareto alphas of the piecewise Pareto distribution
#' \item \code{truncation} Numerical. If \code{truncation} is not \code{NULL} and \code{truncation > max(t)}, then the distribution is truncated at \code{truncation}.
#' \item \code{truncation_type} Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#' \item \code{dispersion} Numerical. Dispersion of the Panjer distribution (i.e. variance to mean ratio).
#' \item \code{Status} Numerical indicator: 0 = success, 1 = some information has been ignored, 2 = no solution found
#' \item \code{Comment} Character. Information on whether the fit was successful
#' }
#'
#' @examples
#' AP <- Example1_AP
#' EL <- Example1_EL
#' PiecewisePareto_Match_Layer_Losses(AP, EL)
#' EL_unlimited <- rev(cumsum(rev(Example1_EL)))
#' PiecewisePareto_Match_Layer_Losses(AP, EL_unlimited, Unlimited_Layers = TRUE)
#' PiecewisePareto_Match_Layer_Losses(AP, EL, FQ_at_lowest_AttPt = 0.5)
#' Example1_FQ <- c(0.3, 0.15, 0.08, 0.02, 0.005)
#' PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = Example1_FQ)
#'
#' @export

PiecewisePareto_Match_Layer_Losses <- function(Attachment_Points, Expected_Layer_Losses, Unlimited_Layers = FALSE, Frequencies = NULL, FQ_at_lowest_AttPt = NULL, FQ_at_highest_AttPt = NULL, TotalLoss_Frequencies = NULL, minimize_ratios = TRUE, Use_unlimited_Layer_for_FQ = TRUE, truncation = NULL, truncation_type = "lp", dispersion = 1, tolerance = 1e-10, alpha_max = 100, merge_tolerance = 1e-6, RoL_tolerance = 1e-6) {

  Results <- PPP_Model()
  Results$Comment <- ""
  if (!is.positive.finite.number(dispersion)) {
    warning("Dispersion must be a positive number.")
    Results$Comment <- "Dispersion must be a positive number."
    Results$Status <- 2
    return(Results)
  } else {
    Results$dispersion <- dispersion
  }


  if (!is.positive.finite.vector(Attachment_Points)) {
    warning("Attachment_Points must be a vector with positive entries..")
    Results$Comment <- "Attachment_Points must be a vector with positive entries."
    Results$Status <- 2
    return(Results)
  }
  k <-length(Attachment_Points)
  if (!is.positive.finite.vector(Expected_Layer_Losses)) {
    warning("Expected_Layer_Losses must be a vector with positive entries.")
    Results$Comment <- "Expected_Layer_Losses must be a vector with positive entries."
    Results$Status <- 2
    return(Results)
  }
  if (length(Expected_Layer_Losses) != k) {
    warning("Attachment_Points and Expected_Layer_Losses must have the same length.")
    Results$Comment <- "Attachment_Points and Expected_Layer_Losses must have the same length."
    Results$Status <- 2
    return(Results)
  }
  if (k > 1 && min(diff(Attachment_Points)) <= 0) {
    warning("Attachment_Points must be increasing.")
    Results$Comment <- "Attachment_Points must be increasing."
    Results$Status <- 2
    return(Results)
  }
  if (!is.TRUEorFALSE(Unlimited_Layers)) {
    warning("Unlimited_Layers must be TRUE or FALSE.")
    Results$Comment <- "Unlimited_Layers must be TRUE or FALSE."
    Results$Status <- 2
    return(Results)
  }
  if (!is.null(Frequencies)) {
    if (!is.positive_or_NA.finite.vector(Frequencies)) {
      warning("Frequencies must be NULL or a vector with positive entries (NAs allowed).")
      Results$Comment <- "Frequencies must be NULL or a vector with positive entries (NAs allowed)."
      Results$Status <- 2
      return(Results)
    }
    if (length(Frequencies) != k) {
      warning("Attachment_Points and Frequencies must have the same length.")
      Results$Comment <- "Attachment_Points and Frequencies must have the same length."
      Results$Status <- 2
      return(Results)
    }
  }
  if (!is.null(TotalLoss_Frequencies)) {
    if (!is.positive.finite.vector(TotalLoss_Frequencies)) {
      warning("TotalLoss_Frequencies must be NULL or a vector with positive entries.")
      Results$Comment <- paste0(Results$Comment, "TotalLoss_Frequencies must be NULL or a vector with positive entries. TotalLoss_Frequencies are ignored! ")
      Results$Status <- 1
      TotalLoss_Frequencies <- NULL
    } else if (length(TotalLoss_Frequencies) != (k-1)) {
      warning("TotalLoss_Frequencies must have length of Frequencies - 1.")
      Results$Comment <- paste0(Results$Comment, "TotalLoss_Frequencies must have length of Frequencies - 1. TotalLoss_Frequencies are ignored! ")
      Results$Status <- 1
      TotalLoss_Frequencies <- NULL
    } else if (is.null(Frequencies)) {
      warning("TotalLoss_Frequencies must be NULL if Frequencies is NULL.")
      Results$Comment <- paste0(Results$Comment, "TotalLoss_Frequencies must be NULL if Frequencies is NULL. TotalLoss_Frequencies are ignored! ")
      Results$Status <- 1
      TotalLoss_Frequencies <- NULL
    } else if (!is.positive.finite.vector(Frequencies)) {
      warning("Frequencies must be positive if TotalLoss_Frequencies are used.")
      Results$Comment <- paste0(Results$Comment, "Frequencies must be positive if TotalLoss_Frequencies are used (no NAs allowed). TotalLoss_Frequencies are ignored! ")
      Results$Status <- 1
      TotalLoss_Frequencies <- NULL
    }
  }
  if (!is.null(FQ_at_lowest_AttPt)) {
    if (!is.positive.finite.number(FQ_at_lowest_AttPt)) {
      warning("FQ_at_lowest_AttPt must be NULL or a positive number. FQ_at_lowest_AttPt is ignored!")
      Results$Comment <- paste0(Results$Comment, "FQ_at_lowest_AttPt must be NULL or a positive number. FQ_at_lowest_AttPt is ignored! ")
      Results$Status <- 1
      FQ_at_lowest_AttPt <- NULL
    }
  }
  if (!is.null(FQ_at_highest_AttPt)) {
    if (!is.positive.finite.number(FQ_at_highest_AttPt)) {
      warning("FQ_at_highest_AttPt must be NULL or a positive number.")
      Results$Comment <- paste0(Results$Comment, "FQ_at_highest_AttPt must be NULL or a positive number. FQ_at_highest_AttPt is ignored! ")
      Results$Status <- 1
      FQ_at_highest_AttPt <- NULL
    }
  }
  if (!is.TRUEorFALSE(minimize_ratios)) {
    minimize_ratios <- TRUE
    warning("minimize_ratios must be TRUE or FALSE. TRUE used.")
    Results$Comment <- paste0(Results$Comment, "minimize_ratios must be TRUE or FALSE. TRUE used. ")
    Results$Status <- 1
  }
  if (!is.TRUEorFALSE(Use_unlimited_Layer_for_FQ)) {
    Use_unlimited_Layer_for_FQ <- TRUE
    warning("Use_unlimited_Layer_for_FQ must be TRUE or FALSE. TRUE used.")
    Results$Comment <- paste0(Results$Comment, "Use_unlimited_Layer_for_FQ must be TRUE or FALSE. TRUE used. ")
    Results$Status <- 1
  }
  if (!is.positive.finite.number(tolerance)) {
    warning("tolerance must be a positive number.")
    Results$Comment <- paste0(Results$Comment, "tolerance must be a positive number.")
    Results$Status <- 2
    return(Results)
  }
  if (!is.positive.finite.number(merge_tolerance)) {
    warning("merge_tolerance must be a positive number.")
    Results$Comment <- paste0(Results$Comment, "merge_tolerance must be a positive number.")
    Results$Status <- 2
    return(Results)
  }
  if (!is.positive.finite.number(RoL_tolerance)) {
    warning("RoL_tolerance must be a positive number.")
    Results$Comment <- paste0(Results$Comment, "RoL_tolerance must be a positive number.")
    Results$Status <- 2
    return(Results)
  }
  if (!is.positive.finite.number(alpha_max)) {
    warning("alpha_max must be a positive number.")
    Results$Comment <- paste0(Results$Comment, "alpha_max must be a positive number.")
    Results$Status <- 2
    return(Results)
  }
  if (!is.atomic(truncation_type) || length(truncation_type) != 1 || !(truncation_type %in% c("lp", "wd"))) {
    warning("truncation_type must be 'lp' or 'wd'.")
    Results$Comment <- paste0(Results$Comment, "truncation_type must be 'lp' or 'wd'.")
    Results$Status <- 2
    return(Results)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      Results$Comment <- paste0(Results$Comment, "truncation must be NULL or a positive number ('Inf' allowed).")
      Results$Status <- 2
      return(Results)
    }
    if (truncation <= max(Attachment_Points)) {
      warning("truncation must be greater than max(Attachment_Points).")
      Results$Comment <- paste0(Results$Comment, "truncation must be greater than max(Attachment_Points).")
      Results$Status <- 2
      return(Results)
    }
    if (is.infinite(truncation)) {
      truncation <- NULL
    }
  }

  if (!is.null(truncation)) {
    Results$truncation <- truncation
    Results$truncation_type <- truncation_type
  }

  if (k == 1) {
    Results$t <- Attachment_Points
    fq <- NULL
    if (!is.null(Frequencies) && !is.na(Frequencies)) {
      fq <- Frequencies
    }
    if (!is.null(FQ_at_highest_AttPt)) {
      fq <- FQ_at_highest_AttPt
    }
    if (!is.null(FQ_at_lowest_AttPt)) {
      fq <- FQ_at_lowest_AttPt
    }
    if (is.null(fq)) {
      alpha <- 2
      Results$alpha <- alpha
      Results$FQ <- Expected_Layer_Losses / Pareto_Layer_Mean(Inf, Attachment_Points, alpha, truncation = truncation)
    } else {
      alpha <- NA
      try(alpha <- Pareto_Find_Alpha_btw_FQ_Layer(Attachment_Points, fq, Inf, Attachment_Points, Expected_Layer_Losses, truncation = truncation), silent = TRUE)
      if (is.na(alpha)) {
        warning("truncation too low.")
        Results$t <- NULL
        Results$alpha <- NULL
        Results$Comment <- paste0(Results$Comment, "truncation too low.")
        Results$Status <- 2
        return(Results)
      } else {
        Results$alpha <- alpha
        Results$FQ <- fq
      }
    }
    if (Results$Comment == "") {Results$Comment <- "OK"}
    return(Results)
  }

  if (Unlimited_Layers) {
    ELL <- Expected_Layer_Losses[1:(k-1)] - Expected_Layer_Losses[2:k]
    ELL <- c(ELL, Expected_Layer_Losses[k])
  } else {
    ELL <- Expected_Layer_Losses
  }
  if (!is.positive.finite.vector(ELL)) {
    warning("Expected losses of layers Attachment_Points[i+1] - Attachment_Points[i] xs Attachment_Points[i] must be positive.")
    Results$Comment <- paste0(Results$Comment, "Expected losses of layers Attachment_Points[i+1] - Attachment_Points[i] xs Attachment_Points[i] must be positive.")
    Results$Status <- 2
    return(Results)
  }

  Limits <- Attachment_Points[2:k] - Attachment_Points[1:(k-1)]
  Limits <- c(Limits, Inf)

  #######################################################################
  # Transform truncation of whole distribution into untruncated problem #
  #######################################################################
  truncation_wd <- FALSE
  if (!is.null(truncation) && !is.infinite(truncation) && truncation_type =="wd") {
    truncation_wd <- TRUE

    k <- k + 1
    Attachment_Points <- c(Attachment_Points, truncation)
    truncation <- NULL
    Limits[k-1] <- Attachment_Points[k] - Attachment_Points[k-1]
    Limits <- c(Limits, Inf)

    if (is.null(Frequencies)) {
      Frequencies <- rep(NaN, k)
    } else {
      Frequencies <- c(Frequencies, NaN)
    }
    if (!is.null(FQ_at_lowest_AttPt)) {
      Frequencies[1] <- FQ_at_lowest_AttPt
      FQ_at_lowest_AttPt <- NULL
    }
    if (!is.null(FQ_at_highest_AttPt)) {
      Frequencies[k-1] <- FQ_at_highest_AttPt
      FQ_at_highest_AttPt <- NULL
    }
    if (is.na(Frequencies[k-1])) {
      alpha_wd <- NA
      try(alpha_wd <- Pareto_Find_Alpha_btw_Layers(Limits[k-2], Attachment_Points[k-2], ELL[k-2], Limits[k-1], Attachment_Points[k-1], ELL[k-1], truncation = Attachment_Points[k]), silent = TRUE)
      if (is.na(alpha_wd)) {
        warning("truncation too low.")
        Results$Comment <- paste0(Results$Comment, "truncation too low.")
        Results$Status <- 2
        return(Results)
      }
      Frequencies[k-1] <- ELL[k-1] / Pareto_Layer_Mean(Limits[k-1], Attachment_Points[k-1], alpha_wd, truncation = Attachment_Points[k])
    }
    alpha_wd <- NA
    try(alpha_wd <- Pareto_Find_Alpha_btw_FQ_Layer(Attachment_Points[k-1], Frequencies[k-1], Limits[k-1], Attachment_Points[k-1], ELL[k-1], truncation = Attachment_Points[k]), silent = TRUE)
    if (is.na(alpha_wd)) {
      warning("truncation too low.")
      Results$Comment <- paste0(Results$Comment, "truncation too low.")
      Results$Status <- 2
      return(Results)
    }
    alpha_wd <- max(alpha_wd, 1.2)

    Frequencies[k] <- 0 # frequency xs truncation point is zero
    f_wd <- (Attachment_Points[k-1] / Attachment_Points[k])^alpha_wd # frequency extrapolation factor without truncation
    # calculate add_fq such that f_wd * (Frequencies[k-1] + add_fq) = Frequencies[k] + add_fq = add_fq:
    add_fq <- Frequencies[k-1] * f_wd / (1 - f_wd) # Frequency xs truncation point if distribution is not truncated
    # redefine Frequencies (from now on frequencies of the corresponding non-truncated distribution):
    Frequencies <- Frequencies + add_fq

    if (!is.null(TotalLoss_Frequencies)) {
      TotalLoss_Frequencies <- c(TotalLoss_Frequencies, 0)
      TotalLoss_Frequencies <- TotalLoss_Frequencies + add_fq
    }
    # redefine ELL (from now on expected layer losses of the corresponding non-truncated distribution):
    ELL <- ELL + Limits[1:(k-1)] * add_fq
    ELL <- c(ELL, NaN)
    ELL[k] <- Frequencies[k] * Pareto_Layer_Mean(Inf, Attachment_Points[k], alpha_wd)
  }

  #########################
  # End of transformation #
  #########################


  RoLs <- ELL / Limits
  Merged_Layer <- rep(FALSE, k)
  if (max(RoLs[2:k] / RoLs[1:(k-1)]) >= 1 - RoL_tolerance) {
    repeat {
      if (k<3) { # should not happen ...
        warning("Layers cannot be merged to obtain strictly decreasing RoLs.")
        Results$Comment <- paste0(Results$Comment, "Layers cannot be merged to obtain strictly decreasing RoLs.")
        Results$Status <- 2
        return(Results)
      }
      pos <- order(RoLs[2:k] / RoLs[1:(k-1)])[k-1]
      ELL[pos] <- ELL[pos] + ELL[pos+1]
      ELL <- ELL[-(pos+1)]
      Attachment_Points <- Attachment_Points[-(pos+1)]
      Limits[pos] <- Limits[pos] + Limits[pos+1]
      Limits <- Limits[-(pos+1)]
      if (!is.null(Frequencies)) {
        Frequencies <- Frequencies[-(pos+1)]
      }
      if (!is.null(TotalLoss_Frequencies)) {
        TotalLoss_Frequencies <- TotalLoss_Frequencies[-pos]
      }
      Merged_Layer <- Merged_Layer[-(pos+1)]
      Merged_Layer[pos] <- TRUE
      k <- k-1
      RoLs <- ELL / Limits
      if (max(RoLs[2:k] / RoLs[1:(k-1)]) < 1 - RoL_tolerance) {break}
    }
    warning("RoLs not strictly decreasing. Layers have been merged.")
    Results$Comment <- paste0(Results$Comment, "RoLs not strictly decreasing. Layers have been merged. ")
    Results$Status <- 1
  }

  if (!is.null(truncation)) {
    RoLs_truncation <- RoLs
    RoLs_truncation[k] <- ELL[k] / (truncation - Attachment_Points[k])
    if (RoLs_truncation[k] / RoLs_truncation[k-1] > 1 - RoL_tolerance) {
      warning("truncation too low!")
      Results$Comment <- paste0(Results$Comment, "truncation too low!")
      Results$Status <- 2
      return(Results)
    }
  }

  if (!is.null(Frequencies)) {
    if (max(RoLs/Frequencies, na.rm = T) >= 1 - RoL_tolerance / 2) {
      Frequencies <- NULL
      TotalLoss_Frequencies <- NULL
      warning("Layer entry frequencies not strictly greater than RoLs.")
      Results$Comment <- paste0(Results$Comment, "Layer entry frequencies not strictly greater than RoLs. Frequencies not used!  ")
      Results$Status <- 1
    }
    if (min(RoLs[1:(k-1)]/Frequencies[2:k], na.rm = T) <= 1 + RoL_tolerance / 2) {
      Frequencies <- NULL
      TotalLoss_Frequencies <- NULL
      warning("Layer exit frequencies not strictly less than RoLs.")
      Results$Comment <- paste0(Results$Comment, "Layer exit frequencies not strictly less than RoLs. Frequencies not used! ")
      Results$Status <- 1
    }
  }
  if (is.null(Frequencies)) {
    Frequencies <- rep(NaN, k)
  }
  alpha_between_layers <- rep(NaN, k-1)
  for (i in 1:(k-1)) {
    if (i < k-1 && is.na(Frequencies[i+1])) {
      try(alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Limits[i+1], Attachment_Points[i+1], ELL[i+1]), silent = TRUE)
      if (is.na(alpha_between_layers[i])) {
        Frequencies[i+1] <- (RoLs[i] + RoLs[i+1]) / 2
      } else {
        Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Limits[i+1], Attachment_Points[i+1], alpha_between_layers[i])
      }
      if (Merged_Layer[i] & !Merged_Layer[i+1]) {
        Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
      } else if (!Merged_Layer[i] & Merged_Layer[i+1]) {
        Frequencies[i+1] <- RoLs[i+1] * (1 + RoL_tolerance / 2)
      }
    } else if (is.na(Frequencies[i+1])) {
      if (is.null(truncation)) {
        if (Use_unlimited_Layer_for_FQ) {
          try(alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Inf, Attachment_Points[i+1], ELL[i+1]), silent = TRUE)
          if (!is.na(alpha_between_layers[i])) {
            Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Inf, Attachment_Points[i+1], alpha_between_layers[i])
          } else {
            Frequencies[i+1] <- RoLs[i] / 2
          }
          if (Merged_Layer[i]) {
            Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
          }
        } else {
          Frequencies[i+1] <- Frequencies[i] * (Attachment_Points[i]/Attachment_Points[i+1])^alpha_between_layers[i-1]
        }
      } else {
        try(alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Inf, Attachment_Points[i+1], ELL[i+1], truncation = truncation), silent = TRUE)
        if (!is.na(alpha_between_layers[i])) {
          Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Inf, Attachment_Points[i+1], alpha_between_layers[i], truncation = truncation)
        } else {
          warning("truncation too low!")
          Results$Comment <- paste0(Results$Comment, "truncation too low!")
          Results$Status <- 2
          return(Results)
        }
        if (Frequencies[i+1] >= RoLs[i] * (1 - RoL_tolerance / 2)) {
          Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
          warning("Option Use_unlimited_Layer_for_FQ not used!")
          Results$Comment <- paste0(Results$Comment, "Option Use_unlimited_Layer_for_FQ not used! ")
          Results$Status <- 1
        }
        if (Merged_Layer[i]) {
          Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
        }
      }
    }
  }
  if (is.na(Frequencies[1])) {
    if (!Merged_Layer[1]) {
      alpha_between_layers[1] <-  Pareto_Find_Alpha_btw_Layers(Limits[1], Attachment_Points[1], ELL[1], Limits[2], Attachment_Points[2], ELL[2])
      Frequencies[1] <- ELL[1] / Pareto_Layer_Mean(Limits[1], Attachment_Points[1], alpha_between_layers[1])
    } else {
      Frequencies[1] <- RoLs[1] * (1 + RoL_tolerance / 2)
    }
  }

  if (!is.null(FQ_at_lowest_AttPt)) {
    if (FQ_at_lowest_AttPt > RoLs[1] * (1 + RoL_tolerance / 2)) {
      Frequencies[1] <- FQ_at_lowest_AttPt
    } else {
      warning("FQ_at_lowest_AttPt too small. Not used!")
      Results$Comment <- paste0(Results$Comment, "FQ_at_lowest_AttPt too small. Not used! ")
      Results$Status <- 1
    }
  }
  if (!is.null(FQ_at_highest_AttPt)) {
    if (FQ_at_highest_AttPt < RoLs[k-1] * (1 - RoL_tolerance / 2)) {
      Frequencies[k] <- FQ_at_highest_AttPt
    } else {
      warning("FQ_at_highest_AttPt too large.")
      Results$Comment <- paste0(Results$Comment, "FQ_at_highest_AttPt too large. Not used! ")
      Results$Status <- 1
    }
  }
  if (!is.null(TotalLoss_Frequencies)) {
    if (max(TotalLoss_Frequencies - RoLs[1:(k-1)]) >= 0) {
      TotalLoss_Frequencies <- NULL
      warning("TotalLoss_Frequencies not strictly less than RoLs.")
      Results$Comment <- paste0(Results$Comment, "TotalLoss_Frequencies not strictly less than RoLs. TotalLoss_Frequencies not used! ")
      Results$Status <- 1
    }
    if (max(TotalLoss_Frequencies - Frequencies[2:k]) < 0) {
      TotalLoss_Frequencies <- NULL
      warning("TotalLoss_Frequencies not greater than or equal to Frequencies.")
      Results$Comment <- paste0(Results$Comment, "TotalLoss_Frequencies not greater than or equal to Frequencies. TotalLoss_Frequencies not used! ")
      Results$Status <- 1
    }
  }

  # Check whether frequency at highest attachment point is large enough in case of truncation
  if (!is.null(truncation)) {
    alpha_test <- NaN
    try(alpha_test <- Pareto_Find_Alpha_btw_FQ_Layer(Attachment_Points[k], Frequencies[k], Inf, Attachment_Points[k], ELL[k], truncation = truncation), silent = TRUE)
    if (is.na(alpha_test)) {
      warning("truncation too low.")
      Results$Comment <- paste0(Results$Comment, "truncation too low. ")
      Results$Status <- 2
      return(Results)

    }
  }

  if (!is.null(TotalLoss_Frequencies)) {
    if (max(TotalLoss_Frequencies - Frequencies[2:k]) > 0) {
      repeat {
        pos <- order(TotalLoss_Frequencies - Frequencies[2:k])[k-1]
        new_AttPoint <- Attachment_Points[pos+1] * (Frequencies[pos+1]/TotalLoss_Frequencies[pos])^(1/alpha_max)
        if (new_AttPoint <= Attachment_Points[pos]) {
          warning("It is not possible to fulfill all TotalLoss_Frequencies.")
          Results$Comment <- paste0(Results$Comment, "It is not possible to fulfill all TotalLoss_Frequencies. ")
          Results$Status <- 1
          TotalLoss_Frequencies[pos] <- Frequencies[pos+1]
        } else {
          Attachment_Points <- c(Attachment_Points[1:pos], new_AttPoint, Attachment_Points[(pos+1):k])
          Frequencies <- c(Frequencies[1:pos], TotalLoss_Frequencies[pos], Frequencies[(pos+1):k])
          if (pos <= k-2) {
            TotalLoss_Frequencies <- c(TotalLoss_Frequencies[1:pos], Frequencies[pos+2], TotalLoss_Frequencies[(pos+1):(k-1)])
          } else {
            TotalLoss_Frequencies <- c(TotalLoss_Frequencies[1:pos], Frequencies[pos+2])
          }
          new_ELL <- Pareto_Layer_Mean(Attachment_Points[pos+2] - Attachment_Points[pos+1], Attachment_Points[pos+1], alpha_max) * Frequencies[pos+1]
          ELL <- c(ELL[1:pos], new_ELL, ELL[(pos+1):k])
          ELL[pos] <- ELL[pos] - new_ELL
          k <- k+1
        }
        if (max(TotalLoss_Frequencies - Frequencies[2:k]) <= 0) {break}
      }
    }
  }


  if (!is.null(truncation)) {
    AvLossLastLayer <- Pareto_Layer_Mean(Inf, Attachment_Points[k], alpha = tolerance, truncation = truncation)
    if (Frequencies[k] * AvLossLastLayer <= ELL[k]) {
      warning("truncation too low!")
      Results$Comment <- paste0(Results$Comment, "truncation too low!")
      Results$Status <- 2
      return(Results)
    }
  }


  s <- Frequencies / Frequencies[1]
  l <- ELL / Frequencies[1]

  Results_Fit_PP <- NULL
  try(Results_Fit_PP <- Fit_PP(Attachment_Points, s, l, truncation = truncation, alpha_max = alpha_max, minimize_ratios = minimize_ratios, merge_tolerance = merge_tolerance), silent = TRUE)
  if (is.null(Results_Fit_PP)) {
    Results$Comment <- paste0(Results$Comment, "Status of Fit_PP not OK")
    Results$Status <- 2
    return(Results)
  }

  if (Results_Fit_PP$Status != "OK") {
    Results$Comment <- paste0(Results$Comment, "Status of Fit_PP not OK")
    Results$Status <- 2
    return(Results)
  } else {
    Results$Status <- 0
  }

  Results$t <- Results_Fit_PP$t
  Results$alpha <- Results_Fit_PP$alpha
  if (truncation_wd) {
    Results$FQ <- Frequencies[1] - Frequencies[k]
  } else {
    Results$FQ <- Frequencies[1]
  }
  if (Results$Comment == "") {Results$Comment <- "OK"}

  if (!is.valid.PPP_Model(Results)) {
    warning("Resulting PPP_Model not valid.")
    Results$Comment <- paste0(Results$Comment, "Resulting PPP_Model not valid.")
    Results$Status <- 2
  }
  return(Results)
}


#' Distribution Function of the Pareto Distribution
#'
#' @description Calculates the cumulative distribution function of a Pareto distribution. This function is deprecated. Use \code{pPareto} instead.
#'
#' @param x Numeric. The function evaluates the CDF at \code{x}.
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param alpha Numeric. Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return Distribution function of the Pareto distribution with parameters \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' x <- 0:10 * 1000
#' pPareto(x, 1000, 2)
#' pPareto(x, 1000, 2, truncation = 5000)
#'
#' @export

Pareto_CDF <- function(x, t, alpha, truncation = NULL) {
  .Deprecated("pPareto")
  pPareto(x, t, alpha, truncation)
}


#' Distribution Function of the Pareto Distribution
#'
#' @description Calculates the cumulative distribution function of a Pareto distribution
#'
#' @param x Numeric. The function evaluates the CDF at \code{x}.
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param alpha Numeric. Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return Distribution function of the Pareto distribution with parameters \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' x <- 0:10 * 1000
#' pPareto(x, 1000, 2)
#' pPareto(x, 1000, 2, truncation = 5000)
#'
#'
#' @export

pPareto <- function(x, t, alpha, truncation = NULL) {
  if (!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE)) {
    warning(valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE, comment = TRUE))
    return(rep(NaN, length(x)))
  }
  sapply(x, FUN = function(x) pPareto_s(x, t, alpha, truncation))
}

pPareto_s <- function(x, t, alpha, truncation = NULL) {
  if (!valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE)) {
    warning(valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE, comment = TRUE))
    return(NaN)
  }
  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }

  if (alpha <= 1e-6) {
    if (is.positive.finite.number(truncation)) {
      if (x <= t) {
        return(0)
      } else if (x < truncation) {
        return(log(x/t) / log(truncation / t))
      } else {
        return(1)
      }
    } else {
      if (is.infinite(x) && is.positive.number(x)) {
        return(1)
      } else {
        return(0)
      }
    }

  }

  if (x <= t) {
    return(0)
  } else if (is.null(truncation)) {
    Result <- 1 - (t / x)^alpha
    return(Result)
  } else if (x >= truncation) {
    return(1)
  } else {
    Result <- (1 - (t / x)^alpha) / (1 - (t / truncation)^alpha)
    return(Result)
  }
}



#' Density of the Pareto Distribution
#'
#' @description Calculates the density function of the Pareto distribution. This function is deprecated. Use \code{dPareto} instead.
#'
#' @param x Numeric. The function evaluates the density at \code{x}
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param alpha Numeric. Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return Density function of the Pareto distribution with parameters \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' x <- 0:10 * 1000
#' dPareto(x, 1000, 2)
#' dPareto(x, 1000, 2, truncation = 5000)
#'
#' @export

Pareto_PDF <- function(x, t, alpha, truncation = NULL) {
  .Deprecated("dPareto")
  dPareto(x, t, alpha, truncation)
}





#' Density of the Pareto Distribution
#'
#' @description Calculates the density function of the Pareto distribution
#'
#' @param x Numeric. The function evaluates the density at x.
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param alpha Numeric. Pareto alpha.
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation.
#'
#' @return Density function of the Pareto distribution with parameters \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' x <- 0:10 * 1000
#' dPareto(x, 1000, 2)
#' dPareto(x, 1000, 2, truncation = 5000)
#'
#' @export

dPareto <- function(x, t, alpha, truncation = NULL) {
  if (!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE)) {
    warning(valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE, comment = TRUE))
    return(rep(NaN, length(x)))
  }
  sapply(x, FUN = function(x) dPareto_s(x, t, alpha, truncation))
}

dPareto_s <- function(x, t, alpha, truncation = NULL) {
  if (!valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE)) {
    warning(valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE, comment = TRUE))
    return(NaN)
  }
  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }

  if (alpha <= 1e-6) {
    if (is.positive.finite.number(truncation)) {
      if (x <= t) {
        return(0)
      } else if (x < truncation) {
        return((1/x) / log(truncation / t))
      } else {
        return(0)
      }
    } else {
      if (is.infinite(x) && is.positive.number(x)) {
        return(Inf)
      } else {
        return(0)
      }
    }
  }

  if (x <= t) {
    return(0)
  } else if (is.null(truncation)) {
    Result <- t^alpha * alpha / x^(alpha + 1)
    return(Result)
  } else if (x >= truncation) {
    return(0)
  } else {
    Result <- t^alpha * alpha / x^(alpha + 1) / pPareto(truncation, t, alpha)
    return(Result)
  }
}




#' Distribution Function of the Piecewise Pareto Distribution
#'
#' @description Calculates the cumulative distribution function of a Piecewise Pareto Distribution. This function is deprecated. Use \code{pPiecewisePareto} instead.
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param x Numeric. The function evaluates the CDF at \code{x}.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return Distribution function of the piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' x <- 0:10 * 1000
#' pPiecewisePareto(x, t, alpha)
#' pPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "lp")
#' pPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

PiecewisePareto_CDF <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  .Deprecated("pPiecewisePareto")
  pPiecewisePareto(x, t, alpha, truncation, truncation_type)
}








#' Distribution Function of the Piecewise Pareto Distribution
#'
#' @description Calculates the cumulative distribution function of a Piecewise Pareto Distribution
#'
#' @references Riegel, U. (2018) Matching tower information with piecewise Pareto. European Actuarial Journal 8(2): 437--460
#'
#' @param x Numeric. The function evaluates the CDF at \code{x}.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return Distribution function of the piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' x <- 0:10 * 1000
#' pPiecewisePareto(x, t, alpha)
#' pPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "lp")
#' pPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

pPiecewisePareto <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(rep(NaN, length(x)))
  }
  sapply(x, FUN = function(x) pPiecewisePareto_s(x, t, alpha, truncation, truncation_type))
}

pPiecewisePareto_s <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }

  n <- length(t)

  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }

  if (n == 1) {
    Result <- pPareto(x, t, alpha, truncation)
    return(Result)
  }


  if (x <= t[1]) {
    return(0)
  }

  if (is.null(truncation)) {
    t <- t[t<x]
    t <- c(t,x)
    n <- length(t)

    factors_t <- t[2:n] / t[1:(n-1)]
    return(1 - prod((1/factors_t)^alpha[1:(n-1)]))
  } else if (truncation_type == "wd") {
    if (x >= truncation) {
      return(1)
    }
    factors_t <- t[2:n] / t[1:(n-1)]
    scaling <- 1 / (1 - prod((1/factors_t)^alpha[1:(n-1)]) * (t[n] / truncation)^alpha[n])

    t <- t[t<x]
    t <- c(t,x)
    n <- length(t)

    factors_t <- t[2:n] / t[1:(n-1)]
    return(scaling * (1 - prod((1/factors_t)^alpha[1:(n-1)])))

  } else {
    if (x >= truncation) {
      return(1)
    } else if (x <= t[n]) {
      t <- t[t<x]
      t <- c(t,x)
      n <- length(t)

      factors_t <- t[2:n] / t[1:(n-1)]
      return(1 - prod((1/factors_t)^alpha[1:(n-1)]))
    } else {
      factors_t <- t[2:n] / t[1:(n-1)]
      excess_prob <- prod((1/factors_t)^alpha[1:(n-1)])
      return(1 - excess_prob * (1 - (1-(t[n]/x)^alpha[n]) / (1-(t[n]/truncation)^alpha[n])))
    }
  }
}


#' Density of the Piecewise Pareto Distribution
#'
#' @description Calculates the density function of the piecewise Pareto distribution. This function is deprecated. Use \code{dPiecewisePareto} instead.
#'
#' @param x Numeric. The function evaluates the density at \code{x}.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return Density function of the piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' x <- 0:10 * 1000
#' dPiecewisePareto(x, t, alpha)
#' dPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "lp")
#' dPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

PiecewisePareto_PDF <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  .Deprecated("dPiecewisePareto")
  dPiecewisePareto(x, t, alpha, truncation, truncation_type)
}





#' Density of the Piecewise Pareto Distribution
#'
#' @description Calculates the density function of the piecewise Pareto distribution
#'
#' @param x Numeric. The function evaluates the density at \code{x}.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return Density function of the piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} evaluated at \code{x}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' x <- 0:10 * 1000
#' dPiecewisePareto(x, t, alpha)
#' dPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "lp")
#' dPiecewisePareto(x, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

dPiecewisePareto <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(rep(NaN, length(x)))
  }
  sapply(x, FUN = function(x) dPiecewisePareto_s(x, t, alpha, truncation, truncation_type))
}

dPiecewisePareto_s <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }

  n <- length(t)

  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }

    if (n == 1) {
    Result <- dPareto(x, t, alpha, truncation)
    return(Result)
  }


  if (x <= t[1]) {
    return(0)
  }


  if (is.null(truncation)) {
    t <- t[t<=x]
    n <- length(t)
    alpha <- alpha[1:n]
    excess_prob <- 1 - pPiecewisePareto(x, t, alpha)
    return(excess_prob * alpha[n] / x)

  } else if (truncation_type == "wd") {
    if (x >= truncation) {
      return(0)
    }
    scaling <- 1 / pPiecewisePareto(truncation, t, alpha)
    t <- t[t<=x]
    n <- length(t)
    alpha <- alpha[1:n]
    excess_prob <- 1 - pPiecewisePareto(x, t, alpha)
    return(excess_prob * alpha[n] / x * scaling)

  } else {
    if (x >= truncation) {
      return(0)
    } else if (x <= t[n]) {
      t <- t[t<=x]
      n <- length(t)
      alpha <- alpha[1:n]
      excess_prob <- 1 - pPiecewisePareto(x, t, alpha)
      return(excess_prob * alpha[n] / x)
    } else {
      n <- length(t)
      excess_prob <- 1 - pPiecewisePareto(t[n], t, alpha)
      excess_prob_trunc <- 1 - pPiecewisePareto(truncation, t, alpha)
      scaling <- excess_prob / (excess_prob - excess_prob_trunc)
      return(excess_prob * scaling * dPareto(x, t[n], alpha[n]))
    }
  }
}




#' Quantile Function of the Piecewise Pareto Distribution
#'
#' @description Calculates the quantile function of a piecewise Pareto distribution
#'
#' @param p Numeric. The function evaluates the quantile function at \code{p}.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#'
#' @return Quantile function of the piecewise Pareto distribution with parameter vectors \code{t} and \code{alpha} evaluated at \code{p}
#'
#' @examples
#' t <- c(1000, 2000, 3000)
#' alpha <- c(1, 1.5, 2)
#' p <- 0:10 * 0.1
#' qPiecewisePareto(p, t, alpha)
#' qPiecewisePareto(p, t, alpha, truncation = 5000, truncation_type = "lp")
#' qPiecewisePareto(p, t, alpha, truncation = 5000, truncation_type = "wd")
#'
#' @export

qPiecewisePareto <- function(p, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.positive.finite.vector(p) || max(p) > 1) {
    warning("p must be a vector with elements in [0,1].")
    return(rep(NaN, length(p)))
  }
  if (!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(rep(NaN, length(p)))
  }
  sapply(p, FUN = function(y) qPiecewisePareto_s(y, t, alpha, truncation, truncation_type))
}

qPiecewisePareto_s <- function(y, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type)) {
    warning(valid.parameters.PiecewisePareto(t, alpha, truncation, truncation_type, comment = TRUE))
    return(NaN)
  }

  n <- length(t)

  if (!is.positive.finite.number(y) || y > 1) {
    warning("y must be a number in [0,1].")
    return(NaN)
  }

  if (n == 1) {
    Result <- qPareto(y, t, alpha, truncation)
    return(Result)
  }

  if (y == 1) {
    if (is.null(truncation)) {
      return(Inf)
    } else {
      return(truncation)
    }
  }


  CDF_untruncated_at_t <- pPiecewisePareto(t, t, alpha)
  if (!is.null(truncation)) {
    if (is.infinite(truncation)) {
      CDF_untruncated_at_truncation <- 1
    } else {
      CDF_untruncated_at_truncation <- pPiecewisePareto(truncation, t, alpha)
    }
    if (truncation_type == "wd") {
      y <- y * CDF_untruncated_at_truncation
    } else {
      if (y > CDF_untruncated_at_t[n]) {
        y <- CDF_untruncated_at_t[n] + (y - CDF_untruncated_at_t[n]) * (CDF_untruncated_at_truncation - CDF_untruncated_at_t[n]) / (1 - CDF_untruncated_at_t[n])
      }
    }
  }

  k_0 <- sum(CDF_untruncated_at_t <= y)
  t_0 <- t[k_0]
  F_0 <- CDF_untruncated_at_t[k_0]


  result <- t_0 / ((1 - y) / (1 - F_0))^(1 / alpha[k_0])
  return(result)
}




#' Quantile Function of the Pareto Distribution
#'
#' @description Calculates the quantile function of a Pareto distribution
#'
#' @param p Numeric. The function evaluates the inverse CDF at \code{p}.
#' @param t Numeric. Threshold of the piecewise Pareto distribution.
#' @param alpha Numeric. Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return Quantile function of the Pareto distribution with parameters \code{t} and \code{alpha}, evaluated at \code{p}
#'
#' @examples
#' p <- 0:10 * 0.1
#' qPareto(p, 1000, 2)
#' qPareto(p, 1000, 2, truncation = 5000)
#'
#' @export

qPareto <- function(p, t, alpha, truncation = NULL) {
  if (!is.nonnegative.finite.vector(p) || max(p) > 1) {
    warning("p must be a vector with elements in [0,1].")
    return(rep(NaN, length(p)))
  }
  if (!valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE)) {
    warning(valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE, comment = TRUE))
    return(rep(NaN, length(p)))
  }
  sapply(p, FUN = function(y) qPareto_s(y, t, alpha, truncation))
}

qPareto_s <- function(y, t, alpha, truncation = NULL) {
  if (!valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE)) {
    warning(valid.parameters.Pareto(t, alpha, truncation, allow.alpha.zero = TRUE, comment = TRUE))
    return(NaN)
  }
  if (!is.nonnegative.finite.number(y) || y > 1) {
    warning("y must be a number in [0,1].")
    return(NaN)
  }

  if (alpha <= 1e-6) {
    if (is.positive.finite.number(truncation)) {
      if (y == 0) {
        return(t)
      } else if (y == 1) {
        return(truncation)
      } else {
        return(t * (truncation / t)^y)
      }
    } else {
      return(Inf)
    }

  }


  if (y == 1) {
    if (is.null(truncation)) {
      return(Inf)
    } else {
      return(truncation)
    }
  }

  if (!is.null(truncation)) {
    if (!is.infinite(truncation)) {
      scale <- pPareto(truncation, t, alpha)
      y <- y * scale
    }
  }
  result <- t / (1 - y)^(1 / alpha)

}


#' Maximum Likelihood Estimation of the Pareto Alpha
#'
#' @description Calculates the maximum likelihood estimator of the parameter alpha of a Pareto distribution
#'
#' @param losses Numeric vector. Losses that are used for the ML estimation.
#' @param t Numeric or numeric vector. Threshold of the Pareto distribution. Alternatively, \code{t} can be a vector of same length as \code{losses}. In this case \code{t[i]} is the reporting threshold of \code{losses[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#' @param tol Numeric. Desired accuracy  (only relevant in the truncated case).
#' @param max_iterations Numeric. Maximum number of iteration in the case \code{truncation < Inf}  (only relevant in the truncated case).
#' @param alpha_min Numeric. Deprecated.
#' @param alpha_max Numeric. Deprecated.
#'
#' @return Maximum likelihood estimator for the parameter \code{alpha} of a Pareto distribution with threshold \code{t} given the observations \code{losses}
#'
#' @examples
#' losses <- rPareto(100, 1000, 2)
#' Pareto_ML_Estimator_Alpha(losses, 1000)
#' losses <- rPareto(100, 1000, 2, truncation = 2000)
#' Pareto_ML_Estimator_Alpha(losses, 1000)
#' Pareto_ML_Estimator_Alpha(losses, 1000, truncation = 2000)
#'
#' t <- rPareto(10000, 100, 2)
#' alpha <- 2
#' losses <- rPareto(10000, t, alpha)
#' Pareto_ML_Estimator_Alpha(losses, t)
#' losses <- rPareto(10000, t, alpha, truncation = 2 * max(t))
#' Pareto_ML_Estimator_Alpha(losses, t, truncation = 2 * max(t))
#' @export

Pareto_ML_Estimator_Alpha <- function(losses, t, truncation = NULL, tol = 1e-7, max_iterations = 1000, alpha_min = 0, alpha_max = Inf) {
  if (!missing(alpha_min) || !missing(alpha_max)) {
    warning("arguments alpha_min and alpha_max are deprecated and are ignored", call. = FALSE)
  }
  if (!is.nonnegative.finite.vector(losses)) {
    warning("losses must be non-negative.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be positive.")
    return(NaN)
  }
  if (length(t) != 1 && length(t) != length(losses)) {
    warning("t must have length 1 or same length as losses.")
    return(NaN)
  }

  losses <- losses[losses > t]
  n <- length(losses)

  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
    }
    if (truncation <= max(t)) {
      warning("truncation must be larger than t")
    }
    if (max(losses) >= truncation) {
      warning("Losses must be < truncation.")
      return(NaN)
    }
  }
  if (is.null(truncation) || is.infinite(truncation)) {
    alpha_hat <- n / sum(log(losses / t))
  } else {
    alpha_hat_iteration <- numeric(max_iterations)
    alpha_hat_iteration[1] <- n / sum(log(losses / t))
    iteration <- function(alpha_hat) {
      if (length(t) == 1) {
        result <- n * (sum(log(losses / t)) - n * log(t / truncation) * (t / truncation)^alpha_hat / (1 - (t / truncation)^alpha_hat))^{-1}
      } else {
        result <- n * (sum(log(losses / t)) - sum(log(t / truncation) * (t / truncation)^alpha_hat / (1 - (t / truncation)^alpha_hat)))^{-1}
      }
      return(result)
    }


    if(!is.positive.finite.number(max_iterations) || max_iterations < 2) {
      warning("max_iterations must be > 1")
      return(NaN)
    } else {
      max_iterations <- ceiling(max_iterations)
    }

    if (!is.positive.finite.number(tol)) {
      warning("tol must be > 1")
      return(NaN)
    }
    check_tol <- function(i) {
      if (abs(alpha_hat_iteration[i] - alpha_hat_iteration[i-1]) < tol) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }

    for (i in 2:max_iterations) {
      alpha_hat_iteration[i] <- iteration(alpha_hat_iteration[i-1])
      if (check_tol(i)) {
        alpha_hat <- alpha_hat_iteration[i]
        return(alpha_hat)
      }
    }
    warning("desired accuracy not achieved; increase max_iterations")
    alpha_hat <- alpha_hat_iteration[max_iterations]
  }
  return(alpha_hat)
}


#' Maximum Likelihood Estimation of the Alphas of the Piecewise Pareto Distribution
#'
#' @description Calculates the maximum likelihood estimator of the parameter vector alpha of a piecewise Pareto distribution
#'
#' @param losses Numeric vector. Losses that are used for the ML estimation.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > max(t)}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Character. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#' @param tol Numeric. Desired accuracy (only relevant in the truncated case).
#' @param max_iterations Numeric. Maximum number of iteration in the case \code{truncation < Inf}  (only relevant in the truncated case).
#'
#' @return Maximum likelihood estimator for the parameter \code{alpha} of a Pareto distribution with threshold \code{t} given the observations \code{losses}
#'
#' @examples
#' losses <- rPiecewisePareto(10000, t = c(100,200,300), alpha = c(1,2,3))
#' PiecewisePareto_ML_Estimator_Alpha(losses, c(100,200,300))
#' losses <- rPiecewisePareto(10000, t = c(100,200,300), alpha = c(1,2,3),
#'                            truncation = 500, truncation_type = "lp")
#' PiecewisePareto_ML_Estimator_Alpha(losses, c(100,200,300))
#' PiecewisePareto_ML_Estimator_Alpha(losses, c(100,200,300),
#'                                    truncation = 500, truncation_type = "lp")
#' losses <- rPiecewisePareto(10000, t = c(100,200,300), alpha = c(1,2,3),
#'                            truncation = 500, truncation_type = "wd")
#' PiecewisePareto_ML_Estimator_Alpha(losses, c(100,200,300))
#' PiecewisePareto_ML_Estimator_Alpha(losses, c(100,200,300),
#'                                    truncation = 500, truncation_type = "wd")
#'
#' @export

PiecewisePareto_ML_Estimator_Alpha <- function(losses, t, truncation = NULL, truncation_type = "lp", tol = 1e-7, max_iterations = 1000) {
  if (!is.nonnegative.finite.vector(losses)) {
    warning("losses must be non-negative.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be positive.")
    return(NaN)
  }
  if (!is.atomic(truncation_type) || length(truncation_type) != 1 || !(truncation_type %in% c("lp", "wd"))) {
    warning("truncation_type must be 'lp' or 'wd'.")
    return(NaN)
  }
  if (!is.positive.finite.number(tol)) {
    warning("tol must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(max_iterations) || max_iterations < 2) {
    warning("max_iterations must be a number >= 2.")
    return(NaN)
  } else {
    max_iterations <- ceiling(max_iterations)
  }


  k <- length(t)



  if (min(diff(t)) <= 0) {
    warning("t must be strictily ascending!")
    return(NaN)
  }
  if (max(losses) <= max(t)) {
    warning("Number of losses > max(t) must be positive.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be a positive number ('Inf' allowed).")
    }
    if (truncation <= t[k]) {
      warning("truncation must be larger than max(t)")
      return(NaN)
    }
    if (max(losses) >= truncation) {
      warning("Losses must be < truncation.")
      return(NaN)
    }
  } else {
    truncation <- Inf
  }

  if (k == 1) {
    Result <- Pareto_ML_Estimator_Alpha(losses, t, truncation = truncation, tol = tol, max_iterations = max_iterations)
    return(Result)
  }


  t <- c(t, truncation)
  losses_xs_t <- lapply(t, function(x) losses[losses >= x])
  count_xs_t <- unlist(lapply(losses_xs_t, "length"))

  alpha_hat <- numeric(k)
  for (i in 1:k) {
    alpha_hat[i] <- (count_xs_t[i] - count_xs_t[i+1]) / (sum(log(pmin(losses_xs_t[[i]], t[i+1]) / t[i])))
  }
  if (!is.infinite(truncation)) {
    if (truncation_type == "lp") {
      alpha_hat[k] <- Pareto_ML_Estimator_Alpha(losses_xs_t[[k]], t[k], truncation = truncation, tol = tol, max_iterations = max_iterations)
    } else {
      iteration <- function(alpha_hat) {
        result <- numeric(k)
        for (i in 1:k) {
          product <- prod((t[-(k+1)] / t[-1])^alpha_hat)
          result[i] <- (count_xs_t[i] - count_xs_t[i+1]) / (sum(log(pmin(losses_xs_t[[i]], t[i+1]) / t[i])) - count_xs_t[1] * log(t[i] / t[i+1]) * product  / (1 - product))
        }
        return(result)
      }

      if(max_iterations < 2) {
        warning("max_iterations must be > 1")
        return(NaN)
      }

      alpha_hat_iteration <- matrix(NaN, nrow = k, ncol = max_iterations)
      alpha_hat_iteration[, 1] <- alpha_hat

      check_tol <- function(i) {
        if (max(abs(alpha_hat_iteration[, i] - alpha_hat_iteration[, i-1])) < tol) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }

      for (i in 2:max_iterations) {
        alpha_hat_iteration[, i] <- iteration(alpha_hat_iteration[, i-1])
        if (check_tol(i)) {
          alpha_hat <- alpha_hat_iteration[, i]
          return(alpha_hat)
        }
      }
      warning("desired accuracy not achieved; increase max_iterations")
      alpha_hat <- alpha_hat_iteration[, max_iterations]
    }
  }
  return(alpha_hat)
}


#' Local Pareto Alpha
#'
#' @description Calculates the local Pareto alpha of the normal, lognormal and gamma distribution
#'
#' @references Riegel, U. (2008) Generalizations of common ILF models. Blaetter der DGVFM 29: 45--71
#'
#' @param x Numeric. Vector of thresholds at which the local Pareto alpha is calculated.
#' @param distribution Character. \itemize{
#' \item \code{'lnorm'} for lognormal distribution (arguments: \code{meanlog}, \code{sdlog})
#' \item \code{'norm'} for normal distribution (arguments: \code{mean}, \code{sd})
#' \item \code{'gamma'} for gamma distribution (arguments: \code{shape}, \code{rate}, \code{scale})
#' }
#' @param ... Arguments for the selected distribution
#'
#' @return Local Pareto alpha of the selected distribution at \code{x}
#'
#' @examples
#' x <- 1:10
#' Local_Pareto_Alpha(x, "norm", mean = 1, sd = 5)
#' x <- 1:10 * 1000000
#' Local_Pareto_Alpha(x, "lnorm", meanlog = 1, sdlog = 5)
#'
#' @export


Local_Pareto_Alpha <- function(x, distribution, ...) {
  if (!is.atomic(x) || !is.numeric(x)) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!(distribution %in% c("lnorm", "norm", "gamma", "weibull", "exp"))) {
    warning("distribution must be 'lnorm', 'norm', 'gamma', 'weibull' or 'exp'.")
    return(rep(NaN, length(x)))
  }
  if (distribution == "lnorm") {
    Result <- x * stats::dlnorm(x, log = FALSE, ...) / (1 - stats::plnorm(x, log.p = FALSE, ...))
  }
  if (distribution == "norm") {
    Result <- x * stats::dnorm(x, log = FALSE, ...) / (1 - stats::pnorm(x, log.p = FALSE, ...))
  }
  if (distribution == "gamma") {
    Result <- x * stats::dgamma(x, log = FALSE, ...) / (1 - stats::pgamma(x, log.p = FALSE, ...))
  }
  if (distribution == "weibull") {
    Result <- x * stats::dweibull(x, log = FALSE, ...) / (1 - stats::pweibull(x, log.p = FALSE, ...))
  }
  if (distribution == "exp") {
    Result <- x * stats::dexp(x, log = FALSE, ...) / (1 - stats::pexp(x, log.p = FALSE, ...))
  }
  Result[is.na(x) | is.infinite(x) | (x < 0)] <- NaN
  return(Result)
}


















#' Distribution Function of the generalized Pareto Distribution
#'
#' @description Calculates the cumulative distribution function of a generalized Pareto distribution
#'
#' @param x Numeric. The function evaluates the CDF at \code{x}.
#' @param t Numeric. Threshold of the generalized Pareto distribution.
#' @param alpha_ini Numeric. Initial Pareto alpha.
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the generalized Pareto distribution is truncated at \code{truncation}.
#'
#' @return Distribution function of the generalized Pareto distribution with parameters \code{t}, \code{alpha_ini} and \code{alpha_tail} evaluated at \code{x}
#'
#' @examples
#' x <- 0:10 * 1000
#' pGenPareto(x, 1000, 1, 3)
#' pGenPareto(x, 1000, 1, 3, truncation = 5000)
#'
#'
#' @export

pGenPareto <- function(x, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation)) {
    warning(valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation, comment = TRUE))
    return(rep(NaN, length(x)))
  }
  sapply(x, FUN = function(x) pGenPareto_s(x, t, alpha_ini, alpha_tail, truncation))
}

pGenPareto_s <- function(x, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation)) {
    warning(valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation, comment = TRUE))
    return(NaN)
  }
  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }


  if (x <= t) {
    return(0)
  } else if (is.null(truncation)) {
    Result <- 1 - (1 + alpha_ini / alpha_tail * (x / t - 1))^(-alpha_tail)
    return(Result)
  } else if (x >= truncation) {
    return(1)
  } else {
    Result <- (1 - (1 + alpha_ini / alpha_tail * (x / t - 1))^(-alpha_tail)) / (1 - (1 + alpha_ini / alpha_tail * (truncation / t - 1))^(-alpha_tail))
    return(Result)
  }
}







#' Density of the generalized Pareto Distribution
#'
#' @description Calculates the density function of the generalized Pareto distribution
#'
#' @param x Numeric. The function evaluates the density at x.
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param alpha_ini Numeric. Initial Pareto alpha.
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the generalized Pareto distribution is truncated at truncation.
#'
#' @return Density function of the Pareto distribution with parameters \code{t}, \code{alpha_ini} and \code{alpha_tail} evaluated at \code{x}
#'
#' @examples
#' x <- 0:10 * 1000
#' dGenPareto(x, 1000, 1, 3)
#' dGenPareto(x, 1000, 1, 3, truncation = 5000)
#'
#' @export

dGenPareto <- function(x, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    warning("x must be a numeric vector.")
    return(rep(NaN, length(x)))
  }
  if (!valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation)) {
    warning(valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation, comment = TRUE))
    return(rep(NaN, length(x)))
  }
  sapply(x, FUN = function(x) dGenPareto_s(x, t, alpha_ini, alpha_tail, truncation))
}

dGenPareto_s <- function(x, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation)) {
    warning(valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation, comment = TRUE))
    return(NaN)
  }
  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }


  if (x <= t) {
    return(0)
  } else if (is.null(truncation)) {
    Result <- alpha_tail * (1 + alpha_ini / alpha_tail * (x / t - 1))^(-alpha_tail - 1) * alpha_ini / alpha_tail / t
    return(Result)
  } else if (x >= truncation) {
    return(0)
  } else {
    Result <- alpha_tail * (1 + alpha_ini / alpha_tail * (x / t - 1))^(-alpha_tail - 1) * alpha_ini / alpha_tail / t / pGenPareto(truncation, t, alpha_ini, alpha_tail)
    return(Result)
  }
}


#' Quantile Function of the generalized Pareto Distribution
#'
#' @description Calculates the quantile function of a generalized Pareto distribution
#'
#' @param p Numeric. The function evaluates the inverse CDF at \code{p}.
#' @param t Numeric. Threshold of the piecewise Pareto distribution.
#' @param alpha_ini Numeric. Initial Pareto alpha.
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the generalized Pareto distribution is truncated at \code{truncation}.
#'
#' @return Quantile function of the Pareto distribution with parameters \code{t}, \code{alpha_ini} and \code{alpha_tail}, evaluated at \code{p}
#'
#' @examples
#' p <- 0:10 * 0.1
#' qGenPareto(p, 1000, 2, 3)
#' qGenPareto(p, 1000, 2, 3, truncation = 5000)
#'
#' @export

qGenPareto <- function(p, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!is.nonnegative.finite.vector(p) || max(p) > 1) {
    warning("p must be a vector with elements in [0,1].")
    return(rep(NaN, length(p)))
  }
  if (!valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation)) {
    warning(valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation, comment = TRUE))
    return(rep(NaN, length(p)))
  }
  sapply(p, FUN = function(y) qGenPareto_s(y, t, alpha_ini, alpha_tail, truncation))
}

qGenPareto_s <- function(y, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation)) {
    warning(valid.parameters.GenPareto(t, alpha_ini, alpha_tail, truncation, comment = TRUE))
    return(NaN)
  }
  if (!is.nonnegative.finite.number(y) || y > 1) {
    warning("y must be a number in [0,1].")
    return(NaN)
  }


  if (y == 1) {
    if (is.null(truncation)) {
      return(Inf)
    } else {
      return(truncation)
    }
  }

  if (!is.null(truncation)) {
    if (!is.infinite(truncation)) {
      scale <- pGenPareto(truncation, t, alpha_ini, alpha_tail)
      y <- y * scale
    }
  }
  result <- t * (1 + alpha_tail / alpha_ini * ((1 - y)^(-1/alpha_tail) - 1))
}







#' Simulation of the generalized Pareto Distribution
#'
#' @description Generates random deviates of a generalized Pareto distribution
#'
#' @param n Numeric. Number of observations.
#' @param t Numeric vector. Thresholds of the generalized Pareto distributions
#' @param alpha_ini Numeric vector. Initial Pareto alphas of the generalized Pareto distributions.
#' @param alpha_tail Numeric vector. Tail Pareto alphas of the generalized Pareto distributions.
#' @param truncation NULL or Numeric vector. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the generalized Pareto distributions are truncated at \code{truncation} (resampled generalized Pareto)
#'
#' @return A vector of \code{n} samples from the (truncated) generalized Pareto distribution with parameters \code{t}, \code{alpha_ini} and \code{alpha_tail}
#'
#' @examples
#' rGenPareto(100, 1000, 2, 3)
#' rGenPareto(100, 1000, 2, 3, truncation = 2000)
#' rGenPareto(100, t = c(1, 10, 100, 1000), alpha_ini = 1, alpha_tail = c(2, 5))
#'
#' @export

rGenPareto <- function(n, t, alpha_ini, alpha_tail, truncation = NULL) {
  if (!is.positive.finite.number(n)) {
    warning("n must be a positive number.")
    return(NaN)
  }
  n <- ceiling(n)
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (n %% length(t) != 0) {
    warning("n is not a multiple of length(t).")
    t <- rep(t, ceiling(n / length(t)))[1:n]
  }
  if (!is.positive.finite.vector(alpha_ini)) {
    warning("alpha_ini must be a positive vector.")
    return(NaN)
  }
  if (n %% length(alpha_ini) != 0) {
    warning("n is not a multiple of length(alpha_ini).")
    alpha_ini <- rep(alpha_ini, ceiling(n / length(alpha_ini)))[1:n]
  }
  if (!is.positive.finite.vector(alpha_tail)) {
    warning("alpha_tail must be a positive vector.")
    return(NaN)
  }
  if (n %% length(alpha_tail) != 0) {
    warning("n is not a multiple of length(alpha_tail).")
    alpha_tail <- rep(alpha_tail, ceiling(n / length(alpha_tail)))[1:n]
  }

  if (!is.null(truncation)) {
    if (!is.positive.vector(truncation)) {
      warning("truncation must be NULL or a positive vector (elements may be 'Inf').")
      return(NaN)
    }
    if (n %% length(truncation) != 0) {
      warning("n is not a multiple of length(truncation).")
      truncation <- rep(truncation, ceiling(n / length(truncation)))[1:n]
    }
    if (sum(truncation <= t) > 0) {
      warning("truncation must be > t.")
      return(NaN)
    }
  }


  FinvGenPareto <- function(x, t, alpha_ini, alpha_tail) {
    return(t * (1 + alpha_tail / alpha_ini * ((1 - x)^(-1 / alpha_tail) - 1)))
  }


  u <- 0
  o <- 1

  if (!is.null(truncation)) {
    o <- 1 - (1 + alpha_ini / alpha_tail * (truncation / t - 1))^(-alpha_tail)
  }

    return(FinvGenPareto(stats::runif(n, u, o), t , alpha_ini, alpha_tail))
}



#' Layer Mean of the generalized Pareto Distribution

#' @description  Calculates the expected loss of a generalized Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha_ini Numeric. Initial Pareto alpha (at \code{t}).
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The expected loss of the (truncated) Pareto distribution with parameters \code{t} and \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' GenPareto_Layer_Mean(4000, 1000, 1, 3)
#' GenPareto_Layer_Mean(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000)
#' GenPareto_Layer_Mean(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 5000)
#' GenPareto_Layer_Mean(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000, truncation = 5000)
#' GenPareto_Layer_Mean(9000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000, truncation = 5000)
#'
#' @export



GenPareto_Layer_Mean <- function(Cover, AttachmentPoint, alpha_ini, alpha_tail, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_ini)) {
    warning("alpha_ini must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_tail)) {
    warning("alpha_tail must be a positive number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point in zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  if (is.infinite(Cover)) {
    if (alpha_tail <= 1) {
      return(Inf)
    } else if (t <= AttachmentPoint) {
      EP <- - t * alpha_tail / (alpha_ini * (1 - alpha_tail)) * (1 + alpha_ini / alpha_tail * (AttachmentPoint / t - 1))^(1 - alpha_tail)
    } else {
      EP <- t - AttachmentPoint
      EP <- EP - t * alpha_tail / (alpha_ini * (1 - alpha_tail))
    }
    return(EP)
  } else {
    # Calculation ignoring truncation
    if (t <= AttachmentPoint) {
      if (alpha_tail == 1) {
        # Pareto: EP <- t * (log(Cover + AttachmentPoint) - log(AttachmentPoint))
        EP <- t * alpha_tail / alpha_ini * (log(1 + alpha_ini / alpha_tail * ((Cover + AttachmentPoint) / t - 1)) - log(1 + alpha_ini / alpha_tail * (AttachmentPoint / t - 1)))
      } else {
        # Pareto: EP <- t / (1 - alpha) * (((Cover + AttachmentPoint) / t)^(1 - alpha) - (AttachmentPoint / t)^(1 - alpha))
        EP <- t * alpha_tail / (alpha_ini * (1 - alpha_tail)) * ((1 + alpha_ini / alpha_tail * ((Cover + AttachmentPoint) / t - 1))^(1 - alpha_tail) - (1 + alpha_ini / alpha_tail * (AttachmentPoint / t - 1))^(1 - alpha_tail))
      }
    } else if (t >= AttachmentPoint + Cover) {
      EP <- Cover
    } else {
      EP <- t - AttachmentPoint
      if (alpha_tail == 1) {
        # Pareto: EP <- EP + t * (log(Cover + AttachmentPoint) - log(t))
        EP <- EP +  t * alpha_tail / alpha_ini * log(1 + alpha_ini / alpha_tail * ((Cover + AttachmentPoint) / t - 1))
      } else {
        # Pareto: EP <- EP + t / (1 - alpha) * (((Cover + AttachmentPoint) / t)^(1 - alpha) - 1)
        EP <- EP + t * alpha_tail / (alpha_ini * (1 - alpha_tail)) * ((1 + alpha_ini / alpha_tail * ((Cover + AttachmentPoint) / t - 1))^(1 - alpha_tail) - 1)
      }
    }

    if (is.positive.finite.number(truncation)) {
      # then Cover + AttachmentPoint <= truncation
      # Adjustment for truncation
      # Pareto: FQ_at_truncation <- (t / truncation)^alpha
      FQ_at_truncation <- 1 - pGenPareto(truncation, t, alpha_ini, alpha_tail)
      EP <- (EP - FQ_at_truncation * Cover) / (1 - FQ_at_truncation)
    }
    return(EP)
  }

}


#' Second Layer Moment of the Generalized Pareto Distribution
#'
#' @description Calculates the second moment of a generalized Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha_ini Numeric. Initial Pareto alpha (at \code{t}).
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return The second moment of the (truncated) generalized Pareto distribution with parameters \code{t}, \code{alpha_ini} and \code{alpha_tail} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' GenPareto_Layer_SM(4000, 1000, 1, 2)
#' GenPareto_Layer_SM(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000)
#' GenPareto_Layer_SM(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 5000)
#' GenPareto_Layer_SM(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000, truncation = 5000)
#' GenPareto_Layer_SM(9000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000, truncation = 5000)
#'
#' @export


GenPareto_Layer_SM <- function(Cover, AttachmentPoint, alpha_ini, alpha_tail, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_ini)) {
    warning("alpha_ini must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_tail)) {
    warning("alpha_tail must be a positive number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point is zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  if (is.infinite(Cover)) {
    if (alpha_tail <= 2) {
      return(Inf)
    }
  }

  G <- function(x) {
    if (alpha_tail != 1) {
      if (x > t) {
        result <- t * alpha_tail / (alpha_ini * ( 1- alpha_tail)) * (1 + alpha_ini / alpha_tail * (x / t - 1))^(1 - alpha_tail)
      } else {
        result <- x - t + t * alpha_tail / (alpha_ini * ( 1 - alpha_tail))
      }
    } else {
      if (x > t) {
        result <- t * alpha_tail / alpha_ini * log(1 + alpha_ini / alpha_tail * (x / t - 1))
      } else {
        result <- x - t
      }
    }
    return(result)
  }
  H <- function(x) {
    if (alpha_tail == 1) {
      if (x > t) {
        result <- t * alpha_tail * x / alpha_ini * log(1 + alpha_ini / alpha_tail * (x / t - 1)) - t^2 * alpha_tail^2 / alpha_ini^2 * ((1 + alpha_ini / alpha_tail * (x / t - 1)) * log(1 + alpha_ini / alpha_tail * (x / t - 1)) - (1 + alpha_ini / alpha_tail * (x / t - 1)))
      } else {
        result <- x^2 / 2 - t^2 / 2 + t^2 * alpha_tail^2 / alpha_ini^2
      }
    } else if (alpha_tail == 2) {
      if (x > t) {
        result <- t * alpha_tail * x / (alpha_ini * (1 - alpha_tail)) * (1 + alpha_ini / alpha_tail * (x / t - 1))^(1 - alpha_tail) - t^2 * alpha_tail^2 / (alpha_ini^2 * (1 - alpha_tail)) * log(1 + alpha_ini / alpha_tail * (x / t - 1))
      } else {
        result <- x^2 / 2 - t^2 / 2 + t^2 * alpha_tail / (alpha_ini * (1 - alpha_tail))
      }
    } else {
      if (x > t) {
        result <- t * alpha_tail * x / (alpha_ini * ( 1 - alpha_tail)) * (1 + alpha_ini / alpha_tail * (x / t - 1))^(1 - alpha_tail) - t^2 * alpha_tail^2 / (alpha_ini^2 * (1 - alpha_tail) * (2 - alpha_tail)) * (1 + alpha_ini / alpha_tail * (x / t - 1))^(2 - alpha_tail)
      } else {
        result <- x^2 / 2 - t^2 / 2 + t^2 * alpha_tail / (alpha_ini * ( 1 - alpha_tail)) - t^2 * alpha_tail^2 / (alpha_ini^2 * (1 - alpha_tail) * (2 - alpha_tail))
      }
    }
    return(result)
  }

  # calculation w/o truncation

  result <- AttachmentPoint^2 * (pGenPareto(Cover + AttachmentPoint, t, alpha_ini, alpha_tail) - pGenPareto(AttachmentPoint, t, alpha_ini, alpha_tail))
  result <- result - 2 * AttachmentPoint * (G(Cover + AttachmentPoint) - G(AttachmentPoint) - (Cover + AttachmentPoint) * (1 - pGenPareto(Cover + AttachmentPoint, t, alpha_ini, alpha_tail)) + AttachmentPoint * (1 - pGenPareto(AttachmentPoint, t, alpha_ini, alpha_tail)))
  result <- result + (2 * (H(Cover + AttachmentPoint) - H(AttachmentPoint)) - (Cover + AttachmentPoint)^2 * (1 - pGenPareto(Cover + AttachmentPoint, t, alpha_ini, alpha_tail)) + AttachmentPoint^2 * (1 - pGenPareto(AttachmentPoint, t, alpha_ini, alpha_tail)))
  result <- result + Cover^2 * (1 - pGenPareto(Cover + AttachmentPoint, t, alpha_ini, alpha_tail))

  # adjustment for truncation

  if (!is.null(truncation)) {
    p <- 1 - pGenPareto(truncation, t, alpha_ini, alpha_tail)
    result <- (result - p * Cover^2) / (1 - p)
  }

  return(result)
}


#' Layer Variance of the Generalized Pareto Distribution
#'
#' @description Calculates the variance of a generalized Pareto distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha_ini Numeric. Initial Pareto alpha (at \code{t}).
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#'
#' @return Variance of the (truncated) generalized Pareto distribution with parameters \code{t}, \code{alpha_ini} and \code{alpha_tail} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' GenPareto_Layer_Var(4000, 1000, 1, 2)
#' GenPareto_Layer_Var(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000)
#' GenPareto_Layer_Var(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 5000)
#' GenPareto_Layer_Var(4000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000, truncation = 5000)
#' GenPareto_Layer_Var(9000, 1000, alpha_ini = 1, alpha_tail = 3, t = 1000, truncation = 5000)
#'
#' @export


GenPareto_Layer_Var <- function(Cover, AttachmentPoint, alpha_ini, alpha_tail, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_ini)) {
    warning("alpha_ini must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_tail)) {
    warning("alpha_tail must be a positive number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point is zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  if (is.infinite(Cover)) {
    if (alpha_tail <= 2) {
      return(Inf)
    }
  }

  result <- GenPareto_Layer_SM(Cover, AttachmentPoint, alpha_ini, alpha_tail, t, truncation) - GenPareto_Layer_Mean(Cover, AttachmentPoint, alpha_ini, alpha_tail, t, truncation)^2

  return(result)
}




