
#' This function the expected loss of the Pareto Distribution Pareto(t, alpha) in the layer Cover xs AttachmentPoint
#' @param Cover Numeric. Cover of the reinsurance layer. Use Inf for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Numeric. Pareto alpha.
#' @param t Numeric.. Threshold of the Pareto distribution. If t = NULL (default) then t <- Attachment Point
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation.
#' @export

Pareto_Layer_Mean <- function(Cover, AttachmentPoint, alpha, t=NULL, truncation = NULL) {
  if (is.null(t)) {
    t <- AttachmentPoint
    if (AttachmentPoint == 0) {
      warning("If Attachment Point in zero, then a t>0 has to be entered.")
      return(NA)
    }
  }
  if (t <= 0) {
    warning("t must be greater than 0.")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (truncation <= AttachmentPoint) {
      warning("truncation must be larger than AttachmentPoint")
      return(NA)
    }
    if (!is.null(t)) {
      if (truncation <= t) {
        warning("truncation must be larger than t")
      }
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  if (Cover == Inf && AttachmentPoint >= 0 && alpha >= 0) {
    if (alpha <= 1) {
      #warning("alpha must be > 1 for unlimited covers!")
      return(Inf)
    } else if (t <= AttachmentPoint) {
      EP <- -(t / AttachmentPoint)^alpha / (1 - alpha) * AttachmentPoint
    } else {
      EP <- t - AttachmentPoint
      EP <- EP - t / (1 - alpha)
    }
    return(EP)
  } else if (Cover >= 0 && AttachmentPoint >= 0 && alpha >= 0) {
    EP <- NaN
    while (is.nan(EP)) {
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
      if (is.nan(EP) || is.infinite(EP)) {
        alpha <- alpha / 2
        EP <- NaN
        warning(paste("alpha reduced to", round(alpha, 2)))
      }
    }
    if (!is.null(truncation)) {
      # then Cover + AttachmentPoint <= truncation
      if (truncation < Inf) {
        FQ_at_truncation <- (t / truncation)^alpha
        EP <- (EP - FQ_at_truncation * Cover) / (1 - FQ_at_truncation)
      }
    }
    return(EP)
  } else {
    warning("Cover, Attachment Point and alpha must be >= 0; NA produced")
    return(NA)
  }

}


#' This function generates random deviates of the Pareto distribution Pareto(t, alpha).
#' @param n Number of observations.
#' @param t Threshold of the Pareto distribution
#' @param alpha Pareto alpha.
#' @param truncation If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation (resampled Pareto)
#' @export

rPareto <- function(n, t, alpha, truncation = NULL) {
  FinvPareto <- function(x,t,alpha) {
    return(t/(1-x)^(1/alpha))
  }
  u <- 0
  o <- 1
  if (!is.null(truncation)) {
    if (is.numeric(truncation)) {
      if (truncation > t) {
        o <- 1 - (t / truncation)^alpha
      }
    }
  }

  return(FinvPareto(runif(n, u, o),t,alpha))
}


#' This function extrapolates from the layer Cover_1 xs AttachmentPoint_1 to the layer Cover_2 xs AttachmentPoint_2 using a Pareto distribution Pareto(t, alpha) with t sufficiently small.
#' @param Cover_1 Numeric. Cover of the layer from which we extrapolate. Use Inf for unlimited layers.
#' @param AttachmentPoint_1 Numeric. Attachment point of the layer from which we extrapolate.
#' @param Cover_2 Numeric. Cover of the layer to which we extrapolate. Use Inf for unlimited layers.
#' @param AttachmentPoint_2 Numeric. Attachment point of the layer to which we extrapolate.
#' @param alpha Numeric. Pareto alpha used for the extrapolation.
#' @param ExpLoss_1 Numeric. Expected loss of the layer from which we extrapolate. If NULL (default) then the function provides only the ratio between the expected losses of the layers.
#' @param truncation Numeric. If truncation is not NULL and truncation > AttachmentPoint_1, then the Pareto distribution is truncated at truncation.
#' @export

Pareto_Extrapolation <- function(Cover_1, AttachmentPoint_1, Cover_2, AttachmentPoint_2, alpha, ExpLoss_1 = NULL, truncation = NULL) {
  if (is.null(ExpLoss_1)) {
    ExpLoss_1 <- 1
  }
  # if (is.null(truncation)) {
  #   if (Cover_1 == Inf|| Cover_2 == Inf) {
  #     if (alpha <= 1) {
  #       warning("alpha must be > 1 for unlimited covers!")
  #       return(NA)
  #     } else if (Cover_1 == Inf && Cover_2 == Inf) {
  #       Result <- ((AttachmentPoint_2)^(1-alpha)) / ((AttachmentPoint_1)^(1-alpha)) * ExpLoss_1
  #     } else if (Cover_2 == Inf) {
  #       Result <- ( - (AttachmentPoint_2)^(1-alpha)) / ((Cover_1 + AttachmentPoint_1)^(1-alpha) - (AttachmentPoint_1)^(1-alpha)) * ExpLoss_1
  #     } else {
  #       Result <- ((Cover_2 + AttachmentPoint_2)^(1-alpha) - (AttachmentPoint_2)^(1-alpha)) / ( - (AttachmentPoint_1)^(1-alpha)) * ExpLoss_1
  #     }
  #     return(Result)
  #   }
  #   if (alpha <= 0) {
  #     Result <- Cover_2 / Cover_1 * ExpLoss_1
  #   } else if (alpha == 1) {
  #     Result <- (log(Cover_2 + AttachmentPoint_2) - log(AttachmentPoint_2)) / (log(Cover_1 + AttachmentPoint_1) - log(AttachmentPoint_1)) * ExpLoss_1
  #   } else {
  #     Result <- ((Cover_2 + AttachmentPoint_2)^(1-alpha) - (AttachmentPoint_2)^(1-alpha)) / ((Cover_1 + AttachmentPoint_1)^(1-alpha) - (AttachmentPoint_1)^(1-alpha)) * ExpLoss_1
  #   }
  #   return(Result)
  # } else {

  if (!is.null(truncation)) {
    if (truncation <= AttachmentPoint_1) {
      warning("truncation must be greater than AttachmentPoint_1")
      return(NA)
    }
  }
  Smaller_AP <- min(AttachmentPoint_1, AttachmentPoint_2)
  Result <- ExpLoss_1 * Pareto_Layer_Mean(Cover_2, AttachmentPoint_2, alpha, t = Smaller_AP, truncation = truncation) / Pareto_Layer_Mean(Cover_1, AttachmentPoint_1, alpha, t = Smaller_AP, truncation = truncation)
  return(Result)
  # }

}



#' This function finds the Pareto alpha between two layers.
#' @param Cover_1 Numeric. Cover of the first layer.
#' @param AttachmentPoint_1 Numeric. Cover of the first layer.
#' @param ExpLoss_1 Numeric. Expected loss of the first layer.
#' @param Cover_2 Numeric. Cover of the second layer.
#' @param AttachmentPoint_2 Numeric. Cover of the second layer.
#' @param ExpLoss_2 Numeric. Expected loss of the second layer.
#' @param max_alpha Numeric. Upper limit for the alpha that is returned.
#' @param tolerance Numeric. Accuracy of the result.
#' @param truncation Numeric. If truncation is not NULL then the Pareto distribution is truncated at truncation.
#' @export

Pareto_Find_Alpha_btw_Layers <- function(Cover_1, AttachmentPoint_1, ExpLoss_1, Cover_2, AttachmentPoint_2, ExpLoss_2, max_alpha = 20, tolerance = 10^(-10), truncation = NULL) {
  if (Cover_1 <= 0 || AttachmentPoint_1 <= 0 || ExpLoss_1 <= 0 || Cover_2 <= 0 || AttachmentPoint_2 <= 0 || ExpLoss_2 <= 0) {
    warning("All input parameters must be positive!")
    return(NA)
  }
  min_alpha <- 0
  if (!is.null(truncation)) {
    if (truncation <= max(AttachmentPoint_1, AttachmentPoint_2)) {
      warning("tuncation must be NULL or greater than both attachment points!")
      return(NA)
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
    try(Result <- uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = T)
    if (AttachmentPoint_1 > AttachmentPoint_2 && AttachmentPoint_1 + Cover_1 >= AttachmentPoint_2 + Cover_2 && f(max_alpha) < 0) {
      Result <- max_alpha
    } else if (AttachmentPoint_1 < AttachmentPoint_2 && AttachmentPoint_1 + Cover_1 <= AttachmentPoint_2 + Cover_2 && f(max_alpha) > 0) {
      Result <- max_alpha
    }
  } else {
    try(Result <- uniroot(f, c(1 + tolerance, max_alpha), tol = tolerance)$root, silent = T)
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
    return(NA)
  }
  alpha <- Result
  return(alpha)

}


#' This function finds the Pareto alpha between an excess frequency and a layer.
#' @param Threshold Numeric. Threshold
#' @param Frequency Numeric. Expected frequency in excess of Thershold
#' @param Cover Numeric. Cover of the second layer.
#' @param AttachmentPoint Numeric. Cover of the layer.
#' @param ExpLoss Numeric. Expected loss of the layer.
#' @param max_alpha Numeric. Upper limit for the alpha that is returned.
#' @param tolerance Numeric. Accuracy of the result.
#' @export

Pareto_Find_Alpha_btw_FQ_Layer <- function(Threshold, Frequency, Cover, AttachmentPoint, ExpLoss, max_alpha = 20, tolerance = 10^(-10)) {
  if (Threshold <= 0 || Frequency <= 0 || Cover <= 0 || AttachmentPoint <= 0 || ExpLoss <= 0) {
    warning("All input parameters must be positive!")
    return(NA)
  }
  if (Threshold > AttachmentPoint && Threshold < AttachmentPoint + Cover) {
    warning("Threshold must be <= AttachmentPoint or >= Cover + AttachmentPoint")
    return(NA)
  }

  f <- function(alpha) {
    Pareto_Layer_Mean(Cover, AttachmentPoint, alpha) * (Threshold / AttachmentPoint)^alpha * Frequency - ExpLoss
  }

  Result <- NA
  if (is.numeric(Cover)) {
    try(Result <- uniroot(f, c(0, max_alpha), tol = tolerance)$root, silent = T)
    if (AttachmentPoint >= Threshold && f(max_alpha) > 0) {
      Result <- max_alpha
    } else if (AttachmentPoint + Cover <= Threshold && f(max_alpha) < 0) {
      Result <- max_alpha
    }
  } else {
    try(Result <- uniroot(f, c(1 + tolerance, max_alpha), tol = tolerance)$root, silent = T)
    if (is.na(Result)) {
      if (Cover == Inf) {
        if (AttachmentPoint >= Threshold && f(max_alpha) < 0) {
          Result <- max_alpha
        }
      }
    }
  }
  if (is.na(Result)) {
    warning("Did not find a solution!")
    return(NA)
  }
  alpha <- Result
  return(alpha)

}





#' This function calculates the expected loss of a Piecewise Pareto Distribution in the layer Cover xs AttachmentPoint
#' @param Cover Numeric. Cover of the reinsurance layer.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. Pareto alpha[i] = Pareto alpha in excess of t[i].
#' @export

PiecewisePareto_Layer_Mean <- function(Cover, AttachmentPoint, t, alpha) {
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (min(t) <= 0) {
    waring("t must have positive elements!")
    return(NA)
  }
  if (min(alpha) < 0) {
    warning("alpha must have non-negative elements!")
    return(NA)
  }
  if (min(diff(t)) <= 0) {
    warning("t must be strictily ascending!")
    return(NA)
  }
  if (length(AttachmentPoint) != 1 || length(Cover) != 1) {
    warning("AttachmentPoint and Cover must have lenght  1!")
    return(NA)
  }
  if (!is.numeric(AttachmentPoint)) {
    warning("AttachmentPoint must be numeric!")
    return(NA)
  }
  if (AttachmentPoint < 0) {
    warning("AttachmentPoint must be non-negative!")
    return(NA)
  }
  if (!is.numeric(Cover)) {
    waring("AttachmentPoint must be numeric!")
    return(NA)
  } else if (Cover < 0) {
    warning("Cover must be non-negative!")
    return(NA)
  }
  if (Cover == 0) {
    return(0)
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
      Result <- Result + Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
    }
  } else if (Cover == Inf) {
    if (k1 < k2) {
      for (i in k1:(k2-1)) {
        Result <- Result + Pareto_Layer_Mean(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
      }
    }
    Result <- Result + Pareto_Layer_Mean(Inf, Att[n], alpha[n], t[n]) * excess_prob[n]
  }

  return(Result)
}





#' This function generates random deviates of the Piecewise Pareto distribution.
#' @param n Numeric. Number of simulations
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. Pareto alpha[i] = Pareto alpha in excess of t[i].
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation.
#' @param truncation_type Charakter. If truncation_type = "wd" then the whole distribution is truncated. If truncation_type = "lp" then a truncated Pareto is used for the last piece.
#' @export

rPiecewisePareto <- function(n, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  k <- length(t)
  if (min(t) <= 0) {
    warning("t must have positive elements!")
    return(NA)
  }
  if (min(alpha) < 0) {
    warning("alpha must have non-negative elements!")
    return(NA)
  }
  if (min(diff(t)) <= 0) {
    warning("t must be strictily ascending!")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (!is.numeric(truncation)) {
      warning("truncation must be NULL or numeric")
      return(NA)
    }
    if (truncation <= t[k]) {
      warning("truncation must be greater than max(t)")
      return(NA)
    }
    if (truncation_type != "wd" && truncation_type != "lp") {
      warning("truncation_type must be wd or lp")
      return(NA)
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
    Result[Simulated_Pieces == i] <- FinvPareto(runif(NumberOfSimulations_for_Pieces[i], min = 0, max = CDF), t[i], alpha[i])
  }

  return(Result)
}




#' This function matches the expected losses of consecutive layers using a piecewise Pareto severity
#' @param Attachment_Points Numeric vector. Vector containing the attachment points of consecutive layers in increasing order
#' @param Expected_Layer_Losses Numeric vector. Vector containing the expected losses of layers xs the attachment points.
#' @param Unlimited_Layers Logical. If true, then Expected_Layer_Losses[i] contains the expected loss of unlimited xs Attachment_Points[i].
#'        If FALSE then Expected_Layer_Losses[i] contains the expected loss of the layers Attachment_Points[i+1] xs Attachment_Points[i]
#' @param Frequencies. Numeric vector. Expected frequencies excess the attachment points. If NULL then the function calculates frequencies.
#' @param FQ_at_lowest_AttPt. Numerical. Expected frequency excess Attachment_Points[1]
#' @param FQ_at_highest_AttPt. Numerical. Expected frequency excess Attachment_Points[k]
#' @param TotalLoss_Frequencies. Numeric vector. TotalLoss_Frequencies[i] is the frequency of total losses to layer i (i.e. Attachment_Points[i+1]-Attachment_Points[i] xs Attachment_Points[i])
#'        TotalLoss_Frequencies[i] is the frequency for losses >= Attachment_Points[i+1], whereas Frequencies[i] is the frequency of losses > Attachment_Points[i].
#'        TotalLoss_Frequencies[i] > Frequencies[i+1] means that there is a point mass of the severity at Attachment_Points[i+1].
#' @param minimize_ratios. Logical. If TRUE then ratios between alphas are minimized.
#' @param Use_unlimited_Layer_for_FQ. Logical. Only relevant if no frequency is provided for the highest attachment point by the user. If TRUE then the frequency is calculated using the Pareto alpha between the last two layers.
#' @param alpha_max. Numerical. Maximum alpha to be used for the matching.
#'
#' @return t. Numeric vector. Vector containing the thresholds for the piecewise Pareto distribution.
#' @return alpha. Numeric vector. Vector containing the Pareto alphas of the piecewise Pareto distribution.
#' @return FQ. Numerical. Frequency in excess of the lowest threshold of the piecewise Pareto distribution.
#' @export

PiecewisePareto_Match_Layer_Losses <- function(Attachment_Points, Expected_Layer_Losses, Unlimited_Layers = FALSE, Frequencies = NULL, FQ_at_lowest_AttPt = NULL, FQ_at_highest_AttPt = NULL, TotalLoss_Frequencies = NULL, minimize_ratios = TRUE, Use_unlimited_Layer_for_FQ = TRUE, tolerance = 10^(-10), alpha_max = 20) {
  if (!is.numeric(Attachment_Points)) {
    stop("Attachment_Points must be numeric.")
  }
  if (!is.numeric(Expected_Layer_Losses)) {
    stop("Expected_Layer_Losses must be numeric.")
  }
  if (min(Expected_Layer_Losses) <= 0) {
    stop("Expected_Layer_Losses must be positive.")
  }
  k <-length(Attachment_Points)
  if (k<2) {
    stop("Attachment_Points must have lenght >= 2.")
  }
  if (length(Expected_Layer_Losses) != k) {
    stop("Attachment_Points and Expected_Layer_Losses must have the same lenght.")
  }
  if (min(diff(Attachment_Points)) <= 0) {
    stop("Attachment_Points must be increasing.")
  }
  if (min(Attachment_Points) <= 0) {
    stop("Attachment_Points must be positive.")
  }
  if (!is.logical(Unlimited_Layers)) {
    stop("Unlimited_Layers must be locigal.")
  }
  if (length(Unlimited_Layers) != 1) {
    stop("Unlimited_Layers must have length 1.")
  }
  if (!is.null(Frequencies)) {
    if (!is.numeric(Frequencies)) {
      stop("Frequencies must be numeric or NULL.")
    }
    if (length(Frequencies) != k) {
      stop("Attachment_Points and Frequencies must have the same lenght.")
    }
  }
  if (!is.null(TotalLoss_Frequencies)) {
    if (!is.numeric(TotalLoss_Frequencies)) {
      stop("TotalLoss_Frequencies must be numeric or NULL.")
    }
    if (length(TotalLoss_Frequencies) != (k-1)) {
      stop("TotalLoss_Frequencies must have lenght of Frequencies - 1.")
    }
    if (is.null(Frequencies)) {
      stop("TotalLoss_Frequencies must be NULL if Frequencies is NULL.")
    }
  }
  if (!is.null(FQ_at_lowest_AttPt)) {
    if (!is.numeric(FQ_at_lowest_AttPt)) {
      stop("FQ_at_lowest_AttPt must be numeric or NULL.")
    }
    if (length(FQ_at_lowest_AttPt) != 1) {
      stop("FQ_at_lowest_AttPt must have lenght 1.")
    }
  }
  if (!is.null(FQ_at_highest_AttPt)) {
    if (!is.numeric(FQ_at_highest_AttPt)) {
      stop("FQ_at_highest_AttPt must be numeric or NULL.")
    }
    if (length(FQ_at_highest_AttPt) != 1) {
      stop("FQ_at_highest_AttPt must have lenght 1.")
    }
  }
  if (!is.logical(minimize_ratios)) {
    stop("minimize_ratios must be locigal.")
  }
  if (length(minimize_ratios) != 1) {
    stop("minimize_ratios must have length 1.")
  }
  if (!is.logical(Use_unlimited_Layer_for_FQ)) {
    stop("Use_unlimited_Layer_for_FQ must be locigal.")
  }
  if (length(Use_unlimited_Layer_for_FQ) != 1) {
    stop("Use_unlimited_Layer_for_FQ must have length 1.")
  }
  if (!is.numeric(tolerance)) {
    stop("tolerance must be numeric.")
  }
  if (length(tolerance) != 1) {
    stop("tolerance must have length 1.")
  }
  if (tolerance <= 0) {
    stop("tolerance must be positive.")
  }
  if (!is.numeric(alpha_max)) {
    stop("alpha_max must be numeric.")
  }
  if (length(alpha_max) != 1) {
    stop("alpha_max must have length 1.")
  }
  if (alpha_max <= 0) {
    stop("alpha_max must be positive.")
  }

  Status <- ""

  if (Unlimited_Layers) {
    ELL <- Expected_Layer_Losses[1:(k-1)] - Expected_Layer_Losses[2:k]
    ELL <- c(ELL, Expected_Layer_Losses[k])
  } else {
    ELL <- Expected_Layer_Losses
  }
  Limits <- Attachment_Points[2:k] - Attachment_Points[1:(k-1)]
  Limits <- c(Limits, Inf)
  RoLs <- ELL / Limits
  if (max(diff(RoLs)) >= 0) {
    repeat {
      if (k<3) {
        stop("Layers cannot be merged to obtain strictly decreasing RoLs.")
      }
      pos <- order(diff(RoLs))[k-1]
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
      k <- k-1
      RoLs <- ELL / Limits
      if (max(diff(RoLs)) < 0) {break}
    }
    Status <- paste0(Status, "RoLs not strictly decreasing. Layers have been merged. ")
  }
  if (!is.null(Frequencies)) {
    if (min(Frequencies - RoLs) <= 0) {
      Frequencies <- NULL
      TotalLoss_Frequencies <- NULL
      Status <- paste0(Status, "Layer entry frequencies not strictly greater than RoLs. Frequencies not used! ")
    }
    if (max(Frequencies[2:k] - RoLs[1:(k-1)]) >= 0) {
      Frequencies <- NULL
      TotalLoss_Frequencies <- NULL
      Status <- paste0(Status, "Layer exit frequencies not strictly less than RoLs. Frequencies not used! ")
    }
  }
  if (is.null(Frequencies)) {
    alpha_between_layers <- numeric(k-1)
    for (i in 1:(k-1)) {
      if (i < k-1) {
        alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Limits[i+1], Attachment_Points[i+1], ELL[i+1])
        Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Limits[i+1], Attachment_Points[i+1], alpha_between_layers[i])
      } else {
        if (Use_unlimited_Layer_for_FQ) {
          alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Inf, Attachment_Points[i+1], ELL[i+1])
          Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Inf, Attachment_Points[i+1], alpha_between_layers[i])
        } else {
          Frequencies[i+1] = Frequencies[i] * (Attachment_Points[i]/Attachment_Points[i+1])^alpha_between_layers[i-1]
        }
      }
    }
    Frequencies[1] <- ELL[1] / Pareto_Layer_Mean(Limits[1], Attachment_Points[1], alpha_between_layers[1])
  }
  if (!is.null(FQ_at_lowest_AttPt)) {
    if (FQ_at_lowest_AttPt > RoLs[1]) {
      Frequencies[1] <- FQ_at_lowest_AttPt
    } else {
      Status <- paste0(Status, "FQ_at_lowest_AttPt too small. Not used! ")
    }
  }
  if (!is.null(FQ_at_highest_AttPt)) {
    if (FQ_at_highest_AttPt < RoLs[k-1]) {
      Frequencies[k] <- FQ_at_highest_AttPt
    } else {
      Status <- paste0(Status, "FQ_at_highest_AttPt too large. Not used! ")
    }
  }
  if (!is.null(TotalLoss_Frequencies)) {
    if (max(TotalLoss_Frequencies - RoLs[1:(k-1)]) >= 0) {
      TotalLoss_Frequencies <- NULL
      Status <- paste0(Status, "TotalLoss_Frequencies not strictly less than RoLs. TotalLoss_Frequencies not used! ")
    }
    if (max(TotalLoss_Frequencies - Frequencies[2:k]) < 0) {
      TotalLoss_Frequencies <- NULL
      Status <- paste0(Status, "TotalLoss_Frequencies not greater than or equal to Frequencies. TotalLoss_Frequencies not used! ")
    }
  }


  if (!is.null(TotalLoss_Frequencies)) {
    if (max(TotalLoss_Frequencies - Frequencies[2:k]) > 0) {
      repeat {
        pos <- order(TotalLoss_Frequencies - Frequencies[2:k])[k-1]
        new_AttPoint <- Attachment_Points[pos+1] * (Frequencies[pos+1]/TotalLoss_Frequencies[pos])^(1/alpha_max)
        if (new_AttPoint <= Attachment_Points[pos]) {
          # not possible
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





  s <- Frequencies / Frequencies[1]
  l <- ELL / Frequencies[1]
  Results <- Fit_PP(Attachment_Points, s, l, alpha_max = alpha_max, minimize_ratios = minimize_ratios)

  Results$FQ <- Frequencies[1]
  if (Status == "") {Status <- "OK."}
  Results$Status <- Status
  return(Results)
}



#' This function calculates the cumulative distribution function of a Pareto Distribution
#' @param x Numeric. The function evaluates the CDF at x.
#' @param t Numeric. Threshold of the piecewise Pareto distribution.
#' @param alpha Numeric. Pareto alpha.
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation.
#' @export

Pareto_CDF <- function(x, t, alpha, truncation = NULL) {
  if (!is.numeric(t) || !is.numeric(alpha) || !is.numeric(x)) {
    waring("x, t and alpha must be numeric.")
    return(NA)
  }
  if (length(t) != 1 || length(alpha) != 1 || length(x) != 1) {
    warning("x, t and alpha must have length 1")
    return(NA)
  }
  if (t <= 0) {
    warning("t must be positive.")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (truncation <= t) {
      warning("truncation must be NULL or greater that t.")
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




#' This function calculates the cumulative distribution function of a Piecewise Pareto Distribution
#' @param x Numeric. The function evaluates the CDF at x.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. Pareto alpha[i] = Pareto alpha in excess of t[i].
#' @param truncation Numeric. If truncation is not NULL and truncation > t, then the Pareto distribution is truncated at truncation.
#' @param truncation_type Charakter. If truncation_type = "wd" then the whole distribution is truncated. If truncation_type = "lp" then a truncated Pareto is used for the last piece.
#' @export

PiecewisePareto_CDF <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.numeric(t) || !is.numeric(alpha)) {
    waring("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (n == 1) {
    Result <- Pareto_CDF(x, t, alpha, truncation)
    return(Result)
  }
  if (min(t) <= 0) {
    warning("t must have positive elements!")
    return(NA)
  }
  if (min(alpha) < 0) {
    warning("alpha must have non-negative elements!")
    return(NA)
  }
  if (min(diff(t)) <= 0) {
    warning("t must be strictily ascending!")
    return(NA)
  }
  if (length(x) != 1) {
    warning("x must have lenght  1!")
    return(NA)
  }
  if (!is.numeric(x)) {
    warning("x must be numeric!")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (!is.numeric(truncation)) {
      warning("truncation must be NULL or numeric")
      return(NA)
    }
    if (truncation <= t[n]) {
      warning("truncation must be greater than max(t)")
      return(NA)
    }
    if (truncation_type != "wd" && truncation_type != "lp") {
      warning("truncation_type must be wd or lp")
      return(NA)
    }
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


