
#' Layer Mean of the Pareto Distribution
#'
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
    if (truncation <= t) {
      warning("truncation must be larger than t")
    }
    if (truncation <= AttachmentPoint) {
      return(0)
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





# #' Calculates the second moment of the xs loss of a Pareto distribution in a reinsurance layer
# #'
# #' Not visible to the user. Used in Pareto_Layer_Var.
# #' @param Cover Numeric. Cover of the reinsurance layer. Use Inf for unlimited layers.
# #' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
# #' @param alpha Numeric. Pareto alpha.


Pareto_Layer_Second_Moment_simple <- function(Cover, AttachmentPoint, alpha) {
  if (AttachmentPoint <= 0) {
    return(NA)
  }
  if (Cover < 0) {
    return(NA)
  }
  if (Cover == 0) {
    return(0)
  }
  if (alpha < 0) {
    return(NA)
  }

  if (Cover == Inf) {
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

Pareto_Layer_Var <- function(Cover, AttachmentPoint, alpha, t=NULL, truncation = NULL) {
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
    if (truncation <= t) {
      warning("truncation must be larger than t")
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
  if (!is.null(truncation)) {
    # probability of truncation if t = AttachmentPoint
    p <- 1- pPareto(truncation, AttachmentPoint, alpha)
    # consider truncation in second moment
    SM <- (SM - p * Cover^2) / (1 - p)
  }
  # consider thresholds t < AttachmentPoint
  p <- 1 - pPareto(AttachmentPoint, t, alpha, truncation = truncation)
  SM <- p * SM

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
#' @param n Number of observations.
#' @param t Threshold of the Pareto distribution
#' @param alpha Pareto alpha.
#' @param truncation If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation} (resampled Pareto)
#'
#' @return A vector of \code{n} samples from the (truncated) Pareto distribution with parameters \code{t} and \code{alpha}
#'
#' @examples
#' rPareto(100, 1000, 2)
#' rPareto(100, 1000, 2, truncation = 2000)
#'
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

  return(FinvPareto(stats::runif(n, u, o),t,alpha))
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
  if (is.null(ExpLoss_1)) {
    ExpLoss_1 <- 1
  }

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
    try(Result <- stats::uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = T)
    if (AttachmentPoint_1 > AttachmentPoint_2 && AttachmentPoint_1 + Cover_1 >= AttachmentPoint_2 + Cover_2 && f(max_alpha) < 0) {
      Result <- max_alpha
    } else if (AttachmentPoint_1 < AttachmentPoint_2 && AttachmentPoint_1 + Cover_1 <= AttachmentPoint_2 + Cover_2 && f(max_alpha) > 0) {
      Result <- max_alpha
    }
  } else {
    try(Result <- stats::uniroot(f, c(1 + tolerance, max_alpha), tol = tolerance)$root, silent = T)
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
  if (Threshold <= 0 || Frequency <= 0 || Cover <= 0 || AttachmentPoint <= 0 || ExpLoss <= 0) {
    warning("All input parameters must be positive!")
    return(NA)
  }
  if (Threshold > AttachmentPoint && Threshold < AttachmentPoint + Cover) {
    warning("Threshold must be <= AttachmentPoint or >= Cover + AttachmentPoint")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (truncation <= AttachmentPoint || truncation <= Threshold) {
      warning("Threshold and AttachmentPoint must be less than truncation")
      return(NA)
    }
    Cover <- min(Cover, truncation - AttachmentPoint)
  }

  f <- function(alpha) {
  #  Pareto_Layer_Mean(Cover, AttachmentPoint, alpha) * (Threshold / AttachmentPoint)^alpha * Frequency - ExpLoss
    if (AttachmentPoint < Threshold) {
      FQ_Factor <- 1 / (1 - pPareto(Threshold, AttachmentPoint, alpha, truncation = truncation))
    } else {
      FQ_Factor <- 1 - pPareto(AttachmentPoint, Threshold, alpha, truncation = truncation)
    }

    Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, truncation = truncation) * FQ_Factor * Frequency - ExpLoss
  }

  Result <- NA
  if (is.infinite(Cover) || !is.null(truncation)) {
    min_alpha <- tolerance
  } else {
    min_alpha <- 0
  }
  try(Result <- stats::uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = T)
  if (AttachmentPoint >= Threshold && f(max_alpha) > 0) {
    Result <- max_alpha
  } else if (AttachmentPoint + Cover <= Threshold && f(max_alpha) < 0) {
    Result <- max_alpha
  }


  if (is.na(Result)) {
    warning("Did not find a solution!")
    return(NA)
  }
  alpha <- Result
  return(alpha)

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
  if (Threshold_1 <= 0 || Frequency_1 <= 0 || Threshold_2 <= 0 || Frequency_2 <= 0) {
    warning("All input parameters must be positive!")
    return(NA)
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
    return(NA)
  }
  if (Frequency_2 > Frequency_1) {
    warning("Frequency of larger threshold must be less than or equal to frequency at lower threshold")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (truncation <= Threshold_2) {
      warning("Thresholds must be less than truncation")
      return(NA)
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
      try(Result <- stats::uniroot(f, c(min_alpha, max_alpha), tol = tolerance)$root, silent = T)
    }
    if (is.na(Result)) {
      warning("Did not find a solution!")
      return(NA)
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (n == 1) {
    return(Pareto_Layer_Mean(Cover, AttachmentPoint, alpha, t, truncation))
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
    warning("AttachmentPoint must be numeric!")
    return(NA)
  } else if (Cover < 0) {
    warning("Cover must be non-negative!")
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
  } else if (Cover == Inf) {
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (n == 1) {
    return(Pareto_Layer_SM(Cover, AttachmentPoint, alpha, t, truncation))
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
    warning("AttachmentPoint must be numeric!")
    return(NA)
  } else if (Cover < 0) {
    warning("Cover must be non-negative!")
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
#' @param scale_pieces Numeric vector. If not \code{NULL} then the density of the i-th Pareto piece (on the Intervall (\code{t[i], t[i+1])}) is scaled with the factor \code{const * scale_pieces[i]} (where \code{const} is a normalization constant)
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
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  k <- length(t)
  if (k == 1) {
    return(rPareto(n, t, alpha, truncation))
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
    if (!is.null(scale_pieces)) {
      warning("either truncation or scale_pieces must be NULL")
      return(NA)
    }
  }
  if (!is.null(scale_pieces)) {
    if (!is.numeric(scale_pieces)) {
      warning("scale_pieces must be NULL or numeric")
      return(NA)
    }
    if (length(scale_pieces) != length(t)) {
      warning("t and scale_pieces must have the same length")
      return(NA)
    }
    if (min(scale_pieces) < 0) {
      warning("all entries of scale_pieces must be non-negative")
      return(NA)
    }
    if (sum(scale_pieces) <= 0) {
      warning("scale_pieces must have a positive entry")
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
#' @param Frequencies Numeric vector. Expected frequencies excess the attachment points. If \code{NULL} then the function calculates frequencies.
#' @param FQ_at_lowest_AttPt Numerical. Expected frequency excess \code{Attachment_Points[1]}
#' @param FQ_at_highest_AttPt Numerical. Expected frequency excess \code{Attachment_Points[k]}
#' @param TotalLoss_Frequencies Numeric vector. \code{TotalLoss_Frequencies[i]} is the frequency of total losses to layer \code{i} (i.e. \code{Attachment_Points[i+1] - Attachment_Points[i]} xs \code{Attachment_Points[i]}).    \code{TotalLoss_Frequencies[i]} is the frequency for losses larger than or equal to \code{Attachment_Points[i+1]}, whereas \code{Frequencies[i]} is the frequency of losses larger than \code{Attachment_Points[i]}.    \code{TotalLoss_Frequencies[i] > Frequencies[i+1]} means that there is a point mass of the severity at \code{Attachment_Points[i+1]}.
#' @param minimize_ratios Logical. If \code{TRUE} then ratios between alphas are minimized.
#' @param Use_unlimited_Layer_for_FQ Logical. Only relevant if no frequency is provided for the highest attachment point by the user. If \code{TRUE} then the frequency is calculated using the Pareto alpha between the last two layers.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > max(Attachment_Points)}, then the last Pareto piece is truncated at \code{truncation} (\code{truncation_type = "lp"}).
#' @param tolerance Numeric. Numerical tolerance.
#' @param alpha_max Numerical. Maximum alpha to be used for the matching.
#' @param merge_tolerance Numerical. Consecutive Pareto pieces are merged if the alphas deviate by less than merge_tolerance.
#' @param RoL_tolerance Numerical. Consecutive layers are merged if RoL decreases less than factor \code{1 - RoL_tolerance}.

#' @return A list containing the following objects: \itemize{
#' \item \code{t} Numeric vector. Vector containing the thresholds for the piecewise Pareto distribution
#' \item \code{alpha} Numeric vector. Vector containing the Pareto alphas of the piecewise Pareto distribution
#' \item \code{Status} Character. Information on whether the fit was succesful
#' \item \code{FQ} Numerical. Frequency in excess of the lowest threshold of the piecewise Pareto distribution
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

PiecewisePareto_Match_Layer_Losses <- function(Attachment_Points, Expected_Layer_Losses, Unlimited_Layers = FALSE, Frequencies = NULL, FQ_at_lowest_AttPt = NULL, FQ_at_highest_AttPt = NULL, TotalLoss_Frequencies = NULL, minimize_ratios = TRUE, Use_unlimited_Layer_for_FQ = TRUE, truncation = NULL, tolerance = 1e-10, alpha_max = 100, merge_tolerance = 1e-6, RoL_tolerance = 1e-6) {
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
  if (!is.null(truncation)) {
    if (!is.numeric(truncation)) {
      warning("truncation must be NULL or numeric")
      return(NA)
    }
    if (truncation <= max(Attachment_Points)) {
      warning("truncation must be greater than max(Attachment_Points)")
      return(NA)
    }
    if (is.infinite(truncation)) {
      truncation <- NULL
    }
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
  Merged_Layer <- rep(FALSE, k)
  if (max(RoLs[2:k] / RoLs[1:(k-1)]) >= 1 - RoL_tolerance) {
    repeat {
      if (k<3) {
        warning("Layers cannot be merged to obtain strictly decreasing RoLs.")
        return(NA)
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
    Status <- paste0(Status, "RoLs not strictly decreasing. Layers have been merged. ")
  }
  if (!is.null(Frequencies)) {
    if (max(RoLs/Frequencies) >= 1 - RoL_tolerance / 2) {
      Frequencies <- NULL
      TotalLoss_Frequencies <- NULL
      Status <- paste0(Status, "Layer entry frequencies not strictly greater than RoLs. Frequencies not used! ")
    }
    if (min(RoLs[1:(k-1)]/Frequencies[2:k]) <= 1 + RoL_tolerance / 2) {
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
        if (Merged_Layer[i] & !Merged_Layer[i+1]) {
          #if (RoLs[i] * (1 - RoL_tolerance) > RoLs[i+1]) {
            Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
          #}
        } else if (!Merged_Layer[i] & Merged_Layer[i+1]) {
          #if (RoLs[i+1] * (1 + RoL_tolerance) < RoLs[i]) {
            Frequencies[i+1] <- RoLs[i+1] * (1 + RoL_tolerance / 2)
          #}
        }
      } else {
        if (Use_unlimited_Layer_for_FQ) {
          if (is.null(truncation)) {
            alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Inf, Attachment_Points[i+1], ELL[i+1])
            Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Inf, Attachment_Points[i+1], alpha_between_layers[i])
            if (Merged_Layer[i]) {
              Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
            }
          } else {
            suppressWarnings(alpha_between_layers[i] <-  Pareto_Find_Alpha_btw_Layers(Limits[i], Attachment_Points[i], ELL[i], Inf, Attachment_Points[i+1], ELL[i+1], truncation = truncation))
            if (!is.na(alpha_between_layers[i])) {
              Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Inf, Attachment_Points[i+1], alpha_between_layers[i], truncation = truncation)
            }
            if (is.na(alpha_between_layers[i]) || Frequencies[i+1] >= RoLs[i]) {
              Frequencies[i+1] <- ELL[i+1] / Pareto_Layer_Mean(Inf, Attachment_Points[i+1], alpha = tolerance, truncation = truncation)
              Frequencies[i+1] <- (Frequencies[i+1] + RoLs[i]) / 2
            }
            if (Frequencies[i+1] >= RoLs[i] * (1 - RoL_tolerance)) {
              Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
              Status <- paste0(Status, "Option Use_unlimited_Layer_for_FQ not used! ")
            }
            if (Merged_Layer[i]) {
              Frequencies[i+1] <- RoLs[i] * (1 - RoL_tolerance / 2)
            }
          }
        } else {
          Frequencies[i+1] <- Frequencies[i] * (Attachment_Points[i]/Attachment_Points[i+1])^alpha_between_layers[i-1]
        }
      }
    }
    if (!Merged_Layer[1]) {
      Frequencies[1] <- ELL[1] / Pareto_Layer_Mean(Limits[1], Attachment_Points[1], alpha_between_layers[1])
    } else {
      Frequencies[1] <- RoLs[1] * (1 + RoL_tolerance / 2)
    }
  }
  if (!is.null(FQ_at_lowest_AttPt)) {
    if (FQ_at_lowest_AttPt > RoLs[1] * (1 + RoL_tolerance / 2)) {
      Frequencies[1] <- FQ_at_lowest_AttPt
    } else {
      Status <- paste0(Status, "FQ_at_lowest_AttPt too small. Not used! ")
    }
  }
  if (!is.null(FQ_at_highest_AttPt)) {
    if (FQ_at_highest_AttPt < RoLs[k-1] * (1 - RoL_tolerance)) {
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


  if (!is.null(truncation)) {
    AvLossLastLayer <- Pareto_Layer_Mean(Inf, Attachment_Points[k], alpha = tolerance, truncation = truncation)
    if (Frequencies[k] * AvLossLastLayer <= ELL[k]) {
      warning("truncation too low!")
      return(NA)
    }
  }


  s <- Frequencies / Frequencies[1]
  l <- ELL / Frequencies[1]

  Results <- Fit_PP(Attachment_Points, s, l, truncation = truncation, alpha_max = alpha_max, minimize_ratios = minimize_ratios, merge_tolerance = merge_tolerance)

  Results$FQ <- Frequencies[1]
  if (Status == "") {Status <- "OK."}
  Results$Status <- Status
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
  sapply(x, FUN = function(x) pPareto_s(x, t, alpha, truncation))
}

pPareto_s <- function(x, t, alpha, truncation = NULL) {
  if (!is.numeric(t) || !is.numeric(alpha) || !is.numeric(x)) {
    warning("x, t and alpha must be numeric.")
    return(NA)
  }
  if (length(t) != 1 || length(alpha) != 1 || length(x) != 1) {
    warning("t and alpha must have length 1")
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
  sapply(x, FUN = function(x) dPareto_s(x, t, alpha, truncation))
}

dPareto_s <- function(x, t, alpha, truncation = NULL) {
  if (!is.numeric(t) || !is.numeric(alpha) || !is.numeric(x)) {
    warning("x, t and alpha must be numeric.")
    return(NA)
  }
  if (length(t) != 1 || length(alpha) != 1 || length(x) != 1) {
    warning("t and alpha must have length 1")
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
  sapply(x, FUN = function(x) pPiecewisePareto_s(x, t, alpha, truncation, truncation_type))
}

pPiecewisePareto_s <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (n == 1) {
    Result <- pPareto(x, t, alpha, truncation)
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


#' Density of the Piecewise Pareto Distribution
#'
#' @description Calculates the density function of the piecewise Pareto distribution. This function is deprecated. Use \code{dPiecewisePareto} instead.
#'
#' @param x Numeric. The function evaluates the density at \code{x}.
#' @param t Numeric vector. Thresholds of the piecewise Pareto distribution.
#' @param alpha Numeric vector. \code{alpha[i]} is the Pareto alpha in excess of \code{t[i]}.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the distribution is truncated at \code{truncation}.
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
  sapply(x, FUN = function(x) dPiecewisePareto_s(x, t, alpha, truncation, truncation_type))
}

dPiecewisePareto_s <- function(x, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (n == 1) {
    Result <- dPareto(x, t, alpha, truncation)
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
#' @param truncation_type Charakter. If \code{truncation_type = "wd"} then the whole distribution is truncated. If \code{truncation_type = "lp"} then a truncated Pareto is used for the last piece.
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
  sapply(p, FUN = function(y) qPiecewisePareto_s(y, t, alpha, truncation, truncation_type))
}

qPiecewisePareto_s <- function(y, t, alpha, truncation = NULL, truncation_type = "lp") {
  if (!is.numeric(t) || !is.numeric(alpha)) {
    warning("alpha and t must be numeric.")
    return(NA)
  }
  if (length(t) != length(alpha)) {
    warning("t and alpha must have the same length")
    return(NA)
  }
  n <- length(t)
  if (n == 1) {
    Result <- qPareto(y, t, alpha, truncation)
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
  if (length(y) != 1) {
    warning("y must have lenght  1!")
    return(NA)
  }
  if (!is.numeric(y)) {
    warning("y must be numeric!")
    return(NA)
  }
  if (y < 0 | y > 1) {
    warning("y must be in the interval [0,1]!")
  }
  if (y == 1) {
    if (is.null(truncation)) {
      return(Inf)
    } else {
      return(truncation)
    }
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

  n <- length(t)

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
  sapply(p, FUN = function(y) qPareto_s(y, t, alpha, truncation))
}

qPareto_s <- function(y, t, alpha, truncation = NULL) {
  if (!is.numeric(t) || !is.numeric(alpha) || !is.numeric(y)) {
    warning("x, t and alpha must be numeric.")
    return(NA)
  }
  if (length(t) != 1 || length(alpha) != 1) {
    warning("t and alpha must have length 1")
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
  if (y < 0 | y > 1) {
    warning("y must be in the interval [0,1]!")
    return(NA)
  } else if (y == 1) {
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
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#' @param alpha_min Numeric. Lower bound for alpha.
#' @param alpha_max Numeric. Upper bound for alpha.
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
#' @export

Pareto_ML_Estimator_Alpha <- function(losses, t, truncation = NULL, alpha_min = 0.001, alpha_max = 10) {
  if (!is.numeric(losses) || !is.numeric(t)) {
    warning("losses and t must be numeric.")
    return(NA)
  }
  if (length(t) != 1) {
    warning("t must have lenght 1.")
    return(NA)
  }
  if (t <= 0) {
    warning("t be positive.")
    return(NA)
  }
  losses <- losses[losses > t]
  n <- length(losses)
  if (n < 1) {
    warning("Number of losses > t must be positive.")
    return(NA)
  }
  if (!is.null(truncation)) {
    if (truncation <= t) {
      warning("truncation must be larger than t")
    }
    if (max(losses) >= truncation) {
      warning("Losses must be < truncation.")
      return(NA)
    }
  }
  if (is.null(truncation) || is.infinite(truncation)) {
    alpha_hat <- n / sum(log(losses / t))
    alpha_hat <- min(alpha_max, max(alpha_min, alpha_hat))
  } else {
    LogLikelihood <- function(alpha) {
      n * (alpha * log(t) + log(alpha) - log(1 - (t/truncation)^alpha)) - (alpha + 1) * sum(log(losses))
    }
    alpha_hat <- stats::optimise(LogLikelihood, c(alpha_min, alpha_max), maximum = T)$maximum
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

  if (distribution == "lnorm") {
    Result <- x * stats::dlnorm(x, log = FALSE, ...) / (1 - stats::plnorm(x, log.p = FALSE, ...))
  }
  if (distribution == "norm") {
    Result <- x * stats::dnorm(x, log = FALSE, ...) / (1 - stats::pnorm(x, log.p = FALSE, ...))
  }
  if (distribution == "gamma") {
    Result <- x * stats::dgamma(x, log = FALSE, ...) / (1 - stats::pgamma(x, log.p = FALSE, ...))
  }

  return(Result)
}
