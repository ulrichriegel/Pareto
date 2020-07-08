#' Expected Loss of a Reinsurance Layer
#'
#' @description  Calculates the expected loss of a reinsurance layer for a collective model
#'
#' @param CollectiveModel A collective model object. Currently only \code{PPP_Models} are handled.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The expected loss of the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
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

Layer_Mean <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) UseMethod("Layer_Mean")



#' Expected Loss of a Reinsurance Layer
#'
#' @description  Calculates the expected loss of a reinsurance layer for a PPP_Model
#'
#' @param CollectiveModel PPP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The expected loss of the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
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

Layer_Mean.PPP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) PPP_Model_Exp_Layer_Loss_v(Cover, AttachmentPoint, CollectiveModel)


#' Expected Loss of a Reinsurance Layer
#'
#' @description  Calculates the expected loss of a reinsurance layer for a PGP_Model
#'
#' @param CollectiveModel PGP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The expected loss of the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Example1_AP
#' Example1_Cov
#' Example1_EL
#' Layer_Mean(PGPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Mean.PGP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) PGP_Model_Exp_Layer_Loss_v(Cover, AttachmentPoint, CollectiveModel)










#' Variance of a Reinsurance Layer
#'
#' @description  Calculates the variance of the loss in a reinsurance layer for a collective model
#'
#' @param CollectiveModel A collective model object. Currently only \code{PPP_Models} are handled.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The variance of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Var(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Var <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) UseMethod("Layer_Var")



#' Variance of a Reinsurance Layer
#'
#' @description  Calculates the variance of the loss in a reinsurance layer for a PPP_model
#'
#' @param CollectiveModel PPP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The variance of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Var(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Var.PPP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) PPP_Model_Layer_Var_v(Cover, AttachmentPoint, CollectiveModel)



#' Variance of a Reinsurance Layer
#'
#' @description  Calculates the variance of the loss in a reinsurance layer for a PGP_model
#'
#' @param CollectiveModel PGP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The variance of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Var(PGPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Var.PGP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) PGP_Model_Layer_Var_v(Cover, AttachmentPoint, CollectiveModel)












#' Standard Deviation of a Reinsurance Layer
#'
#' @description  Calculates the standard deviation of the loss in a reinsurance layer for a collective model
#'
#' @param CollectiveModel A collective model object. Currently only \code{PPP_Models} are handled.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The standard deviation of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Sd(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Sd <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) UseMethod("Layer_Sd")



#' Standard Deviation of a Reinsurance Layer
#'
#' @description  Calculates the standard deviation of the loss in a reinsurance layer for a PPP_model
#'
#' @param CollectiveModel PPP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The standard deviation of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Sd(PPPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Sd.PPP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) PPP_Model_Layer_Sd_v(Cover, AttachmentPoint, CollectiveModel)



#' Standard Deviation of a Reinsurance Layer
#'
#' @description  Calculates the standard deviation of the loss in a reinsurance layer for a PGP_model
#'
#' @param CollectiveModel PGP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The standard deviation of the loss in the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' Example1_Cov <- c(diff(Example1_AP), Inf)
#' Layer_Sd(PGPM, Example1_Cov, Example1_AP)
#'
#' @export

Layer_Sd.PGP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) PGP_Model_Layer_Sd_v(Cover, AttachmentPoint, CollectiveModel)














#' Expected Frequency in Excess of a Threshold
#'
#' @description  Calculates the expected frequency in excess of a threshold for a collective model
#'
#' @param CollectiveModel A collective model object. Currently only \code{PPP_Models} are handled.
#' @param x Numeric. Threshold.
#'
#' @return The expected frequency in excess of \code{x} for the given \code{CollectiveModel}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Excess_Frequency(PPPM, c(-Inf, 0, 1000, 2000, 3000, Inf))
#'
#' @export

Excess_Frequency <- function(CollectiveModel, x = 0) UseMethod("Excess_Frequency")



#' Expected Frequency in Excess of a Threshold
#'
#' @description  Calculates the expected frequency in excess of a threshold for a PPP_model
#'
#' @param CollectiveModel PPP_Model object.
#' @param x Numeric. Threshold.
#'
#' @return The expected frequency in excess of \code{x} for the given \code{CollectiveModel}
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(Example1_AP, Example1_EL)
#' PPPM
#' Excess_Frequency(PPPM, c(-Inf, 0, 1000, 2000, 3000, Inf))
#'
#' @export

Excess_Frequency.PPP_Model <- function(CollectiveModel, x = 0) PPP_Model_Excess_Frequency_v(x, CollectiveModel)


#' Expected Frequency in Excess of a Threshold
#'
#' @description  Calculates the expected frequency in excess of a threshold for a PGP_model
#'
#' @param CollectiveModel PGP_Model object.
#' @param x Numeric. Threshold.
#'
#' @return The expected frequency in excess of \code{x} for the given \code{CollectiveModel}
#'
#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' Excess_Frequency(PGPM, c(-Inf, 0, 1000, 2000, 3000, Inf))
#'
#' @export

Excess_Frequency.PGP_Model <- function(CollectiveModel, x = 0) PGP_Model_Excess_Frequency_v(x, CollectiveModel)






#' Simulate Losses with a Collective Model
#'
#' @description  Simulates losses with a collective model
#'
#' @param CollectiveModel A collective model object. Currently only \code{PPP_Models} are handled.
#' @param nyears Integer. Number of simulated years.
#'
#' @return A matrix where row k contains the simulated losses of the kth simulated year.
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(c(1000, 2000, 3000), c(2000, 1000, 500),
#'                                            truncation = 10000, truncation_type = "wd")
#' PPPM
#' Simulate_Losses(PPPM, 100)
#'
#' @export

Simulate_Losses <- function(CollectiveModel, nyears = 1) UseMethod("Simulate_Losses")


#' Simulate Losses with a PPP_Model
#'
#' @description  Simulates losses with a PPP_Model
#'
#' @param CollectiveModel PPP_Model object.
#' @param nyears Integer. Number of simulated years.
#'
#' @return A matrix where row k contains the simulated losses of the kth simulated year.
#'
#' @examples
#' PPPM <- PiecewisePareto_Match_Layer_Losses(c(1000, 2000, 3000), c(2000, 1000, 500),
#'                                            truncation = 10000, truncation_type = "wd")
#' PPPM
#' Simulate_Losses(PPPM, 100)
#'
#' @export

Simulate_Losses.PPP_Model <- function(CollectiveModel, nyears = 1) {
  if (!is.valid.PPP_Model(CollectiveModel)) {
    warning(is.valid.PPP_Model(CollectiveModel, comment = TRUE))
    return(NaN)
  }
  if (!is.positive.finite.number(nyears)) {
    warning("nyears must be a positive number.")
    return(NaN)
  } else {
    nyears <- ceiling(nyears)
  }
  claim_count <- rPanjer(nyears, CollectiveModel$FQ, CollectiveModel$dispersion)
  claims <- rPiecewisePareto(sum(claim_count), CollectiveModel$t, CollectiveModel$alpha, CollectiveModel$truncation, CollectiveModel$truncation_type)
  result <- matrix(NaN, nrow = nyears, ncol = max(claim_count))
  result[col(result) <= claim_count] <- claims
  return(result)
}



#' Simulate Losses with a PGP_Model
#'
#' @description  Simulates losses with a PGP_Model
#'
#' @param CollectiveModel PGP_Model object.
#' @param nyears Integer. Number of simulated years.
#'
#' @return A matrix where row k contains the simulated losses of the kth simulated year.
#'
#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' Simulate_Losses(PGPM, 100)
#'
#' @export

Simulate_Losses.PGP_Model <- function(CollectiveModel, nyears = 1) {
  if (!is.valid.PGP_Model(CollectiveModel)) {
    warning(is.valid.PGP_Model(CollectiveModel, comment = TRUE))
    return(NaN)
  }
  if (!is.positive.finite.number(nyears)) {
    warning("nyears must be a positive number.")
    return(NaN)
  } else {
    nyears <- ceiling(nyears)
  }
  claim_count <- rPanjer(nyears, CollectiveModel$FQ, CollectiveModel$dispersion)
  claims <- rGenPareto(sum(claim_count), alpha_ini = CollectiveModel$alpha_ini, alpha_tail = CollectiveModel$alpha_tail, t = CollectiveModel$t, truncation = CollectiveModel$truncation)
  result <- matrix(NaN, nrow = nyears, ncol = max(claim_count))
  result[col(result) <= claim_count] <- claims
  return(result)
}



