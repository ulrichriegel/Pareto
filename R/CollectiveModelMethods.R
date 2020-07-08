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

