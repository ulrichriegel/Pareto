#' PGP_Model (Collective Panjer & Generalized Pareto Model) Object
#'
#' @description Constructor function for the PGP_Model object
#'
#' @param FQ Numerical. Expected claim count of the collective model.
#' @param t Numeric. Threshold of the Pareto distribution. If \code{t} is \code{NULL} (default) then \code{t <- Attachment Point} is used
#' @param alpha_ini Numeric. Initial Pareto alpha (at \code{t}).
#' @param alpha_tail Numeric. Tail Pareto alpha.
#' @param truncation Numeric. If \code{truncation} is not \code{NULL} and \code{truncation > t}, then the Pareto distribution is truncated at \code{truncation}.
#' @param dispersion Numerical. Dispersion of the Panjer distribution (i.e. variance to mean ratio).
#' @param Status Numerical indicator if a function returns a PGP_Model object: 0 = success, 1 = some information has been ignored, 2 = no solution found
#' @param Comment Charakter. An optional comment.

#' @examples
#' PGPM <- PGP_Model(2, t = 1000, alpha_ini = 1, alpha_tail = 2 , dispersion = 2)
#' PGPM
#'
#' @export

PGP_Model <- function(FQ = NULL, t = NULL, alpha_ini = NULL, alpha_tail = NULL, truncation = NULL, dispersion = 1, Status = 0, Comment = "OK") {
  obj <- list(FQ = FQ, t = t, alpha_ini = alpha_ini, alpha_tail = alpha_tail, truncation = truncation, dispersion = dispersion, Status = Status, Comment = Comment)
  class(obj) <- "PGP_Model"

  if (!is.valid.PGP_Model(obj)) {
    obj <- list(FQ = NULL, t = NULL, alpha_ini = NULL, alpha_tail = NULL, truncation = NULL, dispersion = 1, Status = 2, Comment = is.valid.PGP_Model(obj, comment = TRUE))
    class(obj) <- "PGP_Model"
  }
  return(obj)

}

#' Print a PGP_Model Object(Collective Panjer & Generalized Pareto Model) Object
#'
#' @description Print method for PGP_Model objects
#'
#' @param x PGP_Model object.
#' @param ... Other arguments, all currently ignored.
#'
#' @export

print.PGP_Model <- function(x, ...) {
  if (!is.positive.finite.number(x$dispersion)) {
    fq_dist <- "Panjer"
  } else if (x$dispersion == 1) {
    fq_dist <- "Poisson"
  } else if (x$dispersion > 1) {
    fq_dist <- "Negative Binomial"
  } else {
    fq_dist <- "Binomial"
  }
  cat("\nPanjer & Generalized Pareto model\n\n")
  cat("Collective model with a ", fq_dist, " distribution for the claim count and a generalized Pareto distributed severity.", sep = "")
  cat("\n\n", fq_dist, " Distribution:", sep = "")

  cat("\nExpected Frequency:   ", x$FQ, sep = "")
  if (is.positive.finite.number(x$dispersion)  && x$dispersion != 1) {
    cat("\nDispersion:           ", x$dispersion, sep = "")
    if (is.positive.finite.number(x$FQ) && x$dispersion > 1) {
      cat(" (i.e. contagion = ", (x$dispersion - 1)/x$FQ, ")", sep = "")
    }
  }
  cat("\nGeneralized Pareto Distribution:")
  cat("\nThreshold:         ", x$t, sep = "   ")
  cat("\nalpha_ini:         ", x$alpha_ini, sep = "   ")
  cat("\nalpha_tail:        ", x$alpha_tail, sep = "   ")
  if (!is.null(x$truncation)) {
    cat("\nTruncation:           ", x$truncation, sep = "")
  } else {
    cat("\nThe distribution is not truncated.")
  }
  cat("\n\nStatus:              ", x$Status)
  cat("\nComments:            ", x$Comment)
  if (!is.valid.PGP_Model(x)) {
    cat("\n\nThe model is not valid.\n")
    cat(is.valid.PGP_Model(x, comment = TRUE))
  }
  cat("\n\n")

}

#' Check if an object is a PGP_Model
#'
#' @description Checks if the class of an object is 'PGP_Model'
#'
#' @param x Object to be checked.

#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' is.valid.PGP_Model(PGPM)
#' is.valid.PGP_Model(PGPM, comment = TRUE)
#'
#' PGPM$alpha_tail <- -2
#' is.PGP_Model(PGPM)
#' is.valid.PGP_Model(PGPM)
#' is.valid.PGP_Model(PGPM, comment = TRUE)
#'
#' @export

is.PGP_Model <- function(x) {
  if (class(x) == "PGP_Model") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' Check if an object is a valid PGP_Model
#'
#' @description Checks if an object is a PGP_Model object and whether it is valid for the use in functions like \code{Layer_Mean}
#'
#' @param x Object to be checked.
#' @param comment If FALSE then the function returns a boolean indicating whether \code{x} is a valid PGP_Model. If TRUE then the function returns a comment instead.

#' @examples
#' PGPM <- PGP_Model(2, 1000, 1, 2, dispersion = 2)
#' PGPM
#' is.valid.PGP_Model(PGPM)
#' is.valid.PGP_Model(PGPM, comment = TRUE)
#'
#' PGPM$alpha_tail <- -2
#' is.valid.PGP_Model(PGPM)
#' is.valid.PGP_Model(PGPM, comment = TRUE)
#'
#' @export

is.valid.PGP_Model <- function(x, comment = FALSE) {
  if (class(x) != "PGP_Model" || typeof(x) != "list") {
    if (!comment) {
      return(FALSE)
    } else {
      return("Object does not have class PGP_Model.")
    }
  }
  required_elements <- c("FQ", "t", "alpha_ini", "alpha_tail", "truncation", "dispersion", "Status", "Comment")
  available <- required_elements %in% names(x)
  if (sum(!available) > 0) {
    if (!comment) {
      return(FALSE)
    } else {
      return(paste("Not all required list elements available. Missing elements:", paste(required_elements[!available], collapse = ", ")))
    }
  }


  if (!is.positive.finite.number(x$FQ)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("FQ must be a positive number.")
    }
  }
  if (!valid.parameters.GenPareto(x$t, x$alpha_ini, x$alpha_tail, x$truncation)) {
    if (!comment) {
      return(FALSE)
    } else {
      return(valid.parameters.GenPareto(x$t, x$alpha_ini, x$alpha_tail, x$truncation, comment = TRUE))
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











PGP_Model_Exp_Layer_Loss_s <- function(Cover, AttachmentPoint, PGP_Model) {
  if (!is.valid.PGP_Model(PGP_Model)) {
    warning(is.valid.PGP_Model(PGP_Model, comment = TRUE))
    return(NaN)
  } else {
    return(PGP_Model$FQ * GenPareto_Layer_Mean(Cover, AttachmentPoint, alpha_ini = PGP_Model$alpha_ini, alpha_tail = PGP_Model$alpha_tail, t = PGP_Model$t, truncation = PGP_Model$truncation))
  }
}




PGP_Model_Exp_Layer_Loss_v <- Vectorize(PGP_Model_Exp_Layer_Loss_s, c("Cover", "AttachmentPoint"))



PGP_Model_Layer_Var_s <- function(Cover, AttachmentPoint, PGP_Model) {
  if (!is.valid.PGP_Model(PGP_Model)) {
    warning(is.valid.PGP_Model(PGP_Model, comment = TRUE))
    return(NaN)
  } else {
    E_N <- PGP_Model$FQ
    Var_N <- E_N * PGP_Model$dispersion
    E_X <- GenPareto_Layer_Mean(Cover, AttachmentPoint, alpha_ini = PGP_Model$alpha_ini, alpha_tail = PGP_Model$alpha_tail, t = PGP_Model$t, truncation = PGP_Model$truncation)
    Var_X <- GenPareto_Layer_Var(Cover, AttachmentPoint, alpha_ini = PGP_Model$alpha_ini, alpha_tail = PGP_Model$alpha_tail, t = PGP_Model$t, truncation = PGP_Model$truncation)
    return(E_N * Var_X + Var_N * E_X^2)
  }
}

PGP_Model_Layer_Var_v <- Vectorize(PGP_Model_Layer_Var_s, c("Cover", "AttachmentPoint"))







PGP_Model_Layer_Sd_s <- function(Cover, AttachmentPoint, PGP_Model) {
  if (!is.valid.PGP_Model(PGP_Model)) {
    warning(is.valid.PGP_Model(PGP_Model, comment = TRUE))
    return(NaN)
  } else {
    return(sqrt(PGP_Model_Layer_Var_v(Cover, AttachmentPoint, PGP_Model)))
  }
}


PGP_Model_Layer_Sd_v <- Vectorize(PGP_Model_Layer_Sd_s, c("Cover", "AttachmentPoint"))




PGP_Model_Excess_Frequency_s <- function(x, PGP_Model) {
  if (!is.valid.PGP_Model(PGP_Model)) {
    warning(is.valid.PGP_Model(PGP_Model, comment = TRUE))
    return(NaN)
  } else if (!is.atomic(x) || !is.numeric(x) || length(x) != 1 || is.na(x)) {
    warning("x must be a number.")
    return(NaN)
  } else {
    return(PGP_Model$FQ * (1 - pGenPareto(x, PGP_Model$alpha_ini, PGP_Model$alpha_tail, t = PGP_Model$t, truncation = PGP_Model$truncation)))
  }
}
PGP_Model_Excess_Frequency_v <- Vectorize(PGP_Model_Excess_Frequency_s, c("x"))






