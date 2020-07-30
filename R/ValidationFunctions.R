is.positive.vector <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    return(FALSE)
  }
  k <- length(x)
  if (sum(x > 0, na.rm = TRUE) < k) {
    return(FALSE)
  }
  return(TRUE)
}


is.positive_or_NA.finite.vector <- function(x) {
  if(!is.atomic(x) || length(x) < 1) {
    return(FALSE)
  }
  if (length(x[!is.na(x)]) > 0 && !is.positive.finite.vector(x[!is.na(x)])) {
    return(FALSE)
  }
  return(TRUE)
}


is.nonnegative_or_NA.finite.vector <- function(x) {
  if(!is.atomic(x) || length(x) < 1) {
    return(FALSE)
  }
  if (length(x[!is.na(x)]) > 0 && !is.nonnegative.finite.vector(x[!is.na(x)])) {
    return(FALSE)
  }
  return(TRUE)
}

is.NA.vector <- function(x) {
  if(!is.atomic(x) || length(x) < 1) {
    return(FALSE)
  }
  if (length(x[!is.na(x)]) > 0) {
    return(FALSE)
  }
  return(TRUE)
}



is.nonnegative.vector <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    return(FALSE)
  }
  k <- length(x)
  if (sum(x >= 0, na.rm = TRUE) < k) {
    return(FALSE)
  }
  return(TRUE)
}


is.positive.finite.vector <- function(x) {
  if(!is.positive.vector(x)) {
    return(FALSE)
  }
  if (sum(is.infinite(x)) > 0) {
    return(FALSE)
  }
  return(TRUE)
}


is.nonnegative.finite.vector <- function(x) {
  if(!is.nonnegative.vector(x)) {
    return(FALSE)
  }
  if (sum(is.infinite(x)) > 0) {
    return(FALSE)
  }
  return(TRUE)
}

is.number <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) != 1 || is.na(x)) {
    return(FALSE)
  }
  return(TRUE)
}

is.positive.number <- function(x) {
  if (!is.positive.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.positive.finite.number <- function(x) {
  if (!is.positive.finite.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.nonnegative.number <- function(x) {
  if (!is.nonnegative.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.nonnegative.finite.number <- function(x) {
  if (!is.nonnegative.finite.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.TRUEorFALSE <- function(x) {
  if (isTRUE(x) || isFALSE(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


valid.parameters.PiecewisePareto <- function(t, alpha, truncation, truncation_type, comment = FALSE) {
  if (!is.positive.finite.vector(t) || !is.nonnegative.finite.vector(alpha)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("t must be a positive vector, alpha must be a non-negative vector.")
    }
  }
  k <- length(t)
  if (length(alpha) != k) {
    if (!comment) {
      return(FALSE)
    } else {
      return("t and alpha must have the same length.")
    }
  }
  if (k > 1 && min(diff(t))<=0) {
    if (!comment) {
      return(FALSE)
    } else {
      return("t must be increasing.")
    }
  }
  if (alpha[k] <= 0) {
    if (!comment) {
      return(FALSE)
    } else {
      return("Last alpha must be positive.")
    }
  }
  if (!is.atomic(truncation_type) || length(truncation_type) != 1 || !(truncation_type %in% c("lp", "wd"))) {
    if (!comment) {
      return(FALSE)
    } else {
      return("truncation_type must be 'lp' or 'wd'.")
    }
  }
  if (!is.null(truncation) && (!is.positive.number(truncation) || truncation <= max(t))) {
    if (!comment) {
      return(FALSE)
    } else {
      return("truncation must be a positive number > max(t).")
    }
  }
  if (!comment) {
    return(TRUE)
  } else {
    return("OK")
  }
}



valid.parameters.Pareto <- function(t, alpha, truncation, allow.alpha.zero = FALSE, comment = FALSE) {
  if (allow.alpha.zero) {
    if (!is.positive.finite.number(t) || !is.nonnegative.finite.number(alpha)) {
      if (!comment) {
        return(FALSE)
      } else {
        return("t must be positive and alpha must be nonnegative.")
      }
    }

  } else {
    if (!is.positive.finite.number(t) || !is.positive.finite.number(alpha)) {
      if (!comment) {
        return(FALSE)
      } else {
        return("t and alpha must be positive numbers.")
      }
    }
  }


  if (!is.null(truncation) && (!is.positive.number(truncation) || truncation <= t)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("truncation must be a positive number > t.")
    }
  }
  if (!comment) {
    return(TRUE)
  } else {
    return("OK")
  }
}


valid.parameters.GenPareto <- function(t, alpha_ini, alpha_tail, truncation, comment = FALSE) {
    if (!is.positive.finite.number(t) || !is.positive.finite.number(alpha_ini) || !is.positive.finite.number(alpha_tail)) {
      if (!comment) {
        return(FALSE)
      } else {
        return("t, alpha_ini and alpha_tail must be positive numbers.")
      }
    }


  if (!is.null(truncation) && (!is.positive.number(truncation) || truncation <= t)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("truncation must be a positive number > t.")
    }
  }
  if (!comment) {
    return(TRUE)
  } else {
    return("OK")
  }
}




