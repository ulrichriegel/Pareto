Fit_PP <- function(a, s, l, truncation, tolerance = 1e-10, alpha_max = 100, minimize_ratios = T, merge_tolerance = 1e-6) {
  # a vector of attachment points
  # s vector of frequencies
  # l[i] exp loss of layer a[i+1] - a[i] xs a[i]

  Result <- list(t = NA, alpha = NA)

  if (length(a) != length(s)) {
    Result$Status <- "a and s must have same lenght!"
    return(Result)
  }
  n <- length(a)
  if (min(diff(a))<=0) {
    Result$Status <- "a must be ascending"
    return(Result)
  }
  if (max(diff(s))>=0) {
    Result$Status <- "s must be descending"
    return(Result)
  }
  if (a[1] <= 0) {
    Result$Status <- "a must be positive."
    return(Result)
  }
  if (s[n] <= 0) {
    Result$Status <- "s must be positive."
    return(Result)
  }

  t <- numeric(2*n-1)
  alpha <- numeric(2*n-1)
  t[2*(1:n)-1] <- a
  # alpha[2*n-1] <- s[n] * a[n] / l[n] + 1 # (Formula without truncation!)
  alpha[2*n-1] <- Pareto_Find_Alpha_btw_FQ_Layer(a[n], s[n], Inf, a[n], l[n], max_alpha = alpha_max, tolerance = tolerance, truncation = truncation)
  for (k in 1:(n-1)) {
    taus <- Calculate_taus(s[k], s[k+1], a[k], a[k+1], l[k], tolerance = tolerance)
    lower_bound <- min(a[k] * (s[k] / s[k+1])^(1/alpha_max), (a[k] + a[k+1]) / 2)
    upper_bound <- max(a[k+1] * (s[k+1] / s[k])^(1/alpha_max), (a[k] + a[k+1]) / 2)
    if (taus[2] < lower_bound) {
      taus <- rep(lower_bound, 2)
    } else if (taus[1] > upper_bound) {
      taus <- rep(upper_bound, 2)
    } else {
      taus[1] <- max(taus[1], lower_bound)
      taus[2] <- min(taus[2], upper_bound)
    }
    t[2*k] <- (taus[1] + taus[2]) / 2
    if (minimize_ratios) {
      if (taus[1]<taus[2]) {
        t_temp <- (taus[1] + taus[2]) / 2
        penalty <- function(t) {
          alphas <- Calculate_alphas(s[k], s[k+1], a[k], a[k+1], l[k], t, tolerance = tolerance, alpha_max = alpha_max)
          #return(abs(alphas[1]-alphas[2]))
          if (alphas[1]<alphas[2]) {
            if (alphas[2]>1000*alphas[1]) {
              return(1000)
            } else {
              return(alphas[2]/alphas[1])
            }
          } else if (alphas[2]<alphas[1]) {
            if (alphas[1]>1000*alphas[2]) {
              return(1000)
            } else {
              return(alphas[1]/alphas[2])
            }
          } else {
            return(1)
          }
        }
        t_temp <- stats::optim(t_temp, fn = penalty, lower = taus[1], upper = taus[2], method = "L-BFGS-B")$par
        t[2*k] <- t_temp
      } else {
        t[2*k] <- taus[1]
      }
    }
    alpha[(2*k-1):(2*k)] <- Calculate_alphas(s[k], s[k+1], a[k], a[k+1], l[k], t[2*k], tolerance = tolerance, alpha_max = alpha_max)
  }

  q <- 2*n-1
  is_equal <- abs(alpha[2:q] - alpha[1:(q-1)]) < merge_tolerance
  is_equal <- c(FALSE, is_equal)
  t <- t[!is_equal]
  alpha <- alpha[!is_equal]

  Result <- list(t = t, alpha = alpha, Status = "OK")
  return(Result)
}



Calculate_taus <- function(s_0, s_1, a_0, a_1, l_0, tolerance = 1e-10) {
  # s_0 = s_k
  # s_1 = s_{k+1}
  # a_0 = a_k
  # a_1 = a_{k+1}
  # l_0 = l_k

  LL <- function(a, b, alpha) {
    if (alpha == 1) {
      return(a * (log(b) - log(a)))
    } else {
      return(a / (1-alpha) * ((b/a)^(1-alpha) - 1))
    }
  }

  lambda <- function(t, alpha) {
    # use s_0, a_0, s_1 and a_1 from environment
    Result <- s_0 * LL(a_0, t, alpha) + s_0 * (a_0 / t)^alpha * LL(t, a_1, (log(s_1/s_0) - alpha*log(a_0/t)) / log(t/a_1) )
    return(Result)
  }

  f <- function(x) {
    lambda(x, log(s_1/s_0) / log(a_0/x)) - l_0
  }

  delta <- tolerance *(a_1 - a_0)
  tau_u <- NA
  try(tau_u <- stats::uniroot(f, interval = c(a_0 + delta, a_1 - delta), tol = tolerance * (a_1 - a_0))$root, silent = T)
  if (is.na(tau_u)) {
    if (f((a_0 + a_1) / 2) < 0) {
      tau_u <- a_1
    } else {
      tau_u <- a_0
    }
  }

  g <- function(x) {
    lambda(x, 0) - l_0
  }

  tau_l <- NA
  try(tau_l <- stats::uniroot(g, interval = c(a_0, a_1), tol = tolerance * (a_1 - a_0))$root, silent = T)
  if (is.na(tau_l)) {
    if (g((a_0 + a_1) / 2) > 0) {
      tau_l <- a_0
    } else {
      tau_l <- a_1
    }
  }
  # if (is.na(tau_l)) {
  #   tau_l <- a_0
  # }

  return(c(tau_l, tau_u))
}




Calculate_alphas <- function(s_0, s_1, a_0, a_1, l_0, t, alpha_max = 100, tolerance = 1e-10) {
  # s_0 = s_k
  # s_1 = s_{k+1}
  # a_0 = a_k
  # a_1 = a_{k+1}
  # l_0 = l_k

  LL <- function(a, b, alpha) {
    if (alpha == 1) {
      return(a * (log(b) - log(a)))
    } else {
      return(a / (1-alpha) * ((b/a)^(1-alpha) - 1))
    }
  }

  lambda <- function(t, alpha) {
    # use s_0, a_0, s_1 and a_1 from environment
    Result <- s_0 * LL(a_0, t, alpha) + s_0 * (a_0 / t)^alpha * LL(t, a_1, (log(s_1/s_0) - alpha*log(a_0/t)) / log(t/a_1) )
    return(Result)
  }

  f <- function(alpha) {
    Result <- lambda(t, alpha) - l_0
    if (!is.nan(Result)) {
      return(Result)
    } else {
      for (i in 1:20) {
        Result <- lambda(t, alpha / 1.1^i) - l_0
        if (!is.nan(Result)) {
          return(Result)
        }
      }
      return(Result)
    }
  }

  alpha_0 <- NA
  try(alpha_0 <- stats::uniroot(f, c(0, alpha_max), tol = tolerance * (a_1 - a_0))$root, silent = T)
  if (is.na(alpha_0)) {
    if (abs(f(0) < abs(f(alpha_max)))) {
      alpha_0 <- 0
    } else {
      alpha_0 <- alpha_max
    }
  }
  alpha_1 <- (log(s_1/s_0) - alpha_0 * log(a_0/t)) / log(t/a_1)
  alpha_1 <- min(alpha_max, alpha_1)

  return(pmax(c(alpha_0, alpha_1), 0))
}


