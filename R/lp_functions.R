
calculate_layer_losses <- function(df_layers, df_thresholds, overlapping, ignore_inconsistent_references) {

  if (overlapping) {
    lp_result <- solve_lp(df_layers, df_thresholds)
    status_all_info <- lp_result$status
  } else {
    if (nrow(df_layers) > 1) {
      df_layers <- df_layers[order(df_layers$attachment_point), ]
      tower <- data.frame(limit = df_layers$limit, att = df_layers$attachment_point, frequency = NA, exp_loss = df_layers$exp_loss, exp_loss_info_avaliable = TRUE)
    } else {
      tower <- data.frame(limit = numeric(0), att = numeric(0), frequency = numeric(0), exp_loss = numeric(0), exp_loss_info_avaliable = logical(0))
    }
    if (nrow(df_thresholds) > 0) {
      df_thresholds <- df_thresholds[order(df_thresholds$threshold), ]
      for (i in 1:nrow(df_thresholds)) {
        index <- df_thresholds$threshold[i] == tower$att
        if (sum(index) > 0) {
          tower$frequency[index] <- df_thresholds$frequency[i]
        } else {
          new_row <- data.frame(limit = NA, att = df_thresholds$threshold[i], frequency = df_thresholds$frequency[i], exp_loss = NA, exp_loss_info_avaliable = FALSE)
          tower <- rbind(tower, new_row)
        }
      }
    }
    tower <- tower[order(tower$att), ]
    for (i in 1:nrow(tower)) {
      if (is.na(tower$limit[i]) && i < nrow(tower)) {
        tower$limit[i] <- tower$att[i+1] - tower$att[i]
      } else if (is.na(tower$limit[i])) {
        tower$limit[i] <- Inf
      }
    }
    # fill gaps
    if (nrow(tower) > 1) {
      new_rows <- NULL
      for (i in 1:(nrow(tower)-1)) {
        if (tower$limit[i] + tower$att[i] < tower$att[i+1]) {
          new_row <- data.frame(limit = tower$att[i+1] - (tower$limit[i] + tower$att[i]), att = tower$limit[i] + tower$att[i], frequency = NA, exp_loss = NA, exp_loss_info_avaliable = FALSE)
          if (is.null(new_rows)) {
            new_rows <- new_row
          } else {
            new_rows <- rbind(new_rows, new_row)
          }
        }
      }
      if (!is.null(new_rows)) tower <- rbind(tower, new_rows)
      tower <- tower[order(tower$att), ]
    }
    if (!is.infinite(tower$limit[nrow(tower)])) {
      new_row <- data.frame(limit = Inf, att = tower$limit[nrow(tower)] + tower$att[nrow(tower)], frequency = NA, exp_loss = NA, exp_loss_info_avaliable = FALSE)
      tower <- rbind(tower, new_row)
    }

    # check that RoLs and frequencies are consistent
    fq_rol <- c(t(matrix(c(tower$frequency, tower$exp_loss / tower$limit), ncol = 2)))
    fq_rol <- fq_rol[!is.na(fq_rol)]
    if (max(diff(fq_rol)) >= 0) {
      if (!requireNamespace("lpSolve", quietly = TRUE)) {
        status_all_info <- 3
        tower$exp_loss_new <- NA
        tower$frequency_new <- NA
        return(list(tower = tower, status = status_all_info))
      } else {
        lp_result <- list(tower = tower, status = 2)
        status_all_info <- 2
      }
    } else {
      lp_result <- list(tower = tower, status = 0)
      status_all_info <- 0
    }

  }

  if (status_all_info != 0) {
    if (!ignore_inconsistent_references) {
      status_all_info <- 4
      tower$exp_loss_new <- NA
      tower$frequency_new <- NA
      return(list(tower = tower, status = status_all_info))
    }
    use_layer <- rep(F, nrow(df_layers))
    use_threshold <- rep(F, nrow(df_thresholds))
    use_layer[1] <- T
    if (nrow(df_layers) >= 2) {
      for (i in 2:nrow(df_layers)) {
        use_layer[i] <- T
        lp_result <- solve_lp(df_layers[use_layer, ], df_thresholds[use_threshold, , drop = F])
        if (lp_result$status > 0) {
          use_layer[i] <- F
        }
      }
    }
    if (nrow(df_thresholds) >= 2) {
      for (i in 1:nrow(df_thresholds)) {
        use_threshold[i] <- T
        lp_result <- solve_lp(df_layers[use_layer, ], df_thresholds[use_threshold, , drop = F])
        if (lp_result$status > 0) {
          use_threshold[i] <- F
        }
      }
    }
    lp_result <- solve_lp(df_layers[use_layer, ], df_thresholds[use_threshold, , drop = F])
  }

  att <- (lp_result$tower)$att
  limit <- (lp_result$tower)$limit
  frequency <- (lp_result$tower)$frequency
  exp_loss <- (lp_result$tower)$exp_loss
  el_info <- (lp_result$tower)$exp_loss_info_avaliable


  tower <- lp_result$tower
  n <- nrow(tower)

  if (n >= 2) {
    info_available <- ifelse(el_info | !is.na(frequency), TRUE, FALSE)
    index_info <- which(info_available)
    for (i in 1:(n-1)) {
      if (!el_info[i]) {
        index1 <- max(index_info[index_info <= i])
        index2 <- min(index_info[index_info > i])
        if (el_info[index1] && !is.na(frequency[index2])) {
          alpha <- Pareto_Find_Alpha_btw_FQ_Layer(att[index2], frequency[index2], limit[index1], att[index1], exp_loss[index1])
          fq1 <- frequency[index2] * (att[index2] / att[index1])^alpha
          exp_loss[i] <- fq1 * Pareto_Layer_Mean(limit[i], att[i], alpha, t=att[index1])
        } else if (!is.na(frequency[index1]) && !is.na(frequency[index2])) {
          alpha <- Pareto_Find_Alpha_btw_FQs(att[index1], frequency[index1], att[index2], frequency[index2])
          exp_loss[i] <- frequency[index1] * Pareto_Layer_Mean(limit[i], att[i], alpha, t=att[index1])
        } else if (el_info[index1] && el_info[index2]) {
          alpha <- Pareto_Find_Alpha_btw_Layers(limit[index1], att[index1], exp_loss[index1], limit[index2], att[index2], exp_loss[index2])
          exp_loss[i] <- Pareto_Extrapolation(limit[index1], att[index1], limit[i], att[i], alpha) * exp_loss[index1]
        } else {
          alpha <- Pareto_Find_Alpha_btw_FQ_Layer(att[index1], frequency[index1], limit[index2], att[index2], exp_loss[index2])
          exp_loss[i] <- frequency[index1] * Pareto_Layer_Mean(limit[i], att[i], alpha, t=att[index1])
        }
      }
    }
  }

  # unlimited layer
  if (!el_info[n]) {
    alpha <- 2
    if (!is.na(frequency[n])) {
      exp_loss[n] <- frequency[n] * Pareto_Layer_Mean(Inf, att[n], alpha)
    } else {
      exp_loss[n] <- Pareto_Extrapolation(limit[n-1], att[n-1], Inf, att[n], alpha) * exp_loss[n-1]
    }
  }

  # Calulate missing Frequencies
  for (i in 1:n) {
    if (is.na(frequency[i])) {
      if (i == 1) {
        alpha = Pareto_Find_Alpha_btw_Layers(limit[1], att[1], exp_loss[1], limit[2], att[2], exp_loss[2])
        frequency[1] <- exp_loss[1] / Pareto_Layer_Mean(limit[1], att[1], alpha)
      } else {
        alpha = Pareto_Find_Alpha_btw_Layers(limit[i-1], att[i-1], exp_loss[i-1], limit[i], att[i], exp_loss[i])
        frequency[i] <- exp_loss[i] / Pareto_Layer_Mean(limit[i], att[i], alpha)
      }

    }
  }

  if (n > 1) {
    set_fq_na <- rep(FALSE, n)
    for (i in 1:(n-1)) {
      if (!is.na(frequency[i]) && !is.na(frequency[i+1]) && frequency[i+1] > frequency[i] * (1 - 1e-5)) {
        set_fq_na[i] <- TRUE
        set_fq_na[i+1] <- TRUE
      }
    }
    frequency[set_fq_na] <- NaN
  }

  tower$exp_loss_new <- exp_loss
  tower$frequency_new <- frequency
  return(list(tower = tower, status = status_all_info))
}

solve_lp <- function(df_layers, df_thresholds) {
  # create increasing list of all attachment points, exit points and thresholds; last entry infinite:
  att <- union(union(df_layers$attachment_point, df_layers$attachment_point + df_layers$limit), df_thresholds$threshold)
  att <- att[order(att)]
  if (!is.infinite(max(att))) att <- c(att, Inf)

  # write tower of layers to df_result
  lim <- diff(att)
  df_result <- data.frame(limit = lim, att = att[-length(att)])

  # number of layers (in RBM input), attachment points (in tower) and thresholds (in RBM input):
  n_layers <- nrow(df_layers)
  n_att <- length(att) - 1
  n_thr <- nrow(df_thresholds)

  # create matrix with conditions for lpSolve (match layer losses)
  lp_con <- matrix(0, nrow = n_layers, ncol = n_att)
  for (i in 1:(length(att) - 1)) {
    lp_con[1:n_layers, i] <- ifelse(df_layers$attachment_point <= att[i] & df_layers$attachment_point + df_layers$limit >=att[i+1], 1, 0)
  }
  # right hand side and direction for lpSolve  (match layer losses)
  lp_rhs <- df_layers$exp_loss
  lp_dir <- rep("=", nrow(df_layers))

  # for which layers in the tower do we have information?
  info_available <- colSums(lp_con) > 0

  # add conditions to ensure decreasing RoLs in the tower
  if (n_att >= 2) {
    temp <- matrix(0, nrow = n_att - 1, ncol = n_att)
    for (i in 1:(n_att - 1)) {
      temp[i, i] <- 1 / lim[i]
      temp[i, i + 1] <- -1 / lim[i + 1]
    }
    lp_con <- rbind(lp_con, temp)
    lp_rhs <- c(lp_rhs, rep(0, n_att - 1))
    lp_dir <- c(lp_dir, rep(">=", n_att - 1))
  }

  # add conditions to ensure that RoLs and excess frequencies are consistent
  freq <- rep(NA, n_att)
  if (n_thr >= 1) {
    for (i in 1:n_thr) {
      index <- match(df_thresholds$threshold[i], att)
      freq[index] <- df_thresholds$frequency[i]
      temp <- matrix(0, nrow = 1, ncol = n_att)
      temp[1, index] <- 1 / lim[index]
      lp_con <- rbind(lp_con, temp)
      lp_dir <- c(lp_dir, "<=")
      lp_rhs <- c(lp_rhs, df_thresholds$frequency[i])
      if (index > 1) {
        temp <- matrix(0, nrow = 1, ncol = n_att)
        temp[1, index - 1] <- 1 / lim[index - 1]
        lp_con <- rbind(lp_con, temp)
        lp_dir <- c(lp_dir, ">=")
        lp_rhs <- c(lp_rhs, df_thresholds$frequency[i])
      }
    }
  }

  # for each layer in the tower with info available, maximize the expected loss of the layer when solving the lp problem
  # then take the average of the solutions; the idea is to obtain a solution which is not at the edge of the condition
  indices <- which(info_available)
  solution <- numeric(ncol(lp_con))
  status <- 0
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    for (i in indices) {
      lp_obj <- rep(0,ncol(lp_con))
      lp_obj[i] <- 1
      lp_result <- lpSolve::lp("max", lp_obj, lp_con, lp_dir, lp_rhs)
      status <- status + lp_result$status
      solution <- solution + lp_result$solution
    }
  }
  solution <- solution / length(indices)
  status <- min(status, 1)

  # add frequencies from RBM input to data frame
  df_result$frequency <- freq
  # add expected losses for the tower to data frame
  if (status == 0) {
    df_result$exp_loss <- solution
  } else {
    df_result$exp_loss <- NA
  }
  # add indicator wheter info is available for a layer in the tower
  df_result$exp_loss_info_avaliable <- info_available

  return(list(tower = df_result, status = status))

}
