library(Pareto)
context("test functions Match_Layer_Losses and PPP_Model")

test_that("PiecewisePareto_Match_Layer_Losses", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
})


test_that("PiecewisePareto_Match_Layer_Losses truncated lp", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, truncation = 10000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
})

test_that("PiecewisePareto_Match_Layer_Losses truncated wd", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, truncation = 10000, truncation_type = 'wd')
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
})


test_that("PiecewisePareto_Match_Layer_Losses truncated lp & frequencies", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
  expect_equal(Excess_Frequency(Fit, c(1000, 2000, 5000)), c(1.1, 0.95, 0.5))

  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, FQ_at_lowest_AttPt = 1.2, FQ_at_highest_AttPt = 0.4, truncation = 10000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
  expect_equal(Excess_Frequency(Fit, c(1000, 2000, 5000)), c(1.2, 0.95, 0.4))

})

test_that("PiecewisePareto_Match_Layer_Losses truncated wd & frequencies", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd')
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
  expect_equal(Excess_Frequency(Fit, c(1000, 2000, 5000)), c(1.1, 0.95, 0.5))

  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, FQ_at_lowest_AttPt = 1.2, FQ_at_highest_AttPt = 0.4, truncation = 10000, truncation_type = 'wd')
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(Layer_Mean(Fit, Cover, AP), EL)
  expect_equal(Excess_Frequency(Fit, c(1000, 2000, 5000)), c(1.2, 0.95, 0.4))

})

test_that("PiecewisePareto_Match_Layer_Losses only two layers", {
  Fit <- PiecewisePareto_Match_Layer_Losses(c(1000,2000), c(100,100))
  expect_equal(Fit$alpha, 2)
  Fit <- PiecewisePareto_Match_Layer_Losses(c(1000,2000), c(100,100), truncation = 10000, truncation_type = "wd")
  expect_equal(Fit$alpha, 1.4364891695020006)

})


test_that("PiecewisePareto_Match_Layer_Losses where exp. loss vector contains zeroes", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 0)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, alpha_max = 1000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(round(Layer_Mean(Fit, Cover, AP)[-5], 5), round(EL[-5], 5))
  expect_equal(Excess_Frequency(Fit, AP[5]) < 3, TRUE)
  expect_equal(round(Excess_Frequency(Fit, c(1000, 2000, 5000)), 5), round(c(1.1, 0.95, 0.5), 5))

  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, alpha_max = 1000, FQ_at_lowest_AttPt = 1.5, FQ_at_highest_AttPt = 0.01)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(round(Layer_Mean(Fit, Cover, AP)[-5], 5), round(EL[-5], 5))
  expect_equal(round(Excess_Frequency(Fit, c(1000, 2000, 5000)), 5), round(c(1.5, 0.95, 0.01), 5))

  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, alpha_max = 1000, FQ_at_lowest_AttPt = 1.5, FQ_at_highest_AttPt = 0.01, truncation = 7000, truncation_type = "wd")
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(round(Layer_Mean(Fit, Cover, AP)[-5], 5), round(EL[-5], 5))
  expect_equal(round(Excess_Frequency(Fit, c(1000, 2000, 5000)), 5), round(c(1.5, 0.95, 0.01), 5))

  EL <- c(1000, 900, 800, 0, 0)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, alpha_max = 1000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(round(Layer_Mean(Fit, Cover, AP)[-(4:5)], 5), round(EL[-(4:5)], 5))
  expect_equal(Excess_Frequency(Fit, AP[4:5]) < 3, rep(TRUE, 2))
  expect_equal(round(Excess_Frequency(Fit, c(1000, 2000)), 5), round(c(1.1, 0.95), 5))


})


test_that("PiecewisePareto_Match_Layer_Losses with TotalLoss_Frequencies", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 0)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.1)
  TLFQs <- c(0.97, NA, NA, 0.2)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, TotalLoss_Frequencies = TLFQs, alpha_max = 1000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(round(Layer_Mean(Fit, Cover, AP)[-5], 5), round(EL[-5], 5))
  expect_equal(round(Excess_Frequency(Fit, c(1000, 2000, 5000)), 5), round(c(1.1, 0.95, 0.2), 5))
  expect_equal(sum(round(Fit$alpha, 5) == 1000), 2)

  EL <- c(1000, 900, 800, 600, 300)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, 0.81, 0.7, 0.1)
  TLFQs <- c(0.97, 0.85, 0.75, 0.2)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, TotalLoss_Frequencies = TLFQs, alpha_max = 1000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(round(Layer_Mean(Fit, Cover, AP), 5), round(EL, 5))
  expect_equal(round(Excess_Frequency(Fit, AP), 5), round(FQs, 5))
  expect_equal(sum(round(Fit$alpha, 5) == 1000), 4)

})


test_that("Fit_References with option PiecewisePareto", {
  covers <- c(1000, 1000, 1000)
  att_points <- c(1000, 2000, 5000)
  exp_losses <- c(100, 50, 10)
  thresholds <- c(4000, 10000)
  fqs <- c(0.04, 0.005)
  Fit <- Fit_References(covers, att_points, exp_losses, thresholds, fqs)
  expect_equal(Layer_Mean(fit, covers, att_points), exp_losses)
  expect_equal(Excess_Frequency(fit, thresholds), fqs)

  if (requireNamespace("lpSolve", quietly = TRUE)) {
    covers <- c(10000, 10000, 10000)
    att_points <- c(1000, 2000, 3000)
    exp_losses <- c(150, 70, 30)
    thresholds <- c(1000, 6000, 10000)
    fqs <- c(0.15, 0.003, 0.002)
    Fit <- Fit_References(covers, att_points, exp_losses, thresholds, fqs)
    expect_equal(Layer_Mean(Fit, covers, att_points), exp_losses)
    expect_equal(Excess_Frequency(Fit, thresholds), fqs)
  }

  if (requireNamespace("lpSolve", quietly = TRUE)) {
    covers <- c(10000, 10000, 10000)
    att_points <- c(1000, 2000, 3000)
    exp_losses <- c(150, 70, 30)
    thresholds <- c(1000, 6000, 10000, 11000)
    fqs <- c(0.15, 0.003, 0.002, 0.03)
    Fit <- Fit_References(covers, att_points, exp_losses, thresholds, fqs, ignore_inconsistent_references = TRUE)
    expect_equal(Layer_Mean(Fit, covers, att_points), exp_losses)
    expect_equal(Excess_Frequency(Fit, thresholds[-4]), fqs[-4])
  }

})




test_that("Panjer distribution", {
  Fit <- PPP_Model(FQ = 10, t = 1000, alpha = 2, dispersion = 1)
  expect_equal(Layer_Mean(Fit, 1, 0), 10)
  expect_equal(Layer_Var(Fit, 1, 0), 10)

  Fit <- PPP_Model(FQ = 10, t = 1000, alpha = 2, dispersion = 0.5)
  expect_equal(Layer_Mean(Fit, 1, 0), 10)
  expect_equal(Layer_Var(Fit, 1, 0), 5)

  Fit <- PPP_Model(FQ = 10, t = 1000, alpha = 2, dispersion = 2)
  expect_equal(Layer_Mean(Fit, 1, 0), 10)
  expect_equal(Layer_Var(Fit, 1, 0), 20)

})


test_that("Simulation, Sd and Var, dispersion = 1", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd')

  set.seed(1972)

  losses <- Simulate_Losses(Fit, 100000)
  layer_losses <- matrix(pmin(1000, pmax(losses - 2000, 0)), nrow = 100000)
  annual_losses <- apply(layer_losses, 1, function(x) sum(x, na.rm = T))

  expect_equal(round(sd(annual_losses), -1), round(Layer_Sd(Fit, 1000, 2000), -1))

  expect_equal(round(Layer_Sd(Fit, 1000, 2000), 3), 939.264)
  expect_equal(round(Layer_Var(Fit, 1000, 2000), 3), 882217.474)
})

test_that("Simulation, Sd and Var, dispersion = 0.63", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd', dispersion = 0.63)

  set.seed(1972)

  losses <- Simulate_Losses(Fit, 100000)
  layer_losses <- matrix(pmin(1000, pmax(losses - 2000, 0)), nrow = 100000)
  annual_losses <- apply(layer_losses, 1, function(x) sum(x, na.rm = T))

  expect_equal(round(sd(annual_losses), -1), round(Layer_Sd(Fit, 1000, 2000), -1))

  expect_equal(round(Layer_Sd(Fit, 1000, 2000), 3), 780.873)
  expect_equal(round(Layer_Var(Fit, 1000, 2000), 3), 609762.928)
})


test_that("Simulation, Sd and Var, dispersion = 2", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd', dispersion = 2)

  set.seed(1972)

  losses <- Simulate_Losses(Fit, 100000)
  layer_losses <- matrix(pmin(1000, pmax(losses - 2000, 0)), nrow = 100000)
  annual_losses <- apply(layer_losses, 1, function(x) sum(x, na.rm = T))

  expect_equal(round(sd(annual_losses) - Layer_Sd(Fit, 1000, 2000), -1), 0)

  expect_equal(round(Layer_Sd(Fit, 1000, 2000), 3), 1272.235)
  expect_equal(round(Layer_Var(Fit, 1000, 2000), 3), 1618581.110)
})
