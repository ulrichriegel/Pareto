library(Pareto)
context("test functions Pareto")

test_that("Pareto_Layer_Mean", {
  expect_equal(Pareto_Layer_Mean(8000, 2000, 2), 1600)
  expect_equal(Pareto_Layer_Mean(8000, 2000, 2, t = 1000), 400)
  expect_equal(Pareto_Layer_Mean(8000, 2000, 2, t = 5000, truncation = 10000), 4666.66666666667)
  expect_equal(Pareto_Layer_Mean(2000,1000,0, truncation = 5000, t = 500), 835.16520831955449)
})

test_that("Pareto_Layer_SM", {
  expect_equal(Pareto_Layer_SM(8000, 2000, 2), 6475503.2994728)
  expect_equal(Pareto_Layer_SM(8000, 2000, 2, t = 1000), 1618875.8248682)
  expect_equal(Pareto_Layer_SM(8000, 2000, 2, t = 5000, truncation = 10000), 23543145.370663)
  expect_equal(Pareto_Layer_SM(2000,1000,0, truncation = 5000, t = 1500), 2584318.901925616)
})

test_that("Pareto_Layer_Var", {
  expect_equal(Pareto_Layer_Var(8000, 2000, 2), 3915503.2994728)
  expect_equal(Pareto_Layer_Var(8000, 2000, 2, t = 1000), 1458875.8248682)
  expect_equal(Pareto_Layer_Var(8000, 2000, 2, t = 5000, truncation = 10000), 1765367.59288524)
  expect_equal(Pareto_Layer_Var(2000,1000,0, truncation = 5000, t = 500), 667015.32799764269)
})


################################
# test rPareto with truncation #
################################

test_that("rPareto with truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- 1000
  Cover <- 8000
  AttPoint <- 2000
  truncation <- 20000
  alpha <- 1.5

  expect_equal(Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t, truncation = truncation), 700.14314962210699)

  losses <- rPareto(NumberOfSimulations, t, alpha, truncation = truncation)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t, truncation = truncation), 2)

  expect_equal(ratio, 1)
})

############################################
# test rPareto without truncation #
############################################

test_that("rPareto without truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- 1000
  Cover <- 8000
  AttPoint <- 2000
  alpha <- 1.5

  expect_equal(Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t), 781.75803033941929)

  losses <- rPareto(NumberOfSimulations, t, alpha)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t), 2)

  expect_equal(ratio, 1)
})


test_that("Pareto_Extrapolation", {
  expect_equal(Pareto_Extrapolation(1000, 1000, 2000, 2000, alpha = 2), 0.5)
  expect_equal(Pareto_Extrapolation(1000, 1000, 2000, 2000, alpha = 2, ExpLoss_1 = 1000), 500)
  expect_equal(Pareto_Extrapolation(1000, 1000, 2000, 2000, alpha = 2, truncation = 3000, ExpLoss_1 = 1000), 142.85714285714292)
  expect_equal(Pareto_Extrapolation(Inf, 1000, Inf, 2000, alpha = 0, truncation = 3000, ExpLoss_1 = 100), 20.975411735345432)
})

test_that("Pareto_Find_Alpha_btw_FQ_Layer", {
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 2000, 2000, 100), 2.9330042139247037)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 1000, 1000, 500), 2)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(2000, 0.25, 1000, 1000, 500), 2)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 1000, 1000, 500, truncation = 5000), 1.8363401129702193)
})

test_that("Pareto_Find_Alpha_btw_Layers", {
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 100), 1)
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 50), 2)
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 50, truncation = 5000), 1.3871313763147217)
})

test_that("Pareto_Find_Alpha_btw_FQs", {
  expect_equal(Pareto_Find_Alpha_btw_FQs(1000, 1, 2000, 0.5), 1)
  expect_equal(Pareto_Find_Alpha_btw_FQs(2000, 0.25, 1000, 1), 2)
  expect_equal(Pareto_Find_Alpha_btw_FQs(2000, 0.25, 1000, 1, truncation = 4000), 1.5849625007211574)
})


##################################
# test Pareto_ML_Estimator_alpha #
##################################

test_that("Pareto_ML_Estimator_alpha", {
  losses <- c(1622.49986584698, 1025.1735923535, 1142.67198754259, 1598.2131674777, 1369.79742768744, 1006.5249344124,
              2019.3663238659, 1007.2758879241, 1377.79293040511, 1605.21438984656)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000), 4), 3.4060)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, truncation = 3000), 4), 2.9601)
})


