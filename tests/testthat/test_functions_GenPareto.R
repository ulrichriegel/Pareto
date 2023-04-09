library(Pareto)
context("test functions GenPareto")

test_that("GenPareto_Layer_Mean", {
  expect_equal(GenPareto_Layer_Mean(8000, 2000, 2000, 2, 2), 1600)
  expect_equal(GenPareto_Layer_Mean(8000, 2000, 1000, 2, 2), 400)
  expect_equal(GenPareto_Layer_Mean(8000, 2000, 5000, 2, 1, truncation = 10000), 4619.7960825054124)
  expect_equal(GenPareto_Layer_Mean(2000, 1000, 500, 0.01, 3, truncation = 5000), 1308.8417333286236)
  expect_equal(GenPareto_Layer_Mean(8000, 2000, c(2000, 1000, 5000), 2, c(2, 2, 1), truncation = c(Inf, Inf, 10000)), c(1600, 400, 4619.7960825054124))
})

test_that("GenPareto_Layer_SM", {
  expect_equal(GenPareto_Layer_SM(8000, 2000, 2000, 2, 2), 6475503.2994728)
  expect_equal(GenPareto_Layer_SM(8000, 2000, 1000, 2.5, 2.5), 705034.42336793197)
  expect_equal(GenPareto_Layer_SM(8000, 2000, 5000, 0.8, 2, truncation = 10000), 27842963.469263397)
  expect_equal(GenPareto_Layer_SM(2000, 1000, 500, 0.1, 0.9, truncation = 1000), 0)
  expect_equal(GenPareto_Layer_SM(c(8000, 2000), c(2000, 1000), c(2000, 500), c(2, 0.1), c(2, 0.9), truncation = c(Inf, 1000)), c(6475503.2994728, 0))
})

test_that("GenPareto_Layer_Var", {
  expect_equal(GenPareto_Layer_Var(8000, 2000, 1000, 2.1, 4), 398133.28976341535)
  expect_equal(GenPareto_Layer_Var(8000, 2000, 1000, 2, 2), 1458875.8248682)
  expect_equal(GenPareto_Layer_Var(8000, 2000, t = 5000, 1.5, 2.7, truncation = 10000), 1869703.2726483792)
  expect_equal(GenPareto_Layer_Var(2000, 1000, 500, 3, 1, truncation = 5000), 196333.67104755351)
  expect_equal(GenPareto_Layer_Var(c(8000, 2000), c(2000, 1000), c(1000, 500), c(2.1, 3), c(4, 1), truncation = c(Inf, 5000)), c(398133.28976341535, 196333.67104755351))
})

test_that("pGenPareto", {
  expect_equal(pGenPareto(c(1:3) * 2000, t = (1:3) * 1000, alpha_ini = c(1, 1.5, 2), alpha_tail = 1:3, truncation = NULL), c(0.5, 0.67346938775510212, 0.78399999999999992))
})

test_that("dGenPareto", {
  expect_equal(dGenPareto(c(1:3) * 2000, t = (1:3) * 1000, alpha_ini = c(1, 1.5, 2), alpha_tail = 1:3, truncation = NULL), c(2.5000000000000001e-04, 1.3994169096209912e-04, 8.6400000000000027e-05))
})




###################################
# test rGenPareto with truncation #
###################################

test_that("rPareto with truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- 1000
  Cover <- 8000
  AttPoint <- 2000
  truncation <- 20000
  alpha_ini <- 1.5
  alpha_tail <- 1

  expect_equal(GenPareto_Layer_Mean(Cover, AttPoint,  t, alpha_ini, alpha_tail, truncation = truncation), 932.32300743380154)

  losses <- rGenPareto(NumberOfSimulations, t, alpha_ini, alpha_tail, truncation = truncation)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / GenPareto_Layer_Mean(Cover, AttPoint,  t, alpha_ini, alpha_tail, truncation = truncation), 2)

  expect_equal(ratio, 1)
})

############################################
# test rGenPareto without truncation       #
############################################

test_that("rPareto without truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- 1000
  Cover <- 8000
  AttPoint <- 2000
  alpha_ini <- 1.5
  alpha_tail <- 3

  expect_equal(GenPareto_Layer_Mean(Cover, AttPoint, t,  alpha_ini, alpha_tail), 411.38659320477501)

  losses <- rGenPareto(NumberOfSimulations, t,  alpha_ini, alpha_tail)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / GenPareto_Layer_Mean(Cover, AttPoint, t,  alpha_ini, alpha_tail), 2)

  expect_equal(ratio, 1)
})



#####################################
# test GenPareto_ML_Estimator_alpha #
#####################################

test_that("GenPareto_ML_Estimator_alpha", {
  losses <- c(1622.49986584698, 1025.1735923535, 1142.67198754259, 1598.2131674777, 1369.79742768744, 1006.5249344124,
              2019.3663238659, 1007.2758879241, 1377.79293040511, 1605.21438984656, 2579.4568112321, 4500.45681)
  expect_equal(GenPareto_ML_Estimator_Alpha(losses, 1000), c(2.1210190911501012, 2.5019159656778251))
  expect_equal(GenPareto_ML_Estimator_Alpha(losses, 1000, truncation = 10000), c(2.3410152490683958, 1.6923583375172637))

  reporting_thresholds <- rep(1000, 12)
  reporting_thresholds[1]<- 1500
  expect_equal(round(GenPareto_ML_Estimator_Alpha(losses, 1000, reporting_thresholds = reporting_thresholds), 5), c(2.54341, 2.15354))
  expect_equal(round(GenPareto_ML_Estimator_Alpha(losses, 1000, truncation = 10000, reporting_thresholds = reporting_thresholds), 5), c(2.91339, 1.41459))

  is.censored <- rep(FALSE, 12)
  is.censored[1:2] <- TRUE
  expect_equal(round(GenPareto_ML_Estimator_Alpha(losses, 1000, reporting_thresholds = reporting_thresholds, is.censored = is.censored), 5), c(1.73377, 2.90203))
  expect_equal(round(GenPareto_ML_Estimator_Alpha(losses, 1000, is.censored = is.censored), 5), c(1.47613, 3.68474))
  expect_equal(round(GenPareto_ML_Estimator_Alpha(losses, 1000, truncation = 10000, reporting_thresholds = reporting_thresholds, is.censored = is.censored), 5), c(1.88626, 1.93255))

  w <- rep(1, length(losses))
  index <- c(1, 4, 6)
  w[index] <- 2
  losses2 <- c(losses, losses[index])
  expect_equal(GenPareto_ML_Estimator_Alpha(losses, 1000, weights = w), GenPareto_ML_Estimator_Alpha(losses2, 1000))
  expect_equal(round(GenPareto_ML_Estimator_Alpha(losses, 1000, weights = w, truncation = 10000), 4), round(GenPareto_ML_Estimator_Alpha(losses2, 1000, truncation = 10000), 4))

})


