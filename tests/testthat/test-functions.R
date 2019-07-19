library(Pareto)

expect_equal(Pareto_Layer_Mean(8000, 2000, 2), 1600)
expect_equal(Pareto_Layer_Mean(8000, 2000, 2, t = 1000), 400)
expect_equal(Pareto_Layer_Mean(8000, 2000, 2, t = 5000, truncation = 10000), 4666.66666666667)

expect_equal(Pareto_Layer_SM(8000, 2000, 2), 6475503.2994728)
expect_equal(Pareto_Layer_SM(8000, 2000, 2, t = 1000), 1618875.8248682)
expect_equal(Pareto_Layer_SM(8000, 2000, 2, t = 5000, truncation = 10000), 23543145.370663)

expect_equal(Pareto_Layer_Var(8000, 2000, 2), 3915503.2994728)
expect_equal(Pareto_Layer_Var(8000, 2000, 2, t = 1000), 1458875.8248682)
expect_equal(Pareto_Layer_Var(8000, 2000, 2, t = 5000, truncation = 10000), 1765367.59288524)

expect_equal(PiecewisePareto_Layer_Mean(8000, 2000, 2, t = 1000), 400)
expect_equal(PiecewisePareto_Layer_Mean(8000, 2000, 2, t = 5000, truncation = 10000), 4666.66666666667)

expect_equal(PiecewisePareto_Layer_SM(8000, 2000, 2, t = 1000), 1618875.8248682)
expect_equal(PiecewisePareto_Layer_SM(8000, 2000, 2, t = 5000, truncation = 10000), 23543145.370663)

expect_equal(PiecewisePareto_Layer_Var(8000, 2000, 2, t = 1000), 1458875.8248682)
expect_equal(PiecewisePareto_Layer_Var(8000, 2000, 2, t = 5000, truncation = 10000), 1765367.59288524)



#####################################################
# test rPiecewisePareto with truncation_type = "wd" #
#####################################################

set.seed(1972)
NumberOfSimulations <- 1000000

t <- c(1000,3000,5000)
Cover <- 8000
AttPoint <- 2000
truncation <- 6000
alpha <- c(0.7,1.5,2)
truncation_type ="wd"

expect_equal(PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 868.729076663466)

losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
mean(xs_losses)

ratio <- round(mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 2)

expect_equal(ratio, 1)

#####################################################
# test rPiecewisePareto with truncation_type = "lp" #
#####################################################

set.seed(1972)
NumberOfSimulations <- 1000000

t <- c(1000,3000,5000)
Cover <- 8000
AttPoint <- 2000
truncation <- 6000
alpha <- c(0.7,1.5,2)
truncation_type ="lp"

expect_equal(PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 1255.52081303827)

losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
mean(xs_losses)

ratio <- round(mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 2)

expect_equal(ratio, 1)

############################################
# test rPiecewisePareto without truncation #
############################################

set.seed(1972)
NumberOfSimulations <- 1000000

t <- c(1000,3000,5000)
Cover <- 8000
AttPoint <- 2000
alpha <- c(0.7,1.5,2)

expect_equal(PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha), 1696.1079667876)

losses <- rPiecewisePareto(NumberOfSimulations, t, alpha)
xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
mean(xs_losses)

ratio <- round(mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha), 2)

expect_equal(ratio, 1)

