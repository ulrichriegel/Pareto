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
