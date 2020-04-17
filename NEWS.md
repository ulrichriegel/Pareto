# Pareto 2.0.0

* PiecewisePareto_Match_Layer_Losses now returns a PPP_Model object. PPP stands for Panjer & Piecewise Pareto. The Panjer class contains the
  Poisson, the Negative Binomial and the Binomial distribution. A PPP_Model object contains the information required to specify a collective model
  with a Panjer distributed claim count and a Piecewise Pareto distributed severity.
* The package provides additional functions for PPP_Model objects:
    * PPP_Model_Exp_Layer_Loss: Calculates the expected loss of a reinsurance layer for a PPP_Model
    * PPP_Model_Layer_Var: Calculates the variance of the loss in a reinsurance layer for a PPP_Model
    * PPP_Model_Layer_Sd: Calculates the standard deviation of the loss in a reinsurance layer for a PPP_Model
    * PPP_Model_Excess_Frequency: Calculates the expected frequency in excess of a threshold for a PPP_Model
    * PPP_Model_Simulate: Simulates losses of a PPP_Model

# Pareto 1.1.5

* PiecewisePareto_Match_Layer_Losses now also works for only one layer
* Improved error handling in PiecewisePareto_Match_Layer_Losses

# Pareto 1.1.3

* Added maximum likelihood estimation of the alphas of a piecewise Pareto distribution.
* Allow for a different reporting threshold for each loss in Pareto_ML_Estimator_Alpha and in rPareto.
* Improved fitting algorithm in Pareto_ML_Estimator_Alpha.
* Better error handling in in Pareto_Find_Alpha_btw_FQ_Layer.

# Pareto 1.1.0

Stable version.
