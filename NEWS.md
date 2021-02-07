# Pareto 2.3.1

* Added distributions in function Local_Pareto_Alpha:
    * Pareto distribution
    * Generalized Pareto distribution
    * Piecewise Pareto distribution


# Pareto 2.3.0

* Vectorization of the following functions:
    * Pareto_Layer_Mean
    * Pareto_Layer_Var
    * Pareto_Layer_SM
    * Pareto_Extrapolation
    * Pareto_Find_Alpha_btw_Layers
    * Pareto_Find_Alpha_btw_FQ_Layer
    * Pareto_Find_Alpha_btw_FQs
    * PiecewisePareto_Layer_Mean (only parameters Cover and AttachmentPoint)
    * PiecewisePareto_Layer_SM (only parameters Cover and AttachmentPoint)
    * PiecewisePareto_Layer_Var (only parameters Cover and AttachmentPoint)
    * pPareto
    * dPareto
    * qPareto
    * pGenPareto
    * dGenPareto
    * qGenPareto
    * GenPareto_Layer_Mean
    * GenPareto_Layer_Var
    * GenPareto_Layer_SM

# Pareto 2.2.2

* Added function Fit_PML_Curve which fits a PPP_Model to a PML curve..

# Pareto 2.2.1

* Added the option to use weights in Pareto_ML_Estimator_Alpha, PiecewisePareto_ML_Estimator_Alpha and GenPareto_ML_Estimator_Alpha.

# Pareto 2.2.0

* Added function Fit_References for the piecewise Pareto distribution. This function fits a PPP model to the expected losses of 
  given reference layers and excess frequencies
* It is now possible to have layers with an expected loss of zero in PiecewisePareto_Match_Layer_Losses
* Improved handling of Frequencies and TotalLoss_Frequencies in PiecewisePareto_Match_Layer_Losses

# Pareto 2.1.0

* Added functions for the generalized Pareto distribution
* Added the class PGP_Model. PGP stands for Panjer & Generalized Pareto. A PGP_Model object contains the information to specify a 
  collective model with a Panjer distributed claim count and a generalized Pareto distributed severity
* The following functions have been replaced by generics for PPP_Models and PGP_Models:
    * PPP_Model_Exp_Layer_Loss has been replaced by Layer_Mean
    * PPP_Model_Layer_Var has been replaced by Layer_Var
    * PPP_Model_Layer_Sd has been replaced by Layer_Sd
    * PPP_Model_Excess_Frequency has been replaced by Excess_Frequency
    * PPP_Model_Simulate has been replaced by Simulate_Losses

# Pareto 2.0.0

* PiecewisePareto_Match_Layer_Losses now returns a PPP_Model object. PPP stands for Panjer & Piecewise Pareto. The Panjer class contains the
  Poisson, the Negative Binomial and the Binomial distribution. A PPP_Model object contains the information required to specify a collective
  model with a Panjer distributed claim count and a Piecewise Pareto distributed severity.
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
