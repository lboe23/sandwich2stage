# sandwich2stage
We provide functions that directly compute estimates of the sandwich variance formed by stacking the stage 1 and stage 2 model estimating equations for the two stage regression model setting of regression calibration. 
Regression calibration is often used when a true exposure variable of interest is difficult to measure. 
The stage 1 model represents a linear calibration model, which uses the observed data, including an error-prone measure of the exposure and other confounding variables, to estimate the true, unknown exposure value. 
The stage 2 model represents an outcome model, which is fit to the estimated exposure from stage 1 as well as other confounders. 
The sandwich variance estimate computed by our function provides an estimate of the variance of the stage 2 model parameters, adjusted for the uncertainty in the estimated exposure, since usual model-based standard errors will be too small.
This code can compute the sandwich variance for a subset of stage 2 models including linear models, generalized linear models (GLM), and Cox proportional hazards models, but this code can be adapted to work more generally for any model for which a sandwich variance estimator exists.  

Functions included in SandwichVar.R are as follows:

SandwichRegCal.lm() - computes the sandwich when the stage 2 model is a linear model (returns a list)

SandwichRegCal.glm() - computes the sandwich when the stage 2 model is a generalized linear model (returns a list)

SandwichRegCal.coxph() - computes the sandwich when the stage 2 model is a Cox PH model (returns a list)

print.SandwichRegCal() - print function for SandwichRegCal

coef.SandwichRegCal() - returns vector of coefficients from stage 1 and stage 2 models 

vcov.SandwichRegCal() - returns sandwich variance matrix from stacked estimating equation approach 

Each of SandwichRegCal.lm, SandwichRegCal.glm, and SandwichRegCal.coxph require the user to provide the fitted stage 1 (stage1.model) and stage 2 model (stage2.model) objects, as well as the variable names for the error-prone exposure variable (xstar) and estimated exposure variable (xhat) in the data. 
