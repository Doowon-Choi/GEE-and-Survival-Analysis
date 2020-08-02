# GEE-and-Survival-Analysis   
This session contains R code for a generalized estimating equation (GEE) and survival analysis.
GEE is used to estimate the parameters of generalized linear model under correlated responses, and survival analysis is to model
time-to-event data.

For survival analysis, this repository covers Kaplan-Meier curve, Cox proportional hazards model, Aalen's additive regression model, and
a Weibull regression tree. 


## 1. Compare two KM curves  
![Two KM curves](https://user-images.githubusercontent.com/69023373/89114662-bfdcba80-d444-11ea-8779-14e57db61f3d.png)

## 2. Compare KM and Cox curves  
![Compare_KM_Cox](https://user-images.githubusercontent.com/69023373/89114681-db47c580-d444-11ea-878f-22f4980c64b0.png)

## 3. Aalen's additive regression model
This model thinks of the cumulatve hazard function for certain subject as a(t) + XB(t), where a(t) is a time-dependent intercept term, X is the covariates for the subject (possibly time-dependent), and B(t) is a time-dependent matrix of coefficients.

![Alen](https://user-images.githubusercontent.com/69023373/89114684-f4507680-d444-11ea-90a5-f2366ec9c13c.png)

## 4. Weibull regression tree   
This tree is to recusrively partition the covariate space and fit a weibull regression for each of the resulting subgroups. The idea is that the parameter(s) of the weibull regression can be different depending on the condition of other covariates, which identifies heterogeneous risks in different hidden subgroups.

In fact, this tree belongs to the class of model-based recursive partitiong (MOB) algorithm that selects split variable via generalized M-fluctuation tests (i.e., detect a change in model parameters) and obtains cut-point of the selected split variable via minimizing the objective function in the node (i.e., negative log-likelihood).

![Weibull_reg](https://user-images.githubusercontent.com/69023373/89116506-5a93c400-d45a-11ea-91f2-4564059aa5f8.png)

# Reference
Aalen, O. O. (1989). A linear regression model for the analysis of life times. *Statistics in medicine*, 8(8), 907-925.

Aalen, O. O. (1993). Further results on the non‐parametric linear regression model in survival analysis. *Statistics in medicine*, 12(17), 1569-1588.

Rusch, T., & Zeileis, A. (2013). Gaining insight with recursive partitioning of generalized linear models. *Journal of Statistical Computation and Simulation*, 83(7).
