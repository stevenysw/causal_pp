# Citation
Steven Siwei Ye, Yanzhen Chen, Oscar Hernan Madrid Padilla. *Non-parametric interpretable score based estimation of heterogeneous treatment effects.* **2021+.**

# Methods
*  our proposed score-based method: low-dimensional data estimation and high-dimensional data estimation
*  endogeneous stratification
*  propensity-score matching (PSM)
*  prognostic-score matching

Note that all algorithms require the usage of nearestneighbour.m by Richard Brown. See details in https://www.mathworks.com/matlabcentral/fileexchange/12574-nearestneighbour-m.

For ADMM, we use parametric max-flow algorithm from "On Total Variation Minimization and Surface Evolution Using Parametric Maximum Flows" by Antonin Chambolle and Jérôme Darbon (https://link.springer.com/article/10.1007/s11263-009-0238-9). Users need to compile "TVexact" first to enable the "graphtv" function.

# Datasets
*  NMES Data
*  
