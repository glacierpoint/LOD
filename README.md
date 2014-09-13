LOD
===

Impute exposure values for measurements below the limit of detection assuming a Generalized Gamma Distribution

In the environmental health sciences, measurements of toxic exposures are often constrained by a lower limit called the limit of detection (LOD), with observations below this limit called non-detects. Although valid inference may be obtained by excluding non-detects in the estimation of exposure effects, this practice can lead to substantial reduction in power to detect a significant effect, depending on the proportion of censoring and the closeness of the effect size to the null value. Therefore, a variety of methods have been commonly used in the environmental science literature to substitute values for the non-detects for the purpose of estimating exposure effects, including ad hoc values such as LOD/2,LOD/sqrt(2) and LOD. Another method substitutes the expected value of the non-detects, i.e., E[X|X ≤ LOD] but this requires that the inference be robust to mild miss-specifications in the distribution of the exposure variable.  In this R package, We use the generalized gamma distribution to estimate imputed values for the non-detects and avoid the risk of distribution miss-specification among the class of distributions represented by the generalized gamma distribution. Multiple imputated datasets are produces which can be analysed using multiple imputation based methods.

Load the package in R as follows:

library(devtools)

install_github("lod","glacierpoint")

Reference:

Arunajadai SG., Rauh VA (2012) Handling covariates subject to limits of detection in regression. Volume 19 issue 3 Pages 369-391.

