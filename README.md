# MetabolicSurv
R package : A biomarker validation approach for predicting survival using metabolic signature, this package develope biomarker signature for metabolic data. It contains a set of functions and cross validation methods  to validate and select biomarkers when the outcome of interest is survival. The package can handle prognostic factors and mainly metabolite matrix as input, the package can served as biomarker validation tool.

## Why use the package
* It can be used with any form of high dimensional/omics data such as: Metabolic data, Gene expression matrix, incase you dont have a data it can simulate hypothetical scinerio of a high dimensional data based on the desired biological parameters
* It developed any form of signature from the high dimensional data to be used for other purpose
* It also employs data reduction techniques such as PCA, PLS and Lasso 
* It classifies subjects based on the signatures into Low and high risk group
* It incorporate the use of subject prognostic information for the to enhance the biomarker for classification
* It gives information about the surival rate of subjects depending on the classificationhjkkjkllkkl



## Installation

You can install the released version of MetabolicSurv from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MetabolicSurv")

```

## A quick Demostration to solve a problem

``` r
"Problem of interest"

"Given a set of subjects with known riskscores and prognostic features how can we use this information to obtain their risk of surving and what group does each respective subject belongs to?"

```

``` r
##  Loading the package
library("MetabolicSurv")

##  Loading one of the inbuilt data
data(DataHR)
names(DataHR)

##  This function does Classification, Survival Estimation and Visualization
Result = EstimateHR(Risk.Scores=DataHR[,1],Data.Survival=DataHR[,2:3]
,Prognostic=DataHR[,4:5],Plots=FALSE,Quantile=0.50)

## Survival information
Result$SurvResult


## Group information
Result$Riskgroup
```
