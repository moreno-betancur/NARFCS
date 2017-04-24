# `mice`: The NARFCS procedure for sensitivity analyses
M. Moreno-Betancur, F. Leacy, D. Tompsett, I. White  

([click here if you don't see the Latex formulas in this document properly](https://rawgit.com/moreno-betancur/NARFCS/master/README.html))

Here we illustrate how to perform sensitivity analyses using the Not At Random Fully Conditional Specification (NARFCS) procedure with the package `mice`. This procedure was first studied by Leacy (2016), and Daniel et al. (in preparation) provide further recommendations about its application in practice. 

Here we focus on the syntax for implementation of NARFCS using `mice`, specifically through dedicated imputation methods. At the moment, methods are available for continuous and binary variables only. These are called `normSens` and `logregSens` but require additional arguments to be passed on to `mice`, as described below.

## Loading and looking at example dataset
To illustrate, we use a simulated dataset `datmis` (available for download [here](https://rawgit.com/moreno-betancur/NARFCS/master/datmis.csv)) with one binary variable $X$ and a continuous variable $Y$, both having missing data, and $M_x$ and $M_y$ being the corresponding missingness indicators (i.e. equal to 1 if variable is missing and 0 otherwise).


```r
datmis<-read.csv("datmis.csv")
head(datmis)
```

```
##    X          Y Mx My
## 1  0   3.693582  0  0
## 2  0 -12.398012  0  0
## 3  0 -17.878958  0  0
## 4  0 -19.245551  0  0
## 5 NA         NA  1  1
## 6  0   9.126956  0  0
```

Let's look at the percentage of missing values in $X$ and $Y$, and the percentage of incomplete cases: 

```r
100 * c(mean(datmis$Mx), mean(datmis$My))
```

```
## [1] 19.0 21.6
```

```r
100 * c(mean(datmis$Mx == 1 | datmis$My == 1))
```

```
## [1] 29.4
```

Assume that the aim is to estimate the coefficient of $X$ in a linear regression of $Y$ on $X$. Before illustrating NARFCS, we provide a recap of a standard FCS analysis one may do in this case. This will be the point of departure for the NARFCS sensitivity analysis and will also help to draw some parallels between the arguments provided to `mice` in the FCS analysis and those additional ones required for NARFCS.

We first need to turn X into a factor for `mice` to recognise this a binary variable:

```r
datmis$X<-as.factor(datmis$X)
```

## FCS analysis

### FCS univariate imputation models
A typical FCS analysis with `mice` would be to use a logistic model to impute $X$ and a linear model to impute $Y$, with $Y$ included as predictor in the imputation model for $X$ and vice-versa. Thus, the univariate imputation model for $X$ would specify the conditional expectation of $X$ given the other variables as follows:
$$\text{logit}\{E(X|Y,M_y,M_x)\}=\beta_{X0}+ \beta_{XY} Y.$$
Similary, the model for $Y$ would specify:
$$ E(Y|X,M_x,M_y)=\beta_{Y0}+ \beta_{YX} X.$$

### Running it with mice
The FCS analysis with `mice` can be done by specifying the linear predictor of each of these models in either of two ways: using the matrix syntax (through the `predictorMatrix` argument) or using the formula syntax (through the `form` argument). 

#### Matrix syntax
With matrix syntax, the analysis is perfomed as follows:

```r
library(mice)

#Set-up predictor matrix
predMatrix<-matrix(rep(0,16),ncol=4)
predMatrix[1,]<-c(0,1,0,0)
predMatrix[2,]<-c(1,0,0,0)

predMatrix
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    0    1    0    0
## [2,]    1    0    0    0
## [3,]    0    0    0    0
## [4,]    0    0    0    0
```

```r
#FCS analysis with matrix syntax
impFCS<-mice(datmis,m=5, method=c("logreg","norm","",""), predictorMatrix=predMatrix,seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impFCS,lm(Y~X)))$qbar
```

```
## (Intercept)          X2 
##   -1.270058   13.169915
```

#### Formula syntax

The formula syntax is sometimes more straightforward, especially when it comes to interactions (though for now we have none - see further below):

```r
#FCS analysis using formula syntax
impFCS<-mice(datmis,m=5, method=c("logreg","norm","",""), form=c("~1+Y","~1+X","",""),seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impFCS,lm(Y~X)))$qbar
```

```
## (Intercept)          X2 
##   -1.270058   13.169915
```

At a given MICE iteration, the $\beta$-parameters are estimated from the observed and current imputed values of all variables, but ultimately depend on the expectations and associations of the variables in the observed data. One can thus ask what would happen if some aspects of the distribution of the missing data differed from that of the observed data. Answering this question is the aim of a sensitvity analysis which we will perform using the NARFCS procedure.

## NARFCS analysis to assess sensitivity

### NARFCS univariate imputation models
NARFCS (Leacy 2016) is a generalisation to the case with multiple incomplete variables of the so-called *delta-adjustment* method. We have previously used this approach to perform sensitivity analyses when there is a single univariate imputation model, both for continuous (Moreno-Betancur and Chavance, 2016) and binary (Leacy et al. 2017; Moreno-Betancur et al. 2015) variables, and various other examples exist in the literature.

A senstivity analysis using the NARFCS procedure aims at shifting the imputations drawn at each iteration by a quantity determined by a linear combination of the variables already included in the imputation model. Since only the imputed values are shifted, and not the observed, this shift only concerns individuals with the corresponding value missing (i.e. with $M_x=1$ if $X$ is being imputed and $M_y=1$ if $Y$ is being imputed). Formally, this comes to specifying a univariate imputation model with an additional term representing an interaction between the missingness indicator and the chosen linear combination for the shift. In our example, this could be as follows:

$$ \text{logit}\{E(X|Y,M_y,M_x)\}=\underbrace{\beta_{X0}+ \beta_{XY} Y}_{\text{IDENTIFIABLE PART}} +\underbrace{M_x*(\delta_{X0}+\delta_{XY}Y)}_\text{UNIDENTIFIABLE PART}.$$
and 
$$ E(Y|X,M_x,M_y)=\underbrace{\beta_{Y0}+ \beta_{YX} X}_{\text{IDENTIFIABLE PART}} +\underbrace{M_y*(\delta_{Y0}+\delta_{YX}X)}_\text{UNIDENTIFIABLE PART}.$$

The $\delta$-parameters that make up the linear combination in the shift are called *sensitvity parameters*. These values pertain to the way in which the distribution of the missing data departs from that of the observed data, and as such are not estimable from the observed data, that is, they are unidentifiable. Thus, the NARFCS model consists of an identfiable part, with $\beta$-parameters estimated from observed data, and an unidentifiable part, with $\delta$- (sensitivity) parameters. 

Typically plausible ranges for the latter would be elicited from experts. For this, it helps to understand their interpretation. Re-writing the equations above, we see that in this example the $\delta$-parameters correspond to shifts in both the intercept and slope of the model (in the corresponding scale) for the missing data:

$$ \text{logit}\{E(X|Y,M_y,M_x)\}= (\beta_{X0}+\delta_{X0}M_x)+ (\beta_{XY} +\delta_{XY}M_x)Y$$
and 
$$ E(Y|X,M_x,M_y)= (\beta_{Y0}+\delta_{Y0}M_y)+ (\beta_{YX} +\delta_{YX}M_y)X$$

Of course, given that each sensitivity parameter will be varied across a range of values, for the sake of parsimony we will often assume that many of the $\delta$-parameters are zero (i.e. exclude the term from the unidentifiable part), and we will vary only a key subset of them across a range of non-zero values according to the context and target parameter of interest. Indeed, each vector of values for the vector of sensitvity parameters will yield one set of results for the analysis. Thus, a key aspect of sensitivity analyses is how to process and report the large volume of results that arises from these analyses. This is out of the scope of this document and we refer the reader to the aforementioned references for some examples. 

Here, we show how to perform the NARFCS analysis with `mice` for one given vector of values for the 
$\delta$-parameters in our example. Specifically, we set $\delta_{XY}=0$ (that is, this term is excluded), $\delta_{Y0}=2$, $\delta_{YX}=5$ and $\delta_{X0}=-4$.

### Running it with mice
The specification of the predictors in the identifiable part of the model is still done the same way as in the FCS analysis, using the `predictorMatrix` or `form` arguments. To specify the unidentifable part, the user must do three changes to the way `mice` was called for the FCS analysis:

1. The predictors in the unidentifiable part of the NARFCS model are specified through either of two additional arguments, `predictorSens` or `formSens`, in a similar way as for the identifiable part. 

2. The user needs to specify the values of the sensitivity parameters through the argument `parmSens`, by passing a list of vectors, with the jth entry being the vector of sensitivity parameters for the imputation model of the jth variable in the dataset.  The size of the jth entry must match the number of predictors speficied for unidentified part of the imputation model of the jth variable through 
`predictorSens` or `formSens`, plus 1 (for the intercept). If one of the predictors is a categorical variable, say with $K$ categories, then only one sensitivity parameter value needs to be passed on for that predictor and this value is used for all $K$ categories, i.e. it is replicated $K$ times. 
If one wants a different value for each category, then one can always add the corresponding 
dummy variables to the dataset prior to calling `mice` and treat these as binary variables, 
with each assigned its own sensitivity parameter value. 

3. The imputation methods for the variables that have sensitivity paramters (in this case both $X$ and $Y$) are no longer `"norm"` and `"logreg"` but rather `"normSens"` and `"logregSens"`. No other imputation methods have been extended for NARFCS at the moment.

#### Matrix syntax
Using the matrix syntax, the NARFCS procedure with the above models and choice of parameters is done as follows:


```r
#Set-up predictor matrix for identifiable part as before:
predMatrix<-matrix(rep(0,16),ncol=4)
predMatrix[1,]<-c(0,1,0,0)
predMatrix[2,]<-c(1,0,0,0)

predMatrix
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    0    1    0    0
## [2,]    1    0    0    0
## [3,]    0    0    0    0
## [4,]    0    0    0    0
```

```r
#Set-up predictor matrix for unidentifiable part:
predSens<-matrix(rep(0,16),ncol=4)
predSens[1,]<-c(0,0,0,0)
predSens[2,]<-c(1,0,0,0)

predSens
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    0    0    0    0
## [2,]    1    0    0    0
## [3,]    0    0    0    0
## [4,]    0    0    0    0
```

```r
# NARFCS analysis using matrix syntax
impNARFCS<-mice(datmis,m=5, method=c("logregSens","normSens","",""), predictorMatrix=predMatrix,
                 predictorSens=predSens, parmSens=list(c(-4),c(2,5)), seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impNARFCS,lm(Y~X)))$qbar
```

```
## (Intercept)          X2 
##   0.7805584  10.8163702
```

Note that no predictors have been specified for $X$ in this part of the imputation model. Hence,
only an intercept shift will be included. If the user wants to exclude this intercept shift
(i.e. exclude the unidentifiable part of the model altogether), then it suffices to set the corresponding sensitivity parameter to zero.

If we extract the predictor matrix used for the unidentifiable part from the object returned by `mice`, we see the names of the colums are preceded by a colon (":").


```r
impNARFCS$predictorSens
```

```
##    :X :Y :Mx :My
## X   0  0   0   0
## Y   1  0   0   0
## Mx  0  0   0   0
## My  0  0   0   0
```
These names highlight that, contrary to the terms in `predictorMatrix`, the terms in `predictorSens` formally come into the model only as an interaction with the missingness indicator (i.e. only for the missing data) as described above.

#### Formula syntax
With the formula syntax, the analysis is performed as follows:

```r
# NARFCS analysis using formula syntax
impNARFCS<-mice(datmis,m=5, method=c("logregSens","normSens","",""),  form=c("~1+Y","~1+X","",""),
                  formSens=c("~1","~1+X","",""), parmSens=list(c(-4),c(2,5)), seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impNARFCS,lm(Y~X)))$qbar
```

```
## (Intercept)          X2 
##   0.7805584  10.8163702
```

#### A note on interactions
When including interactions between predictors in the unidentifiable part of the model, similar issues arise as in the typical FCS analysis. Specifically, using the matrix syntax, the user must use passive imputation and carefully set up the visiting sequence as described by [van Buuren and Groothuis-Oudshoorn (2011)](https://www.jstatsoft.org/article/view/v045i03/v45i03.pdf). The passive imputation requires in particular adding the variable(s) representing the interaction term (e.g. the product in the case of two continuous variables) to the dataset prior to calling `mice`. 

With the formula syntax none of this is required as the interaction term can be specified directly in the formula. In this case the two approaches will however yield slighlty different results (even when using the same seed) because the random values used to initialise missing values of the interaction variable (used by the matrix syntax) will be different to the product of the random values used to initialise missing values of the interacting variables (used by the formula syntax).

## References

Leacy FP. Multiple imputation under missing not at random assumptions via fully conditional 
specification 2016 (PhD thesis).

Tompsett D, Leacy FP, Moreno-Betancur M, White IR. On the use of the not at random fully conditional
specification procedure (NARFCS) in practice (in preparation).

Moreno-Betancur M, Chavance M. Sensitivity analysis of incomplete longitudinal data departing from the missing at random assumption: Methodology and application in a clinical trial with drop-outs. *Statistical Methods in Medical Research* 2016; 25(4) 1471–1489.

Leacy FP, Floyd S, Yates TA, White IR. Analyses of Sensitivity to the Missing-at-Random Assumption Using Multiple Imputation With Delta Adjustment: Application to a Tuberculosis/HIV Prevalence
Survey With Incomplete HIV-Status Data. *American Journal of Epidemiology*, 185(4):304-315.

Moreno-Betancur M, Rey G, Latouche A. Direct likelihood inference and sensitivity analysis for competing risks regression with missing causes of failure. *Biometrics* 2015; 71(2):498-507.

#### Cite this document as:

Moreno-Betancur M, Leacy FP, Tompsett D, White I. "mice: The NARFCS procedure for sensitivity analyses" (Available at: https://rawgit.com/moreno-betancur/NARFCS/master/README.html)
