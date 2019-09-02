`mice`: The NARFCS procedure for sensitivity analyses
================
M. Moreno-Betancur, F. Leacy, D. Tompsett, I. White

This document illustrates how to perform sensitivity analyses using the Not At Random Fully Conditional Specification (NARFCS) procedure implemented in an extension of the R package `mice`.

Background
----------

Standard multiple imputation procedures to deal with missing data typically rely on the marginal and conditional distributions of variables in the observed data to impute the missing values. This is for instance the basis of the Fully Conditional Specification (FCS) approach to handle multiple incomplete variables, which is generally implemented using Multiple Imputation by Chained Equations (MICE) as in the R package `mice` ([van Buuren S, Groothuis-Oudshoorn K., 2011](https://www.jstatsoft.org/article/view/v045i03/v45i03.pdf)).

A *sensitivity analysis* in this setting aims at assessing how our results would change if the distribution of the missing data differed in a systematic way from that of the observed data, which is something that we cannot of course assess using the latter.

The NARFCS procedure was first formally studied by Leacy (2016), and generalises the so-called *delta-adjustment* sensitivity analysis method to the case with multiple incomplete variables within the FCS framework. In practical terms, the NARFCS procedure consists in shifting the imputations drawn at each iteration of MICE by a user-specified quantity that can vary across subjects, to reflect systematic departures of the missing data from the observed data distribution.

Tompsett et al. (in preparation) provide recommendations for the application of NARFCS in practice, particularly regarding the interpretation and elicitation of the sensitivity parameters required by the procedure. Here we present the syntax for applying NARFCS to a dataset using an extension of the package `mice` that includes dedicated imputation methods. At the moment, methods are available for continuous and binary variables only, called `normSens` and `logregSens` respectively. Their use requires additional arguments to be passed on to `mice`, as described below.

Installation of the extension
-----------------------------

The currently recommended version of the NARFCS extension of the package `mice` can be installed directly from GitHub using the `devtools` package, which must be first installed from CRAN by the user to proceed. Then it suffices to execute the following commands to install the NARFCS extension of `mice`:

``` r
library(devtools)
install_github("moreno-betancur/mice")
```

Loading and looking at example dataset
--------------------------------------

([click here if you see raw rather than compiled LaTeX formulas below](https://raw.githack.com/moreno-betancur/NARFCS/master/README.html))

To illustrate, we use a simulated dataset `datmis` (available for download [here](https://rawgit.com/moreno-betancur/NARFCS/master/datmis.csv)) which consists of:

-   a binary variable *X*, taking the values "No" (the reference level) and "Yes",
-   a continuous variable *Y*, and
-   a categorical variable *Z*, with 3 categories: "Cat0" (the reference level), "Cat1" and "Cat2".

``` r
datmis<-read.csv("datmis.csv")
head(datmis)
```

    ##      X          Y    Z
    ## 1   No   3.693582 Cat1
    ## 2   No -12.398012 Cat2
    ## 3   No -17.878958 Cat2
    ## 4   No -19.245551 Cat2
    ## 5 <NA>         NA Cat0
    ## 6   No   9.126956 Cat1

``` r
summary(datmis)
```

    ##     X             Y              Z      
    ##  No  :264   Min.   :-32.700   Cat0:167  
    ##  Yes :141   1st Qu.: -5.713   Cat1:166  
    ##  NA's: 95   Median :  3.491   Cat2:167  
    ##             Mean   :  3.524             
    ##             3rd Qu.: 12.249             
    ##             Max.   : 49.663             
    ##             NA's   :108

The variables *X* and *Y* are incomplete and *Z* is complete, and it is easily verified that a *complete case* analysis would discard around 30% of the records:

``` r
100 * c(sum(is.na(datmis$X)), sum(is.na(datmis$Y)), sum(is.na(datmis$Z)))/nrow(datmis)
```

    ## [1] 19.0 21.6  0.0

``` r
100 * sum(!is.na(datmis$X) & !is.na(datmis$Y) & !is.na(datmis$Z))/nrow(datmis)
```

    ## [1] 70.6

Suppose that the aim is to estimate the coefficient of *X* in a linear regression of *Y* on *X* and *Z*. Before illustrating NARFCS, we provide a recap of a standard FCS analysis one may do in this case. This will be the point of departure for the NARFCS sensitivity analysis and will also help to draw some parallels between the arguments provided to `mice` in the FCS analysis and those additional ones required for NARFCS.

FCS analysis
------------

### FCS univariate imputation models

A typical FCS analysis with `mice` would be to use a logistic model to impute *X* and a linear model to impute *Y*, both including *Z* as predictor, and with *Y* included as predictor in the imputation model for *X* and vice-versa. A usual, for mathematical purposes, *X* is taken to be the binary 0/1 indicator for the "Yes" category, and including *Z* as predictor means that we include the binary indicators *Z*<sub>1</sub> and *Z*<sub>2</sub> for categories "Cat1" and "Cat2" in the models. Thus, the univariate imputation model for *X* would specify the conditional expectation of *X* given the other variables as follows:
logit{*E*(*X*|*Y*, *Z*)} = *β*<sub>*X*0</sub> + *β*<sub>*X**Y*</sub>*Y* + *β*<sub>*X**Z*<sub>1</sub></sub>*Z*<sub>1</sub> + *β*<sub>*X**Z*<sub>2</sub></sub>*Z*<sub>2</sub>.
 Similarly, the model for *Y* would specify:
*E*(*Y*|*X*, *Z*)=*β*<sub>*Y*0</sub> + *β*<sub>*Y**X*</sub>*X* + *β*<sub>*Y**Z*<sub>1</sub></sub>*Z*<sub>1</sub> + *β*<sub>*Y**Z*<sub>2</sub></sub>*Z*<sub>2</sub>.
 and Var(*Y*|*X*, *Z*)=*σ*<sup>2</sup>, assuming constant variance.

### Running FCS with mice

This FCS analysis can be performed with `mice` by specifying the linear predictor of each of these models in either of two ways: using the matrix syntax (through the `predictorMatrix` argument) or using the formula syntax (through the `form` argument).

#### Matrix syntax for FCS

With the matrix syntax, the analysis requires specifying the predictors of each imputation model through a matrix of dimension equal to the number of variables in the dataset (in this case 3). A cell with a 1 indicates that the column variable is included as predictor in the imputation model of the row variable, with both rows and columns representing the variables in the order in which they appear in the dataset. Other cells are set to zero:

``` r
library(mice)
```

    ## Loading required package: lattice

``` r
#Set-up predictor matrix
predMatrix<-matrix(rep(0,9),ncol=3)
colnames(predMatrix)<-rownames(predMatrix)<-names(datmis)
predMatrix["X",]<-c(0,1,1)
predMatrix["Y",]<-c(1,0,1)

predMatrix
```

    ##   X Y Z
    ## X 0 1 1
    ## Y 1 0 1
    ## Z 0 0 0

``` r
#FCS analysis with matrix syntax
impFCS<-mice(datmis,m=5, method=c("logreg","norm",""), predictorMatrix=predMatrix,seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impFCS,lm(Y~X+Z)))$qbar
```

    ## (Intercept)          X2          Z2          Z3 
    ##   15.958084    6.594639  -13.450028  -26.004113

**Note 1**: For categorical variables, the function `pool` names the coefficients in the output according to the rank of the category. Thus "X2" corresponds to the second category of *X*, and "Z2" and "Z3" to the second and third category of *Z*, respectively. It is thus important to note that, if not instructed otherwise, R automatically orders the levels of a factor with the reference level first and the rest following in alphabetical order. We can check the order of the levels of a factor variable using the `levels` function:

``` r
levels(datmis$X)
```

    ## [1] "No"  "Yes"

``` r
levels(datmis$Z)
```

    ## [1] "Cat0" "Cat1" "Cat2"

To set the desired reference level, the user can use the `relevel` function (see `?relevel`).

#### Formula syntax for FCS

The formula syntax is sometimes more straightforward to use, especially when it comes to interactions (though for now we have none - see **Note 2** below). Each imputation model is specified using a formula in the same manner used for specifying models in R's standard model-fitting functions:

``` r
#FCS analysis using formula syntax
impFCS<-mice(datmis,m=5, method=c("logreg","norm",""), form=c("~1+Y+Z","~1+X+Z",""),seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impFCS,lm(Y~X+Z)))$qbar
```

    ## (Intercept)          X2          Z2          Z3 
    ##    15.98365     6.58963   -13.48023   -26.04045

At a given MICE iteration, the *β*-parameters in the imputation models above are estimated from the observed and current imputed values of all variables, but ultimately depend on the expectations and associations of the variables in the observed data. One can thus ask what would happen if some aspects of the distribution of the missing data differed from what can be inferred from the observed data. Answering this question is the aim of a sensitivity analysis, which we will perform using the NARFCS procedure.

NARFCS analysis to assess sensitivity
-------------------------------------

### NARFCS univariate imputation models

NARFCS (Leacy 2016) is a generalisation to the case with multiple incomplete variables of the so-called *delta-adjustment* method for sensitivity analyses. In the case of one incomplete variable, we have previously used the method both for continuous ([Moreno-Betancur and Chavance, 2016](http://journals.sagepub.com/doi/abs/10.1177/0962280213490014)) and binary ([Leacy et al. 2017](https://academic.oup.com/aje/article-lookup/doi/10.1093/aje/kww107); [Moreno-Betancur et al. 2015](http://onlinelibrary.wiley.com/doi/10.1111/biom.12295/full)) variables, and various other examples exist in the literature.

A sensitivity analysis using the NARFCS procedure aims at shifting the imputations drawn at each iteration by a constant or, more generally, by a quantity determined by a linear combination of the variables already included in the imputation model. When imputing a binary variable, it is the linear predictor of the imputation model that is shifted prior to drawing the imputed value.

Since only the imputed values are shifted, and not the observed, this shift only concerns individuals with the corresponding value missing. That is, if *M*<sub>*X*</sub> and *M*<sub>*Y*</sub> denote the missingness indicators of *X* and *Y* (i.e. *M*<sub>*X*</sub> = 1 if *X* is missing and *M*<sub>*X*</sub> = 0 otherwise, and similarly for *M*<sub>*Y*</sub>) then the shift concerns only individuals with *M*<sub>*X*</sub> = 1 if *X* is being imputed, and *M*<sub>*Y*</sub> = 1 if *Y* is being imputed. Formally, this amounts to specifying, for each incomplete variable, a univariate imputation model with an additional term representing an interaction between the missingness indicator and the chosen linear combination for the shift. In our example, this could be as follows:

and

$$ E(Y|X,M\_Y)=\\underbrace{\\beta\_{Y0}+ \\beta\_{YX} X + \\beta\_{YZ\_{1}} Z\_1+\\beta\_{YZ\_{2}} Z\_2}\_{\\text{IDENTIFIABLE PART}} +\\underbrace{M\_Y\*(\\delta\_{Y0}+\\delta\_{YX}X+\\delta\_{YZ\_1}Z\_1+\\delta\_{YZ\_2}Z\_2)}\_\\text{UNIDENTIFIABLE PART}.$$

The *δ*-parameters that make up the linear combination in the shifts are called *sensitivity parameters*. These values describe the way in which the distribution of the missing data departs from that of the observed data, and as such are not estimable from the observed data, that is, they are unidentifiable. Thus, the NARFCS model consists of an identifiable part, with *β*-parameters estimated from observed data, and an unidentifiable part, with *δ*- (sensitivity) parameters.

Typically, plausible ranges for the latter would be elicited from experts. Tompsett et al. (in preparation) provide a detailed study of the issue of sensitivity parameter elicitation in NARFCS, including methodology and guidelines. It is important to note here that sensitivity parameters must be specified on the log-odds scale when imputing a binary variable.

Of course, given that each sensitivity parameter will be varied across a range of values, for the sake of parsimony we will often assume that many of the *δ*-parameters are zero (i.e. exclude the term from the unidentifiable part), and we will vary only a key subset of them across a range of non-zero values according to the context and target parameter of interest. Indeed, each set of values for the set of sensitivity parameters will yield one set of results for the analysis. Thus, a key aspect of sensitivity analyses is how to process and report the large volume of results that arises from these analyses. This is beyond the scope of this document and we refer the reader to the aforementioned references for some examples.

Here, we show how to perform the NARFCS analysis with `mice` for one given set of values for the *δ*-parameters in our example. Specifically, in the model for *X*, we set *δ*<sub>*X*0</sub> = −4 and *δ*<sub>*X**Y*</sub> = *δ*<sub>*X**Z*<sub>1</sub></sub> = *δ*<sub>*X**Z*<sub>2</sub></sub> = 0; that is, all terms except the intercept are excluded, so the shift when imputing *X* is just a constant. For *Y*, we set *δ*<sub>*Y*0</sub> = 2, *δ*<sub>*Y**X*</sub> = 0 (i.e. the *X*-term is excluded), *δ*<sub>*Y**Z*<sub>1</sub></sub> = 1 and *δ*<sub>*Y**Z*<sub>2</sub></sub> = −3.

### Running NARFCS with mice

The specification of the predictors in the identifiable part of the model is still done the same way as in the FCS analysis, using the `predictorMatrix` or `form` arguments. To specify the unidentifiable part, the user must make three changes to the way `mice` was called for the FCS analysis:

1.  The user must specify special imputation methods for the incomplete variables that are to be imputed using NARFCS, which are not necessarily (nor usually) all. In our example we will apply NARFCS to *X* and *Y*. Instead of `"norm"` and `"logreg"`, the user must use `"normSens"` and `"logregSens"` according to the type (continuous or binary, respectively) of the incomplete variable. No other imputation methods have been extended for NARFCS at the moment.

2.  The predictors in the unidentifiable part of the NARFCS model are specified through either of two additional arguments, `predictorSens` or `formSens`, in a similar way as for the identifiable part. By default only an intercept term is included.

3.  The values of the sensitivity parameters have to be specified through the argument `parmSens`, which is a list object with one element per variable in the dataset in the order in which they appear there. Its elements are lists of vectors, specified as follows:

-   For variables in the dataset that are complete or to be imputed via standard FCS, the corresponding element is a list with an empty string (`""`) as unique element, i.e. `list("")` (this is the default for all variables).

-   For variables to be imputed via NARFCS, it is a list of length equal to 1 (for the intercept) plus the number of variables included as predictors in the unidentifiable part of the imputation model, as specified through `predictorSens` or `formSens`. The first component is the sensitivity parameter corresponding to the intercept term, which is always included, and is thus a scalar (i.e. a vector of length 1). The following components are sensitivity parameter vectors, of length 1 except for factors (see below), corresponding to each of the predictor variables. If using the matrix syntax, these must be provided in the order in which they appear in the dataset (which is the order in which they appear in the matrix). If using the formula syntax, they are provided in the order in which they appear in the corresponding formula. For a factor predictor variable with *K* categories, the element is a vector of length *K* − 1 containing the sensitivity parameters for the *K* − 1 indicator variables of the non-reference categories, following the order of the levels of that factor variable (see related **Note 1** above). These parameters are interpreted relative to the reference category, as usual. As a reminder, sensitivity parameters must be specified in the log-odds scale when imputing a binary variable (`"logregSens"` method).

Since the specification of the `parmSens` argument is delicate, the program runs various checks to verify that the lengths of the lists and vectors are consistent with what was provided through `predictorSens` or `formSens` arguments, and provides informative error messages if this is not the case. When the procedure runs successfully and `printFlag=TRUE` (the default), the sensitivity parameters used for each predictor variable/factor level in each imputation model are printed out on the console so that the user can check that these correspond to what was intended. Carefully checking this output is highly recommended.

We will now illustrate how, together, these three steps allow us to apply NARFCS in our example using each possible syntax.

#### Matrix syntax for NARFCS

Using the matrix syntax, the NARFCS procedure with the above models and choice of parameters is performed as follows:

``` r
#Set-up predictor matrix for identifiable part as before:
predMatrix<-matrix(rep(0,9),ncol=3)
colnames(predMatrix)<-rownames(predMatrix)<-names(datmis)
predMatrix["X",]<-c(0,1,1)
predMatrix["Y",]<-c(1,0,1)

predMatrix
```

    ##   X Y Z
    ## X 0 1 1
    ## Y 1 0 1
    ## Z 0 0 0

``` r
#Set-up predictor matrix for unidentifiable part:
predSens<-matrix(rep(0,9),ncol=3)
colnames(predSens)<-paste(":",names(datmis),sep="")
rownames(predSens)<-names(datmis)

predSens["X",]<-c(0,0,0)
predSens["Y",]<-c(0,0,1)

predSens
```

    ##   :X :Y :Z
    ## X  0  0  0
    ## Y  0  0  1
    ## Z  0  0  0

``` r
#Set-up list with sensitivity parameter values

pSens<-rep(list(list("")), ncol(datmis))
names(pSens)<-names(datmis)
pSens[["X"]]<-list(-4)
pSens[["Y"]]<-list(2,c(1,-3))

pSens
```

    ## $X
    ## $X[[1]]
    ## [1] -4
    ## 
    ## 
    ## $Y
    ## $Y[[1]]
    ## [1] 2
    ## 
    ## $Y[[2]]
    ## [1]  1 -3
    ## 
    ## 
    ## $Z
    ## $Z[[1]]
    ## [1] ""

``` r
#NARFCS analysis using matrix syntax
impNARFCS<-mice(datmis,m=5, method=c("logregSens","normSens",""), predictorMatrix=predMatrix,
                 predictorSens=predSens, parmSens=pSens, seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impNARFCS,lm(Y~X+Z)))$qbar
```

    ## (Intercept)          X2          Z2          Z3 
    ##   19.040240    3.314205  -14.497868  -28.276516

Note that no predictors have been specified for *X* in this part of the imputation model. Hence, only an intercept shift will be included.

If we extract the `predictorSens` matrix from the object returned by `mice`, we see that the names of the columns are preceded by a colon (":"), and this will be the case even if we do not pre-specify these names ourselves as in the code above.

``` r
impNARFCS$predictorSens
```

    ##   :X :Y :Z
    ## X  0  0  0
    ## Y  0  0  1
    ## Z  0  0  0

These names highlight that, contrary to the terms in `predictorMatrix`, the terms in `predictorSens` formally come into the model only as an interaction with the missingness indicator (i.e. only for the missing data) as described previously.

#### Formula syntax for NARFCS

With the formula syntax, the analysis is performed as follows:

``` r
#NARFCS analysis using formula syntax
impNARFCS<-mice(datmis,m=5, method=c("logregSens","normSens",""),  form=c("~1+Y+Z","~1+X+Z",""),
                  formSens=c("~1","~1+Z",""), parmSens=pSens, seed=234235,print=F)

#Pool results from linear regression from imputed datasets:
pool(with(impNARFCS,lm(Y~X+Z)))$qbar
```

    ## (Intercept)          X2          Z2          Z3 
    ##   19.357351    3.457046  -15.519064  -28.400211

As with the matrix syntax, only an intercept shift will be included for *X* in the unidentifiable part of the model since no other predictors have been specified. This is however made more explicit with the formula syntax by the specification of the non-empty formula `"~1"`.

**Note 2:** When including interactions between predictors in the unidentifiable part of the model, similar issues arise as in the typical FCS analysis. Specifically, using the matrix syntax, the user must use passive imputation and carefully specify the visiting sequence as described by [van Buuren and Groothuis-Oudshoorn (2011)](https://www.jstatsoft.org/article/view/v045i03/v45i03.pdf). The passive imputation requires in particular adding the variable(s) representing the interaction term (e.g. the product in the case of two continuous variables) to the dataset prior to calling `mice`. With the formula syntax none of this is required as the interaction term can be specified directly in the formula. In this case the two approaches will yield slightly different results (even when using the same seed) because the random values used to initialise missing values of the interaction variable (used by the matrix syntax) will be different to the product of the random values used to initialise missing values of the interacting variables (used by the formula syntax).

References
----------

Leacy FP. *Multiple imputation under missing not at random assumptions via fully conditional specification* (PhD thesis). University of Cambridge, 2016.

[Tompsett DM, Leacy FP, Moreno-Betancur M, Heron J, White IR.](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.7643) On the use of the not at random fully conditional specification procedure (NARFCS) in practice. *Statistics in Medicine* 2018; Epub ahead of print 2 April 2018.

[van Buuren S, Groothuis-Oudshoorn K.](https://www.jstatsoft.org/article/view/v045i03/v45i03.pdf) mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software* 2011; 45(3), 1-67.

[Moreno-Betancur M, Chavance M.](http://journals.sagepub.com/doi/abs/10.1177/0962280213490014) Sensitivity analysis of incomplete longitudinal data departing from the missing at random assumption: Methodology and application in a clinical trial with drop-outs. *Statistical Methods in Medical Research* 2016; 25(4) 1471-1489.

[Leacy FP, Floyd S, Yates TA, White IR.](https://academic.oup.com/aje/article-lookup/doi/10.1093/aje/kww107) Analyses of Sensitivity to the Missing-at-Random Assumption Using Multiple Imputation With Delta Adjustment: Application to a Tuberculosis/HIV Prevalence Survey With Incomplete HIV-Status Data. *American Journal of Epidemiology* 2017, 185(4):304-315.

[Moreno-Betancur M, Rey G, Latouche A.](http://onlinelibrary.wiley.com/doi/10.1111/biom.12295/full) Direct likelihood inference and sensitivity analysis for competing risks regression with missing causes of failure. *Biometrics* 2015; 71(2):498-507.

#### Cite this document as:

Moreno-Betancur M, Leacy FP, Tompsett D, White I. "mice: The NARFCS procedure for sensitivity analyses" (Available at: <https://github.com/moreno-betancur/NARFCS>)
