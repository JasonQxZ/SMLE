
Usage
=====

A demo code for SMLE-screening
------------------------------

First we show how to use `SMLE` to conduct feature screening and
post-screening selection via a simulated example. We generate a dataset
with *n* = 400 observations and *p* = 1000 features. We generate the
feature matrix *X* from a multivariate normal distribution with an
auto-regressive structure, where the adjacent features have a high
correlation of *ρ* = 0.9. The response variable *Y* is generated based
on the following logistic model with success rate *π* and linear
predictor:
logit(π) = 2*x*<sub>1</sub> + 3*x*<sub>3</sub> − 3*x*<sub>5</sub> + 3*x*<sub>7</sub> − 4*x*<sub>9</sub>.

    library(SMLE)

    set.seed(1)

    Data_eg <- Gen_Data(n = 400, p = 1000, family = "binomial",
    correlation = "AR", rho = 0.9, pos_truecoef = c(1,3,5,7,9),
    effect_truecoef = c(2,3,-3,3,-4))

    Data_eg

    ## Call:
    ##  Gen_Data(n = 400, p = 1000, pos_truecoef = c(1, 3, 5, 7, 9), 
    ##     effect_truecoef = c(2, 3, -3, 3, -4), correlation = "AR", 
    ##     rho = 0.9, family = "binomial")
    ##  
    ## An object of class sdata
    ##  
    ## Simulated Dataset Properties:
    ##  Length of response: 400
    ##  Dim of features: 400 x 1000
    ##  Correlation: auto regressive
    ##  Rho: 0.9
    ##  Index of Causal Features: 1,3,5,7,9
    ##  Model Type: binomial

In this setup, the feature matrix contains only five features that are
causally-related to the response, as indicated in the model. Some
features have marginal effects to response due to the correlation
structure. From the true model, we know that *x*<sub>2</sub> is not
causally-related to the response. Yet, we can see that the marginal
effect of *x*<sub>2</sub> appears to be pretty high; thus, this
irrelevant feature is likely to be retained in the model if the
screening is done based on marginal effects only.

    coef(summary(glm(Data_eg$Y ~ Data_eg$X[,2], family = "binomial")))

    ##                  Estimate Std. Error   z value     Pr(>|z|)
    ## (Intercept)    0.02440072  0.1125979 0.2167067 8.284369e-01
    ## Data_eg$X[, 2] 1.10465766  0.1337063 8.2618250 1.434619e-16

The following code shows the simplest function call to `SMLE()`, where
we aim to retain only *k* = 10 important features out of *p* = 1000.

    fit1 <- SMLE(Data_eg$Y, Data_eg$X, k = 10, family = "binomial")

    summary(fit1)

    ## Call:
    ##   SMLE(formula = Data_eg$Y, X = Data_eg$X, k = 10, family = "binomial")
    ##  
    ## An object of class summary.smle
    ##  
    ## Summary:
    ## 
    ##   Length of response: 400
    ##   Dim of features: 400 x 1000
    ##   Model type: binomial
    ##   Retained model size: 10
    ##   Retained features: 1,3,5,7,9,68,430,536,661,709
    ##   Coefficients estimated by IHT:  1.572  3.127 -1.893  2.534 -4.499 -0.348 -0.316  0.579  0.717  0.444
    ##   Number of IHT iteration steps: 59

The function returns a `'smle'` object and `summary()` function confirms
that a refined set of 10 features is selected after 59 IHT iterations.
We can see that all 5 causal features used to generate the response are
retained in the refined set. This indicates that screening is
successful; the dimensionality of the feature space is reduced from
*p* = 1000 down to *k* = 10 without losing any important information. In
this example, `SMLE()` accurately removes
*x*<sub>2</sub>, *x*<sub>4</sub>, *x*<sub>6</sub>, *x*<sub>8</sub>, as
its screening naturally incorporates the joint effects among features.

Further selection after screeening
----------------------------------

Note that the refined set returned in the model still contains some
irrelevant features; this is to be expected (the `k` always chosen to be
larger than the actual number of casual features), as the goal of
feature screening is merely to remove most irrelevant features before
conducting an in-depth analysis. One may conduct an elaborate selection
on the refined set to further identify the causal features.

As can be seen below, `smle_select()` returns a `'selection'` object
`"fit1_s"`, which exactly identifies the five features in the true data
generating model.

    fit1_s <- smle_select(fit1, criterion = "ebic")

    summary(fit1_s)

    ## Call:
    ##   Select(object = fit1, criterion = "ebic")
    ##  
    ## An object of class summary.selection
    ##  
    ## Summary:
    ##  
    ##   Length of response: 400
    ##   Dim of features: 400 x 1000
    ##   Model type: binomial
    ##   Selected model size: 5
    ##   Selected features: 1,3,5,7,9
    ##   Selection criterion: ebic
    ##   Gamma for ebic: 0.5

An example for categorical features.
------------------------------------

Categorical features fed in the package will be convert to `'factor'`
and dummy coded during the iterations. In this example, we generate a
dataset with causal categorical features and separate it into training
and testing groups in order to perform a prediction task.

    set.seed(1)
    Data_sim2 <- Gen_Data(n = 420, p = 1000, family = "gaussian", num_ctgidx = 5, 
                          pos_ctgidx = c(1,3,5,7,9), effect_truecoef= c(1,2,3,-4,-5),
                          pos_truecoef = c(1,3,5,7,8), level_ctgidx = c(3,3,3,4,5))

    train_X <- Data_sim2$X[1:400,]; test_X <- Data_sim2$X[401:420,]
    train_Y <- Data_sim2$Y[1:400]; test_Y <- Data_sim2$Y[401:420]

    test_X[1:5,1:10]

    ##     C1          N2 C3         N4 C5         N6 C7         N8 C9         N10
    ## 401  C -0.17405549  B  0.3337833  B -1.8054836  B  0.9696772  C -0.88066391
    ## 402  C  0.96129056  C -0.2113226  A -0.6780407  B -2.1994065  C -0.48558301
    ## 403  B  0.29382666  A -0.5510979  B -0.4733581  C  1.9480938  B  0.22743281
    ## 404  B  0.08099936  B  0.2583611  A  1.0274171  B  0.1798532  B -0.06646135
    ## 405  B  0.18366184  A -1.3752104  B -0.5973876  C  0.4150568  B  0.35161359

Users may specify whether to treat those dummy covariates as a single
group feature or as individual features, and which type of dummy coding
is used by arguments: `gourp` and `codyingtype`. Note that the number of
features retained in the model may less than the `k` specified when
`group` is `FALSE` since one categorical feature may be chose several
times by its covariates. More details see the package manual.

    fit_1 <- SMLE(train_Y, train_X, family = "gaussian", group = TRUE, codingtype = "standard", k = 10)
    fit_1

    ## Call:
    ##   SMLE(formula = train_Y, X = train_X, k = 10, family = "gaussian", 
    ##     group = TRUE, codingtype = "standard")
    ##  
    ## An object of class smle
    ##  
    ## Subset:
    ##   Model size: 10
    ##   Feature Name: C1,C3,C5,C7,N8,N297,N327,N671,N727,N992
    ##   Feature Index: 1,3,5,7,8,297,327,671,727,992

    predict(fit_1, newdata = test_X)

    ##         401         402         403         404         405         406 
    ##  -0.9015259  13.0801313 -14.7328909  -2.6187980  -6.8932628  -9.4658403 
    ##         407         408         409         410         411         412 
    ##   9.3106969   3.5823446   4.8032411  -3.3804086  -5.4039198   5.4538378 
    ##         413         414         415         416         417         418 
    ##  -6.2625712   8.2882945  -5.1134303   0.5188218   8.0068460  -3.4739318 
    ##         419         420 
    ##  -9.5672038   6.3228009

    fit_2 <- SMLE(train_Y, train_X, family = "gaussian", group = FALSE, codingtype = "all", k = 10)
    fit_2

    ## Call:
    ##   SMLE(formula = train_Y, X = train_X, k = 10, family = "gaussian", 
    ##     group = FALSE, codingtype = "all")
    ##  
    ## An object of class smle
    ##  
    ## Subset:
    ##   Model size: 5
    ##   Feature Name: C1,C3,C5,C7,N8
    ##   Feature Index: 1,3,5,7,8

    predict(fit_2, newdata = test_X)

    ##         401         402         403         404         405         406 
    ##  -0.9243726  13.1374968 -15.0696386  -3.3249395  -7.3601888 -10.1484939 
    ##         407         408         409         410         411         412 
    ##   8.7424408   3.6110151   4.5296146  -3.6479939  -5.6094304   7.0635413 
    ##         413         414         415         416         417         418 
    ##  -7.6369911   8.4131558  -3.9716724   0.5468016   8.0184504  -3.8099618 
    ##         419         420 
    ##  -8.5262404   4.3285765

Formula interface
-----------------

`SMLE` always works in low dimensional as a selection method. Although
it is not recommend, interface to `'formula'` object provides user a
better understanding to the package in high dimension.

    library(datasets)
    data("attitude")
    SMLE(rating+complaints ~ complaints + complaints:privileges + learning + raises*advance, data = attitude)

    ## Call:
    ##   SMLE(X = x, Y = y, categorical = FALSE)
    ##  
    ## An object of class smle
    ##  
    ## Subset:
    ##   Model size: 5
    ##   Feature Name: complaints,learning,raises,advance,complaints:privileges
    ##   Feature Index: 1,2,3,4,5

Xu, Chen, and Jiahua Chen. 2014. “The Sparse Mle for
Ultrahigh-Dimensional Feature Screening.” *Journal of the American
Statistical Association* 109 (507): 1257–69.
