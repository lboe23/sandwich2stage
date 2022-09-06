
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sandwich2stage

<!-- badges: start -->
<!-- badges: end -->

`sandwich2stage` introduces a function for computing an estimate of the
sandwich variance for the two-stage regression model setting of
regression calibration. The sandwich is one approach for obtaining
standard errors in two-stage regression settings that account for the
extra uncertainty added by the calibration model step. Our function
computes an estimate of the sandwich variance obtained by stacking the
stage 1 and stage 2 estimating equation contributions.

## Installation

`sandwich2stage` may be installed from github as follows:

``` r
library(devtools)
#> Loading required package: usethis
install_github("lboe23/sandwich2stage", subdir="pkg")
#> Skipping install of 'sandwich2stage' from a github remote, the SHA1 (87f44922) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Example

After we have loaded our package, we will want to load in the sample
data assumed to be from a simple random sample, which is contained in
our package.

``` r
library(sandwich2stage)
data("sandwichdata_SRS")
```

Next, we need to load the survey package and create a survey design
object. For the case of the data `sandwichdata_SRS` from a simple random
sample, we will want to specify a simple random sampling design.

``` r
library(survey)
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
sampdesign <- svydesign(id=~1, data=sandwichdata_SRS)
#> Warning in svydesign.default(id = ~1, data = sandwichdata_SRS): No weights or
#> probabilities supplied, assuming equal probability
```

We may then fit the stage 1 and stage 2 models, saving the estimated
nuisance parameters from the stage 1 model and using them to obtain an
estimated of the unknown exposure (xhat). This estimate exposure, xhat,
will be used as a covariate in the stage 2 model. Note that the stage 1
model is only fit to the subset which contains validation data
(i.e. where v=1).

``` r
stage1.model<-survey::svyglm(xstarstar~xstar+z,design=sampdesign,family=gaussian(),subset=v==1)
alphas.stage1<-coef(stage1.model)
sampdesign <- update(sampdesign,xhat =predict(stage1.model,newdata=sampdesign$variables) )
```

We will then fit the stage 2 model with the estimated exposure (xhat) as
a covariate.

``` r
stage2.model<-  survey::svyglm(y ~ xhat+z,design=sampdesign,family=binomial())
```

Finally, we can obtain an estimate of the sandwich variance using our
function`sandwich2stage()`. The sandwich variance matrix is saved below
in the object `sandwichvar`.

``` r
sandwich.object<-sandwich2stage(stage1.model,stage2.model,xstar="xstar",xhat="xhat",Stage1ID="ID",Stage2ID="ID")
sandwichvar<-vcov(sandwich.object) 
```
