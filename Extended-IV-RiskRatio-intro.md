Tutorial on G-estimation of multiplicative structural mean models with potentially invalid instrumental variables
=================================================================================================================

Setting up
----------

First, source the code for the structural mean model (SMM) moments used
in this tutorial. The function SMM.moments() includes the asymptotically
efficient moments described in Equations 4 and 6 (and derived in
Appendix B) of Shardell and Ferrucci (2016).

    source("SMM.moments.R")

Make sure that the gmm package is installed and call it.

    if("gmm" %in% rownames(installed.packages()) == FALSE) {install.packages("gmm")}
    library(gmm)

Example data
------------

The following produces the data in Table 3 of Shardell and Ferrucci
(2016) and converts it from wide format to long format:

    ### data from Table 3 in Shardell and Ferrucci (2016)
    ### Y = outcome, X = exposure, Z = invalid IV, C = covariate
    var.combos <- list(Y = 0:1, X = 0:1, Z = 0:1, C = 0:1)
    wide.data <- expand.grid(var.combos)
    wide.data$combo.freq <- c(180, 129, 211, 235, 220, 85, 504, 374, 226, 62, 94, 26, 282, 70, 334, 115)

    ### converting data from wide to long
    the.data = wide.data[rep(1:nrow(wide.data), wide.data$combo.freq), 
                                 -grep("combo.freq", names(wide.data))]

Implement gmm() using SMM.moments() for asymptotically efficient estimation
---------------------------------------------------------------------------

G-estimation of multiplicative structural mean models is carried out as
a special case of generalized method of moments (GMM). GMM can be
sensitive to starting values; for this reason, generate estimates using
gmm() with starting coefficient value of 0, then update the values by
implementing gmm() again. Results are the coefficients and standard
errors for the exposure (X) and the potentially invalid instrument (Z).
Results are used to derive causal risk ratios and 95% confidence
intervals reported in the “EIV, Optimal” section of Table 2 in Shardell
and Ferrucci (2016).

    ### using starting coefficiant values of 0 and 
    ### assuming iid for estimating variance covariance 
    SMM.pre <- gmm(SMM.moments, x=the.data, t0=c(0,0), vcov="iid")

    ###update results
    SMM <- gmm(SMM.moments, x=the.data, t0=c(coef(SMM.pre)), vcov="iid")

    ###psi.z = psi[1]; psi.x = psi[2]
    psi <- coef(SMM)
    print(psi)

    ##   Theta[1]   Theta[2] 
    ## -0.4784239  1.5304685

    ###SE(psi.z) = SE[1]; SE(psi.x) = SE[2]
    SE <- sqrt(diag(SMM$vcov))
    print(SE)

    ##   Theta[1]   Theta[2] 
    ## 0.03406614 0.17917564

    ##print Causal Risk Ratio and 95% CI
    lcl <- exp(psi - 1.96*SE)
    ucl <- exp(psi + 1.96*SE)

    print(cbind(exp(psi), lcl, ucl))

    ##                          lcl       ucl
    ## Theta[1] 0.6197594 0.5797296 0.6625533
    ## Theta[2] 4.6203409 3.2520457 6.5643451
