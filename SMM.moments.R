
#######################################################################################
# Date: 03/16/2016
# Author: Michelle Shardell, PhD
# This function provides the aymptotically efficient moments for G-estimation of 
# multiplicative structural mean models with potentially invalid IVs described in 
#
#  Shardell M, Ferrucci L. Instrumental variable analysis of multiplicative models with 
#      potentially invalid instruments. Stat Med. 2016 Dec 20;35(29):5430-5447. 
#      doi: 10.1002/sim.7069. Epub 2016 Aug 16. PMID: 27527517; PMCID: PMC5118169. 
######################################################################################

###moment functions (score equations) for G-estimation of structural mean models (SMMs)
SMM.moments <- function(psi,x){
             C <- x[,“C”] ###covariate
             Y <- x[,“Y”] ###outcome
             X <- x[,“X”] ###exposure
             Z <- x[,“Z”] ###(possibly invalid) instrument
             ###predicted treatment-free counterfactual
             H <- Y*exp(-psi[1]*Z-psi[2]*X)
             ###q = E[H | C]
             q <- lm(H∼C)$fitted.values
             dH.dpsi.1 <- -Y*Z*exp(-psi[1]*Z-psi[2]*X)
             ###D1 = E[dH.dpsi.1 | Z, C]
             D1 <- lm(dH.dpsi.1∼Z*C)$fitted.values
             dH.dpsi.2 <- -Y*X*exp(-psi[1]*Z-psi[2]*X)
             ###D2 = E[dH.dpsi.2 | Z, C]
             D2 <- lm(dH.dpsi.2∼Z*C)$fitted.values
             sq.dev <- (H-q)ˆ2
             ###w = 1/E[(H-q)ˆ2 | Z, C]
             w <- 1/(lm(sq.dev∼Z*C)$fitted.values)
             ###Ew = E[w | C]
             Ew <- lm(w∼C)$fitted.values
             wD1 <- w*D1
             ###EwD1 = E[wD1 | C]
             EwD1 <- lm(wD1∼C)$fitted.values
             wD2 <- w*D2
             ###EwD2 = E[wD2 | C]
             EwD2 <- lm(wD2∼C)$fitted.values
             ###Optimal weights for the estimating equations
             dopt.1 <- w*(D1 - EwD1/Ew)
             dopt.2 <- w*(D2 - EwD2/Ew)
             ###moments
             m1 <- dopt.1*(H - q)
             m2 <- dopt.2*(H - q)
             return(cbind(m1,m2))
}
