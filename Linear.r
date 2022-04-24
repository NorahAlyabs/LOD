#################################################################### 
#                                                                  #
#                  Simulation for linear model                     #
#                                                                  #
####################################################################

#required packages for the simulation.
library(mice)
#The mice package is used to impute the missing values by the predictive mean matching (PMM) method.

#'Define variables:
#'
#' param n is the sample size n = 50, 100, 500
#' @param mr missing proportion for X1 and X2: we used mr= 20%, 40%, and 60%
#' @param lod  LOD type where 0:lower, 1:upper, 2:interval
#' @param X1 covariate subjected to LOD, X1 ~ Weibul(1,1/3)
#' @param X2 covariate subjected to LOD, X2 ~ Normal(0,0.8^2)
#' @param Z fully observed covariate , Z ~ Normal(0,1)
#' @param V1 missing indicator for X1
#' @param V2 missing indicator for X2
#' @param Y is the outcome variable Y= -2 + X1 - X2 +  Z + e
#' @param e error term e ~ Normal(0,1)


# 1. Generating data:

gen.data.lod.linear <- function (n, mr, lod) {
  X1 <- rweibull(n, 1, 1/3)   
  X2 <- rnorm(n, 0, 0.8)
  Z  <- rnorm(n, 0, 1)
  if (!(mr %in% seq(0, 1, 0.01)))    stop("Undefinded censoring rate")
    if (lod == 0){ # Lower LOD
    L1 <- qweibull(mr, 1, 1/3) ; U1 <- qweibull(1, 1, 1/3)
    L2 <- qnorm(mr, 0, 0.8)   ; U2 <- qnorm(1, 0, 0.8)}
    if (lod == 1){ #Upper LOD
      L1 <- qweibull(0, 1, 1/3) ; U1 <- qweibull(1 - mr, 1, 1/3)
      L2 <- qnorm(0, 0, 0.8 )   ; U2 <- qnorm(1 - mr, 0, 0.8)}
    if (lod == 2){ #Interval LOD
      L1 <- qweibull(mr/4, 1, 1/3) ; U1 <- qweibull(1 - 3 * mr/4, 1, 1/3)
      L2 <- qnorm(mr/4, 0, 0.8)   ; U2 <- qnorm(1 - 3 * mr/4, 0, 0.8)}
    if (!(lod %in% c(0,1,2)))    stop("Undefinded LOD type")
    V1 <- 1 * (X1 >= L1 & X1 <= U1)  # 1 => observed, 0 => censored
    V2 <- 1 * (X2 >= L2 & X2 <= U2)   # 1 => observed, 0 => censored
    Y <- -2 + X1 - X2 +  Z + rnorm(n, 0, 1)
    return(data.frame( Y = Y,  V1, X1 = pmin( pmax(X1, L1), U1),
    V2 , X2 = pmin( pmax( X2, L2), U2), Z))
  }
# 2- Estimate the 3 coefficients and their associated variances for the 12 methods for one sample:

One.simulation.linear <- function(n, mr, lod ) { 
  # Generating data:
  data <- gen.data.lod.linear(n, mr, lod )
  # M1-  CC: Complete case
  if (with(data, sum(V1 * V2)) > 2)
  {CCA <-lm(Y ~ X1 + X2 + Z, data = data, subset = V1 * V2 > 0 )
  if (sum(is.na(vcov(CCA))) == 0) {CC1 <- c(coef(CCA)[2:4], diag(vcov(CCA))[2:4])
  }  else{ CC1 <- rep(NA, 6)}} else { CC1 <- rep(NA, 6) } 

  # M2- Subtitution with L/2 or 2U:
  dat1 <- data
  dat1$X1[dat1$V1 == 0 & data$X1 == min(data$X1)] <- min(data$X1)/2
  dat1$X2[dat1$V2 == 0 & data$X2 == min(data$X2)] <- min(data$X2)/2
  dat1$X1[dat1$V1 == 0 & data$X1 == max(data$X1)] <- max(data$X1)*2
  dat1$X2[dat1$V2 == 0 & data$X2 == max(data$X2)] <- max(data$X2)*2
  sub1 <- lm(Y ~ X1 + X2 + Z , data = dat1 )
  if (dim(vcov(sub1))[1] == 4)
  {Sub1 <- c(coef(sub1)[2:4], diag(vcov(sub1))[2:4])
  } else {Sub1 <- c(coef(sub1)[2:4], rep(NA, 3))}

   #M3- Subtitution with L/sqrt(2) or sqrt(2)U:
  dat2 <- data
  dat2$X1[dat2$V1 == 0 & data$X1 == min(data$X1)] <- min(data$X1)/sqrt(2)
  dat2$X2[dat2$V2 == 0 & data$X2 == min(data$X2)] <- min(data$X2)/sqrt(2)
  dat2$X1[dat2$V1 == 0 & data$X1 == max(data$X1)] <- max(data$X1)*sqrt(2)
  dat2$X2[dat2$V2 == 0 & data$X2 == max(data$X2)] <- max(data$X2)*sqrt(2)
  sub2 <- lm(Y ~ X1 + X2 + Z , data = dat2 )
  if (dim(vcov(sub2))[1] == 4)
  {Sub2 <- c(coef(sub2)[2:4], diag(vcov(sub2))[2:4])
  } else {Sub2 <- c(coef(sub2)[2:4], rep(NA, 3))}

  # M4- Single Imputation with conditional expected value E(X|X<L) & E(X|X>U)
  # function to calculate the conditional expectations:
  dat3 <- data                    # Data with NA instead of LODs
  dat3$X1[data$V1 == 0] <- NA
  dat3$X2[data$V2 == 0] <- NA
  dat4 <- dat3
  #The limits of detection:
  L1 <- min(data$X1) 
  L2 <- min(data$X2)
  U1 <- max(data$X1) 
  U2 <- max(data$X2)
  # MLE function for incomplete covariates X1 & X2:
    truncNor <- function(para, X, V, l, u) {
    Xo <- X[V == 1] ;    Xm <- X[V == 0]
    r1 <- sum(Xm == min(X)) ;    r2 <- sum(Xm == max(X))
    mu <- para[1] ;   s1 <- para[2]
    # Likelihood function
    - sum(dnorm(Xo, mu, s1, log = T)) - r1 * log(pnorm(l, mu, s1)) - r2 *
    log(1- pnorm(u, mu, s1))
  }
# Intitial MLE of mu & sigma for X1 & X2 subjected LOD
    para <- with(dat4, c(mean(X1, na.rm = T), sd(X1, na.rm = T),
                         mean(X2, na.rm = T), sd(X2, na.rm = T)))
    # The MLE of mu & sigma for X1 & X2 subjected LOD
    MLE1 <- optim(par = para[1:2], fn = truncNor, X = data$X1, V = data$V1, 
    l = L1, u = U1)[["par"]]
    MLE2 <- optim(par = para[3:4], fn = truncNor, X = data$X2, V = data$V2, 
    l = L2, u = U2)[["par"]]
    phi <- function(x, mle){dnorm((x- mle[1])/ mle[2], log = T)}
    Phi <- function(x, mle){pnorm((x- mle[1])/ mle[2])}
    # Impute missing values with the the conditional expectations
  dat4$X1[data$V1 == 0 & data$X1 == L1 ] <- MLE1[1] - MLE1[2] *
  exp(phi(L1, MLE1))/ Phi(L1, MLE1)
  dat4$X2[data$V2 == 0 & data$X2 == L2 ] <- MLE2[1] - MLE2[2] *
  exp(phi(L2, MLE2))/ Phi(L2, MLE2)
  dat4$X1[data$V1 == 0 & data$X1 == U1 ] <- MLE1[1] + MLE1[2] *
  exp(phi(U1, MLE1))/ (1 - Phi(U1, MLE1))
  dat4$X2[data$V2 == 0 & data$X2 == U2 ] <- MLE2[1] + MLE2[2] *
  exp(phi(U2, MLE2))/(1 - Phi(U2, MLE2))
  #Fitting linear model and remove divergence cases if any.
  Si <- lm(Y ~ X1 + X2 + Z , data = dat4 )
  if (dim(vcov(Si))[1] == 4)
  {SI <- c(coef(Si)[2:4], diag(vcov(Si))[2:4])
  } else {SI <- c(coef(Si)[2:4], rep(NA, 3))}

# M5- mice::PMM
  Data1 <- data[,1:6]
  Data1$X1[data$V1 == 0] <- Data1$X2[data$V2 == 0] <- NA
  predictormat <- 1 - diag(1, ncol(Data1))
  predictormat[,c(2,4)] <- 0
  imp <- mice(Data1, meth=c("","","pmm","","pmm",""), seed=2021,
  predictorMatrix=predictormat, remove_collinear = FALSE, printFlag = FALSE)
  PMM <- mice::pool(with(imp, lm(Y  ~ X1 + X2 + Z)))
  PMM1 <- rep(NA,6)
  if (anyNA(summary(PMM)$std.error)) {PMM1[4:6] <- NA}
  if (length(summary(PMM)$estimate) < 4 )  {PMM1[1:3] <- NA
  } else {PMM1 <- c(summary(PMM)[2:4,2], summary(PMM)[2:4, 3]^2)}

   # M6 - MI from imputing random variable from conditional density:
  # number of missing values due to lower and upper LOD for X1 & X2:
  r1l <- sum(data$V1 == 0 & data$X1 == L1)
  r2l <- sum(data$V2 == 0 & data$X2 == L2)
  r1u <- sum(data$V1 == 0 & data$X1 == U1)
  r2u <- sum(data$V2 == 0 & data$X2 == U2)

  # Imputing & estimating function:
  Sub4 <- function(i, dat = dat5, V){
    # 1- X subject to lower LOD:
    dat$X1[data$V1 == 0 & data$X1 == L1] <- MLE1[1] +  MLE1[2] *
    qnorm(runif(r1l, 0, 1) *  Phi(L1, MLE1))
    dat$X2[data$V2 == 0 & data$X2 == L2] <- MLE2[1] +  MLE2[2] *
    qnorm(runif(r2l, 0, 1) *  Phi(L2, MLE2))
    # 2- X subject to upper LOD:
    dat$X1[data$V1 == 0 & data$X1 == U1] <- MLE1[1] + MLE1[2] *
    qnorm(Phi(U1, MLE1) + runif(r1u, 0, 1) *  (1 - Phi(U1, MLE1)))
    dat$X2[data$V2 == 0 & data$X2 == U2] <- MLE2[1] + MLE2[2] *
    qnorm(Phi(U2, MLE2) + runif(r2u, 0, 1) *  (1 - Phi(U2, MLE2)))
    # linear model:
    if (V == 0){  #M6 without extending
  sub4 <- lm(Y ~ X1 + X2 + Z, data = dat )
    }
    if (V == 2){  # Extending with MID
      sub4 <- lm(Y ~ X1 + X2 + Z + I(1 - V1) + I(1 - V2) + I(Z*X1) +I(Z*X2), data = dat )
    }
    if (V == 3){ # Extending with EMID
      sub4 <- lm(Y ~ X1 + X2 + Z + V1 + V2 + I(X1 * (1- V1)) +I(X2 * (1- V2)) + I(X2 * Z)  + I(X1*Z),  data = dat )
    }
    if (dim(vcov(sub4))[1] == 4 + 2 * V  )
    {Sub4 <- c(coef(sub4)[2:4], diag(vcov(sub4))[2:4])
    } else {Sub4 <- c(coef(sub4)[2:4], rep(NA, 3))}
    Sub4
  } 
  dat5 <- dat3
  # Multiple imputation for 5 times:
  MultiSub4 <-  sapply(1:5, Sub4, V=0)
  # Pooled parameters using Rubin`s Rule:
 MIC <- rep(NA,6)
  MIC[1:3] <- rowMeans(MultiSub4[1:3, ], na.rm = T)
  MIC[4:6] <- sqrt(rowMeans(MultiSub4[4:6, ]^2, na.rm = T) + (1+1/5) * apply(MultiSub4[1:3, ], 1, sd))

  # M7- Missing indicator approach (MDI) :
  MDI <- lm(Y  ~ I(X1 * V1) + I(X2 * V2) + Z + I(1 - V1) + I(1 - V2), data = data )
  if (dim(vcov(MDI))[1] == 6)
  {MDI1 <- c(coef(MDI)[2:4], diag(vcov(MDI))[2:4])
  } else {MDI1 <- c(coef(MDI)[2:4], rep(NA, 3))}

  # M8- Extended missing indicator approach (EMDI):
  EMDI <- lm(Y ~ I(X1 * V1) + I(X2 * V2) + Z + I(1 - V1) + I(1 - V2) +
  I(Z * (1 - V1)) + I(Z * (1 - V2)), data = data)
  if (dim(vcov(EMDI))[1] == 8)
  {EMDI1 <- c(coef(EMDI)[2:4], diag(vcov(EMDI))[2:4])
  } else {EMDI1 <- c(coef(EMDI)[2:4], rep(NA, 3))}

  # M9- combine MDI(m7) and MI(M5):
  MIM <- mice::pool(with(imp, lm(Y ~ X1 + X2 + Z + I(1 - V1) + I(1 - V2) + 
  I(Z*X1) +I(Z*X2) )))
  MIM1 <- rep(NA,6)
  if (!anyNA(summary(MIM)$std.error)) {MIM1[4:6] <- summary(MIM)[2:4, 3]^2}
  if (length(summary(MIM)$estimate) == 8 )  {MIM1[1:3] <- summary(MIM)[2:4,2]}

  # M10 -   combine MDI(M7) & MIC(M6)
  # Multiple imputation for 5 times:
  MultiSub5 <-  sapply(1:5, Sub4, V=2)
  # Pooled parameters using Rubin`s Rule:
  MIMC <- rep(NA,6)
  MIMC[1:3] <- rowMeans(MultiSub5[1:3, ], na.rm = T)
  MIMC[4:6] <- sqrt(rowMeans(MultiSub5[4:6, ]^2, na.rm = T) + (1+1/5) *
  apply(MultiSub5[1:3, ], 1, sd))

  # M11- combine EMDI(M8) and MI(M5):
  EMIM <- mice::pool(with(imp, lm(Y ~ X1 + X2 + Z + V1 + V2 + I(X1 * (1- V1)) +
  I(X2 * (1- V2)) + I(X2 * Z)  + I(X1 * Z) )))
  EMIM1 <- rep(NA,6)
  if (!anyNA(summary(EMIM)$std.error)) {EMIM1[4:6] <- summary(EMIM)[2:4, 3]^2}
  if (length(summary(EMIM)$estimate) == 10 )  {EMIM1[1:3] <- summary(EMIM)[2:4,2]}

  # M12 - combine EMDI(M8) & MIC(M6):
  # Multiple imputation for 5 times:
  MultiSubE <-  sapply(1:5, Sub4, V=3)
  # Pooled parameters using Rubin`s Rule:
  EMIMC <- rep(NA,6)
  EMIMC[1:3] <- rowMeans(MultiSubE[1:3, ], na.rm = T)
  EMIMC[4:6] <- sqrt(rowMeans(MultiSubE[4:6, ]^2, na.rm = T) + (1+1/5) *
  apply(MultiSubE[1:3, ], 1, sd))

 result <- c(CC1, Sub1, Sub2, SI, PMM1, MIC, MDI1, EMDI1, MIM1, MIMC, EMIM1, EMIMC )
  Result <- matrix(result, ncol = 6, byrow = T)
  colnames(Result) <- c("b1", "b2", "b3", "V(b1)", "V(b2)", "V(b3)")
  row.names(Result) <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12")
  Result
}

set.seed(0123)
One.simulation.linear(100, 0.2, 0)
