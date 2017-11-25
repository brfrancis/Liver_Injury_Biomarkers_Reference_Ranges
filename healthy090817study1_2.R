data<-read.csv("C:\\Users/brfra/Dropbox/CDSS/Dan/Healthy/healthy090817study1_2.csv",sep=",",header=TRUE)
data<-read.csv("/home/ben/Dropbox/CDSS/Dan/Healthy/healthy090817study1_2.csv",sep=",",header=TRUE)
data$Eth=relevel(data$Ethnicity,ref="White")
resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "\n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "\n")
  cat("Mean of the outliers:", round(mo, 2), "\n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "\n")
  cat("Mean if we remove outliers:", round(m2, 2), "\n")
  par(resetPar())
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "\n")
    return(invisible(dt))
}

library(lme4)
library(car)
library(MASS)
library(xts)
library(RLRsim)
library(outliers)
library(quantreg)
library(jmuOutlier)
citation()
citation("lme4")
citation("car")
citation("MASS")
citation("xts")
citation("RLRsim")
citation("outliers")
citation("quantreg")
# library(lqmm)

##########################################################################################################################
# ALT 
##########################################################################################################################

bio=data$ALT
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - sqrt
###############################################################

#Trans Data
tr="s"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 0
###############################################################

#Test Data - fixed and random effects
res1=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)


###############################################################
# Univariate Fixed Effects - Sex, Height, Weight, BMI
###############################################################

summary(lmer(tbio ~ Sex + Weight + (1 | ID), data = data, REML = FALSE))
Anova(lmer(tbio ~ Sex + Height + (1 | ID), data = data, REML = FALSE))
Anova(lmer(tbio ~ Sex + BMI + (1 | ID), data = data, REML = FALSE))#

###############################################################
# Multivariate Fixed Effects - Sex, BMI
###############################################################

lmm<-lmer(tbio ~ Sex + BMI + (1 | ID), data = data, REML = FALSE)
summary(lmm);slmm=summary(lmm)
lmm0 <- lm(tbio ~ Sex + BMI, data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? Yes 
###############################################################

bio=data$tbio
if(tr=="l"){exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      cbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}
if(tr=="s"){round((cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2)}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

#Quantile Regression 
srq<-summary(rq(tbio~Sex+BMI,data=data[data$TP==1,],tau = 0.5,ci=FALSE),se="ker")$coefficients
srqCI<-cbind(srq[,1],srq[,1]-1.65*srq[,2],srq[,1]+1.65*srq[,2])
#Male
round((rbind(srqCI[1,]+srqCI[2,]+17*srqCI[3,],
srqCI[1,]+srqCI[2,]+21.7*srqCI[3,],
srqCI[1,]+srqCI[2,]+27.5*srqCI[3,],
srqCI[1,]+srqCI[2,]+35*srqCI[3,]))^2)
#Female
round((rbind(srqCI[1,]+17*srqCI[3,],
             srqCI[1,]+21.7*srqCI[3,],
             srqCI[1,]+27.5*srqCI[3,],
             srqCI[1,]+35*srqCI[3,]))^2)

srq<-summary(rq(tbio~Sex+BMI,data=data[data$TP==1,],tau = 0.025,ci=FALSE),se="ker")$coefficients
srqCI<-cbind(srq[,1],srq[,1]-1.65*srq[,2],srq[,1]+1.65*srq[,2])
#Male
round((rbind(srqCI[1,]+srqCI[2,]+17*srqCI[3,],
             srqCI[1,]+srqCI[2,]+21.7*srqCI[3,],
             srqCI[1,]+srqCI[2,]+27.5*srqCI[3,],
             srqCI[1,]+srqCI[2,]+35*srqCI[3,]))^2)
#Female
round((rbind(srqCI[1,]+17*srqCI[3,],
             srqCI[1,]+21.7*srqCI[3,],
             srqCI[1,]+27.5*srqCI[3,],
             srqCI[1,]+35*srqCI[3,]))^2)

srq<-summary(rq(tbio~Sex+BMI,data=data[data$TP==1,],tau = 0.975,ci=FALSE),se="ker")$coefficients
srqCI<-cbind(srq[,1],srq[,1]-1.65*srq[,2],srq[,1]+1.65*srq[,2])
#Male
round((rbind(srqCI[1,]+srqCI[2,]+17*srqCI[3,],
             srqCI[1,]+srqCI[2,]+21.7*srqCI[3,],
             srqCI[1,]+srqCI[2,]+27.5*srqCI[3,],
             srqCI[1,]+srqCI[2,]+35*srqCI[3,]))^2)
#Female
round((rbind(srqCI[1,]+17*srqCI[3,],
             srqCI[1,]+21.7*srqCI[3,],
             srqCI[1,]+27.5*srqCI[3,],
             srqCI[1,]+35*srqCI[3,]))^2)







###############################################################
# Reference Range
#          lower upper
# 2.5%  11     9    12
# 50%   23    22    24
# 97.5% 44    41    50
#
# Reference Range - Female Underweight
#           lower upper
# Only 1 subject
#
# Reference Range - Female Healthy
#           lower upper
# 2.5%  11   Inf    11
# 50%   19    17    23
# 97.5% 33    31   Inf
#
# Reference Range - Female Overweight
#           lower upper
# 2.5%  13   Inf    15
# 50%   22    20    25
# 97.5% 38    32   Inf
#
# Reference Range - Female Obese
#           lower upper
# 2.5%  16   Inf    17
# 50%   25    23    32
# 97.5% 42    38   Inf
#
# Reference Range - Male Under
#           lower upper
# Only 2 patients
#
# Reference Range - Male Healthy
#           lower upper
# 2.5%  13   Inf    15
# 50%   24    20    28
# 97.5% 46    44   Inf
#
# Reference Range - Male Overweight
#           lower upper
# 2.5%  15   Inf    17
# 50%   30    22    37
# 97.5% 48    42   Inf
#
# Reference Range - Male Obese 
#           lower upper
# 2.5%  19   Inf    25
# 50%   31    25    40
# 97.5% 41    40   Inf
#
#
###############################################################

##########################################################################################################################
# miR.122let 
##########################################################################################################################

bio=data$miR.122let
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - log
###############################################################

#Trans Data
tr="l"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 1
###############################################################

#Test Data - fixed and random effects
res2=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)

###############################################################
# Univariate Fixed Effects - Weight, BMI
###############################################################

Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))#

###############################################################
# Multivariate Fixed Effects - BMI
###############################################################

lmm<-lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ BMI, data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? Yes 
###############################################################

#Reference Ranges
bio=data$tbio
if(tr=="l"){round(exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9)))))}
if(tr=="s"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

#Quantile Regression 
srq<-summary(rq(tbio~BMI,data=data[data$TP==1,],tau = 0.025,ci=FALSE),se="ker")$coefficients
srqCI<-cbind(srq[,1],srq[,1]-1.65*srq[,2],srq[,1]+1.65*srq[,2])
round(exp(rbind(srqCI[1,]+17*srqCI[2,],
                srqCI[1,]+21.7*srqCI[2,],
                srqCI[1,]+27.5*srqCI[2,],
                srqCI[1,]+35*srqCI[2,])),digits = 2)

srq<-summary(rq(tbio~BMI,data=data[data$TP==1,],tau = 0.5,ci=FALSE),se="ker")$coefficients
srqCI<-cbind(srq[,1],srq[,1]-1.65*srq[,2],srq[,1]+1.65*srq[,2])
round(exp(rbind(srqCI[1,]+17*srqCI[2,],
             srqCI[1,]+21.7*srqCI[2,],
             srqCI[1,]+27.5*srqCI[2,],
             srqCI[1,]+35*srqCI[2,])),digits = 2)

#srq<-summary(rq(tbio~Sex+BMI,data=data[data$TP==1,],tau = 0.975,ci=FALSE),se="ker")$coefficients
#srqCI<-cbind(srq[,1],srq[,1]-1.65*srq[,2],srq[,1]+1.65*srq[,2])
#round(exp(rbind(srqCI[1,]+17*srqCI[2,], srqCI[1,]+21.7*srqCI[2,],srqCI[1,]+27.5*srqCI[2,],srqCI[1,]+35*srqCI[2,])),digits = 2)
rbind(c(6.1,4,Inf),c(6.4,4.32,Inf),c(6.77,5.92,Inf),c(8.04,5.58,Inf))

###############################################################
# Reference Range
#           lower upper
# 2.5%  0     0     0
# 50%   1     1     1
# 97.5% 8     6    11
#
#
# Reference Range - Underweight
#           lower upper
# Only 3 subjects
#
# Reference Range - Healthy
#           lower upper
# 2.5%  0.17  0.00  0.22
# 50%   0.95  0.71  1.21
# 97.5% 6.40  4.32   Inf
#
# Reference Range - Overweight
#             lower upper
# 2.5%  0.27  0.00  0.32
# 50%   1.47  1.04  1.85
# 97.5% 6.17  5.92   Inf
#
# Reference Range - Obese
#           lower upper
# 2.5%  0.21  0.00  0.28
# 50%   1.50  0.75  3.18
# 97.5% 8.04  5.58   Inf
#
#
###############################################################

##########################################################################################################################
# miR.122copies
##########################################################################################################################

bio=data$miR.122copies
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())
###############################################################
# Transform Chosen - log
###############################################################

#Trans Data
tr="l"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)
###############################################################
# Outliers Removed - 0
###############################################################

#Test Data - fixed and random effects
Anova(lm(tbio ~ Sex, data = data))
Anova(lm(tbio ~ Age, data = data))
Anova(lm(tbio ~ Height, data = data))
Anova(lm(tbio ~ Weight, data = data))
Anova(lm(tbio ~ BMI, data = data))
Anova(lm(tbio ~ Eth, data = data))

res3=round(rbind(c(summary(lm(tbio ~ Sex, data = data))$coefficients[2,1:2],
                   Anova(lm(tbio ~ Sex, data = data))$`F value`[1],
                   Anova(lm(tbio ~ Sex, data = data))$`Pr(>F)`[1]),
                 c(summary(lm(tbio ~ Age, data = data))$coefficients[2,1:2],
                   Anova(lm(tbio ~ Age, data = data))$`F value`[1],
                   Anova(lm(tbio ~ Age, data = data))$`Pr(>F)`[1]),
                 c(summary(lm(tbio ~ Height, data = data))$coefficients[2,1:2],
                   Anova(lm(tbio ~ Height, data = data))$`F value`[1],
                   Anova(lm(tbio ~ Height, data = data))$`Pr(>F)`[1]),
                 c(summary(lm(tbio ~ Weight, data = data))$coefficients[2,1:2],
                   Anova(lm(tbio ~ Weight, data = data))$`F value`[1],
                   Anova(lm(tbio ~ Weight, data = data))$`Pr(>F)`[1]),
                 c(summary(lm(tbio ~ BMI, data = data))$coefficients[2,1:2],
                   Anova(lm(tbio ~ BMI, data = data))$`F value`[1],
                   Anova(lm(tbio ~ BMI, data = data))$`Pr(>F)`[1]),
                 cbind(summary(lm(tbio ~ Eth, data = data))$coefficients[2:4,1:2],
                       c(Anova(lm(tbio ~ Eth, data = data))$`F value`[1],NA,NA),
                       c(Anova(lm(tbio ~ Eth, data = data))$`Pr(>F)`[1],NA,NA))),digits=6)


###############################################################
# Univariate Fixed Effects - 
###############################################################


###############################################################
# Multivariate Fixed Effects - 
###############################################################

# lmm<-lmer(tbio ~ Sex + BMI + (1 | ID), data = data, REML = FALSE)
# summary(lmm)
# lmm0 <- lm(tbio ~ Sex + BMI, data = data)
# Anova(lmm)
# exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? NA
###############################################################

#Reference Ranges
bio=data$tbio
if(tr=="l"){round(exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9)))))}
if(tr=="s"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

###############################################################
# Reference Range
#                  lower  upper
# 2.5%   198.0410  122.7  215.5
# 50%    668.5235  587.2  736.9
# 97.5% 3547.9841 2912.1 4321.2
#
###############################################################

##########################################################################################################################
# HMGB1
##########################################################################################################################

bio=data$HMGB1
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - norm
###############################################################

#Trans Data
tr="n"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 2
###############################################################

#Test Data - fixed and random effects
res4=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects - 
###############################################################

lmm<-lmer(tbio ~ (1 | ID), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? Yes 
###############################################################

bio=data$tbio
if(tr=="l"){exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}
if(tr=="s"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2}
if(tr=="n"){round(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))),digits=2)}

###############################################################
# Reference Range
#             lower upper
# 2.5%  0.22  0.17  0.32
# 50%   1.24  1.16  1.29
# 97.5% 2.34  2.23  2.42
#
###############################################################

##########################################################################################################################
# ccK18
##########################################################################################################################

bio=data$ccK18
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - sqrt
###############################################################

#Trans Data
tr="s"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 0
###############################################################

#Test Data - fixed and random effects
res5=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)
###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects - 
###############################################################

lmm<-lmer(tbio ~ (1 | ID), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? Yes 
###############################################################

bio=data$tbio
if(tr=="l"){exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}
if(tr=="s"){round((cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2)}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

###############################################################
# Reference Range
#                 lower upper
# 2.5%   57.27741  53.0  59.9
# 50%   132.39992 122.3 141.9
# 97.5% 271.89194 256.2 290.9
#
###############################################################

##########################################################################################################################
# FL.K18
##########################################################################################################################

bio=data$FL.K18
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - sqrt
###############################################################

#Trans Data
tr="s"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 0
###############################################################

#Test Data - fixed and random effects
res6=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)
###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects -
###############################################################

lmm<-lmer(tbio ~ (1 | ID), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? No
###############################################################

bio=data$tbio
if(tr=="l"){exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}
if(tr=="s"){round((cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2)}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

###############################################################
# Reference Range
#                 lower upper
# 2.5%  114.4749 102.4 126.5
# 50%   248.0494 224.7 267.7
# 97.5% 474.5125 455.7 488.4
#
###############################################################

##########################################################################################################################
# GLDH
##########################################################################################################################

bio=data$GLDH
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - sqrt
###############################################################

#Trans Data
tr="s"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 12
###############################################################

#Test Data - fixed and random effects
res7=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects -
###############################################################

lmm<-lmer(tbio ~ (1 | ID), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? No
###############################################################

bio=data$tbio
if(tr=="l"){exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}
if(tr=="s"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

###############################################################
# Reference Range
#           lower upper
# 2.5%   5     5     6
# 50%   13    13    14
# 97.5% 27    26    30
#
###############################################################

##########################################################################################################################
# CSF.1
##########################################################################################################################

bio=data$CSF.1
par(mfrow=c(3, 1))
qqp(bio,"norm")
qqp(bio,"lnorm")
qqp(sqrt(bio),"norm")
par(resetPar())

###############################################################
# Transform Chosen - sqrt
###############################################################

#Trans Data
tr="s"
if(tr=="l"){data$tbio=log(bio)};if(tr=="s"){data$tbio=sqrt(bio)};if(tr=="n"){data$tbio=bio}
outlierKD(data, tbio)

###############################################################
# Outliers Removed - 0
###############################################################

#Test Data - fixed and random effects
res8=round(rbind(c(summary(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Sex + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Age + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Height + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ Weight + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 c(summary(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$coefficients[2,1:2],
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$Chisq,
                   Anova(lmer(tbio ~ BMI + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`),
                 cbind(summary(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$coefficients[2:4,1:2],
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$Chisq,NA,NA),
                       c(Anova(lmer(tbio ~ Eth + (1 | ID), data = data, REML = FALSE))$`Pr(>Chisq)`,NA,NA))),digits=6)
###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects - 
###############################################################

lmm<-lmer(tbio ~ (1 | ID), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Subject Random Effect? No
###############################################################

bio=data$tbio
if(tr=="l"){exp(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}
if(tr=="s"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))^2}
if(tr=="n"){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==1&!is.na(bio)], probs = 0.975, conf.level = 0.9))))}

###############################################################
# Reference Range
#             lower upper
# 2.5%  0.46  0.30  0.56
# 50%   1.40  1.30  1.46
# 97.5% 2.42  2.33  2.88
#
###############################################################

res=rbind(res1,res2,res3,res4,res5,res6,res7,res8)
write.csv(res,file="/home/ben/Dropbox/CDSS/Dan/Healthy/healthy090817study1_2res.csv",quote=FALSE,sep=",")
