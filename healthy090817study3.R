install.packages("outliers")
install.packages("car")
install.packages("xts")
install.packages("lme4")
install.packages("RLRsim")

library(lme4)
library(car)
library(MASS)
library(xts)
library(RLRsim)
data<-read.csv("/home/ben/Dropbox/CDSS/Dan/Healthy/healthy090817study3.csv",sep=",",header=TRUE)
data<-read.csv("C:\\Users/brfra/Dropbox/CDSS/Dan/Healthy/healthy090817study3.csv",sep=",",header=TRUE)

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
# Outliers Removed - 6
###############################################################

#Test Data - fixed and random effects
Anova(lmer(tbio ~ Sex + (1 | TP), data = data, REML = FALSE))#

###############################################################
# Univariate Fixed Effects - Sex
###############################################################

Anova(lmer(tbio ~ Sex + (1 | TP), data = data, REML = FALSE))#

###############################################################
# Multivariate Fixed Effects - Sex
###############################################################

lmm<-lmer(tbio ~ Sex + (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ Sex, data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? No 
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
#                 lower upper
# 2.5%   9.00000   Inf    10
# 50%   17.49643    15    23
# 97.5% 41.84418    39   Inf
#
###############################################################

##########################################################################################################################
# miR.122let 
##########################################################################################################################

bio=data$miR.122
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
# Outliers Removed - 3
###############################################################

#Test Data - fixed and random effects

###############################################################
# Univariate Fixed Effects -
###############################################################

###############################################################
# Multivariate Fixed Effects - 
###############################################################

lmm<-lmer(tbio ~ (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? No 
###############################################################

#Reference Ranges
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
#                 lower upper
# 2.5%  0.2543371  0.00  0.33
# 50%   1.2359612  0.53  1.82
# 97.5% 3.5041277  2.39   Inf
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

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects - 
###############################################################

lmm<-lmer(tbio ~ (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? No
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
#                  lower upper
# 2.5%  0.8912447   Inf  1.18
# 50%   1.8449695  1.43  2.22
# 97.5% 2.6942477  2.46   Inf
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

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects -
###############################################################

lmm<-lmer(tbio ~ (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? No 
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
#                 lower upper
# 2.5%  117.2528    Inf 155.14
# 50%   222.2402 188.11 260.10
# 97.5% 336.3171 308.84    Inf
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

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects -
###############################################################

lmm<-lmer(tbio ~ (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? No
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
#                 lower upper
#2.5%  226.1811    Inf 267.3
#50%   364.6701 291.14 433.8
#97.5% 557.8520 539.00   Inf
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
# Outliers Removed - 0
###############################################################

#Test Data - fixed and random effects

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects -
###############################################################

lmm<-lmer(tbio ~ (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? Yes
###############################################################

bio=data$tbio
if(tr=="l"){lapply(seq(1,9),function(x){exp(cbind(c(quantile(bio[data$TP==x&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                      rbind(quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                            quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                            quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.975, conf.level = 0.9))))})}
if(tr=="s"){lapply(seq(1,9),function(x){(cbind(c(quantile(bio[data$TP==x&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.13, conf.level = 0.9),
                         quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.87, conf.level = 0.9))))^2})}
if(tr=="n"){lapply(seq(1,9),function(x){(cbind(c(quantile(bio[data$TP==1&!is.na(bio)], c(0.025,0.5,0.975), type = 7)),
                   rbind(quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.025, conf.level = 0.9),
                         quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.5, conf.level = 0.9),
                         quantileCI(bio[data$TP==x&!is.na(bio)], probs = 0.975, conf.level = 0.9))))})}
png("C:\\Users/brfra/Dropbox/CDSS/Dan/Healthy/GLDH_Diurnal_Plot.png")
boxplot((tbio)^2~TP,data=data[data$TP>=1&data$TP<9,],ylab="GLDH (ng/mg Cr)", xlab="Measurement Time Point",xaxt='n')
axis(1, at=1:8, labels=c("7am","10am","1pm","4pm","7pm","10pm","Midnight","4am"))
dev.off()
###############################################################
# Reference Range
# 7am
#               lower    upper
# 2.5%   51.57  45.56  68.06
# 50%    87.89  68.06 126.56
# 97.5% 648.34 182.25 992.25
# 
# 10am
#               lower    upper
# 2.5%   56.25  56.25  95.06
# 50%   126.56  95.06 203.06
# 97.5% 491.18 225.00 612.56
# 
# 1pm
#               lower     upper
# 2.5%   56.25  56.25   81.00
# 50%   110.25  81.00  126.56
# 97.5% 931.77 272.25 1040.06
# 
# 4pm
#               lower     upper
# 2.5%   45.56  45.56  144.00
# 50%   272.25 196.00  410.06
# 97.5% 985.18 410.06 1040.06
# 
# 7pm
#               lower  upper
# 2.5%   51.57  45.56  81.00
# 50%   110.25  81.00 144.00
# 97.5% 740.18 162.56 812.25
# 
# 10pm
#                 lower  upper
# 2.5%   25.00    25.00  36.00
# 50%    62.02    36.00 144.00
# 97.5% 763.83   144.00 812.25
# 
# 12am
#                 lower   upper
# 2.5%   56.25   56.25     81.00
# 50%   115.56   81.00    182.25
# 97.5% 834.85  203.06   1190.25 
# 
# 4am
#               lower     upper
# 2.5%   53.29  25.00  162.56
# 50%   260.02 169.00  400.00
# 97.5% 911.29 410.06 1190.25
# 
# 8am
#               lower     upper
# 2.5%   47.09  36.00  100.00
# 50%   126.56 110.25  248.06
# 97.5% 866.57 272.25 1139.06
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

###############################################################
# Univariate Fixed Effects - 
###############################################################

###############################################################
# Multivariate Fixed Effects - 
###############################################################

lmm<-lmer(tbio ~ (1 | TP), data = data, REML = FALSE)
summary(lmm)
lmm0 <- lm(tbio ~ 1, data = data)
#Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###############################################################
# Time Point Random Effect? No
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
#                  lower upper
# 2.5%  0.3402661   Inf  0.39
# 50%   0.9749936  0.83  1.39
# 97.5% 1.9780156  1.67   Inf
#
###############################################################
