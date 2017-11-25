data2<-read.csv("/home/ben/Dropbox/CDSS/Dan/Healthy/healthy090817study2l.csv",sep=",",header=TRUE)


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

##ALT
qqp(data2$ALT,"norm")
qqp(data2$ALT,"lnorm")
qqp(sqrt(data2$ALT),"norm")
data2$sqALT=sqrt(data2$ALT)

outlierKD(data2, sqALT)#0

lmm <- lmer(log(sqALT) ~ Sex + BMI + (1 | ID), data = data2,
            REML = FALSE)
lmm0 <- lm(log(ALT) ~ Sex + BMI , data = data2)
summary(lmm)
summary(lmm0)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

hist(aggregate(ALT~ID,data=data2,FUN=mean)[,2]/aggregate(ALT~ID,data=data2,FUN=sd)[,2])
