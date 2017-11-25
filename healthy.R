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
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
data<-read.csv("/home/ben/Dropbox/CDSS/Dan/Healthy/healthy090817study3.csv",sep=",",header=TRUE)

#ALT
qqp(data$ALT,"norm")
qqp(data$ALT,"lnorm")
lmm <- lmer(log(ALT) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm(log(ALT) ~ Sex + Age + Ethnicity , data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

#CSF-1
qqp(data$CSF.1,"norm")
qqp(data$CSF.1,"lnorm")
qqp(sqrt(data$CSF.1),"norm")
lmm <- lmer((CSF.1) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm((CSF.1) ~ Sex + Age + Ethnicity, data = data)

lmm <- lmer((CSF.1) ~ Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm((CSF.1) ~ Ethnicity, data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

##GLDH
qqp(data$GLDH,"norm")
qqp(data$GLDH,"lnorm")#
qqp(sqrt(data$GLDH),"norm")
lmm <- lmer(log(GLDH) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm(log(GLDH) ~ Sex + Age + Ethnicity , data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)#*

##FL.K18
qqp(data$FL.K18,"norm")
qqp(data$FL.K18,"lnorm")#
qqp(sqrt(data$FL.K18),"norm")
lmm <- lmer(sqrt(FL.K18) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm(sqrt(FL.K18) ~ Sex + Age + Ethnicity , data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)#age
exactLRT(m=lmm,m0=lmm0)

##ccK18
qqp(data$ccK18,"norm")
qqp(data$ccK18,"lnorm")
qqp(sqrt(data$ccK18),"norm")
lmm <- lmer(sqrt(ccK18) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm(sqrt(ccK18) ~ Sex + Age + Ethnicity , data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)#age
exactLRT(m=lmm,m0=lmm0)

##HMGB1
qqp(data$HMGB1,"norm")
qqp(data$HMGB1,"lnorm")
qqp(sqrt(data$HMGB1),"norm")
lmm <- lmer(sqrt(HMGB1) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm(sqrt(HMGB1) ~ Sex + Age + Ethnicity , data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)#age
exactLRT(m=lmm,m0=lmm0)

##miR.122
qqp(data$miR.122,"norm")
qqp(data$miR.122,"lnorm")
qqp(sqrt(data$miR.122),"norm")
lmm <- lmer(log(HMGB1) ~ Sex + Age + Ethnicity + (1 | TP), data = data,
            REML = FALSE)
lmm0 <- lm(log(HMGB1) ~ Sex + Age + Ethnicity , data = data)
summary(lmm)
summary(lmm0)
Anova(lmm)#age
exactLRT(m=lmm,m0=lmm0)



T=data$Time
Boxplot(data$ALT..U.l.,T,labels=data$Volunteer.ID,main="ALT")
Boxplot(data$CT.miR.122.U6,T,labels=data$Volunteer.ID,main="miR-122 U6")
Boxplot(data$CT.miR.122.let.7d,T,labels=data$Volunteer.ID,main="miR-122 7d")
Boxplot(data$miR.122.copy.number.uL,T,labels=data$Volunteer.ID,main="miR-122 copy")
Boxplot(data$HMGB1..ng.ml.,T,labels=data$Volunteer.ID,main="HMGB1")
Boxplot(data$GLDH..U.l.,T,labels=data$Volunteer.ID,main="GLDH")
Boxplot(data$cc.K18..U.l.,T,labels=data$Volunteer.ID,main="cc-K18")
Boxplot(data$FL.K18..U.l.,T,labels=data$Volunteer.ID,main="FL-K18")

boxplot(data$CT.miR.122.U6~T,outline=F,main="miR-122 U6")
boxplot(data$CT.miR.122.let.7d~T,outline=F,main="miR-122 7d")
boxplot(data$miR.122.copy.number.uL~T,outline=F,main="miR-122 copy")
boxplot(data$HMGB1..ng.ml.~T,outline=F,main="HMGB1")
boxplot(data$GLDH..U.l.~T,outline=F,main="GLDH")
boxplot(data$cc.K18..U.l.~T,outline=F,main="cc-K18")
boxplot(data$FL.K18..U.l.~T,outline=F,main="FL-K18")

T.levels=unique(T)

split1<-data[data$Time==T.levels[1],]
split2<-data[data$Time==T.levels[2],]
split3<-data[data$Time==T.levels[3],]
split4<-data[data$Time==T.levels[4],]
split5<-data[data$Time==T.levels[5],]
split6<-data[data$Time==T.levels[6],]
split7<-data[data$Time==T.levels[7],]
split8<-data[data$Time==T.levels[8],]
split9<-data[data$Time==T.levels[9],]

##ALT
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  hist(split$ALT..U.l.,main=T.levels[i])
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$CT.miR.122.U6),main=T.levels[i])
  hist(log(split$CT.miR.122.U6),main=T.levels[i])
  hist(sqrt(split$CT.miR.122.U6),main=T.levels[i])
}
ALT.split=NULL
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(split1$ALT..U.l.)))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(log(split2$ALT..U.l.))))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(split3$ALT..U.l.)))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(split4$ALT..U.l.)))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(split5$ALT..U.l.)))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(log(split6$ALT..U.l.))))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(split7$ALT..U.l.)))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(log(split8$ALT..U.l.))))
ALT.split<-cbind(ALT.split,remove_outliers(remove_outliers(sqrt(split9$ALT..U.l.))))
                            
ALTref=NULL
for (i in 1:9) {
  a<-quantile(ALT.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(ALT.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(ALT.split[,i], probs=c(.5),na.rm=TRUE)
  ALTref<-cbind(ALTref,c(a,b,c))
  Boxplot(ALT.split[,i],labels=rownames(split),main=T.levels[i])
}

ALTref[,1]<-ALTref[,1]
ALTref[,2]<-exp(ALTref[,2])
ALTref[,3]<-ALTref[,3]
ALTref[,4]<-ALTref[,4]
ALTref[,5]<-ALTref[,5]
ALTref[,6]<-exp(ALTref[,6])
ALTref[,7]<-ALTref[,7]
ALTref[,8]<-exp(ALTref[,8])
ALTref[,9]<-(ALTref[,9])^2


##miR-122.U6
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$CT.miR.122.U6),main=T.levels[i])
  hist(log(split$CT.miR.122.U6),main=T.levels[i])
  hist(sqrt(split$CT.miR.122.U6),main=T.levels[i])
}

U6.split=NULL
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(log(split1$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(log(split2$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(sqrt(split3$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(sqrt(split4$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(log(split5$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(log(split6$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(log(split7$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(log(split8$CT.miR.122.U6))))
U6.split<-cbind(U6.split,remove_outliers(remove_outliers(sqrt(split9$CT.miR.122.U6))))

U6ref=NULL
for (i in 1:9) {
  a<-quantile(U6.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(U6.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(U6.split[,i], probs=c(.5),na.rm=TRUE)
  U6ref<-cbind(U6ref,c(a,b,c))
  Boxplot(U6.split[,i],labels=rownames(split),main=T.levels[i])
}

U6ref[,1]<-exp(U6ref[,1])
U6ref[,2]<-exp(U6ref[,2])
U6ref[,3]<-(U6ref[,3])^2
U6ref[,4]<-(U6ref[,4])^2
U6ref[,5]<-exp(U6ref[,5])
U6ref[,6]<-exp(U6ref[,6])
U6ref[,7]<-exp(U6ref[,7])
U6ref[,8]<-exp(U6ref[,8])
U6ref[,9]<-(U6ref[,9])^2

##miR-122.7d
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$CT.miR.122.let.7d),main=T.levels[i])
  hist(log(split$CT.miR.122.let.7d),main=T.levels[i])
  hist(sqrt(split$CT.miR.122.let.7d)^,main=T.levels[i])
}

miR7d.split=NULL
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split1$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split2$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split3$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split4$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split5$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split6$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(sqrt(split7$CT.miR.122.let.7d))))
miR7d.split<-cbind(miR7d.split,split8$CT.miR.122.let.7d)
miR7d.split<-cbind(miR7d.split,remove_outliers(remove_outliers(log(split9$CT.miR.122.let.7d))))



miR7dref=NULL
for (i in 1:9) {
  a<-quantile(miR7d.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(miR7d.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(miR7d.split[,i], probs=c(.5),na.rm=TRUE)
  miR7dref<-cbind(miR7dref,c(a,b,c))
  Boxplot(miR7d.split[,i],labels=rownames(split),main=T.levels[i])
}

miR7dref[,1]<-exp(miR7dref[,1])
miR7dref[,2]<-exp(miR7dref[,2])
miR7dref[,3]<-exp(miR7dref[,3])
miR7dref[,4]<-exp(miR7dref[,4])
miR7dref[,5]<-exp(miR7dref[,5])
miR7dref[,6]<-exp(miR7dref[,6])
miR7dref[,7]<-(miR7dref[,7])^2
miR7dref[,8]<-miR7dref[,8]
miR7dref[,9]<-exp(miR7dref[,9])

##miR-122.copy
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$miR.122.copy.number.uL),main=T.levels[i])
  hist(log(split$miR.122.copy.number.uL),main=T.levels[i])
  hist(sqrt(split$miR.122.copy.number.uL),main=T.levels[i])
}

copy.split=NULL
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(sqrt(split1$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(log(split2$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(sqrt(split3$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(sqrt(split4$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(log(split5$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(sqrt(split6$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(log(split7$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(sqrt(split8$miR.122.copy.number.uL))))
copy.split<-cbind(copy.split,remove_outliers(remove_outliers(sqrt(split9$miR.122.copy.number.uL))))

copyref=NULL
for (i in 1:9) {
  a<-quantile(copy.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(copy.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(copy.split[,i], probs=c(.5),na.rm=TRUE)
  copyref<-cbind(copyref,c(a,b,c))
  Boxplot(copy.split[,i],labels=rownames(split),main=T.levels[i])
}

copyref[,1]<-(copyref[,1])^2
copyref[,2]<-exp(copyref[,2])
copyref[,3]<-(copyref[,3])^2
copyref[,4]<-(copyref[,4])^2
copyref[,5]<-exp(copyref[,5])
copyref[,6]<-(copyref[,6])^2
copyref[,7]<-exp(copyref[,7])
copyref[,8]<-(copyref[,8])^2
copyref[,9]<-(copyref[,9])^2


##HMGB1
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$HMGB1..ng.ml.),main=T.levels[i])
  hist(log(split$HMGB1..ng.ml.),main=T.levels[i])
  hist(sqrt(split$HMGB1..ng.ml.),main=T.levels[i])
}


HMG.split=NULL
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(split1$HMGB1..ng.ml.)))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(split2$HMGB1..ng.ml.)))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(log(split3$HMGB1..ng.ml.))))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(sqrt(split4$HMGB1..ng.ml.))))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(split5$HMGB1..ng.ml.)))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(split6$HMGB1..ng.ml.)))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(split7$HMGB1..ng.ml.)))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(sqrt(split8$HMGB1..ng.ml.))))
HMG.split<-cbind(HMG.split,remove_outliers(remove_outliers(split9$HMGB1..ng.ml.)))

HMGref=NULL
for (i in 1:9) {
  a<-quantile(HMG.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(HMG.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(HMG.split[,i], probs=c(.5),na.rm=TRUE)
  HMGref<-cbind(HMGref,c(a,b,c))
  Boxplot(HMG.split[,i],labels=rownames(split),main=T.levels[i])
}

HMGref[,1]<-HMGref[,1]
HMGref[,2]<-HMGref[,2]
HMGref[,3]<-exp(HMGref[,3])
HMGref[,4]<-(HMGref[,4])^2
HMGref[,5]<-HMGref[,5]
HMGref[,6]<-HMGref[,6]
HMGref[,7]<-HMGref[,7]
HMGref[,8]<-(HMGref[,8])^2
HMGref[,9]<-HMGref[,9]

##GLDH
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$GLDH..U.l.),main=T.levels[i])
  hist(log(split$GLDH..U.l.),main=T.levels[i])
  hist(sqrt(split$GLDH..U.l.),main=T.levels[i])
}
GLDH.split=NULL
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(sqrt(split1$GLDH..U.l.))))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split2$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split3$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split4$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split5$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split6$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split7$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split8$GLDH..U.l.)))
GLDH.split<-cbind(GLDH.split,remove_outliers(remove_outliers(split9$GLDH..U.l.)))

GLDHref=NULL
for (i in 1:9) {
  a<-quantile(GLDH.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(GLDH.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(GLDH.split[,i], probs=c(.5),na.rm=TRUE)
  GLDHref<-cbind(GLDHref,c(a,b,c))
  Boxplot(GLDH.split[,i],labels=rownames(split),main=T.levels[i])
}

GLDHref[,1]<-(GLDHref[,1])^2
GLDHref[,2]<-GLDHref[,2]
GLDHref[,3]<-GLDHref[,3]
GLDHref[,4]<-GLDHref[,4]
GLDHref[,5]<-GLDHref[,5]
GLDHref[,6]<-GLDHref[,6]
GLDHref[,7]<-GLDHref[,7]
GLDHref[,8]<-GLDHref[,8]
GLDHref[,9]<-GLDHref[,9]

##cc.K18
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$cc.K18..U.l.),main=T.levels[i])
  hist(log(split$cc.K18..U.l.),main=T.levels[i])
  hist(sqrt(split$cc.K18..U.l.),main=T.levels[i])
}
cc.split=NULL
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(split1$cc.K18..U.l.)))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(log(split2$cc.K18..U.l.))))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(log(split3$cc.K18..U.l.))))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(split4$cc.K18..U.l.)))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(split5$cc.K18..U.l.)))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(split6$cc.K18..U.l.)))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(split7$cc.K18..U.l.)))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(sqrt(split8$cc.K18..U.l.))))
cc.split<-cbind(cc.split,remove_outliers(remove_outliers(split9$cc.K18..U.l.)))

ccref=NULL
for (i in 1:9) {
  a<-quantile(cc.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(cc.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(cc.split[,i], probs=c(.5),na.rm=TRUE)
  ccref<-cbind(ccref,c(a,b,c))
  Boxplot(cc.split[,i],labels=rownames(split),main=T.levels[i])
}

ccref[,1]<-ccref[,1]
ccref[,2]<-exp(ccref[,2])
ccref[,3]<-exp(ccref[,3])
ccref[,4]<-ccref[,4]
ccref[,5]<-ccref[,5]
ccref[,6]<-ccref[,6]
ccref[,7]<-ccref[,7]
ccref[,8]<-ccref[,8]
ccref[,9]<-(ccref[,9])^2

##FL.K18
for (i in 1:9) {
  split<-data[data$Time==T.levels[i],]
  par(mfrow=c(3,1))
  hist((split$FL.K18..U.l.),main=T.levels[i])
  hist(log(split$FL.K18..U.l.),main=T.levels[i])
  hist(sqrt(split$FL.K18..U.l.),main=T.levels[i])
}
FL.split=NULL
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(log(split1$FL.K18..U.l.))))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(log(split2$FL.K18..U.l.))))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(split3$FL.K18..U.l.)))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(sqrt(split4$FL.K18..U.l.))))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(split5$FL.K18..U.l.)))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(sqrt(split6$FL.K18..U.l.))))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(sqrt(split7$FL.K18..U.l.))))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(split8$FL.K18..U.l.)))
FL.split<-cbind(FL.split,remove_outliers(remove_outliers(split9$FL.K18..U.l.)))

FLref=NULL
for (i in 1:9) {
  a<-quantile(FL.split[,i], probs=c(.025),na.rm=TRUE)
  b<-quantile(FL.split[,i], probs=c(.975),na.rm=TRUE)
  c<-quantile(FL.split[,i], probs=c(.5),na.rm=TRUE)
  FLref<-cbind(FLref,c(a,b,c))
  Boxplot(FL.split[,i],labels=rownames(split),main=T.levels[i])
}

FLref[,1]<-exp(FLref[,1])
FLref[,2]<-exp(FLref[,2])
FLref[,3]<-FLref[,3]
FLref[,4]<-sqrt(FLref[,4])
FLref[,5]<-FLref[,5]
FLref[,6]<-sqrt(FLref[,6])
FLref[,7]<-sqrt(FLref[,7])
FLref[,8]<-FLref[,8]
FLref[,9]<-FLref[,9]

table<-rbind(ALTref,U6ref,miR7dref,copyref,HMGref,GLDHref,ccref,FLref)

write.table(table,file="/home/ben/Dropbox/CDSS/healthy_results.csv",quote=FALSE,sep=",")


