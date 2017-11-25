library(joineR)
library(lattice)
library(ggplot2)
library(methods)
library(lme4)

data<-read.csv("/home/ben/Dropbox/CDSS/healthy.csv",sep=",",header=TRUE)
which(is.na(data)==T)

healthy.surv <- UniqueVariables(data,
                              var.col = c("ST","Surv"),
                              id.col = "Volunteer.ID")
healthy.long <- data[, c("Volunteer.ID", "hours","cc.K18..U.l.","GLDH..U.l.","HMGB1..ng.ml.","miR.122.copy.number.uL","CT.miR.122.let.7d","CT.miR.122.U6", "FL.K18..U.l.", "ALT")]
healthy.cov <- UniqueVariables(data,
                             c("age"),
                             id.col = "Volunteer.ID")
healthy.jd <- jointdata(longitudinal = healthy.long,
                            baseline = healthy.cov,
                            survival = healthy.surv,
                            id.col = "Volunteer.ID",
                            time.col = "hours")
summary(healthy.jd)
ALT.fit <- joint(data = healthy.jd,
             long.formula = ALT ~ 1 + hours,
             surv.formula = Surv(ST,Surv) ~ 1,
             model = "intslope")
summary(ALT.fit)

FL.fit <- joint(data = healthy.jd,
                 long.formula = FL.K18..U.l. ~ 1 + hours,
                 surv.formula = Surv(ST,Surv) ~ 1,
                 model = "intslope")
summary(FL.fit)

cc.fit <- joint(data = healthy.jd,
                long.formula = cc.K18..U.l. ~ 1 + hours,
                surv.formula = Surv(ST,Surv) ~ 1,
                model = "intslope")
summary(cc.fit)

GLDH.fit <- joint(data = healthy.jd,
                long.formula = GLDH..U.l. ~ 1 + hours,
                surv.formula = Surv(ST,Surv) ~ 1,
                model = "intslope")
summary(GLDH.fit)

HMGB1.fit <- joint(data = healthy.jd,
                  long.formula = HMGB1..ng.ml. ~ 1 + hours,
                  surv.formula = Surv(ST,Surv) ~ 1,
                  model = "intslope")
summary(HMGB1.fit)

miRcopy.fit <- joint(data = healthy.jd,
                  long.formula = miR.122.copy.number.uL ~ 1 + hours,
                  surv.formula = Surv(ST,Surv) ~ 1,
                  model = "intslope")
summary(miRcopy.fit)

miR7d.fit <- joint(data = healthy.jd,
                     long.formula = CT.miR.122.let.7d ~ 1 + hours,
                     surv.formula = Surv(ST,Surv) ~ 1,
                     model = "intslope")
summary(miR7d.fit)

miRU6.fit <- joint(data = healthy.jd,
                   long.formula = CT.miR.122.U6 ~ 1 + hours,
                   surv.formula = Surv(ST,Surv) ~ 1,
                   model = "intslope")
summary(miRU6.fit)

#No survival outcome needed so
#use LME model instead.

fmALT<-lmer(ALT ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmALT.nullr<-lmer(ALT ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmALT.nullf<-lmer(ALT ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmALT,fmALT.nullr) # test if random is sig
anova(fmALT,fmALT.nullf) # test if fixed is sig

fmFL<-lmer(FL.K18..U.l. ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmFL.nullr<-lmer(FL.K18..U.l. ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmFL.nullf<-lmer(FL.K18..U.l. ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmFL,fmFL.nullr)
anova(fmFL,fmFL.nullf)

fmcc<-lmer(cc.K18..U.l. ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmcc.nullr<-lmer(cc.K18..U.l. ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmcc.nullf<-lmer(cc.K18..U.l. ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmcc,fmcc.nullr)
anova(fmcc,fmcc.nullf)

fmGLDH<-lmer(GLDH..U.l. ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmGLDH.nullr<-lmer(GLDH..U.l. ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmGLDH.nullf<-lmer(GLDH..U.l. ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmGLDH,fmGLDH.nullr)
anova(fmGLDH,fmGLDH.nullf)

fmHMG<-lmer(HMGB1..ng.ml. ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmHMG.nullr<-lmer(HMGB1..ng.ml. ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmHMG.nullf<-lmer(HMGB1..ng.ml. ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmHMG,fmHMG.nullr)
anova(fmHMG,fmHMG.nullf)

fmcopy<-lmer(miR.122.copy.number.uL ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmcopy.nullr<-lmer(miR.122.copy.number.uL ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmcopy.nullf<-lmer(miR.122.copy.number.uL ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmcopy,fmcopy.nullr)
anova(fmcopy,fmcopy.nullf)

fm7d<-lmer(log(CT.miR.122.let.7d) ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fm7d.nullr<-lmer(log(CT.miR.122.let.7d) ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fm7d.nullf<-lmer(log(CT.miR.122.let.7d) ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fm7d,fm7d.nullr)
anova(fm7d,fm7d.nullf)

fmU6<-lmer(log(CT.miR.122.U6) ~ hours + (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
fmU6.nullr<-lmer(log(CT.miR.122.U6) ~ hours + (1|Volunteer.ID),data=data,REML=FALSE)
fmU6.nullf<-lmer(log(CT.miR.122.U6) ~ (1|hours) + (1|Volunteer.ID),data=data,REML=FALSE)
anova(fmU6,fmU6.nullr)
anova(fmU6,fmU6.nullf)





## Plotting

plot(ALT.jd,Y.col="ALT..U.l.")
plot(ALT.jd,Y.col="CT.miR.122.U6")
plot(ALT.jd,Y.col="CT.miR.122.let.7d")
plot(ALT.jd,Y.col="miR.122.copy.number.uL")
plot(ALT.jd,Y.col="HMGB1..ng.ml.")
plot(ALT.jd,Y.col="cc.K18..U.l.")
plot(ALT.jd,Y.col="FL.K18..U.l.")

pat<-data[data$Volunteer.ID==1,]

## Per marker
res.lme <- lme(ALT ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
ALTgg <-ggplot(data = data, aes(x = hours, y = ALT, group = Volunteer.ID))
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/ALT.png", width=800, height=1200,res=120)
ALTgg + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 4)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(log(CT.miR.122.U6) ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/miRU6.png", width=800, height=1200,res=120)
U6 <- ggplot(data = data, aes(x = hours, y = log(CT.miR.122.U6), group = Volunteer.ID))
U6 + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 4)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(log(CT.miR.122.let.7d) ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
data$resid <- with(data, fitted - CT.miR.122.let.7d)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/miR7d.png", width=800, height=1200,res=120)
miR7d <- ggplot(data = data, aes(x = hours, y = log(CT.miR.122.let.7d), group = Volunteer.ID))
miR7d + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(log(miR.122.copy.number.uL) ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/miRcopy.png", width=800, height=1200,res=120)
miR <- ggplot(data = data, aes(x = hours, y = log(miR.122.copy.number.uL), group = Volunteer.ID))
miR + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(HMGB1..ng.ml. ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/HMG.png", width=800, height=1200,res=120)
HMG <- ggplot(data = data, aes(x = hours, y = HMGB1..ng.ml., group = Volunteer.ID))
HMG + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 4)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(GLDH..U.l. ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/GLDH.png", width=800, height=1200,res=120)
GLDH <- ggplot(data = data, aes(x = hours, y = data$GLDH..U.l., group = Volunteer.ID))
GLDH + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 3)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(cc.K18..U.l. ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/ccK18.png", width=800, height=1200,res=120)
cc <- ggplot(data = data, aes(x = hours, y = cc.K18..U.l., group = Volunteer.ID))
cc + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 4)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

res.lme <- lme(FL.K18..U.l. ~ 1, data = data, random = ~ 1|Volunteer.ID)
data$fitted<-predict(res.lme)
png("/home/ben/Dropbox/CDSS/PlotsPerMarker/FLK18.png", width=800, height=1200,res=120)
FL <- ggplot(data = data, aes(x = hours, y = FL.K18..U.l., group = Volunteer.ID))
FL + geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 4)) + geom_line(aes(colour=Volunteer.ID)) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

#### Per patient

ALT <-ggplot(data = pat, aes(x = T, y = ALT..U.l., group = Volunteer.ID))

ALT + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)

U6 <- ggplot(data = pat, aes(x = T, y = log(CT.miR.122.U6), group = Volunteer.ID))
U6 + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)


png("/home/ben/Dropbox/CDSS/PlotsPerMarker/miRcopy.png", width=1200, height=800,res=120)
miR <- ggplot(data = data, aes(x = T, y = log(miR.122.copy.number.uL), group = Volunteer.ID))
miR + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerMarker/HMG.png", width=1200, height=800,res=120)
HMG <- ggplot(data = data, aes(x = T, y = HMGB1..ng.ml., group = Volunteer.ID))
HMG + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerMarker/GLDH.png", width=1200, height=800,res=120)
GLDH <- ggplot(data = data, aes(x = T, y = data$GLDH..U.l., group = Volunteer.ID))
GLDH + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerMarker/ccK18.png", width=1200, height=800,res=120)
cc <- ggplot(data = data, aes(x = T, y = cc.K18..U.l., group = Volunteer.ID))
cc + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerMarker/FLK18.png", width=1200, height=800,res=120)
FL <- ggplot(data = data, aes(x = T, y = FL.K18..U.l., group = Volunteer.ID))
FL + geom_line(aes(colour=Volunteer.ID)) + stat_smooth(aes(group = 1), method = "lm", se = TRUE) +  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3)
dev.off()

time<-factor(data$T)
vol<-factor(data$Volunteer.ID)

## Repeated Measures ANOVA

ALT.out=aov(data$ALT..U.l. ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID) )
summary(ALT.out)
U6.out=aov(log(data$CT.miR.122.U6) ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID) )
summary(U6.out)
mIR.7D.out=aov(data$CT.miR.122.let.7d ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID) )
summary(mIR.7D.out)
copy.out=aov(log(data$miR.122.copy.number.uL) ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID))
summary(copy.out)
HMGB1.out=aov(data$HMGB1..ng.ml. ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID))
summary(HMGB1.out)
GLDH.out=aov(data$GLDH..U.l. ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID))
summary(GLDH.out)
cc.out=aov(data$cc.K18..U.l. ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID))
summary(cc.out)
FL.out=aov(data$FL.K18..U.l. ~ data$T * data$Volunteer.ID + Error(data$Volunteer.ID))
summary(FL.out)

ALT1.out=aov(data$ALT..U.l. ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID) )
summary(ALT1.out)
U61.out=aov(data$CT.miR.122.U6 ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID) )
summary(U61.out)
miR.7D1.out=aov(data$CT.miR.122.let.7d ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID) )
summary(miR.7D1.out)
copy1.out=aov(data$miR.122.copy.number.uL ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID))
summary(copy1.out)
HMGB11.out=aov(data$HMGB1..ng.ml. ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID))
summary(HMGB11.out)
GLDH1.out=aov(data$GLDH..U.l. ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID))
summary(GLDH1.out)
cc1.out=aov(data$cc.K18..U.l. ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID))
summary(cc1.out)
FL1.out=aov(data$FL.K18..U.l. ~ data$Volunteer.ID + Error(data$T/data$Volunteer.ID))
summary(FL1.out)

## Seperate Linear Models
png("/home/ben/Dropbox/CDSS/PlotsPerID/ALT.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,data$ALT..U.l., xlab="Volunteer", ylab=colnames(data[4]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(data$ALT..U.l. ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab=colnames(data[4]),legend=FALSE)
dev.off()
ints <- by(data, data$Volunteer.ID, 
           function(data) coefficients(lm(data$ALT..U.l. ~ data$T, data=data))[[1]])  
ints1 <- unlist(ints)
names(ints1) <- NULL
mean(ints1)
sqrt(var(ints1))

slopes <- by(data, data$Volunteer.ID,
             function(data) coefficients(lm(data$Volunteer.ID ~ data$T, data=data))[[2]])  
slopes1 <- unlist(slopes)
names(slopes1) <- NULL
mean(slopes1)
sqrt(var(slopes1))

cor( ints1, slopes1)

png("/home/ben/Dropbox/CDSS/PlotsPerID/U6.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,log(data$CT.miR.122.U6), xlab="Volunteer", ylab=colnames(data[5]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(log(data$CT.miR.122.U6) ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab="log(CT.miR.122.U6)",legend=FALSE)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerID/7D.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,log(data$CT.miR.122.let.7d), xlab="Volunteer", ylab=colnames(data[6]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(log(data$CT.miR.122.let.7d) ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab="log(CT.miR.122.let.7d)",legend=FALSE)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerID/copy.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,data$miR.122.copy.number.uL, xlab="Volunteer", ylab=colnames(data[7]),legend=FALSE,ylim=c(0,150000))
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(data$miR.122.copy.number.uL ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab=colnames(data[7]),ylim=c(0,150000),legend=FALSE)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerID/HMGB1.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,data$HMGB1..ng.ml., xlab="Volunteer", ylab=colnames(data[8]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(data$HMGB1..ng.ml. ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab=colnames(data[8]),legend=FALSE)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerID/GLDH.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,data$GLDH..U.l., xlab="Volunteer", ylab=colnames(data[9]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(data$GLDH..U.l. ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab=colnames(data[9]),legend=FALSE)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerID/ccK18.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,data$cc.K18..U.l., xlab="Volunteer", ylab=colnames(data[10]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(data$cc.K18..U.l. ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab=colnames(data[10]),legend=FALSE)
dev.off()

png("/home/ben/Dropbox/CDSS/PlotsPerID/FLK18.png", width=1200, height=800,res=120)
par(mfrow=c(1,2))
interaction.plot(data$T,data$Volunteer.ID,data$FL.K18..U.l., xlab="Volunteer", ylab=colnames(data[11]),legend=FALSE)
fit <- by(data, data$Volunteer.ID, function(data) fitted.values(lm(data$FL.K18..U.l. ~ data$T, data=data))) 
fit <- unlist(fit)
interaction.plot(data$T, data$Volunteer.ID, fit, xlab="Volunteer", ylab=colnames(data[11]),legend=FALSE)
dev.off()
