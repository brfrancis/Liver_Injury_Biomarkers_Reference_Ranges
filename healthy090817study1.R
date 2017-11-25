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
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  par(resetPar())
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}

library(lme4)
library(car)
library(MASS)
library(xts)
library(RLRsim)

##ALT 
qqp(data$ALT,"norm")
qqp(data$ALT,"lnorm")
qqp(sqrt(data$ALT),"norm")#
data$sqALT=sqrt(data$ALT)

outlierKD(data, sqALT)#0

lmm <- lmer(sqALT ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(sqALT ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)



##miR-122let
qqp(data$miR.122let,"norm")
qqp(data$miR.122let,"lnorm")#
qqp(sqrt(data$miR.122let),"norm")
data$lmiR.122let=log(data$miR.122let)

outlierKD(data, lmiR.122let)#1

lmm <- lmer(lmiR.122let ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(lmiR.122let ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

##miR-122copies
qqp(data$miR.122copies,"norm")
qqp(data$miR.122copies,"lnorm")#
qqp(sqrt(data$miR.122copies),"norm")

data$lmiR.122copies=log(data$miR.122copies)
outlierKD(data, lmiR.122copies)#0

summary(lm(log(miR.122copies)~Sex,data))
summary(lm(log(miR.122copies)~Age,data))
summary(lm(log(miR.122copies)~Height,data))
summary(lm(log(miR.122copies)~Weight,data))
summary(lm(log(miR.122copies)~BMI,data))
summary(lm(log(miR.122copies)~Eth,data))

##HMGB1
qqp(data$HMGB1,"norm")#
qqp(data$HMGB1,"lnorm")
qqp(sqrt(data$HMGB1),"norm")

outlierKD(data, HMGB1)#1

lmm <- lmer(HMGB1 ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(HMGB1 ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

##ccK18
qqp(data$ccK18,"norm")
qqp(data$ccK18,"lnorm")
qqp(sqrt(data$ccK18),"norm")#

data$sqccK18=sqrt(data$ccK18)
outlierKD(data, sqccK18)#0

lmm <- lmer(sqccK18 ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(sqccK18 ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)
##FL-K18
qqp(data$FL.K18,"norm")
qqp(data$FL.K18,"lnorm")
qqp(sqrt(data$FL.K18),"norm")

data$sqFL.K18=sqrt(data$FL.K18)
outlierKD(data, sqFL.K18)#0

lmm <- lmer(sqFL.K18 ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(sqFL.K18 ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

##GLDH
qqp(data$GLDH,"norm")
qqp(data$GLDH,"lnorm")#
qqp(sqrt(data$GLDH),"norm")

data$lGLDH=log(data$GLDH)
outlierKD(data, lGLDH)#11

lmm <- lmer(lGLDH ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(lGLDH ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

##CSF-1
qqp(data$CSF.1,"norm")
qqp(data$CSF.1,"lnorm")
qqp(sqrt(data$CSF.1),"norm")#

data$sqCSF.1=sqrt(data$CSF.1)
outlierKD(data, sqCSF.1)#0

lmm <- lmer(sqCSF.1 ~ Sex + Age + Height + Weight + BMI + Eth + (1 | ID), data = data,
            REML = FALSE)
summary(lmm)
lmm0 <- lm(sqCSF.1 ~ Sex + Age + Height + Weight + BMI + Eth , data = data)
Anova(lmm)
exactLRT(m=lmm,m0=lmm0)

###Quantile Regression 
install.packages("quantreg")
install.packages("jmuOutlier")
library(quantreg)
library(jmuOutlier)


x <- seq(0,100,length.out = 100)
qs <- 1:9/10
qr2 <- rq(sqALT ~ Age, data=data, tau = qs)
plot(summary(qr2), parm="Age")

rbind(c(quantile(data$sqALT, 0.025, type = 7),q.CI(data$sqALT,0.025,0.9)),c(quantile(data$sqALT, 0.5, type = 7),q.CI(data$sqALT,0.5,0.9)),c(quantile(data$sqALT, 0.975, type = 7),q.CI(data$sqALT,0.975,0.9)))

cbind(c(quantile(data$sqALT, c(0.025,0.5,0.975), type = 7)),
rbind(quantileCI(data$sqALT, probs = 0.025, conf.level = 0.9),
quantileCI(data$sqALT, probs = 0.5, conf.level = 0.9),
quantileCI(data$sqALT, probs = 0.975, conf.level = 0.9)))
