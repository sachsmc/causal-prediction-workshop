library(lava)
library(Publish)
source("functions/synthesize-td1.R")
library(survival)
library(riskRegression)

beta_cs <- .6
beta_ct <- 1.2
beta_ls <- 2
beta_st <- -0.3

set.seed(70825)
data <- synthesize_td1(1000, beta_cs, beta_ct, beta_ls, beta_st)

fittrue <- CSC(Hist(time_cvd, status_cvd) ~ hidden_Comorbidity + Statin + sex + age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+ I(log2(eGFR)) * I(age > 40) +value_Smoking+value_Motion, data = data)


fitomit <- CSC(Hist(time_cvd, status_cvd) ~ Statin + sex + age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+I(log2(eGFR)) * I(age > 40) +value_Smoking+value_Motion, data = data)


check <- coef(fittrue$models[[1]])["Statin"] < 0 & coef(fitomit$models[[1]])["Statin"] > 0
check

set.seed(123)
datanew <- synthesize_td1(1, beta_cs, beta_ct, beta_ls, beta_st)
datanew <- rbind(datanew, datanew)
datanew$Statin[1] <- 0

riskRegression::predictRisk(fittrue, newdata = datanew, times = 5, cause = 1)
riskRegression::predictRisk(fitomit, newdata = datanew, times = 5, cause = 1)


glmtrue <- glm(cvd_5year ~ hidden_Comorbidity + Statin + sex + age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+ I(log2(eGFR)) * I(age > 40) +value_Smoking+value_Motion, data = data, 
               family = "binomial")

glmomit <- glm(cvd_5year ~  Statin + sex + age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+ I(log2(eGFR)) * I(age > 40) +value_Smoking+value_Motion, data = data, 
               family = "binomial")


predict(glmtrue, newdata = datanew, type = "response")
predict(glmomit, newdata = datanew, type = "response")

predictRisk(glmtrue, newdata = datanew)
predictRisk(glmomit, newdata = datanew)


library(eventglm)

pfittrue <- cumincglm(Surv(time_cvd, factor(status_cvd)) ~ hidden_Comorbidity + Statin + sex + age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+ I(log2(eGFR)) * I(age > 40) +value_Smoking+value_Motion, data = data, time = 5, cause = "1", link = "logit")


pfitomit <- cumincglm(Surv(time_cvd, factor(status_cvd)) ~ Statin + sex + age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+ I(log2(eGFR)) * I(age > 40) +value_Smoking+value_Motion, data = data, time = 5, cause = "1", link = "logit")


plogis(predict(pfittrue, newdata = datanew))
plogis(predict(pfitomit, newdata = datanew))

