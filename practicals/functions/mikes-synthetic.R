library(igraph)
library(causaleffect)
library(TruncatedNormal)
library(riskRegression)
library(plotROC)
library(ids)



pa <- function(dag, v) {
  neighbors(dag, v, mode = "in")
}

ch <- function(dag, v) {
  neighbors(dag, v, mode = "out")
}

exog <- function(dag) {
  lapply(adjacent_vertices(dag, V(dag), mode = "in"), \(d) length(d) == 0)
}

#exog(dag)

#names(V(dag))[unlist(exog(dag))]

dag2 <- graph_from_literal(
  hidden_comorbidity,
  sex_male,
  age +- sex_male,
  diabetes_duration +- sex_male:age,
  time_cens +- sex_male:age:diabetes_duration,
  motion +- sex_male:age:diabetes_duration,
  smoking +- sex_male:age:diabetes_duration,
  HBA1C_time1 +- diabetes_duration:age:sex_male,
  HBA1C_time2 +- HBA1C_time1:motion:smoking:sex_male,
  urine_albumin_time1 +- age:sex_male:diabetes_duration:HBA1C_time1,
  urine_albumin_time2 +- urine_albumin_time1:age:sex_male:hidden_comorbidity,
  LDL_time1 +- age:sex_male:motion:hidden_comorbidity,
  statin +- hidden_comorbidity:LDL_time1:motion,
  LDL_time2 +- LDL_time1:statin:sex_male:motion,
  SBP_time1 +- age:sex_male:motion:smoking:hidden_comorbidity:diabetes_duration:HBA1C_time1,
  SBP_time2 +- age:sex_male:motion:smoking:hidden_comorbidity:SBP_time1:
    statin:HBA1C_time2:HBA1C_time1,
  eGFR_time1 +- SBP_time1:age:sex_male:
    diabetes_duration:HBA1C_time1,
  eGFR_time2 +- SBP_time2:urine_albumin_time1:age:sex_male:
    HBA1C_time1:HBA1C_time2,
  cvd_5year +- LDL_time2:age:diabetes_duration:sex_male:
    smoking:motion:hidden_comorbidity:SBP_time2,
  time_esrd +- urine_albumin_time1:urine_albumin_time2:eGFR_time1:
    eGFR_time2:HBA1C_time1:HBA1C_time2:sex_male

)

plot(dag2)
exog(dag2)

# cat(sapply(V(dag2)$name, \(vs) {
#   pas <- pa(dag2, vs)$name
#   paste(vs, " = c(", paste(pas, collapse = "=, "), "), \n")
# }), sep = "")

coefficients <- list(hidden_comorbidity  = c(  ),
                     sex_male  = c(  ),
                     age  = c( sex_male = -.1),
                     diabetes_duration  = c( sex_male=.5, age = .8),
                     time_cens  = c( sex_male=-.3, age=-.03, diabetes_duration = -.05),
                     motion  = c( sex_male=.3, age=-.1, diabetes_duration=-.5 ),
                     smoking  = c( sex_male=.1, age=-.3, diabetes_duration=-.7 ),
                     HBA1C_time1  = c( sex_male=.4, age=.1, diabetes_duration =.8),
                     HBA1C_time2  = c( sex_male=.4, motion=-1.2, smoking=.8, HBA1C_time1 = 1.25),
                     urine_albumin_time1  = c( sex_male=.4, age=.1, diabetes_duration=.3, HBA1C_time1=.9 ),
                     urine_albumin_time2  = c( hidden_comorbidity=.2, sex_male=.3, age=.1, urine_albumin_time1 =1.1),
                     LDL_time1  = c( hidden_comorbidity=.25, sex_male=1.25, age=.15, motion =-1),
                     statin  = c( hidden_comorbidity=.8, motion=.2, LDL_time1 =1.3),
                     LDL_time2  = c( sex_male=1.25, motion=-1, LDL_time1=1, statin=-2 ),
                     SBP_time1  = c( hidden_comorbidity=.1, sex_male=1, age=.01,
                                     diabetes_duration=.05, motion=-2, smoking=1.5, HBA1C_time1=.05 ),
                     SBP_time2  = c( hidden_comorbidity=.1, sex_male=1.5, age=.02, motion=-2,
                                     smoking=1.25, HBA1C_time1=.2, HBA1C_time2=.5, statin=-4, SBP_time1 =1.02),
                     eGFR_time1  = c( sex_male=.5, age=-.01, diabetes_duration=-.08, HBA1C_time1=-.8, SBP_time1=-1.5 ),
                     eGFR_time2  = c( sex_male=1.2, age=-.05, HBA1C_time1=-.7, HBA1C_time2=-1.2,
                                      urine_albumin_time1=-.005, SBP_time2=-.1 ),
                     cvd_5year  = c( hidden_comorbidity=.8, sex_male=.4, age=.1, diabetes_duration=.1,
                                    motion=-.3, smoking=.2, LDL_time2=.6, SBP_time2 = .2 ),
                     time_esrd  = c( sex_male=.1, HBA1C_time1=.05, HBA1C_time2=.05,
                                     urine_albumin_time1=.4, urine_albumin_time2=.6, eGFR_time1=-1.5, eGFR_time2=-2.25 )
                     )

generate_data <- function(n = 1000, coefficients, dag, intervene = NULL) {

  dag_order <- topo_sort(dag, mode = "out")
  thisdata <- data.frame(dummy = rep(1, n))

  generator <- list(
    hidden_comorbidity = function(thisdata) {
      rnorm(n)

    },
    sex_male = function(thisdata) {
      rbinom(n, 1, .45)
    },
    age = function(thisdata) {

      lpred <- c(thisdata[,pa(dag, "age")$name] * coefficients[["age"]][pa(dag, "age")$name])
      sapply(lpred, \(lp) rtnorm(1, lp + 42, 7.3, 18, 100))

    },
    diabetes_duration = function(thisdata) {

      pas <- pa(dag, "diabetes_duration")$name
      lpred <- -5 + c(as.matrix(thisdata[, pas]) %*% coefficients$diabetes_duration[pas] )
      res <- sapply(1:nrow(thisdata), \(i) {
        rtnorm(1, lpred[i], 3.2, 0, thisdata$age[i])
      })
      res

    },
    time_cens = function(thisdata) {
      pas <- pa(dag, "time_cens")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$time_cens[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(exp(gamma + alpha)) - 15
      }, interval = c(-50, 50))$root

      rexp(n, 1 / exp(alpha0 + gamma))

    },

    motion = function(thisdata) {
      pas <- pa(dag, "motion")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$motion[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .65
      }, interval = c(-500, 500))$root
      rbinom(n, 1, plogis(alpha0 + gamma))

    },
    smoking = function(thisdata) {
      pas <- pa(dag, "smoking")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$smoking[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .065
      }, interval = c(-500, 500))$root
      rbinom(n, 1, plogis(alpha0 + gamma))

    },
    HBA1C_time1 = function(thisdata) {
      pas <- pa(dag, "HBA1C_time1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$HBA1C_time1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 52
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 15, 0, 250))

    },
    LDL_time1 = function(thisdata) {
      pas <- pa(dag, "LDL_time1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$LDL_time1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 2.5
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 1.25, 0, 15))

    },
    HBA1C_time2 = function(thisdata) {
      pas <- pa(dag, "HBA1C_time1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$HBA1C_time1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 58
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 20, 0, 250))

    },
    urine_albumin_time1 = function(thisdata) {
      pas <- pa(dag, "urine_albumin_time1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$urine_albumin_time1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 42
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, abs(lp), 0, 500))

    },
    statin = function(thisdata) {
      pas <- pa(dag, "statin")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$statin[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .25
      }, interval = c(-500, 500))$root
      rbinom(n, 1, plogis(alpha0 + gamma))

    },
    LDL_time2 = function(thisdata) {
      pas <- pa(dag, "LDL_time2")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$LDL_time2[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 3.15
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 1.25, 0, 15))

    },
    urine_albumin_time2 = function(thisdata) {
      pas <- pa(dag, "urine_albumin_time2")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$urine_albumin_time2[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 62
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, abs(lp), 0, 500))
    },
    SBP_time1 = function(thisdata) {
      pas <- pa(dag, "SBP_time1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$SBP_time1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 115
      }, interval = c(-5000, 5000))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 6, 85, 180))
    },
    SBP_time2 = function(thisdata) {
      pas <- pa(dag, "SBP_time2")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$SBP_time2[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 120
      }, interval = c(-5000, 5000))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 6, 85, 190))

    },
    eGFR_time1 = function(thisdata) {
      pas <- pa(dag, "eGFR_time1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$eGFR_time1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 100
      }, interval = c(-5000, 5000))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 10, 0, 180))
    },
    eGFR_time2 = function(thisdata) {
      pas <- pa(dag, "eGFR_time2")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$eGFR_time2[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 100
      }, interval = c(-5000, 5000))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 10, 0, 180))
    },

    cvd_5year = function(thisdata) {
      pas <- pa(dag, "cvd_5year")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$cvd_5year[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .05
      }, interval = c(-5000, 5000))$root

      rbinom(n, 1, plogis(alpha0 + gamma + rnorm(n, sd = 3)))

    },
    time_esrd = function(thisdata) {
      pas <- pa(dag, "time_esrd")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$time_esrd[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(exp(gamma + alpha)) - 17
      }, interval = c(-5000, 5000))$root

      rexp(n, 1 / exp(alpha0 + gamma))

    }

  )

  for(vs in dag_order$name) {
    if(!is.null(intervene)) {
      if(vs %in% names(intervene)) {
        thisdata[[vs]] <- intervene[[vs]]
        next
      }
    }

    thisdata[[vs]] <- generator[[vs]](thisdata)
  }

  thisdata$dummy <- NULL
  thisdata

}



## tests
set.seed(120825)
data2 <- generate_data(1000, coefficients = coefficients, dag = dag2)

newdata <- data2[c(1:1000, 1:1000),]
newdata$statin <- (rep(c(0, 1), each = 1000))

causal.effect("cvd_5year", "statin", G = dag2)
dagitty::adjustmentSets(SEMgraph::graph2dagitty(dag2),
                        exposure = "statin", outcome = "cvd_5year", type = "canonical",
                        max.results = 5)

glmtrue <- glm(cvd_5year ~ statin + HBA1C_time1 + HBA1C_time2 + LDL_time1 + SBP_time1 + age +
                 diabetes_duration + hidden_comorbidity + motion + sex_male + smoking,
               data = data2, family = "binomial")

glmomit <- glm(cvd_5year ~ statin + HBA1C_time1 + HBA1C_time2 + LDL_time1 + SBP_time1 + age +
                 diabetes_duration + motion + sex_male + smoking, data = data2,
               family = "binomial")

glmmed <- glm(cvd_5year ~ statin + HBA1C_time1 + HBA1C_time2 + LDL_time1 + SBP_time1 + age +
                diabetes_duration + motion + sex_male + smoking + SBP_time2 + LDL_time2 , data = data2,
               family = "binomial")


preds <- list(pt = predict(glmtrue, newdata = newdata, type = "response"),
po = predict(glmomit, newdata = newdata, type = "response"),
pm = predict(glmmed, newdata = newdata, type = "response"))

lapply(preds, \(pt) {
  c(mean(pt[1:1000]), mean(pt[1001:2000]))
})

ggplot(data2, aes(m = predict(glmmed, type = "response"), d = cvd_5year)) + geom_roc()

Score(list(glmomit), data = data2, formula = cvd_5year ~ 1,
      summary = "risk", split.method = "cv10") |> summary()



data2$pid <- paste0("train-", random_id(n = nrow(data2), bytes = 4))

datatrain <- data2[, c("pid", "age", "sex_male", "diabetes_duration", "smoking", "motion", "HBA1C_time1",
                       "urine_albumin_time1", "LDL_time1", "SBP_time1", "eGFR_time1", "statin",
                       "HBA1C_time2",
                       "urine_albumin_time2", "LDL_time2", "SBP_time2", "eGFR_time2", "cvd_5year")]

saveRDS(datatrain, file = "../type1-diabetes-train.rds")

## validation data

set.seed(8256)
datav <- generate_data(400, coefficients = coefficients, dag = dag2)
datav$pid <- paste0("valid-", random_id(n = nrow(datav), bytes = 4), "d")
datav$statin <- 0

datav2 <- datav
datav2$statin <- 1
datav2$pid <- paste0(substr(datav2$pid, 1, nchar(datav2$pid[1]) - 1), "p")

valid <- rbind(datav, datav2)

roc <- data.frame(m = c(predict(glmtrue, newdata = valid, type = "response"),
                        predict(glmomit, newdata = valid, type = "response"),
                        predict(glmmed, newdata = valid, type = "response")),
                  d = rep(valid$cvd_5year, 3), model = rep(c("true", "real", "mediate"), each = nrow(valid))) |>
  ggplot(aes(m = m, d = d, color = model)) + geom_roc()
roc

valid <- valid[sample(1:nrow(valid)),c("pid", "age", "sex_male", "hidden_comorbidity", "diabetes_duration", "smoking", "motion", "HBA1C_time1",
                                       "urine_albumin_time1", "LDL_time1", "SBP_time1", "eGFR_time1", "statin",
                                       "HBA1C_time2",
                                       "urine_albumin_time2", "LDL_time2", "SBP_time2", "eGFR_time2", "cvd_5year")]
saveRDS(valid, file = "../type1-diabetes-valid.rds")

## try under intervention


set.seed(8256)
dataint <- generate_data(
  200,
  coefficients = coefficients,
  dag = dag2,
  interven = list(statin = rep(c(0, 1), each = 100))
)[, c(
  "age",
  "sex_male",
  "diabetes_duration",
  "smoking",
  "motion",
  "HBA1C_time1",
  "urine_albumin_time1",
  "LDL_time1",
  "SBP_time1",
  "eGFR_time1",
  "statin",
  "HBA1C_time2",
  "urine_albumin_time2",
  "LDL_time2",
  "SBP_time2",
  "eGFR_time2",
  "cvd_5year"
)]

by(dataint$cvd_5year, dataint$statin, mean)

saveRDS(dataint, file = "../type1-diabetes-statin-intervene.rds")




