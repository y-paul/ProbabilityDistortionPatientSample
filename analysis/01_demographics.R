library(tidyverse)



control <- read.csv("data/processed/control/overall/quest_dat.csv")
patient <- read.csv("data/processed/patient/quest_dat.csv")
control$group <- "control"
patient$group <- rep("Patient", times=60)





dat <- rbind(control, patient)

table(dat$diagnosis)

unique(dat$education)



dat$ed_sum <- dat$education +  dat$occupation







tapply(patient$BDI_GESAMT, patient$diagnosis, FUN=mean)
tapply(patient$GAD_GESAMT, patient$diagnosis, FUN=mean)
tapply(patient$PHQ_GESAMT, patient$diagnosis, FUN=mean)
tapply(patient$PSWQ_GESAMT, patient$diagnosis, FUN=mean, na.rm=T)
tapply(patient$PSWQ_GESAMT, patient$diagnosis, FUN=mean, na.rm=T)



pre_data <- readxl::read_xlsx("data/raw/patient/questionnaire/pre_data.xlsx")
colMeans(pre_data, na.rm=T)

dat$bdi_pre <- c(rep(NA, times=60), pre_data$BDI)
dat$phq_pre <- c(rep(NA, times=60), pre_data$PHQ_Dep)
dat$gad_pre <- c(rep(NA, times=60), pre_data$PHQ_GAS)
patient$bdi_pre <- pre_data$BDI
patient$phq_pre <- pre_data$PHQ_Dep
patient$gad_pre <- pre_data$PHQ_GAS
library(vtable)


sumtable(dat[,c("age", "sex",
              "ed_sum",
              "ERQ_NEUBEWERTUNG_GESAMT",
              "ERQ_SUPRESSION_GESAMT",
              "BDI_GESAMT",
              "PHQ_GESAMT",
              "GAD_GESAMT",
              "UIS_GESAMT",
              "UIS_A",
              "UIS_B",
              "UIS_C",
              "PSWQ_GESAMT",
              "group")],
         add.median=T, group="group",
         group.test=T, digits=2,
         numformat = function(x) {
           if (is.numeric(x)) formatC(x, format = "f", digits = 2) else x
         })



ttestBF(control$BDI_GESAMT, patient$BDI_GESAMT)
ttestBF(control$GAD_GESAMT, patient$GAD_GESAMT)
ttestBF(control$PHQ_GESAMT, patient$PHQ_GESAMT)
ttestBF(control$UIS_GESAMT, patient$UIS_GESAMT)
ttestBF(control$ERQ_NEUBEWERTUNG_GESAMT, patient$ERQ_NEUBEWERTUNG_GESAMT)
ttestBF(control$ERQ_SUPRESSION_GESAMT, patient$ERQ_SUPRESSION_GESAMT)
ttestBF(control$PSWQ_GESAMT, patient$PSWQ_GESAMT)
ttestBF(control$UIS_A, patient$UIS_A)
ttestBF(control$UIS_B, patient$UIS_B)
ttestBF(control$UIS_C, patient$UIS_C)

ttestBF(control$BDI_GESAMT, patient$bdi_pre)
ttestBF(control$GAD_GESAMT, patient$gad_pre)
ttestBF(control$PHQ_GESAMT, patient$phq_pre)




dat_parameters <- read.csv("data/fit/tk_parameters.csv")


vtable::st(dat_parameters[,c( "gamma_de", "gamma_ex", "theta_de", "theta_ex", "iperf_de", "iperf_ex", "group")],
   group="group", group.test=T)



dat_parameters$performance_de_g_level_pat[1]
dat_parameters$performance_ex_g_level_pat[1]
dat_parameters$performance_de_g_level_con[1]
dat_parameters$performance_ex_g_level_con[1]
dat_parameters$gap <- dat_parameters$gamma_de- dat_parameters$gamma_ex

t.test(data=dat_parameters, gamma_ex~group)
t.test(data=dat_parameters, gap~group)


# age demographics


ttestBF(control$age, patient$age)



tab <- table(dat$group, dat$sex)
contingencyTableBF(as.matrix(tab), sampleType="indepMulti", fixedMargin="rows")
chisq.test(tab)


ttestBF(data=dat, formula = ed_sum ~ group, na.rm=T)






