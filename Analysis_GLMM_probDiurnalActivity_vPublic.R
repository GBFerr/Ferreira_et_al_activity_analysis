#####################################################################################

# GLM(M) to estimate probability of diurnal activity

# Paper: Ferreira et al, Limited temporal response of Cerrado mammals to anthropogenic pressure 
# in areas under contrasting levels of protection

# Fig 4

#####################################################################################
setwd("C:/Users/#####/CodeAndData")

# Approach:
# 1) Run full GLMM model for spp recorded >30 sites
# 2) If random effect different from 0 fit GLMMs with all possible combination of variables (no interaction),
# if not, fit GLMs with all possible combination of variables (no interaction).
# 3) Get model selection table - not provided in the paper
# 4) Get model avg estimates (conditional and full) for model set for summed AICw <=0.95 # full avg used in Fig 4
# 5) Get AIC sum for variables - not provided in the paper

# GLMM matrix
d <- read.csv("GLM_Matx_dayANDcoreAct.csv")
library(lme4)
library(MuMIn)
options(na.action = "na.fail")
#### Probability of diurnal activity
#### grey brock deer; Mgou - GLMM ####
# response var
dayrecs.Mgou <- d[d$sppCode == "Mgou", "day_records"] 
nightrecs.Mgou <- d[d$sppCode == "Mgou", "night_records"] 

# explanatory var
PA_type.Mgou <- d[d$sppCode == "Mgou", "pa_type"]
Dist_humset.Mgou <- d[d$sppCode == "Mgou", "Dist_humset"]
NDVI.Mgou <- d[d$sppCode == "Mgou", "NDVImean_500m"]
Dist_wat.Mgou <- d[d$sppCode == "Mgou", "Dist_water"]
# random var
CT_array.Mgou <-as.vector(d[d$sppCode == "Mgou", "survey_area"])

# full model
mFull.mix.Mgou <- glmer(cbind(dayrecs.Mgou, nightrecs.Mgou) ~ 
                          PA_type.Mgou + Dist_humset.Mgou + NDVI.Mgou + Dist_wat.Mgou
                        + (1 | CT_array.Mgou),
                        family = binomial)
summary(mFull.mix.Mgou)

drg.Mgou <- dredge(mFull.mix.Mgou, rank = AICc) # creates the 16 possible models
avg.Mgou <-model.avg(drg.Mgou, subset = cumsum(weight) <= .95)
summary(avg.Mgou)

# importance of variables
import.Mgou <- as.data.frame(importance(drg.Mgou))
import.Mgou$sppCode <- rep("Mgou",4)

#### crab-eating fox; Ctho - GLMM ####
# response var
dayrecs.Ctho <- d[d$sppCode == "Ctho", "day_records"] 
nightrecs.Ctho <- d[d$sppCode == "Ctho", "night_records"] 

# explanatory var
PA_type.Ctho <- as.factor(d[d$sppCode == "Ctho", "pa_type"])
Dist_humset.Ctho <- d[d$sppCode == "Ctho", "Dist_humset"]
NDVI.Ctho <- d[d$sppCode == "Ctho", "NDVImean_500m"]
Dist_wat.Ctho <- d[d$sppCode == "Ctho", "Dist_water"]
# random var
CT_array.Ctho <-as.vector(d[d$sppCode == "Ctho", "survey_area"])

# full model
mFull.mix.Ctho <- glmer(cbind(dayrecs.Ctho, nightrecs.Ctho) ~ 
                          PA_type.Ctho + Dist_humset.Ctho + NDVI.Ctho + Dist_wat.Ctho
                        + (1 | CT_array.Ctho),
                        family = binomial)
summary(mFull.mix.Ctho)
drg.Ctho <- dredge(mFull.mix.Ctho, rank = AICc)
avg.Ctho <-model.avg(drg.Ctho, subset = cumsum(weight) <= .95)
summary(avg.Ctho)

# importance of variables
import.Ctho <- as.data.frame(importance(drg.Ctho))
import.Ctho$sppCode <- rep("Ctho",4)

#### Dasypus armadillos; Dnov - GLM  ####
# response var
dayrecs.Dnov <- d[d$sppCode == "Dnov", "day_records"] 
nightrecs.Dnov <- d[d$sppCode == "Dnov", "night_records"] 

# explanatory var
PA_type.Dnov <- d[d$sppCode == "Dnov", "pa_type"]
Dist_humset.Dnov <- d[d$sppCode == "Dnov", "Dist_humset"]
NDVI.Dnov <- d[d$sppCode == "Dnov", "NDVImean_500m"]
Dist_wat.Dnov <- d[d$sppCode == "Dnov", "Dist_water"]
# random var
CT_array.Dnov <-as.vector(d[d$sppCode == "Dnov", "survey_area"])

# full model
mFull.Dnov <- glm(cbind(dayrecs.Dnov, nightrecs.Dnov) ~ 
                          PA_type.Dnov + Dist_humset.Dnov + NDVI.Dnov + Dist_wat.Dnov,
                        family = binomial)
summary(mFull.Dnov)

drg.Dnov <- dredge(mFull.Dnov, rank = AICc)
avg.Dnov <-model.avg(drg.Dnov, subset = cumsum(weight) <= .95)
summary(avg.Dnov)

# importance of variables
import.Dnov <- as.data.frame(importance(drg.Dnov))
import.Dnov$sppCode <- rep("Dnov",4)

#### puma; Pcon - GLM  ####
# response var
dayrecs.Pcon <- d[d$sppCode == "Pcon", "day_records"] 
nightrecs.Pcon <- d[d$sppCode == "Pcon", "night_records"] 

# explanatory var
PA_type.Pcon <- d[d$sppCode == "Pcon", "pa_type"]
Dist_humset.Pcon <- d[d$sppCode == "Pcon", "Dist_humset"]
NDVI.Pcon <- d[d$sppCode == "Pcon", "NDVImean_500m"]
Dist_wat.Pcon <- d[d$sppCode == "Pcon", "Dist_water"]
# random var
CT_array.Pcon <-as.vector(d[d$sppCode == "Pcon", "survey_area"])

# full model
mFull.Pcon <- glm(cbind(dayrecs.Pcon, nightrecs.Pcon) ~ 
                          PA_type.Pcon + Dist_humset.Pcon + NDVI.Pcon + Dist_wat.Pcon,
                        family = binomial)
summary(mFull.Pcon)

drg.Pcon <- dredge(mFull.Pcon, rank = AICc)
avg.Pcon <-model.avg(drg.Pcon, subset = cumsum(weight) <= .95)
summary(avg.Pcon)

# importance of variables
import.Pcon <- as.data.frame(importance(drg.Pcon))
import.Pcon$sppCode <- rep("Pcon",4)

#### hoary fox; Lvet - GLMM ####
# response var
dayrecs.Lvet <- d[d$sppCode == "Lvet", "day_records"] 
nightrecs.Lvet <- d[d$sppCode == "Lvet", "night_records"] 

# explanatory var
PA_type.Lvet <- d[d$sppCode == "Lvet", "pa_type"]
Dist_humset.Lvet <- d[d$sppCode == "Lvet", "Dist_humset"]
NDVI.Lvet <- d[d$sppCode == "Lvet", "NDVImean_500m"]
Dist_wat.Lvet <- d[d$sppCode == "Lvet", "Dist_water"]
# random var
CT_array.Lvet <-as.vector(d[d$sppCode == "Lvet", "survey_area"])

# full model
mFull.mix.Lvet <- glmer(cbind(dayrecs.Lvet, nightrecs.Lvet) ~ 
                          Dist_humset.Lvet + NDVI.Lvet + Dist_wat.Lvet
                        + (1 | CT_array.Lvet),
                        family = binomial)
summary(mFull.mix.Lvet)

drg.Lvet <- dredge(mFull.mix.Lvet, rank = AICc)
avg.Lvet <-model.avg(drg.Lvet, subset = cumsum(weight) <= .95)
summary(avg.Lvet)

# importance of variables
import.Lvet <- as.data.frame(importance(drg.Lvet))
import.Lvet$sppCode <- rep("Lvet",3)

#### agouti; Daza - GLM ####
# response var
dayrecs.Daza <- d[d$sppCode == "Daza", "day_records"] 
nightrecs.Daza <- d[d$sppCode == "Daza", "night_records"] 

# explanatory var
PA_type.Daza <- d[d$sppCode == "Daza", "pa_type"]
Dist_humset.Daza <- d[d$sppCode == "Daza", "Dist_humset"]
NDVI.Daza <- d[d$sppCode == "Daza", "NDVImean_500m"]
Dist_wat.Daza <- d[d$sppCode == "Daza", "Dist_water"]
# random var
CT_array.Daza <-as.vector(d[d$sppCode == "Daza", "survey_area"])

# full model
mFull.Daza <- glm(cbind(dayrecs.Daza, nightrecs.Daza) ~ 
                          PA_type.Daza + Dist_humset.Daza + NDVI.Daza + Dist_wat.Daza,
                        family = binomial)
summary(mFull.Daza)

drg.Daza <- dredge(mFull.Daza, rank = AICc)
avg.Daza <-model.avg(drg.Daza, subset = cumsum(weight) <= .95)
summary(avg.Daza)

# importance of variables
import.Daza <- as.data.frame(importance(drg.Daza))
import.Daza$sppCode <- rep("Daza",4)

#### oncilla; Ltig - GLM ####
# response var
dayrecs.Ltig <- d[d$sppCode == "Ltig", "day_records"] 
nightrecs.Ltig <- d[d$sppCode == "Ltig", "night_records"] 

# explanatory var
PA_type.Ltig <- d[d$sppCode == "Ltig", "pa_type"]
Dist_humset.Ltig <- d[d$sppCode == "Ltig", "Dist_humset"]
NDVI.Ltig <- d[d$sppCode == "Ltig", "NDVImean_500m"]
Dist_wat.Ltig <- d[d$sppCode == "Ltig", "Dist_water"]
# random var
CT_array.Ltig <-as.vector(d[d$sppCode == "Ltig", "survey_area"])

# full model
mFull.Ltig <- glm(cbind(dayrecs.Ltig, nightrecs.Ltig) ~ 
                          PA_type.Ltig + Dist_humset.Ltig + NDVI.Ltig + Dist_wat.Ltig,
                        family = binomial)
summary(mFull.Ltig)

drg.Ltig <- dredge(mFull.Ltig, rank = AICc)
avg.Ltig <-model.avg(drg.Ltig, subset = cumsum(weight) <= .95)
summary(avg.Ltig)

# importance of variables
import.Ltig <- as.data.frame(importance(drg.Ltig))
import.Ltig$sppCode <- rep("Ltig",4)

#### ocelot; Lpar - GLM ####
# response var
dayrecs.Lpar <- d[d$sppCode == "Lpar", "day_records"] 
nightrecs.Lpar <- d[d$sppCode == "Lpar", "night_records"] 

# explanatory var
PA_type.Lpar <- d[d$sppCode == "Lpar", "pa_type"]
Dist_humset.Lpar <- d[d$sppCode == "Lpar", "Dist_humset"]
NDVI.Lpar <- d[d$sppCode == "Lpar", "NDVImean_500m"]
Dist_wat.Lpar <- d[d$sppCode == "Lpar", "Dist_water"]
# random var
CT_array.Lpar <-as.vector(d[d$sppCode == "Lpar", "survey_area"])

# full model
mFull.Lpar <- glm(cbind(dayrecs.Lpar, nightrecs.Lpar) ~ 
                    PA_type.Lpar + Dist_humset.Lpar + NDVI.Lpar + Dist_wat.Lpar,
                  family = binomial)
summary(mFull.Lpar)

drg.Lpar <- dredge(mFull.Lpar, rank = AICc)
avg.Lpar <-model.avg(drg.Lpar, subset = cumsum(weight) <= .95)
summary(avg.Lpar)

# importance of variables
import.Lpar <- as.data.frame(importance(drg.Lpar))
import.Lpar$sppCode <- rep("Lpar",4)

#### collared peccary; Ptaj - GLM ####
# response var
dayrecs.Ptaj <- d[d$sppCode == "Ptaj", "day_records"] 
nightrecs.Ptaj <- d[d$sppCode == "Ptaj", "night_records"] 

# explanatory var
PA_type.Ptaj <- d[d$sppCode == "Ptaj", "pa_type"]
Dist_humset.Ptaj <- d[d$sppCode == "Ptaj", "Dist_humset"]
NDVI.Ptaj <- d[d$sppCode == "Ptaj", "NDVImean_500m"]
Dist_wat.Ptaj <- d[d$sppCode == "Ptaj", "Dist_water"]
# random var
CT_array.Ptaj <-as.vector(d[d$sppCode == "Ptaj", "survey_area"])

# full model
mFull.Ptaj <- glm(cbind(dayrecs.Ptaj, nightrecs.Ptaj) ~ 
                    PA_type.Ptaj + Dist_humset.Ptaj + NDVI.Ptaj + Dist_wat.Ptaj,
                  family = binomial)
summary(mFull.Ptaj)

drg.Ptaj <- dredge(mFull.Ptaj, rank = AICc)
avg.Ptaj <-model.avg(drg.Ptaj, subset = cumsum(weight) <= .95)
summary(avg.Ptaj)

# importance of variables
import.Ptaj <- as.data.frame(importance(drg.Ptaj))
import.Ptaj$sppCode <- rep("Ptaj",4)

#### giant anteater; Mtri - GLM ####
# response var
dayrecs.Mtri <- d[d$sppCode == "Mtri", "day_records"] 
nightrecs.Mtri <- d[d$sppCode == "Mtri", "night_records"] 

# explanatory var
PA_type.Mtri <- d[d$sppCode == "Mtri", "pa_type"]
Dist_humset.Mtri <- d[d$sppCode == "Mtri", "Dist_humset"]
NDVI.Mtri <- d[d$sppCode == "Mtri", "NDVImean_500m"]
Dist_wat.Mtri <- d[d$sppCode == "Mtri", "Dist_water"]
# random var
CT_array.Mtri <-as.vector(d[d$sppCode == "Mtri", "survey_area"])

# full model
mFull.Mtri <- glm(cbind(dayrecs.Mtri, nightrecs.Mtri) ~ 
                    PA_type.Mtri + Dist_humset.Mtri + NDVI.Mtri + Dist_wat.Mtri,
                  family = binomial)
summary(mFull.Mtri)

drg.Mtri <- dredge(mFull.Mtri, rank = AICc)
avg.Mtri <-model.avg(drg.Mtri, subset = cumsum(weight) <= .95)
summary(avg.Mtri)

# importance of variables
import.Mtri <- as.data.frame(importance(drg.Mtri))
import.Mtri$sppCode <- rep("Mtri",4)

#### tapir; Tter - GLM ####
# response var
dayrecs.Tter <- d[d$sppCode == "Tter", "day_records"] 
nightrecs.Tter <- d[d$sppCode == "Tter", "night_records"] 

# explanatory var
PA_type.Tter <- d[d$sppCode == "Tter", "pa_type"]
Dist_humset.Tter <- d[d$sppCode == "Tter", "Dist_humset"]
NDVI.Tter <- d[d$sppCode == "Tter", "NDVImean_500m"]
Dist_wat.Tter <- d[d$sppCode == "Tter", "Dist_water"]
# random var
CT_array.Tter <-as.vector(d[d$sppCode == "Tter", "survey_area"])

# full model
mFull.Tter <- glm(cbind(dayrecs.Tter, nightrecs.Tter) ~ 
                   Dist_humset.Tter + NDVI.Tter + Dist_wat.Tter,  # no PA variable
                  family = binomial)
summary(mFull.Tter)

drg.Tter <- dredge(mFull.Tter, rank = AICc)
avg.Tter <-model.avg(drg.Tter, subset = cumsum(weight) <= .95)
summary(avg.Tter)

# importance of variables
import.Tter <- as.data.frame(importance(drg.Tter))
import.Tter$sppCode <- rep("Tter",3)

#### skunk; Csem - GLM ####
dayrecs.Csem <- d[d$sppCode == "Csem", "day_records"] 
nightrecs.Csem <- d[d$sppCode == "Csem", "night_records"] 

# explanatory var
PA_type.Csem <- d[d$sppCode == "Csem", "pa_type"]
Dist_humset.Csem <- d[d$sppCode == "Csem", "Dist_humset"]
NDVI.Csem <- d[d$sppCode == "Csem", "NDVImean_500m"]
Dist_wat.Csem <- d[d$sppCode == "Csem", "Dist_water"]
# random var
CT_array.Csem <-as.vector(d[d$sppCode == "Csem", "survey_area"])

# full model

mFull.Csem <- glm(cbind(dayrecs.Csem, nightrecs.Csem) ~ 
                          PA_type.Csem + Dist_humset.Csem + NDVI.Csem + Dist_wat.Csem,
                        family = binomial)
summary(mFull.Csem)

drg.Csem <- dredge(mFull.Csem, rank = AICc)
avg.Csem <-model.avg(drg.Csem, subset = cumsum(weight) <= .95)
summary(avg.Csem)

# importance of variables
import.Csem <- as.data.frame(importance(drg.Csem))
import.Csem$sppCode <- rep("Csem",4)

#### yellow armadillo; Esex - GLM ####
# response var
dayrecs.Esex <- d[d$sppCode == "Esex", "day_records"] 
nightrecs.Esex <- d[d$sppCode == "Esex", "night_records"] 

# explanatory var
PA_type.Esex <- d[d$sppCode == "Esex", "pa_type"]
Dist_humset.Esex <- d[d$sppCode == "Esex", "Dist_humset"]
NDVI.Esex <- d[d$sppCode == "Esex", "NDVImean_500m"]
Dist_wat.Esex <- d[d$sppCode == "Esex", "Dist_water"]
# random var
CT_array.Esex <-as.vector(d[d$sppCode == "Esex", "survey_area"])

# warning message with Hessian problems and huge SE for PA_type
# full model

mFull.Esex <- glm(cbind(dayrecs.Esex, nightrecs.Esex) ~ 
                          Dist_humset.Esex + NDVI.Esex + Dist_wat.Esex,
                        family = binomial)
summary(mFull.Esex)

drg.Esex <- dredge(mFull.Esex, rank = AICc)
avg.Esex <-model.avg(drg.Esex, subset = cumsum(weight) <= .95)
summary(avg.Esex)

# importance of variables
import.Esex <- as.data.frame(importance(drg.Esex))
import.Esex$sppCode <- rep("Esex",3)

#### maned wolf; Cbra - GLM ####
# response var
dayrecs.Cbra <- d[d$sppCode == "Cbra", "day_records"] 
nightrecs.Cbra <- d[d$sppCode == "Cbra", "night_records"] 

# explanatory var
PA_type.Cbra <- d[d$sppCode == "Cbra", "pa_type"]
Dist_humset.Cbra <- d[d$sppCode == "Cbra", "Dist_humset"]
NDVI.Cbra <- d[d$sppCode == "Cbra", "NDVImean_500m"]
Dist_wat.Cbra <- d[d$sppCode == "Cbra", "Dist_water"]
# random var
CT_array.Cbra <-as.vector(d[d$sppCode == "Cbra", "survey_area"])

# full model
mFull.Cbra <- glm(cbind(dayrecs.Cbra, nightrecs.Cbra) ~ 
                          PA_type.Cbra + Dist_humset.Cbra + NDVI.Cbra + Dist_wat.Cbra,
                        family = binomial)
summary(mFull.Cbra)

drg.Cbra <- dredge(mFull.Cbra, rank = AICc)
avg.Cbra <-model.avg(drg.Cbra, subset = cumsum(weight) <= .95)
summary(avg.Cbra)
str(avg.Cbra)
# importance of variables
import.Cbra <- as.data.frame(importance(drg.Cbra))
import.Cbra$sppCode <- rep("Cbra",4)

#### tamandua; Ttet - GLM ####
# response var
dayrecs.Ttet <- d[d$sppCode == "Ttet", "day_records"] 
nightrecs.Ttet <- d[d$sppCode == "Ttet", "night_records"] 

# explanatory var
PA_type.Ttet <- d[d$sppCode == "Ttet", "pa_type"]
Dist_humset.Ttet <- d[d$sppCode == "Ttet", "Dist_humset"]
NDVI.Ttet <- d[d$sppCode == "Ttet", "NDVImean_500m"]
Dist_wat.Ttet <- d[d$sppCode == "Ttet", "Dist_water"]
# random var
CT_array.Ttet <-as.vector(d[d$sppCode == "Ttet", "survey_area"])

# full model
mFull.Ttet <- glm(cbind(dayrecs.Ttet, nightrecs.Ttet) ~ 
                          Dist_humset.Ttet + NDVI.Ttet + Dist_wat.Ttet,
                        family = binomial)
summary(mFull.Ttet)

drg.Ttet <- dredge(mFull.Ttet, rank = AICc)
avg.Ttet <-model.avg(drg.Ttet, subset = cumsum(weight) <= .95)
summary(avg.Ttet)

# importance of variables
import.Ttet <- as.data.frame(importance(drg.Ttet))
import.Ttet$sppCode <- rep("Ttet",3)

save.image("C:/Users/Guilherme/Google Drive/My papers/ActivityPattern_vs_Pressure_SVP_Chap4PhD/CodeAndData/GLMM_probDiurnalActivity_ENV.RData")

