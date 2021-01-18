#####################################################################################

# Code for estimating proportion of time active and overlap in daily activity patterns

# Paper: Ferreira et al, Limited effect of anthropogenic pressure on Cerrado mammals' 
# activity patterns in areas under contrasting levels of protection

# Figs 2 and 3

#####################################################################################
# setwd("#/Users/######/CodeAndData")
# loading data with independent records
indRecs <- read.csv("independentRecs_LargeNativeMammals.csv")
n <- as.data.frame(table(indRecs$Species, indRecs$pa_type)) # number of independent recs/species
names(indRecs)
#### estimating proportion of time active (activity level) at each PAtype ####
# restricted to spp with > 25 rec in each PAtype
library(activity)

# list of spp w/ >25 recs in each PAtype 
spp.index <- c("Mgou","Ctho", "Ltig","Csem", "Dnov", "Esex", "Ptaj", "Daza")
resultsListAct_strict <- list()  # creating placeholder for results, output is a list
resultsListAct_multi <- list()
for(i in spp.index){
  act.strict<- fitact(indRecs[indRecs$sppCode == i & indRecs$pa_type == "strict" , "Time_rad"]
                      ,reps=10000,sample="data",show=TRUE)
    act.multi<- fitact(indRecs[indRecs$sppCode == i & indRecs$pa_type == "multiple-use" , "Time_rad"]
                     ,reps=10000,sample="data",show=TRUE)
  resultsListAct_strict[[i]] <- act.strict
  resultsListAct_multi[[i]] <- act.multi
}

Mgou.strict <- resultsListAct_strict$Mgou@act
Mgou.multi <- resultsListAct_multi$Mgou@act
Ctho.strict <- resultsListAct_strict$Ctho@act
Ctho.multi <- resultsListAct_multi$Ctho@act
Ltig.strict <- resultsListAct_strict$Ltig@act
Ltig.multi <- resultsListAct_multi$Ltig@act
Csem.strict <- resultsListAct_strict$Csem@act
Csem.multi <- resultsListAct_multi$Csem@act
Dnov.strict <- resultsListAct_strict$Dnov@act
Dnov.multi <- resultsListAct_multi$Dnov@act
Esex.strict <- resultsListAct_strict$Esex@act
Esex.multi <- resultsListAct_multi$Esex@act
Ptaj.strict <- resultsListAct_strict$Ptaj@act
Ptaj.multi <- resultsListAct_multi$Ptaj@act
Daza.strict <- resultsListAct_strict$Daza@act
Daza.multi <- resultsListAct_multi$Daza@act


prop.ActResult <- as.data.frame(rbind(Mgou.strict, Mgou.multi,
                    Ctho.strict, Ctho.multi,
                    Ltig.strict, Ltig.multi,
                    Csem.strict, Csem.multi,
                    Dnov.strict, Dnov.multi,
                    Esex.strict, Esex.multi,
                    Ptaj.strict, Ptaj.multi,
                    Daza.strict, Daza.multi))

# Wald test to compare proportion of time active between PA_type
Mgou <- compareAct(list(resultsListAct_strict$Mgou, resultsListAct_multi$Mgou))
Ctho <- compareAct(list(resultsListAct_strict$Ctho, resultsListAct_multi$Ctho))
Ltig <- compareAct(list(resultsListAct_strict$Ltig, resultsListAct_multi$Ltig))
Csem <- compareAct(list(resultsListAct_strict$Csem, resultsListAct_multi$Csem))
Dnov <- compareAct(list(resultsListAct_strict$Dnov, resultsListAct_multi$Dnov))
Esex <- compareAct(list(resultsListAct_strict$Esex, resultsListAct_multi$Esex))
Ptaj <- compareAct(list(resultsListAct_strict$Ptaj, resultsListAct_multi$Ptaj))
Daza <- compareAct(list(resultsListAct_strict$Daza, resultsListAct_multi$Daza))

WaldTestResult <- as.data.frame(rbind(Mgou,Mgou,Ctho,Ctho,Ltig,Ltig,Csem,Csem,Dnov,Dnov,
                        Esex,Esex,Ptaj,Ptaj,Daza,Daza)) 

WaldTestResult$spp <- c("Mgou","Mgou","Ctho","Ctho","Ltig","Ltig","Csem","Csem","Dnov","Dnov",
  "Esex","Esex","Ptaj","Ptaj","Daza","Daza")
prop.ActResult <- cbind(prop.ActResult,WaldTestResult)
prop.ActResult
write.csv(prop.ActResult,"prop_Act_PAtype_8spp.csv")

#### Overlap analysis ####
library(overlap)
# estimates of activity overlap between strict and multi-use PAs ####
# restricted to spp w/ > 25 recs in each PAtype

#### Ptaj - collared peccary ####
# get Dhats
(DhatsPtaj<-overlapEst(indRecs[indRecs$Species == "Pecari tajacu" & indRecs$pa_type == "strict", "Time_rad"], 
                       indRecs[indRecs$Species == "Pecari tajacu" & indRecs$pa_type == "multiple-use", "Time_rad"]))


# smoothed bootstrap for each PAtype; usually use 10,000 reps:
bsPtaj.strict<-resample(indRecs[indRecs$Species == "Pecari tajacu" & indRecs$pa_type == "strict", "Time_rad"],10000) 
bsPtaj.multi<-resample(indRecs[indRecs$Species == "Pecari tajacu" & indRecs$pa_type == "multiple-use", "Time_rad"],10000)

#bootEst
# Following Ridout and Linkie (2009) and Meredith & Ridout 2018  (overview of overlap)) recommendations 
# When smaller sample <50, Dhat1 performed best, while Dhat4 was better when smaller sample >75 (Meredith & Ridout 2018; 
# here I used Dhat 4 also for > 50 recs
bsPtaj<-bootEst(bsPtaj.strict,bsPtaj.multi,adjust=c(NA,1,NA)) #Dhat 4
Ptaj.mean <-mean(bsPtaj[,2])

#Convert column with Dhat to vector and get CIs
bsPtajVec<-as.vector(bsPtaj[,2])#2 because it's Dhat 4
Ptaj.CI <- as.vector(bootCI(DhatsPtaj[2],bsPtaj)['basic0',]) #2 because it's Dhat 4
Ptaj.overlap <- rbind(Ptaj.mean,Ptaj.CI)

# Watson test to compare activities 
library(circular)
Ptaj.Watson <- watson.two.test(indRecs[indRecs$Species == "Pecari tajacu" & indRecs$pa_type == "strict", "Time_rad"],
                indRecs[indRecs$Species == "Pecari tajacu" & indRecs$pa_type == "multiple-use", "Time_rad"])

# Test Statistic: 0.4482 
# P-value < 0.001 

#### Mgou - grey brocket deer ####
# get Dhats
(DhatsMgou<-overlapEst(indRecs[indRecs$Species == "Mazama gouazoubira" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Mazama gouazoubira" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col
# smoothed bootstrap 
bsMgou.strict<-resample(indRecs[indRecs$Species == "Mazama gouazoubira" & indRecs$pa_type == "strict", 18],10000) 
bsMgou.multi<-resample(indRecs[indRecs$Species == "Mazama gouazoubira" & indRecs$pa_type == "multiple-use", 18],10000)

# bootEst: 
bsMgou<-bootEst(bsMgou.strict,bsMgou.multi,adjust=c(NA,1,NA)) 
Mgou.mean <- mean(bsMgou[,2])
#Convert column with Dhat to vector and get CIs
bsMgouVec<-as.vector(bsMgou[,2])#2 because it's Dhat 4
Mgou.CI <- as.vector(bootCI(DhatsMgou[2],bsMgouVec)['basic0',]) #2 because it's Dhat 4
Mgou.overlap <- rbind(Mgou.mean,Mgou.CI)

# Watson test to compare activities 
Mgou.Watson <- watson.two.test(indRecs[indRecs$Species == "Mazama gouazoubira" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Mazama gouazoubira" & indRecs$pa_type == "multiple-use", "Time_rad"])

#Test Statistic: 0.1237 
#P-value > 0.10 

#### Ctho - crab-eating fox ####
# get Dhats
(DhatsCtho<-overlapEst(indRecs[indRecs$Species == "Cerdocyon thous" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Cerdocyon thous" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col

# smoothed bootstrap for eat PAtype #usually use 10,000 reps:
bsCtho.strict<-resample(indRecs[indRecs$Species == "Cerdocyon thous" & indRecs$pa_type == "strict", 18],10000) 
bsCtho.multi<-resample(indRecs[indRecs$Species == "Cerdocyon thous" & indRecs$pa_type == "multiple-use", 18],10000)

#bootEst: 
bsCtho<-bootEst(bsCtho.strict,bsCtho.multi,adjust=c(NA,1,NA)) 
Ctho.mean <-mean(bsCtho[,2])
#Convert column with Dhat to vector and get CIs
bsCthoVec<-as.vector(bsCtho[,2])#2 because it's Dhat 4
Ctho.CI <- as.vector(bootCI(DhatsCtho[2],bsCthoVec)['basic0',]) #2 because it's Dhat 4
Ctho.overlap <- rbind(Ctho.mean,Ctho.CI)

# Watson test to compare activities 
Ctho.Watson <- watson.two.test(indRecs[indRecs$Species == "Cerdocyon thous" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Cerdocyon thous" & indRecs$pa_type == "multiple-use", "Time_rad"])
# Test Statistic: 0.1541 
# 0.05 < P-value < 0.10 

#### Ltig - oncilla ####
# get Dhats
(DhatsLtig<-overlapEst(indRecs[indRecs$Species == "Leopardus tigrinus" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Leopardus tigrinus" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col

# smoothed bootstrap for eat PAtype #usually use 10,000 reps:
bsLtig.strict<-resample(indRecs[indRecs$Species == "Leopardus tigrinus" & indRecs$pa_type == "strict", 18],10000) 
bsLtig.multi<-resample(indRecs[indRecs$Species == "Leopardus tigrinus" & indRecs$pa_type == "multiple-use", 18],10000)

# bootEst: 
bsLtig<-bootEst(bsLtig.strict,bsLtig.multi,adjust=c(NA,1,NA))
Ltig.mean <-mean(bsLtig[,2])
#Convert column with Dhat to vector and get CIs
bsLtigVec<-as.vector(bsLtig[,2])#2 because it's Dhat 4
#bootCI(DhatsLtig[2],bsLtigVec)['norm0',] #2 because it's Dhat 4
Ltig.CI <- as.vector(bootCI(DhatsLtig[2],bsLtigVec)['basic0',]) #2 because it's Dhat 4
Ltig.overlap <- rbind(Ltig.mean,Ltig.CI)

# Watson test to compare activities 
Ltig.Watson <- watson.two.test(indRecs[indRecs$Species == "Leopardus tigrinus" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Leopardus tigrinus" & indRecs$pa_type == "multiple-use", "Time_rad"])
#Test Statistic: 0.121 
#P-value > 0.10

#### Csem - skunk ####
# get Dhats
(DhatsCsem<-overlapEst(indRecs[indRecs$Species == "Conepatus semistriatus" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Conepatus semistriatus" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col

# smoothed bootstrap for eat PAtype #usually use 10,000 reps:
bsCsem.strict<-resample(indRecs[indRecs$Species == "Conepatus semistriatus" & indRecs$pa_type == "strict", 18],10000) 
bsCsem.multi<-resample(indRecs[indRecs$Species == "Conepatus semistriatus" & indRecs$pa_type == "multiple-use", 18],10000)

#Analyse with bootEst: 
bsCsem<-bootEst(bsCsem.strict,bsCsem.multi,adjust=c(NA,1,NA)) 
Csem.mean <-mean(bsCsem[,2])

#Convert column with Dhat to vector and get CIs
bsCsemVec<-as.vector(bsCsem[,2])#2 because it's Dhat 4
Csem.CI <- as.vector(bootCI(DhatsCsem[2],bsCsemVec)['basic0',]) #2 because it's Dhat 4
Csem.overlap <- rbind(Csem.mean,Csem.CI)

Csem.Watson <- watson.two.test(indRecs[indRecs$Species == "Conepatus semistriatus" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Conepatus semistriatus" & indRecs$pa_type == "multiple-use", "Time_rad"])
#Test Statistic: 0.1566 
#0.05 < P-value < 0.10 

#### Dnov - Dasypus armadillo ####
# get Dhats
(DhatsDnov<-overlapEst(indRecs[indRecs$Species == "Dasypus sp" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Dasypus sp" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col

# smoothed bootstrap for eat PAtype #usually use 10,000 reps:
bsDnov.strict<-resample(indRecs[indRecs$Species == "Dasypus sp" & indRecs$pa_type == "strict", 18],10000) 
bsDnov.multi<-resample(indRecs[indRecs$Species == "Dasypus sp" & indRecs$pa_type == "multiple-use", 18],10000)

#Analyse with bootEst: 
bsDnov<-bootEst(bsDnov.strict,bsDnov.multi,adjust=c(0.8,NA,NA)) 
Dnov.mean <-mean(bsDnov[,1])
#Convert column with Dhat to vector and get CIs
bsDnovVec<-as.vector(bsDnov[,1])#1 because it's Dhat1
Dnov.CI <- as.vector(bootCI(DhatsDnov[1],bsDnovVec)['basic0',]) #1 because it's Dhat 1
Dnov.overlap <- rbind(Dnov.mean,Dnov.CI)

Dnov.Watson <- watson.two.test(indRecs[indRecs$Species == "Dasypus sp" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Dasypus sp" & indRecs$pa_type == "multiple-use", "Time_rad"])
#Test Statistic: 0.2497 
#0.01 < P-value < 0.05

#### Daza ####
# get Dhats
(DhatsDaza<-overlapEst(indRecs[indRecs$Species == "Dasyprocta azarae" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Dasyprocta azarae" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col

#smoothed bootstrap for eat PAtype #usually use 10,000 reps:
bsDaza.strict<-resample(indRecs[indRecs$Species == "Dasyprocta azarae" & indRecs$pa_type == "strict", 18],10000) 
bsDaza.multi<-resample(indRecs[indRecs$Species == "Dasyprocta azarae" & indRecs$pa_type == "multiple-use", 18],10000)

#Analyse with bootEst: 
bsDaza<-bootEst(bsDaza.strict,bsDaza.multi,adjust=c(0.8,NA,NA)) 
Daza.mean <-mean(bsDaza[,1])
#Convert column with Dhat to vector and get CIs
bsDazaVec<-as.vector(bsDaza[,1])#1 because it's Dhat1
Daza.CI <- as.vector(bootCI(DhatsDaza[1],bsDazaVec)['basic0',]) 
Daza.overlap <- rbind(Daza.mean,Daza.CI)

Daza.Watson <- watson.two.test(indRecs[indRecs$Species == "Dasyprocta azarae" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Dasyprocta azarae" & indRecs$pa_type == "multiple-use", "Time_rad"])

#Test Statistic: 0.114 
#P-value > 0.10 

#### Esex ####
# get Dhats
(DhatsEsex<-overlapEst(indRecs[indRecs$Species == "Euphractus sexcinctus" & indRecs$pa_type == "strict", 18], 
                       indRecs[indRecs$Species == "Euphractus sexcinctus" & indRecs$pa_type == "multiple-use", 18]))
#18 is the Time_rad col

#Do smoothed bootstrap for eat PAtype #usually use 10,000 reps:
bsEsex.strict<-resample(indRecs[indRecs$Species == "Euphractus sexcinctus" & indRecs$pa_type == "strict", 18],10000) 
bsEsex.multi<-resample(indRecs[indRecs$Species == "Euphractus sexcinctus" & indRecs$pa_type == "multiple-use", 18],10000)

#Analyse with bootEst: 
bsEsex<-bootEst(bsEsex.strict,bsEsex.multi,adjust=c(0.8,NA,NA)) 
Esex.mean <-mean(bsEsex[,1])
#Convert column with Dhat to vector and get CIs
bsEsexVec<-as.vector(bsEsex[,1])#1 because it's Dhat1
Esex.CI <- as.vector(bootCI(DhatsEsex[1],bsEsexVec)['basic0',]) 
Esex.overlap <- rbind(Esex.mean,Esex.CI)

Esex.Watson <- watson.two.test(indRecs[indRecs$Species == "Euphractus sexcinctus" & indRecs$pa_type == "strict", "Time_rad"],
                               indRecs[indRecs$Species == "Euphractus sexcinctus" & indRecs$pa_type == "multiple-use", "Time_rad"])

#Test Statistic: 0.2394 
#0.01 < P-value < 0.05 

#### joining all overlap results ####
overlapResults <- as.data.frame(rbind(Mgou.overlap,Ptaj.overlap, Dnov.overlap, Esex.overlap,
                        Ltig.overlap,Ctho.overlap,Csem.overlap,Daza.overlap)) # joining mean and CIs for each spp
overlapResults$species <- c("Mgou","Mgou","Ptaj","Ptaj","Dnov","Dnov",
                            "Esex","Esex", "Ltig","Ltig","Ctho","Ctho","Csem","Csem","Daza","Daza")
write.csv(overlapResults,"overlapResults_8spp.csv")

