########################################################################

# Estimating activity level and activity pattern from camera trap data #
# Guilherme Ferreira - 14/July/2020

########################################################################
# Key papers:
# 1) Rowcliffe et al. Methods in Ecology and Evolution 5.11 (2014):1170-1179.
# 2) Ridout & Linkie. Journal of Agricultural, Biological, and Environmental Statistics 14.3 (2009):322-337.

##### Function to transform time #####
# gettime function from Marcus Rowcliffe's Github - activity
#Converts character, POSIXct or POSIXlt time of day data to numeric
#ARGUMENTS
# x: vector of character, POSIXct or POSIXlt time data to convert
# format: used only if x is character, see strptime
# scale: scale on which to return times (see below)
#VALUE
#A vector of numeric times of day in units defined by scale:
# radian: [0,2*pi]
# hours: [0,24]
# proportion: [0,1]
# example: ptime <- gettime(BCItime$date, "%d/%m/%Y %H:%M", "proportion")
# SOLAR TIME FUNCTIONS FROM ROWCLIFFE ####
install.packages("insol")
require(insol)

gettime <- function(x, format="%Y-%m-%d %H:%M:%S", scale=c("radian","hour","proportion")){
  if(class(x)[1]=="character") x <- strptime(x, format, "UTC") else
    if(class(x)[1]=="POSIXct") x <- as.POSIXlt(x) else
      if(class(x)[1]!="POSIXlt") stop("x must be character or POSIXt class")
    scale <- match.arg(scale)
    res <- x$hour + x$min/60 + x$sec/3600
    if(scale=="radian") res <- res*pi/12
    if(scale=="proportion") res <- res/24
    if(all(res==0, na.rm=T)) warning("All times are 0: may be just strptime default?")
    res
}


##### loading camera trap data ####
setwd("Z:/biome_health_project_files/country_files/nepal/tagging_photos")
d <- read.csv("Latest_species_meta.csv") # file with CT image metadata and species tags
names(d)
d2 <- d[,c("site_cam.x", "site_id", "CommonName", "date_fixed", "Time", "Image_nam" )] # keeping only cols needed
head(d2)
d2$protection <- d2$site_id # col indicating National Park or outside NP, 2 categories i.o. 3
d2[d2$protection== "BZ", "protection"] <- "outside" 
head(d2)
d2[d2$protection== "OBZ", "protection"] <- "outside"

d2$Time_rad <- gettime(d2$Time, "%H:%M:%S")  # transforming time in radians
d2$Time_prop <- gettime(d2$Time, "%H:%M:%S", "proportion") # transforming time in proportion
range(d2$Time_rad) # just checking - should be between ca. 0 and 6.28
range(d2$Time_prop) # just checking - should be between 0 and 1

d2[is.na(d2$CommonName), "CommonName"] <- "blank" # replacing NAs for blank on species names
which(is.na(d2$CommonName)) # checking

table(d2$CommonName, d2$protection) # number of images per spp
# wild boar and barking deer have a decent number of images in both protection levels (65/243; 130/171, respectively )

##### Wild Boar activity levels ####
library(activity)
# Wild boar
# just checking if selection is ok
length(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"])
length(d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"])

boar_overall<- fitact(d2[d2$CommonName == "Wild Boar", "Time_rad"],
                      reps=1000,sample="model",show=TRUE)

boar_NP<- fitact(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"],
                    reps=1000,sample="model",show=TRUE)

boar_outside<- fitact(d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"],
                         reps=1000,sample="model",show=TRUE)

boar_overall@act
boar_NP@act
boar_outside@act

# Wald test to compare activity level between management types
compareAct(list(boar_NP, boar_outside))

#### No need to run - testing sample= 'data'; no difference observed  ####
boar_overallD<- fitact(d2[d2$CommonName == "Wild Boar", "Time_rad"],
                      reps=10000,sample="data",show=TRUE) # with 10K reps parameters were very similar as "model"

boar_NPD<- fitact(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"],
                 reps=10000,sample="data",show=TRUE)

boar_outsideD<- fitact(d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"],
                      reps=1000,sample="data",show=TRUE)

boar_overallD@act
boar_NPD@act
boar_outsideD@act
#####
# activity pattern graph
plot(boar_overall, centre= "night",  tline=list(lty=1,lwd=2, col="blue"), 
     cline=list(lty=2,lwd=1.5, col="black"), dline=list(col="dark gray"),main= "Wild boar", cex.main=1)

plot(boar_NP, centre= "night",  tline=list(lty=1,lwd=2, col="blue"), 
     cline=list(lty=2,lwd=1.5, col="black"), dline=list(col="dark gray"),main= "Wild boar_NP", cex.main=1)

plot(boar_outside, centre= "night",  tline=list(lty=1,lwd=2, col="blue"), 
     cline=list(lty=2,lwd=1.5, col="black"), dline=list(col="dark gray"),main= "Wild boar_outside", cex.main=1)

#### Wild Boar overlap estimates ####
library(overlap)
# estimating Dhats - overlap metric
(DhatsWB<-overlapEst(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"],
                     d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"]))

#Do smoothed bootstrap for each protection level 
bsWB.NP<-resample(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"],10000) 
bsWB.outside<-resample(d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"],10000)

#Analyse with bootEst: 
# Ridout and Linkie (2009) recommend using adjust=0.8 to estimate Dhat1, adjust=1 for Dhat 4, an adjust=4 for Dhat 5.
# When smaller sample <50, Dhat1 performed best, while Dhat4 was better when smaller sample >75 (Meredith & Ridout 2018; overview of overlap)
bsWB<-bootEst(bsWB.NP,bsWB.outside,adjust=c(NA,1,NA)) #Dhat 4
bsWB.mean <-as.vector(colMeans(bsWB))
bsWB.mean <-bsWB.mean[2]
#Convert column with Dhat4 to vector and get CIs
bsWBVec<-as.vector(bsWB[,2])#2 because it's Dhat 4
#bootCI(DhatsPtaj[2],bsPtajVec)['norm0',] #2 because it's Dhat 4
WB.CI <- as.vector(bootCI(DhatsWB[2],bsWBVec)['basic0',]) #2 because it's Dhat 4
WB.overlap <- cbind(bsWB.mean,WB.CI)

#plotting overlapping activities
overlapPlot(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"],
            d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"], 
            xcenter= "midnight",
            rug=TRUE,
            main="Wild Boar")
legend('topright',c("NP","outside"),lty=c(1,2),col=c("black","blue"),bty='n')

#Watson two Test to compare whether 2 circular data came from the same distribution
# i.e. is WB activity patterns different in the two protection levels?
# Refs:Jammalamadaka & SenGupta 2001; Oliveira-Santos et al. 2013
library(circular)
watson.two.test(d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"],
                d2[d2$CommonName == "Wild Boar" & d2$protection == "outside" , "Time_rad"])


##### Barking Deer activity level####
length(d2[d2$CommonName == "Barking Deer" & d2$protection == "NP" , "Time_rad"])
length(d2[d2$CommonName == "Barking Deer" & d2$protection == "outside" , "Time_rad"])

deer_overall<- fitact(d2[d2$CommonName == "Barking Deer", "Time_rad"],
                      reps=1000,sample="model",show=TRUE)

deer_NP<- fitact(d2[d2$CommonName == "Barking Deer" & d2$protection == "NP" , "Time_rad"],
                 reps=1000,sample="model",show=TRUE)

deer_outside<- fitact(d2[d2$CommonName == "Barking Deer" & d2$protection == "outside" , "Time_rad"],
                      reps=1000,sample="model",show=TRUE)

deer_overall@act
deer_NP@act
deer_outside@act

# Wald test to compare activity level between management types
compareAct(list(deer_NP, deer_outside))

# activity pattern graph
plot(deer_overall, centre= "night",  tline=list(lty=1,lwd=2, col="blue"), 
     cline=list(lty=2,lwd=1.5, col="black"), dline=list(col="dark gray"),main= "Barking Deer", cex.main=1)

plot(deer_NP, centre= "night",  tline=list(lty=1,lwd=2, col="blue"), 
     cline=list(lty=2,lwd=1.5, col="black"), dline=list(col="dark gray"),main= "Barking Deer_NP", cex.main=1)

plot(deer_outside, centre= "night",  tline=list(lty=1,lwd=2, col="blue"), 
     cline=list(lty=2,lwd=1.5, col="black"), dline=list(col="dark gray"),main= "Barking Deer_outside", cex.main=1)

#plotting overlap graph
overlapPlot(d2[d2$CommonName == "Barking Deer" & d2$protection == "NP" , "Time_rad"],
            d2[d2$CommonName == "Barking Deer" & d2$protection == "outside" , "Time_rad"], 
            xcenter= "midnight",
            rug=TRUE,
            #abline(v=c(1.57,4.71),lty=3),
            main="Barking Deer")
legend('topright',c("NP","outside"),lty=c(1,2),col=c("black","blue"),bty='n')

#Watson two Test to compare whether 2 circular data came from the same distribution
library(circular)
watson.two.test(d2[d2$CommonName == "Barking Deer" & d2$protection == "NP" , "Time_rad"],
                d2[d2$CommonName == "Barking Deer" & d2$protection == "outside" , "Time_rad"])


#### other plots
overlapPlot(d2[d2$CommonName == "Tiger" & d2$protection == "NP" , "Time_rad"],
            d2[d2$CommonName == "Wild Boar" & d2$protection == "NP" , "Time_rad"], 
            xcenter= "midnight",
            rug=TRUE,
            #abline(v=c(1.57,4.71),lty=3),
            main="Tiger vs Wild Boar")
legend('topright',c("Tg","WB"),lty=c(1,2),col=c("black","blue"),bty='n')


overlapPlot(d2[d2$CommonName == "Nilgai" & d2$protection == "outside" , "Time_rad"],
            d2[d2$CommonName == "Chital" & d2$protection == "outside" , "Time_rad"], 
            #xcenter= "midnight",
            rug=TRUE,
            #abline(v=c(1.57,4.71),lty=3),
            main="Niglai vs Chital")
legend('topright',c("Ng","Ch"),lty=c(1,2),col=c("black","blue"),bty='n')