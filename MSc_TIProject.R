####MSc Project: Techincal Interactions Data work
####Author: Braden Charles DeMattei
####Created: April 20th, 2020
####Last Modified: July 26th, 2020
#############################################
####READ ME: Shorthands and Identifier Guide######

##Identifying specifc Area+Gear Combinations
#Example: For an object that relates to gear type TR1 which was used in area 27.4
  ##Object name will have '1274' in name to identify it
#A '2' is TR2 and an 'O' is TRother

#If something is named '.melt' it is an object that is the result of a matrix or similar being run through the 'melt()' function

#References to six species (e.g., the spp6 modifier) refers to when six species were 
##being analyzed. One species was removed due to lack of observations. The modifier was kept 
###for consistency.

###Libraries Used######
library(vegan)
library(FactoMineR)
library(factoextra)
library(pvclust)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(grid)
library(abind)
library(rgl)
library(plyr)
library(magrittr)
library(zoo)
library(dplyr)
library(tidyr)
library(DataCombine)
library(bi5009) #custom package for my usage
library(Rmisc)
library(gridExtra)
library(cowplot)
library(ggtext)

############################################
getwd() #making sure working directory is set correctly
load("TechIntTripData.rdata") #loading in data set
trip <- tripdat #Making sure I don't alter original data by creating new data set

head(trip) #checking variables
str(trip)


trip2 <- trip #Redundant data set as I'm going to be adding to it
trip3$spp_id <- as.numeric(as.factor(trip2$spp)) #adding column for numerical codes that correspond to each species
trip3$gear_id <- as.numeric(as.factor(trip2$foCatNat)) #numerical codes that correspond to each gear type
trip3$area_id <- as.numeric(as.factor(trip2$area)) #numerical codes that correspond to each ICES Area
head(trip3)


unique(trip3$area_id) #8 areas
unique(trip3$spp_id) #243 species
unique(trip3$quarter)#4 quarters
unique(trip3$gear_id) 
#5 gear types

spp_guide <- data.frame() #creating spp_id reference guide
spp_guide[1:243,1] <- unique(trip3$spp) #adding one of each scientific names
spp_guide[1:243,2] <- unique(trip3$spp_id) #adding corresponding numerical code
View(spp_guide)
write.csv(spp_guide, "Species_Guide.csv")

##############ALL STRATA####################

ti_as <- data.frame() #Dataframe for the all strata loop that will store the technical interaction strength percentages

for(i in 1:243){#This will loop through each of the 243 species
   sp1 <- subset(trip3, spp_id == i) #Creating dataframe of target species i, with i referring to numerical species code
     for(j in 1:243){#This will be a secondary loop so that sp1 can be compared to each of the other species in the data set
       x1 <- data.frame()#Creating dataframe to store weights for when bycatch species is present and 0s when they are not
       sp2 <- subset(trip3, spp_id == j) #Creating dataframe of bycatch species j, with j referring to numerical species code
       if(i == j){ti_as[i,j] <- 100} #This will add a 100% TI strength for when the species i and j are the same
       else{for(k in 1:nrow(sp1)){ #For loop that will cycle through each observation of the target species
         if((sp1$area[k] %in% sp2$area)& #Conditional statement: Observation k of species i must have same gear, area, and quarter as any row of species j 
            (sp1$quarter[k] %in% sp2$quarter)&(sp1$foCatNat[k] %in% sp2$foCatNat)){
          x1[k,1] <- sp1$wt[k]#If conditional statement met, weight of observation k will be stored in dataframe
         } else{x1[k,1] <- 0}}#If conditional statement not met, a 0 will be placed
     ti_as[i,j] <- ((sum(x1)/sum(sp1$wt))*100)}} #The sum of the stored weights and 0s is divided by the total weight of target species i
}                                                  ###then multiplied by 100 to create percentage. 
                                                          ###This is the Technical Interaction Strength
                                                   #Row i of TI dataframe corresponds to target species i
                                                   #Column j of TI dataframe corresponds to bycatch species j

write.csv(as_mat, "TI_AS.csv") #writing CSV file for All Strata

anyNA(ti_as) #checking for errors
as_mat <- data.matrix(ti_as) #creating matrix for heatmap
ras_mat <- apply(as_mat, 2, rev) #rotating the matrix so the heatmap corresponds with the original matrix's rows and columns
image(ras_mat, axes = F, main = "All Species and All Strata")


################JUST GEAR########################

ti_gear <- data.frame() #dataframe for the all gear loop
  #Refer to annotations above for specific aspects. 
    ##Looking at TI in entire dataset, 
      ##where target species is caught in just the same gear type as bycatch species
for(i in 1:243){
  spg1 <- subset(trip3, spp_id == i)
  for(j in 1:243){
    g1 <- data.frame()
    spg2 <- subset(trip3, spp_id == j)
    if(i == j){ti_gear[i,j] <- 100}
    else{for(k in 1:nrow(spg1)){
      if(spg1$foCatNat[k] %in% spg2$foCatNat){
        g1[k,1] <- spg1$wt[k]}
        else{g1[k,1] <- 0}}
      ti_gear[i,j] <- (sum(g1)/sum(spg1$wt))*100}} 
}      
      
      
anyNA(ti_gear) #checking for errors
gear_mat <- data.matrix(ti_gear) #creating matrix for heat map
rgear_mat <- apply(gear_mat, 2, rev)
image(rgear_mat, axes = F, main = "All Species, Just Gear") #creating heat map
write.csv(gear_mat, "TI_AllGear.csv") #creating CSV file

############JUST AREA##########################
ti_area <- data.frame()
  #Looking at TI in entire dataset, 
    ##where target species was caught in just the same area as bycatch species
for(i in 1:243){
  spa1 <- subset(trip3, spp_id == i)
  for(j in 1:243){
    a1 <- data.frame()
    spa2 <- subset(trip3, spp_id == j)
    if(i == j){ti_area[i,j] <- 100}
    else{for(k in 1:nrow(spa1)){
      if(spa1$area[k] %in% spa2$area){
        a1[k,1] <- spa1$wt[k]}
        else{a1[k,1] <- 0}}
     ti_area[i,j] <- (sum(a1)/sum(spa1$wt))*100}} 
}   

anyNA(ti_area)
area_mat <- data.matrix(ti_area)
rarea_mat <- apply(area_mat, 2, rev)
image(rarea_mat, axes = F, main = "All Species, Just Area")
write.csv(area_mat, "TI_AllArea.csv")
#############JUST QUARTER####################
ti_qtr <- data.frame()
  #Looking at TI in entire dataset, 
    ##where target species is caught in just the same quarter as bycatch species
for(i in 1:243){
  spq1 <- subset(trip3, spp_id == i)
  for(j in 1:243){
    q1 <- data.frame()
    spq2 <- subset(trip3, spp_id == j)
    if(i == j){ti_qtr[i,j] <- 100}
    else{for(k in 1:nrow(spq1)){
      if(spq1$quarter[k] %in% spq2$quarter){
        q1[k,1] <- spq1$wt[k]}
      else{q1[k,1] <- 0}}
      ti_qtr[i,j] <- (sum(q1)/sum(spq1$wt))*100}} 
}   

anyNA(ti_qtr)
qtr_mat <- data.matrix(ti_qtr)
rqtr_mat <- apply(qtr_mat, 2, rev)
image(rqtr_mat, axes = F, main = "All Species, Just Quarter")
write.csv(qtr_mat, "TI_AllQtr.csv")

##########Specific Gear Strata###############

unique(trip3$gear_id)

tripgear1 <- subset(trip3, gear_id == 1) #FPO, just Nephrops so will be skipped

#Creating Unique Species Code for TM Gear (2)
tripgear2 <- subset(trip3, gear_id == 2) #TM
tripgear2$spp2_id <- as.numeric(as.factor(tripgear2$spp)) #creating unique numerical codes for this dataframe
  ##allows for cleaner matrix image later. Refer to guide produced
unique(tripgear2$spp2_id) #2 species
tg2_guide <- data.frame() #creating spp_id reference guide
tg2_guide[1:2,1] <- unique(tripgear2$spp) #adding one of each scientific names
tg2_guide[1:2,2] <- unique(tripgear2$spp_id) #adding corresponding original numerical code
tg2_guide[1:2,3] <- unique(tripgear2$spp2_id)
View(tg2_guide)
write.csv(tg2_guide, "Trip Gear 2 Guide.csv")

#Creating USC for tr1 (3)
tripgear3 <- subset(trip3, gear_id == 3) #TR1
tripgear3$spp2_id <- as.numeric(as.factor(tripgear3$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripgear3$spp2_id) # 220 species
tg3_guide <- data.frame() #creating spp_id reference guide
tg3_guide[1:220,1] <- unique(tripgear3$spp) #adding one of each scientific names
tg3_guide[1:220,2] <- unique(tripgear3$spp_id) #adding corresponding original numerical code
tg3_guide[1:220,3] <- unique(tripgear3$spp2_id)
View(tg3_guide)
write.csv(tg3_guide, "Trip Gear 3 Guide.csv")

#USC for TR2 (4)
tripgear4 <- subset(trip3, gear_id == 4) #TR2
tripgear4$spp2_id <- as.numeric(as.factor(tripgear4$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripgear4$spp2_id) # 141 species
tg4_guide <- data.frame() #creating spp_id reference guide
tg4_guide[1:141,1] <- unique(tripgear4$spp) #adding one of each scientific names
tg4_guide[1:141,2] <- unique(tripgear4$spp_id) #adding corresponding original numerical code
tg4_guide[1:141,3] <- unique(tripgear4$spp2_id)
View(tg4_guide)
write.csv(tg4_guide, "Trip Gear 4 Guide.csv")

#USC for TRother (5)
tripgear5 <- subset(trip3, gear_id == 5) #TRother
tripgear5$spp2_id <- as.numeric(as.factor(tripgear5$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripgear5$spp2_id) # 65 species
tg5_guide <- data.frame() #creating spp_id reference guide
tg5_guide[1:65,1] <- unique(tripgear5$spp) #adding one of each scientific names
tg5_guide[1:65,2] <- unique(tripgear5$spp_id) #adding corresponding original numerical code
tg5_guide[1:65,3] <- unique(tripgear5$spp2_id)
View(tg5_guide)
write.csv(tg5_guide, "Trip Gear 5 Guide.csv")


spp_g2 <- unique(tripgear2$spp2_id)#creating list for 'for' loop so that the numbers that are cycled through are just those relevant to the TM gear type
spp_g2 <- sort(spp_g2)
spp_g3 <- unique(tripgear3$spp2_id)#Just new species numerical codes for the TR1 gear type
spp_g3 <- sort(spp_g3)
spp_g4 <- unique(tripgear4$spp2_id)#Just new species numerical codes for the TR2 gear type
spp_g4 <- sort(spp_g4)
spp_g5 <- unique(tripgear5$spp2_id)#Just new species numerical codes for the TRother gear type
spp_g5 <- sort(spp_g5)

##Disregarding gear ID 1 (FPO) as only Nephrops are caught in this gear type

##Gear ID 2

ti_g2 <- data.frame() 

  #Just looking at TI in gear type TM
for(i in spp_g2){
  sp1 <- subset(tripgear2, spp2_id == i)
  for(j in spp_g2){
    g <- data.frame()
    sp2 <- subset(tripgear2, spp2_id ==j)
    if(i == j){ti_g2[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$area[k] %in% sp2$area)){
        g[k,1] <- sp1$wt[k]}
      else{g[k,1] <- 0}}
      ti_g2[i,j] <- (sum(g)/sum(sp1$wt))*100}} 
}   

anyNA(ti_g2)
g2_mat <- data.matrix(ti_g2)
rg2_mat <- apply(g2_mat, 2, rev)
image(rg2_mat, axes = F, main = "Species Caught in Gear Type TM")

##Gear ID 3
ti_g3 <- data.frame()
  #Looking at TI in just the TR1 gear type dataframe

for(i in spp_g3){
  sp1 <- subset(tripgear3, spp2_id == i)
  for(j in spp_g3){
    g <- data.frame()
    sp2 <- subset(tripgear3, spp2_id ==j)
    if(i == j){ti_g3[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$area[k] %in% sp2$area)){
        g[k,1] <- sp1$wt[k]}
      else{g[k,1] <- 0}}
      ti_g3[i,j] <- (sum(g)/sum(sp1$wt))*100}} 
} 

 
g3_mat <- data.matrix(ti_g3)
rg3_mat <- apply(g3_mat, 2, rev)
image(rg3_mat, axes = F, main = "Species Caught in Gear Type TR1")

##Gear ID 4
ti_g4 <- data.frame()

  #Just looking at TI in the TR2 gear type dataframe

for(i in spp_g4){
  sp1 <- subset(tripgear4, spp2_id == i)
  for(j in spp_g4){
    g <- data.frame()
    sp2 <- subset(tripgear4, spp2_id ==j)
    if(i == j){ti_g4[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$area[k] %in% sp2$area)){
        g[k,1] <- sp1$wt[k]}
      else{g[k,1] <- 0}}
      ti_g4[i,j] <- (sum(g)/sum(sp1$wt))*100}} 
} 

g4_mat <- data.matrix(ti_g4)
rg4_mat <- apply(g4_mat, 2, rev)
image(rg4_mat, axes = F, main = "Species Caught in Gear Type TR2")

##Gear ID 5
ti_g5 <- data.frame()

  #Just looking at TI in the TRother gear type dataframe

for(i in spp_g5){
  sp1 <- subset(tripgear5, spp2_id == i)
  for(j in spp_g5){
    g <- data.frame()
    sp2 <- subset(tripgear5, spp2_id ==j)
    if(i == j){ti_g5[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$area[k] %in% sp2$area)){
        g[k,1] <- sp1$wt[k]}
      else{g[k,1] <- 0}}
      ti_g5[i,j] <- (sum(g)/sum(sp1$wt))*100}} 
} 

g5_mat <- data.matrix(ti_g5)
rg5_mat <- apply(g5_mat, 2, rev)
image(rg5_mat, axes = F, main = "Species Caught in Gear Type TRother")

##############Specific Area Strata###############
unique(trip3$area_id)

#Only Areas 27.6.a, 27.6.b, and 27.4 are of interest
#Discard
# triparea1 <- subset(trip3, area_id == 1) #Area 27.2.a
# triparea1$spp2_id <- as.numeric(as.factor(triparea1$spp)) #creating unique numerical codes for this dataframe
# ##allows for cleaner matrix image later. Refer to guide produced
# unique(triparea1$spp2_id) # 200 species
# ta1_guide <- data.frame() #creating spp_id reference guide
# ta1_guide[1:200,1] <- unique(triparea1$spp) #adding one of each scientific names
# ta1_guide[1:200,2] <- unique(triparea1$spp_id) #adding corresponding original numerical code
# ta1_guide[1:200,3] <- unique(triparea1$spp2_id)
# View(ta1_guide)
# write.csv(ta1_guide, "Trip Area 1 Guide.csv")


triparea2 <- subset(trip3, area_id == 2) #Area 27.4
triparea2$spp2_id <- as.numeric(as.factor(triparea2$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(triparea2$spp2_id) # 173 species
ta2_guide <- data.frame() #creating spp_id reference guide
ta2_guide[1:173,1] <- unique(triparea2$spp) #adding one of each scientific names
ta2_guide[1:173,2] <- unique(triparea2$spp_id) #adding corresponding original numerical code
ta2_guide[1:173,3] <- unique(triparea2$spp2_id)
View(ta2_guide)
write.csv(ta2_guide, "Trip Area 2 Guide.csv")
unique(triparea2$gear_id) #TM/TR1/TR2/TRother
table(triparea2$gear_id) #TM: 38, TR1: 15623, TR2: 4212, TRother 135

#Discard
# triparea3 <- subset(trip3, area_id == 3) #Area 27.7.b, just one observation

#Discard
# triparea4 <- subset(trip3, area_id == 4) #Area 27.5.b.1
# triparea4$spp2_id <- as.numeric(as.factor(triparea4$spp)) #creating unique numerical codes for this dataframe
# ##allows for cleaner matrix image later. Refer to guide produced
# unique(triparea4$spp2_id) # 101 species
# ta4_guide <- data.frame() #creating spp_id reference guide
# ta4_guide[1:101,1] <- unique(triparea4$spp) #adding one of each scientific names
# ta4_guide[1:101,2] <- unique(triparea4$spp_id) #adding corresponding original numerical code
# ta4_guide[1:101,3] <- unique(triparea4$spp2_id)
# View(ta4_guide)
# write.csv(ta4_guide, "Trip Area 4 Guide.csv")

#Discard
# unique(triparea5$area)
# triparea5 <- subset(trip3, area_id == 5) #Area 27.5.b.2
# triparea5 <- transform(triparea5, spp2_id=match(spp, unique(spp))) #creating unique numerical codes for this dataframe
# ##allows for cleaner matrix image later. Refer to guide produced
# unique(triparea5$spp2_id) # 57 species
# ta5_guide <- data.frame() #creating spp_id reference guide
# ta5_guide[1:57,1] <- unique(triparea5$spp) #adding one of each scientific names
# ta5_guide[1:57,2] <- unique(triparea5$spp_id) #adding corresponding original numerical code
# ta5_guide[1:57,3] <- unique(triparea5$spp2_id)
# View(ta5_guide)
# write.csv(ta5_guide, "Trip Area 5 Guide.csv")
# 
# 


triparea6 <- subset(trip3, area_id == 6) #Area 27.6.a
triparea6$spp2_id <- as.numeric(as.factor(triparea6$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
length(unique(triparea6$spp2_id)) # 200 species
ta6_guide <- data.frame() #creating spp_id reference guide
ta6_guide[1:200,1] <- unique(triparea6$spp) #adding one of each scientific names
ta6_guide[1:200,2] <- unique(triparea6$spp_id) #adding corresponding original numerical code
ta6_guide[1:200,3] <- unique(triparea6$spp2_id)
View(ta6_guide)
write.csv(ta6_guide, "Trip Area 6 Guide.csv")
unique(triparea6$gear_id) #All
table(triparea6$gear_id) #TR1 and 2 >3000, FPO/TM/


triparea7 <- subset(trip3, area_id == 7) #Area 27.6.b
triparea7$spp2_id <- as.numeric(as.factor(triparea7$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
length(unique(triparea7$spp2_id)) # 101 species
ta7_guide <- data.frame() #creating spp_id reference guide
ta7_guide[1:101,1] <- unique(triparea7$spp) #adding one of each scientific names
ta7_guide[1:101,2] <- unique(triparea7$spp_id) #adding corresponding original numerical code
ta7_guide[1:101,3] <- unique(triparea7$spp2_id)
View(ta7_guide)
write.csv(ta7_guide, "Trip Area 7 Guide.csv")
unique(triparea7$gear_id) #only TR1 and TRother
table(triparea7$gear_id) #TR1: 875, TRother:193
nrow(triparea7)

#Discard
# triparea8 <- subset(trip3, area_id == 8) #Area 27.5.b.2, all in same quarter with same gear
# triparea8 <- transform(triparea8, spp2_id=match(spp, unique(spp))) #creating unique numerical codes for this dataframe
# ##allows for cleaner matrix image later. Refer to guide produced
# unique(triparea8$spp2_id) # 5 species
# ta8_guide <- data.frame() #creating spp_id reference guide
# ta8_guide[1:5,1] <- unique(triparea8$spp) #adding one of each scientific names
# ta8_guide[1:5,2] <- unique(triparea8$spp_id) #adding corresponding original numerical code
# ta8_guide[1:5,3] <- unique(triparea8$spp2_id)
# View(ta8_guide)
# write.csv(ta8_guide, "Trip Area 8 Guide.csv")


# spp_a1 <- unique(triparea1$spp2_id)
# spp_a1 <- sort(spp_a1)
spp_a2 <- unique(triparea2$spp2_id)#Species numerical codes in Area 27.4 (2)
spp_a2 <- sort(spp_a2)
# spp_a4 <- unique(triparea4$spp2_id)
# spp_a4 <- sort(spp_a4)
# spp_a5 <- unique(triparea5$spp2_id)#Species numerical codes in Area 27.5.b (5)
spp_a6 <- unique(triparea6$spp2_id)#Species numerical codes in Area 27.6.a (6)
spp_a6 <- sort(spp_a6)
spp_a7 <- unique(triparea7$spp2_id)#Species numerical codes in Area 27.6.b(7)
spp_a7 <- sort(spp_a7)
# spp_a8 <- unique(triparea8$spp2_id)#Species numerical codes in Area 25.5.b.2 (8)

#Discard
##Area (1)
# ti_a1 <- data.frame()
#   #Just looking at TI in Area 1 dataframe, 
#     ##where target species is caught in same gear type and quarter as bycatch species
# for(i in spp_a1){
#   sp1 <- subset(triparea1, spp2_id == i)
#   for(j in spp_a1){
#     a <- data.frame()
#     sp2 <- subset(triparea1, spp2_id == j)
#     if(i == j){ti_a1[i,j] <- 100}
#     else{for(k in 1:nrow(sp1)){
#       if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
#         a[k,1] <- sp1$wt[k]}
#       else{a[k,1] <- 0}}
#       ti_a1[i,j] <- (sum(a)/sum(sp1$wt))*100}} 
# }   
# a1_mat <- data.matrix(ti_a1)
# ra1_mat <- apply(a1_mat, 2, rev)
# image(ra1_mat, axes = F, main = "Species Caught in Area 27.6.a")
table(triparea6$foCatNat)

#Area 27.4 (2)
ti_a2 <- data.frame()
  #Just looking at TI in Area 27.4 dataframe
for(i in spp_a2){
  sp1 <- subset(triparea2, spp2_id == i)
  for(j in spp_a2){
    a <- data.frame()
    sp2 <- subset(triparea2, spp2_id == j)
    if(i == j){ti_a2[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
        a[k,1] <- sp1$wt[k]}
      else{a[k,1] <- 0}}
      ti_a2[i,j] <- (sum(a)/sum(sp1$wt))*100}} 
}   
a2_mat <- data.matrix(ti_a2)
ra2_mat <- apply(a2_mat, 2, rev)
image(ra2_mat, axes = F, main = "Species Caught in Area 27.4")

#Disregarded Area 27.7.b as there is nothing to compare in the dataframe with only one observation

#Discard
#Area (4)
# ti_a4 <- data.frame()
#   #Just looking at TI in Area 27.6.b dataframe
# for(i in spp_a4){
#   sp1 <- subset(triparea4, spp2_id == i)
#   for(j in spp_a4){
#     a <- data.frame()
#     sp2 <- subset(triparea4, spp2_id == j)
#     if(i == j){ti_a4[i,j] <- 100}
#     else{for(k in 1:nrow(sp1)){
#       if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
#         a[k,1] <- sp1$wt[k]}
#       else{a[k,1] <- 0}}
#       ti_a4[i,j] <- (sum(a)/sum(sp1$wt))*100}} 
# }   
# a4_mat <- data.matrix(ti_a4)
# ra4_mat <- apply(a4_mat, 2, rev)
# image(ra4_mat, axes = F, main = "Species Caught in Area 27.6.b")

#Discard
# #Area  (5)
# ti_a5 <- data.frame()
#   #Just looking at TI in Area 27.5.b dataframe
# for(i in spp_a5){
#   sp1 <- subset(triparea5, spp2_id == i)
#   for(j in spp_a5){
#     a <- data.frame()
#     sp2 <- subset(triparea5, spp2_id == j)
#     if(i == j){ti_a5[i,j] <- 100}
#     else{for(k in 1:nrow(sp1)){
#       if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
#         a[k,1] <- sp1$wt[k]}
#       else{a[k,1] <- 0}}
#       ti_a5[i,j] <- (sum(a)/sum(sp1$wt))*100}} 
# }   
# a5_mat <- data.matrix(ti_a5)
# ra5_mat <- apply(a5_mat, 2, rev)
# image(ra5_mat, axes = F, main = "Species Caught in Area 27.5.b")
# 

#Area 27.6.a (6)
ti_a6 <- data.frame()
  #Just looking at TI in Area 27.6.a dataframe
for(i in spp_a6){
  sp1 <- subset(triparea6, spp2_id == i)
  for(j in spp_a6){
    a <- data.frame()
    sp2 <- subset(triparea6, spp2_id == j)
    if(i == j){ti_a6[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
        a[k,1] <- sp1$wt[k]}
      else{a[k,1] <- 0}}
      ti_a6[i,j] <- (sum(a)/sum(sp1$wt))*100}}
}
a6_mat <- data.matrix(ti_a6)
ra6_mat <- apply(a6_mat, 2, rev)
image(ra6_mat, axes = F, main = "Species Caught in Area 27.6.a")

#Area 27.6.b (7) All caught in same gear in same quarter
ti_a7 <- data.frame()
  #Just looking at TI in Area 27.6.b dataframe
for(i in spp_a7){
  sp1 <- subset(triparea7, spp2_id == i)
  for(j in spp_a7){
    a <- data.frame()
    sp2 <- subset(triparea7, spp2_id == j)
    if(i == j){ti_a7[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
        a[k,1] <- sp1$wt[k]}
      else{a[k,1] <- 0}}
      ti_a7[i,j] <- (sum(a)/sum(sp1$wt))*100}}
}
a7_mat <- data.matrix(ti_a7)
ra7_mat <- apply(a7_mat, 2, rev)
image(ra7_mat, axes = F, main = "Species Caught in Area 27.6.b")

# #Area (8) All caught in same gear in same quarter
# ti_a8 <- data.frame()
#   #Just looking at TI in Area 25.5.b.2 dataframe
# for(i in spp_a8){
#   sp1 <- subset(triparea8, spp2_id == i)
#   for(j in spp_a8){
#     a <- data.frame()
#     sp2 <- subset(triparea8, spp2_id == j)
#     if(i == j){ti_a8[i,j] <- 100}
#     else{for(k in 1:nrow(sp1)){
#       if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$quarter[k] %in% sp2$quarter)){
#         a[k,1] <- sp1$wt[k]}
#       else{a[k,1] <- 0}}
#       ti_a8[i,j] <- (sum(a)/sum(sp1$wt))*100}} 
# }   
# a8_mat <- data.matrix(ti_a8)
# ra8_mat <- apply(a8_mat, 2, rev)
# image(ra8_mat, axes = F, main = "Species Caught in Area 25.5.b.2")

##############Specific Quarter Strata############

tripqtr1 <- subset(trip3, quarter == 1) #Quarter 1
tripqtr1$spp2_id <- as.numeric(as.factor(tripqtr1$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripqtr1$spp2_id) # 175 species
tq1_guide <- data.frame() #creating spp_id reference guide
tq1_guide[1:175,1] <- unique(tripqtr1$spp) #adding one of each scientific names
tq1_guide[1:175,2] <- unique(tripqtr1$spp_id) #adding corresponding original numerical code
tq1_guide[1:175,3] <- unique(tripqtr1$spp2_id)
View(tq1_guide)
write.csv(tq1_guide, "Trip Quarter 1 Guide.csv")


tripqtr2 <- subset(trip3, quarter == 2) #Quarter 2
tripqtr2$spp2_id <- as.numeric(as.factor(tripqtr2$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripqtr2$spp2_id) # 192 species
tq2_guide <- data.frame() #creating spp_id reference guide
tq2_guide[1:192,1] <- unique(tripqtr2$spp) #adding one of each scientific names
tq2_guide[1:192,2] <- unique(tripqtr2$spp_id) #adding corresponding original numerical code
tq2_guide[1:192,3] <- unique(tripqtr2$spp2_id)
View(tq2_guide)
write.csv(tq2_guide, "Trip Quarter 2 Guide.csv")


tripqtr3 <- subset(trip3, quarter == 3) #Quarter 3
tripqtr3$spp2_id <- as.numeric(as.factor(tripqtr3$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripqtr3$spp2_id) # 162 species
tq3_guide <- data.frame() #creating spp_id reference guide
tq3_guide[1:162,1] <- unique(tripqtr3$spp) #adding one of each scientific names
tq3_guide[1:162,2] <- unique(tripqtr3$spp_id) #adding corresponding original numerical code
tq3_guide[1:162,3] <- unique(tripqtr3$spp2_id)
View(tq3_guide)
write.csv(tq3_guide, "Trip Quarter 3 Guide.csv")


tripqtr4 <- subset(trip3, quarter == 4) #Quarter 4
tripqtr4$spp2_id <- as.numeric(as.factor(tripqtr4$spp)) #creating unique numerical codes for this dataframe
##allows for cleaner matrix image later. Refer to guide produced
unique(tripqtr4$spp2_id) # 137 species
tq4_guide <- data.frame() #creating spp_id reference guide
tq4_guide[1:137,1] <- unique(tripqtr4$spp) #adding one of each scientific names
tq4_guide[1:137,2] <- unique(tripqtr4$spp_id) #adding corresponding original numerical code
tq4_guide[1:137,3] <- unique(tripqtr4$spp2_id)
View(tq4_guide)
write.csv(tq4_guide, "Trip Quarter 4 Guide.csv")


spp_qtr1 <- unique(tripqtr1$spp2_id)#Species numerical codes for Quarter 1
spp_qtr1 <- sort(spp_qtr1)
spp_qtr2 <- unique(tripqtr2$spp2_id) #Species numerical codes for Quarter 2
spp_qtr2 <- sort(spp_qtr2)
spp_qtr3 <- unique(tripqtr3$spp2_id) #Species numerical codes for Quarter 3
spp_qtr3 <- sort(spp_qtr3)
spp_qtr4 <- unique(tripqtr4$spp2_id) #Species numerical codes for Quarter 4
spp_qtr4 <- sort(spp_qtr4)

##Quarter 1
ti_q1 <- data.frame()
  #Just looking at TI in Quarter 1 dataframe, 
    ##Where target species is caught in same gear type and area as bycatch
for(i in spp_qtr1){
  sp1 <- subset(tripqtr1, spp2_id == i)
  for(j in spp_qtr1){
    q <- data.frame()
    sp2 <- subset(tripqtr1, spp2_id == j)
    if(i == j){ti_q1[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$area[k] %in% sp2$area)){
        q[k,1] <- sp1$wt[k]}
      else{q[k,1] <- 0}}
      ti_q1[i,j] <- (sum(q)/sum(sp1$wt))*100}} 
}   
q1_mat <- data.matrix(ti_q1)
rq1_mat <- apply(q1_mat, 2, rev)
image(rq1_mat, axes = F, main = "Species Caught in Quarter 1")

##Quarter 2
ti_q2 <- data.frame()
  #Just looking at TI in Quarter 2 dataframe
for(i in spp_qtr2){
  sp1 <- subset(tripqtr2, spp2_id == i)
  for(j in spp_qtr2){
    q <- data.frame()
    sp2 <- subset(tripqtr2, spp2_id == j)
    if(i == j){ti_q2[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$area[k] %in% sp2$area)){
        q[k,1] <- sp1$wt[k]}
      else{q[k,1] <- 0}}
      ti_q2[i,j] <- (sum(q)/sum(sp1$wt))*100}} 
}   
q2_mat <- data.matrix(ti_q2)
rq2_mat <- apply(q2_mat, 2, rev)
image(rq2_mat, axes = F, main = "Species Caught in Quarter 2")

##Quarter 3
ti_q3 <- data.frame()
  #Just looking at TI in Quarter 3 dataframe
for(i in spp_qtr3){
  sp1 <- subset(tripqtr3, spp2_id == i)
  for(j in spp_qtr3){
    q <- data.frame()
    sp2 <- subset(tripqtr3, spp2_id == j)
    if(i == j){ti_q3[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$area[k] %in% sp2$area)){
        q[k,1] <- sp1$wt[k]}
      else{q[k,1] <- 0}}
      ti_q3[i,j] <- (sum(q)/sum(sp1$wt))*100}} 
}   
q3_mat <- data.matrix(ti_q3)
rq3_mat <- apply(q3_mat, 2, rev)
image(rq3_mat, axes = F, main = "Species Caught in Quarter 3")

##Quarter 4
ti_q4 <- data.frame()
  #Just looking at TI in Quarter 4 dataframe
for(i in spp_qtr4){
  sp1 <- subset(tripqtr4, spp2_id == i)
  for(j in spp_qtr4){
    q <- data.frame()
    sp2 <- subset(tripqtr4, spp2_id == j)
    if(i == j){ti_q4[i,j] <- 100}
    else{for(k in 1:nrow(sp1)){
      if((sp1$foCatNat[k] %in% sp2$foCatNat)&(sp1$area[k] %in% sp2$area)){
        q[k,1] <- sp1$wt[k]}
      else{q[k,1] <- 0}}
      ti_q4[i,j] <- (sum(q)/sum(sp1$wt))*100}} 
}   
q4_mat <- data.matrix(ti_q4)
rq4_mat <- apply(q4_mat, 2, rev)
image(rq4_mat, axes = F, main = "Species Caught in Quarter 4")

##############Important Commercial and Conservation Species######################

spp_38 <- c(142, 199, 86, 198, 161, 141, 56, 129, 128, 143, 177, 147, 146, 229, 179, 150, 
            119, 92, 11, 122, 69, 121, 219, 48, 45, 90, 71, 73, 70, 208, 195, 106, 148, 62, 104, 44, 207, 180) 
length(spp_38)
#These species were determined using the Scottish Sea Fisheries Statistics 2018, Table 3
  ##As well as the IUCN Red List

trip382 <- subset(trip3, spp_id %in% spp_38)
v243 <- c(1:243)
g_other <- setdiff(v243, spp_38)

trip38_other <- subset(trip3, spp_id %in% g_other)
trip38_other$spp <- "Other"
trip382 <- rbind(trip382, trip38_other)

head(trip382)
test38 <- unique(trip382$spp_id)
setdiff(test38, spp_38)


trip382$spp2_id <- as.numeric(as.factor(trip382$spp)) #creating unique species numerical codes for matrix
head(trip382)
spp_382 <- unique(trip382$spp2_id)#using the new numerical codes
spp_382
spp_382 <- sort(spp_382)

spp38_guide <- data.frame() #creating spp_id reference guide
spp38_guide[1:39,1] <- unique(trip382$spp) #adding one of each scientific names
spp38_guide[1:39,2] <- unique(trip382$spp2_id) #adding corresponding new numerical code
View(spp38_guide)
write.csv(spp38_guide, "39_Species_Guide.csv")

length(spp_382)

#########38 Species with All Strata###############
ti_38 <- data.frame() 

for(i in spp_382){
  sp1 <- subset(trip382, spp2_id == i) 
  for(j in spp_382){
    t <- data.frame()
    sp2 <- subset(trip382, spp2_id == j) 
    if(i == j){ti_38[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$area[k] %in% sp2$area)& 
         (sp1$quarter[k] %in% sp2$quarter)&(sp1$foCatNat[k] %in% sp2$foCatNat)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_38[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_38)

row.names(ti_38) <- unique(trip382$spp)
colnames(ti_38) <- unique(trip382$spp)
View(ti_38)

mat_38 <- data.matrix(ti_38)
rmat_38 <- apply(mat_38, 2, rev)
image(rmat_38, axes = F, main = "38 Species and All Strata")



###Overall 38+Other HC##############

dist_38 <- dist(mat_38, method = "euclidean")
hc_38<- hclust(dist_38, method = "ward.D2")
plot(hc_38)

coph_38 <- cophenetic(hc_38)
cor(dist_38, coph_38) #Greater than 0.75, therefore considered good

wss <- (nrow(test_38)-1)*sum(apply(test_38,2,var))
for (i in 2:5) wss[i] <- sum(kmeans(test_38, 
                                    centers=i)$withinss)
plot(1:5, wst, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") #2 Clusters



# Ward Hierarchical Clustering
d <- dist(mat_38, method = "euclidean") # distance matrix
fit1 <- hclust(d, method="ward.D2") 
plot(fit1) # display dendogram
groups <- cutree(fit1, k=4) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(fit1, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2 <- pvclust(t(mat_38), method.hclust="ward.D2",
               method.dist="euclidean")
plot(fit2)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2, alpha=.95, pv = "au", max.only = F) 

?ggdendrogram

dendro.38 <- as.dendrogram(hc_38)
dendro.plot.38 <- ggdendrogram(data = dendro.38, rotate = TRUE, k = 2)
dendro.plot.38 <- dendro.plot.38 + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.38)

melt.38 <- melt(mat_38)
head(melt.38)
colnames(melt.38) <- c("Targeted", "Bycatch", "TI")
View(mat_38)
melt.38$TI_Ranges <- cut(melt.38$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.38 <- order.dendrogram(dendro.38)
melt.38$Targeted <- factor(x = melt.38$Targeted,
                           levels = melt.38$Targeted[order.38],
                           ordered = T)


heatmap.38 <- ggplot(data=melt.38, 
                     aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.38)

grid.newpage()
print(heatmap.38, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.38, vp = viewport(x = 0.89, y = 0.591, width = 0.2, height = 0.6995))

#########Looking at Number of Trips in each Area by year and quarter##########
###Year
?tapply
trips6_year <- tapply(triparea6$trpCode, triparea6$year, length)#Area 27.6.a
print(trips6_year)

trips2_year <- tapply(triparea2$trpCode, triparea2$year, length)#Area 27.4
print(trips2_year)

trips7_year <- tapply(triparea7$trpCode, list(triparea7$year, triparea7$foCatNat), length)#Area 27.6.b 
print(trips7_year)#No trips in 2010

###Quarter
trips6_qtr <- tapply(triparea6$trpCode, triparea6$quarter, length)
print(trips6_qtr)

trips2_qtr <- tapply(triparea2$trpCode, triparea2$quarter, length)
print(trips2_qtr)

trips7_qtr <- tapply(triparea7$trpCode, triparea7$quarter, length)
print(trips7_qtr)

###Year and Quarter

trips6_yq <- tapply(triparea6$trpCode, list(triparea6$year, triparea6$quarter), length)
print(trips6_yq)

trips2_yq <- tapply(triparea2$trpCode, list(triparea2$year, triparea2$quarter), length)
print(trips2_yq)

trips7_yq <- tapply(triparea7$trpCode, list(triparea7$year, triparea7$quarter), length)
print(trips7_yq)

######################################################################
#######################START SUBSETS OF TRIP382######################
######################################################################
##AREA 27.4, TR1 Subsets######################################
tr1_274 <- subset(trip382, gear_id == 3 & area_id == 2) #creating TR1, Area 27.4 subset of 38 species data fraem
tr1_274

unique(tr1_274$year)
unique(tr1_274$foCatNat)
unique(tr1_274$area)
tr1_274$spp_id <- tr1_274$spp2_id #Creating unique species code for for-loop
tr1_274$spp2_id <- as.numeric(as.factor(tr1_274$spp))
spp_tr1274 <- unique(tr1_274$spp2_id)
spp_tr1274 <- sort(spp_tr1274)

ti_tr1274 <- data.frame() 

for(i in spp_tr1274){
  sp1 <- subset(tr1_274, spp2_id == i) 
  for(j in spp_tr1274){
    t <- data.frame()
    sp2 <- subset(tr1_274, spp2_id == j) 
    if(i == j){ti_tr1274[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$year[k] %in% sp2$year)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1274[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1274)

row.names(ti_tr1274) <- sort(unique(tr1_274$spp))
colnames(ti_tr1274) <- sort(unique(tr1_274$spp))
View(ti_tr1274)

mat_tr1274 <- as.matrix(ti_tr1274)
##HClust with TR1 27.4
dist_tr1274 <- dist(mat_tr1274, method = "euclidean") #distance matrix
hc_tr1274<- hclust(dist_tr1274, method = "ward.D2")
plot(hc_tr1274)

coph_1274 <- cophenetic(hc_tr1274)
cor(dist_tr1274, coph_1274) #0.83, Greater than 0.75, therefore considered good


fviz_nbclust(mat_tr1274, FUN = hcut, method = "silhouette") #optimal number of clusters is 2
?fviz_nbclust


# Ward Hierarchical Clustering
plot(hc_tr1274) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr1274, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2 <- pvclust(t(mat_tr1274), method.hclust="ward.D2",
                method.dist="euclidean")
plot(fit2)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2, alpha=.95)

dendro.1274 <- as.dendrogram(hc_tr1274)
dendro.plot.1274 <- ggdendrogram(data = dendro.1274, rotate = TRUE, k = 3)
dendro.plot.1274 <- dendro.plot.1274 + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.1274)

melt.1274 <- melt(mat_tr1274)
head(melt.1274)
colnames(melt.1274) <- c("Targeted", "Bycatch", "TI")
melt.1274$TI_Ranges <- cut(melt.1274$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1274 <- order.dendrogram(dendro.1274)
melt.1274$Targeted <- factor(x = melt.1274$Targeted,
                           levels = melt.1274$Targeted[order.1274],
                           ordered = T)


heatmap.1274 <- ggplot(data=melt.1274, 
                     aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1274)

grid.newpage()
print(heatmap.1274, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1274, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

###
###TR1 27.4 2009########
tr1_27409 <- subset(tr1_274, year == 2009)

unique(tr1_27409$year)
unique(tr1_27409$foCatNat)
unique(tr1_27409$area)
tr1_27409$spp2_id <- as.numeric(as.factor(tr1_27409$spp))
spp_tr127409 <- unique(tr1_27409$spp2_id)
spp_tr127409 <- sort(spp_tr127409)

ti_tr127409 <- data.frame() 

for(i in spp_tr127409){
  sp1 <- subset(tr1_27409, spp2_id == i) 
  for(j in spp_tr127409){
    t <- data.frame()
    sp2 <- subset(tr1_27409, spp2_id == j) 
    if(i == j){ti_tr127409[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127409[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127409)

row.names(ti_tr127409) <- sort(unique(tr1_27409$spp))
colnames(ti_tr127409) <- sort(unique(tr1_27409$spp))
View(ti_tr127409)

mat_tr127409 <- as.matrix(ti_tr127409)
dim(mat_tr127409)
##HClust with TR1 27.4 2009
dist_tr127409 <- dist(mat_tr127409, method = "euclidean") #distance matrix
hc_tr127409<- hclust(dist_tr127409, method = "ward.D2")
plot(hc_tr127409)

coph_127409 <- cophenetic(hc_tr127409)
cor(dist_tr127409, coph_127409) #Greater than 0.75, therefore considered good


fviz_nbclust(mat_tr127409, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127409) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127409, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127409 <- pvclust(t(mat_tr127409), method.hclust="ward.D2",
                method.dist="euclidean")
plot(fit127409)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127409, alpha=.95)

dendro.127409 <- as.dendrogram(hc_tr127409, k = 2)
dendro.plot.127409 <- ggdendrogram(data = dendro.127409, rotate = TRUE, k = 2)
dendro.plot.127409 <- dendro.plot.127409 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127409)

melt.127409 <- melt(mat_tr127409)
head(melt.127409)
colnames(melt.127409) <- c("Targeted", "Bycatch", "TI")
melt.127409$TI_Ranges <- cut(melt.127409$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127409 <- order.dendrogram(dendro.127409)
melt.127409$Targeted <- factor(x = melt.127409$Targeted,
                             levels = melt.127409$Targeted[order.127409],
                             ordered = T)

heatmap.127409 <- ggplot(data=melt.127409, 
                       aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127409)

grid.newpage()
print(heatmap.127409, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127409, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

##
###TR1 27.4 2010##########
tr1_27410 <- subset(tr1_274, year == 2010)

unique(tr1_27410$year)
unique(tr1_27410$foCatNat)
unique(tr1_27410$area)
tr1_27410$spp2_id <- as.numeric(as.factor(tr1_27410$spp))
spp_tr127410 <- unique(tr1_27410$spp2_id)
spp_tr127410 <- sort(spp_tr127410)

ti_tr127410 <- data.frame() 

for(i in spp_tr127410){
  sp1 <- subset(tr1_27410, spp2_id == i) 
  for(j in spp_tr127410){
    t <- data.frame()
    sp2 <- subset(tr1_27410, spp2_id == j) 
    if(i == j){ti_tr127410[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127410[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127410)

row.names(ti_tr127410) <- sort(unique(tr1_27410$spp))
colnames(ti_tr127410) <- sort(unique(tr1_27410$spp))
View(ti_tr127410)

mat_tr127410 <- as.matrix(ti_tr127410)
dim(mat_tr127410)


##HClust with TR1 27.4 2010
dist_tr127410 <- dist(mat_tr127410, method = "euclidean") #distance matrix
hc_tr127410 <- hclust(dist_tr127410, method = "ward.D2")
plot(hc_tr127410)

coph_127410 <- cophenetic(hc_tr127410)
cor(dist_tr127410, coph_127410) #0.65, Less than 0.75, suboptimal

fviz_nbclust(mat_tr127410, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


wss127410 <- (nrow(mat_tr127410)-1)*sum(apply(mat_tr127410,2,var))
for (i in 2:5) wss127410[i] <- sum(kmeans(mat_tr127410, 
                                          centers=i)$withinss)
plot(1:5, wss127410, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") #2 Clusters


# Ward Hierarchical Clustering
plot(hc_tr127410) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127410, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127410 <- pvclust(t(mat_tr127410), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127410)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127410, alpha=.95)

dendro.127410 <- as.dendrogram(hc_tr127410, k = 2)
dendro.plot.127410 <- ggdendrogram(data = dendro.127410, rotate = TRUE, k = 2)
dendro.plot.127410 <- dendro.plot.127410 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127410)

melt.127410 <- melt(mat_tr127410)
head(melt.127410)
colnames(melt.127410) <- c("Targeted", "Bycatch", "TI")
melt.127410$TI_Ranges <- cut(melt.127410$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127410 <- order.dendrogram(dendro.127410)
melt.127410$Targeted <- factor(x = melt.127410$Targeted,
                               levels = melt.127410$Targeted[order.127410],
                               ordered = T)

heatmap.127410 <- ggplot(data=melt.127410, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127410)

grid.newpage()
print(heatmap.127410, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127410, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))



##

###TR1 27.4 2011#######

tr1_27411 <- subset(tr1_274, year == 2011)

unique(tr1_27411$year)
unique(tr1_27411$foCatNat)
unique(tr1_27411$area)
tr1_27411$spp2_id <- as.numeric(as.factor(tr1_27411$spp))
spp_tr127411 <- sort(unique(tr1_27411$spp2_id))


ti_tr127411 <- data.frame() 

for(i in spp_tr127411){
  sp1 <- subset(tr1_27411, spp2_id == i) 
  for(j in spp_tr127411){
    t <- data.frame()
    sp2 <- subset(tr1_27411, spp2_id == j) 
    if(i == j){ti_tr127411[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127411[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127411)

row.names(ti_tr127411) <- sort(unique(tr1_27411$spp))
colnames(ti_tr127411) <- sort(unique(tr1_27411$spp))
View(ti_tr127411)

mat_tr127411 <- as.matrix(ti_tr127411)
dim(mat_tr127411)

##HClust with TR1 27.4 2011
dist_tr127411 <- dist(mat_tr127411, method = "euclidean") #distance matrix
hc_tr127411 <- hclust(dist_tr127411, method = "ward.D2")
plot(hc_tr127411)

coph_127411 <- cophenetic(hc_tr127411)
cor(dist_tr127411, coph_127411) #0.77, greater than 0.75 so that's good

library(FactoMineR)
library(factoextra)

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127411, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127411) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127411, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit127411 <- pvclust(t(mat_tr127411), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127411)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127411, alpha=.95)

dendro.127411 <- as.dendrogram(hc_tr127411, k = 2)
dendro.plot.127411 <- ggdendrogram(data = dendro.127411, rotate = TRUE, k = 2)
dendro.plot.127411 <- dendro.plot.127411 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127411)

melt.127411 <- melt(mat_tr127411)
head(melt.127411)
colnames(melt.127411) <- c("Targeted", "Bycatch", "TI")
melt.127411$TI_Ranges <- cut(melt.127411$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127411 <- order.dendrogram(dendro.127411)
melt.127411$Targeted <- factor(x = melt.127411$Targeted,
                               levels = melt.127411$Targeted[order.127411],
                               ordered = T)

heatmap.127411 <- ggplot(data=melt.127411, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127411)

grid.newpage()
print(heatmap.127411, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127411, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR1 27.4 2012######

tr1_27412 <- subset(tr1_274, year == 2012)

unique(tr1_27412$year)
unique(tr1_27412$foCatNat)
unique(tr1_27412$area)
tr1_27412$spp2_id <- as.numeric(as.factor(tr1_27412$spp))
spp_tr127412 <- sort(unique(tr1_27412$spp2_id))


ti_tr127412 <- data.frame() 

for(i in spp_tr127412){
  sp1 <- subset(tr1_27412, spp2_id == i) 
  for(j in spp_tr127412){
    t <- data.frame()
    sp2 <- subset(tr1_27412, spp2_id == j) 
    if(i == j){ti_tr127412[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127412[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127412)

row.names(ti_tr127412) <- sort(unique(tr1_27412$spp))
colnames(ti_tr127412) <- sort(unique(tr1_27412$spp))
View(ti_tr127412)

mat_tr127412 <- as.matrix(ti_tr127412)
dim(mat_tr127412)

##HClust with TR1 27.4 2012
dist_tr127412 <- dist(mat_tr127412, method = "euclidean") #distance matrix
hc_tr127412 <- hclust(dist_tr127412, method = "ward.D2")
plot(hc_tr127412)

coph_127412 <- cophenetic(hc_tr127412)
cor(dist_tr127413, coph_127412) #0.77, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127412, FUN = hcut, method = "silhouette") #optimal number of clusters is 3


# Ward Hierarchical Clustering
plot(hc_tr127412) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127412, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127412 <- pvclust(t(mat_tr127412), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127412)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127412, alpha=.95)

dendro.127412 <- as.dendrogram(hc_tr127412, k = 2)
dendro.plot.127412 <- ggdendrogram(data = dendro.127412, rotate = TRUE, k = 3)
dendro.plot.127412 <- dendro.plot.127412 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127412)

melt.127412 <- melt(mat_tr127412)
head(melt.127412)
colnames(melt.127412) <- c("Targeted", "Bycatch", "TI")
melt.127412$TI_Ranges <- cut(melt.127412$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127412 <- order.dendrogram(dendro.127412)
melt.127412$Targeted <- factor(x = melt.127412$Targeted,
                               levels = melt.127412$Targeted[order.127412],
                               ordered = T)

heatmap.127412 <- ggplot(data=melt.127412, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127412)

grid.newpage()
print(heatmap.127412, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127412, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.4 2013#####

tr1_27413 <- subset(tr1_274, year == 2013)

unique(tr1_27413$year)
unique(tr1_27413$foCatNat)
unique(tr1_27413$area)
tr1_27413$spp2_id <- as.numeric(as.factor(tr1_27413$spp))
spp_tr127413 <- sort(unique(tr1_27413$spp2_id))


ti_tr127413 <- data.frame() 

for(i in spp_tr127413){
  sp1 <- subset(tr1_27413, spp2_id == i) 
  for(j in spp_tr127413){
    t <- data.frame()
    sp2 <- subset(tr1_27413, spp2_id == j) 
    if(i == j){ti_tr127413[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127413[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127413)

row.names(ti_tr127413) <- sort(unique(tr1_27413$spp))
colnames(ti_tr127413) <- sort(unique(tr1_27413$spp))
View(ti_tr127413)

mat_tr127413 <- as.matrix(ti_tr127413)
dim(mat_tr127413)

##HClust with TR1 27.4 2013
dist_tr127413 <- dist(mat_tr127413, method = "euclidean") #distance matrix
hc_tr127413 <- hclust(dist_tr127413, method = "ward.D2")
plot(hc_tr127413)

coph_127413 <- cophenetic(hc_tr127413)
cor(dist_tr127413, coph_127413) #0.78, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127413, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127413) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127413, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127413 <- pvclust(t(mat_tr127413), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127413)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127413, alpha=.95)

dendro.127413 <- as.dendrogram(hc_tr127413, k = 2)
dendro.plot.127413 <- ggdendrogram(data = dendro.127413, rotate = TRUE, k = 2)
dendro.plot.127413 <- dendro.plot.127413 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127413)

melt.127413 <- melt(mat_tr127413)
head(melt.127413)
colnames(melt.127413) <- c("Targeted", "Bycatch", "TI")
melt.127413$TI_Ranges <- cut(melt.127413$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127413 <- order.dendrogram(dendro.127413)
melt.127413$Targeted <- factor(x = melt.127413$Targeted,
                               levels = melt.127413$Targeted[order.127413],
                               ordered = T)

heatmap.127413 <- ggplot(data=melt.127413, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127413)

grid.newpage()
print(heatmap.127413, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127413, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.4 2014#####

tr1_27414 <- subset(tr1_274, year == 2014)

unique(tr1_27414$year)
unique(tr1_27414$foCatNat)
unique(tr1_27414$area)
tr1_27414$spp2_id <- as.numeric(as.factor(tr1_27414$spp))
spp_tr127414 <- sort(unique(tr1_27414$spp2_id))


ti_tr127414 <- data.frame() 

for(i in spp_tr127414){
  sp1 <- subset(tr1_27414, spp2_id == i) 
  for(j in spp_tr127414){
    t <- data.frame()
    sp2 <- subset(tr1_27414, spp2_id == j) 
    if(i == j){ti_tr127414[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127414[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127414)

row.names(ti_tr127414) <- sort(unique(tr1_27414$spp))
colnames(ti_tr127414) <- sort(unique(tr1_27414$spp))
View(ti_tr127414)

mat_tr127414 <- as.matrix(ti_tr127414)
dim(mat_tr127414)

##HClust with TR1 27.4 2014
dist_tr127414 <- dist(mat_tr127414, method = "euclidean") #distance matrix
hc_tr127414 <- hclust(dist_tr127414, method = "ward.D2")
plot(hc_tr127414)

coph_127414 <- cophenetic(hc_tr127414)
cor(dist_tr127414, coph_127414) #0.73, less than 0.75 so less than optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127414, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127414) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127414, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127414 <- pvclust(t(mat_tr127414), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127414)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127414, alpha=.95)

dendro.127414 <- as.dendrogram(hc_tr127414, k = 2)
dendro.plot.127414 <- ggdendrogram(data = dendro.127414, rotate = TRUE, k = 2)
dendro.plot.127414 <- dendro.plot.127414 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127414)

melt.127414 <- melt(mat_tr127414)
head(melt.127414)
colnames(melt.127414) <- c("Targeted", "Bycatch", "TI")
melt.127414$TI_Ranges <- cut(melt.127414$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127414 <- order.dendrogram(dendro.127414)
melt.127414$Targeted <- factor(x = melt.127414$Targeted,
                               levels = melt.127414$Targeted[order.127414],
                               ordered = T)

heatmap.127414 <- ggplot(data=melt.127414, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127414)

grid.newpage()
print(heatmap.127414, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127414, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.4 2015####

tr1_27415 <- subset(tr1_274, year == 2015)

unique(tr1_27415$year)
unique(tr1_27415$foCatNat)
unique(tr1_27415$area)
tr1_27415$spp2_id <- as.numeric(as.factor(tr1_27415$spp))
spp_tr127415 <- sort(unique(tr1_27415$spp2_id))


ti_tr127415 <- data.frame() 

for(i in spp_tr127415){
  sp1 <- subset(tr1_27415, spp2_id == i) 
  for(j in spp_tr127415){
    t <- data.frame()
    sp2 <- subset(tr1_27415, spp2_id == j) 
    if(i == j){ti_tr127415[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127415[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127415)

row.names(ti_tr127415) <- sort(unique(tr1_27415$spp))
colnames(ti_tr127415) <- sort(unique(tr1_27415$spp))
View(ti_tr127415)

mat_tr127415 <- as.matrix(ti_tr127415)
dim(mat_tr127415)

##HClust with TR1 27.4 2015
dist_tr127415 <- dist(mat_tr127415, method = "euclidean") #distance matrix
hc_tr127415 <- hclust(dist_tr127415, method = "ward.D2")
plot(hc_tr127415)

coph_127415 <- cophenetic(hc_tr127415)
cor(dist_tr127415, coph_127415) #0.80, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127415, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127415) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127415, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127415 <- pvclust(t(mat_tr127415), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127415)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127415, alpha=.95)

dendro.127415 <- as.dendrogram(hc_tr127415, k = 2)
dendro.plot.127415 <- ggdendrogram(data = dendro.127415, rotate = TRUE, k = 2)
dendro.plot.127415 <- dendro.plot.127415 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127415)

melt.127415 <- melt(mat_tr127415)
head(melt.127415)
colnames(melt.127415) <- c("Targeted", "Bycatch", "TI")
melt.127415$TI_Ranges <- cut(melt.127415$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127415 <- order.dendrogram(dendro.127415)
melt.127415$Targeted <- factor(x = melt.127415$Targeted,
                               levels = melt.127415$Targeted[order.127415],
                               ordered = T)

heatmap.127415 <- ggplot(data=melt.127415, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127415)

grid.newpage()
print(heatmap.127415, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127415, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.4 2016####

tr1_27416 <- subset(tr1_274, year == 2016)

unique(tr1_27416$year)
unique(tr1_27416$foCatNat)
unique(tr1_27416$area)
tr1_27416$spp2_id <- as.numeric(as.factor(tr1_27416$spp))
spp_tr127416 <- sort(unique(tr1_27416$spp2_id))


ti_tr127416 <- data.frame() 

for(i in spp_tr127416){
  sp1 <- subset(tr1_27416, spp2_id == i) 
  for(j in spp_tr127416){
    t <- data.frame()
    sp2 <- subset(tr1_27416, spp2_id == j) 
    if(i == j){ti_tr127416[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127416[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127416)

row.names(ti_tr127416) <- sort(unique(tr1_27416$spp))
colnames(ti_tr127416) <- sort(unique(tr1_27416$spp))
View(ti_tr127416)

mat_tr127416 <- as.matrix(ti_tr127416)
dim(mat_tr127416)

##HClust with TR1 27.4 2016
dist_tr127416 <- dist(mat_tr127416, method = "euclidean") #distance matrix
hc_tr127416 <- hclust(dist_tr127416, method = "ward.D2")
plot(hc_tr127416)

coph_127416 <- cophenetic(hc_tr127416)
cor(dist_tr127416, coph_127416) #0.89, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127416, FUN = hcut, method = "silhouette") #optimal number of clusters is 3


# Ward Hierarchical Clustering
plot(hc_tr127416) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127416, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127416 <- pvclust(t(mat_tr127416), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127416)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127416, alpha=.95)

dendro.127416 <- as.dendrogram(hc_tr127416, k = 3)
dendro.plot.127416 <- ggdendrogram(data = dendro.127416, rotate = TRUE, k = 2)
dendro.plot.127416 <- dendro.plot.127416 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127416)

melt.127416 <- melt(mat_tr127416)
head(melt.127416)
colnames(melt.127416) <- c("Targeted", "Bycatch", "TI")
melt.127416$TI_Ranges <- cut(melt.127416$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127416 <- order.dendrogram(dendro.127416)
melt.127416$Targeted <- factor(x = melt.127416$Targeted,
                               levels = melt.127416$Targeted[order.127416],
                               ordered = T)

heatmap.127416 <- ggplot(data=melt.127416, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127416)

grid.newpage()
print(heatmap.127416, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127416, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.4 2017####

tr1_27417 <- subset(tr1_274, year == 2017)

unique(tr1_27417$year)
unique(tr1_27417$foCatNat)
unique(tr1_27417$area)
tr1_27417$spp2_id <- as.numeric(as.factor(tr1_27417$spp))
spp_tr127417 <- sort(unique(tr1_27417$spp2_id))


ti_tr127417 <- data.frame() 

for(i in spp_tr127417){
  sp1 <- subset(tr1_27417, spp2_id == i) 
  for(j in spp_tr127417){
    t <- data.frame()
    sp2 <- subset(tr1_27417, spp2_id == j) 
    if(i == j){ti_tr127417[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127417[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127417)

row.names(ti_tr127417) <- sort(unique(tr1_27417$spp))
colnames(ti_tr127417) <- sort(unique(tr1_27417$spp))
View(ti_tr127417)

mat_tr127417 <- as.matrix(ti_tr127417)
dim(mat_tr127417)

##HClust with TR1 27.4 2017
dist_tr127417 <- dist(mat_tr127417, method = "euclidean") #distance matrix
hc_tr127417 <- hclust(dist_tr127417, method = "ward.D2")
plot(hc_tr127417)

coph_127417 <- cophenetic(hc_tr127417)
cor(dist_tr127417, coph_127417) #0.89, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127417, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127417) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127417, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127417 <- pvclust(t(mat_tr127417), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127417)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127417, alpha=.95)

dendro.127417 <- as.dendrogram(hc_tr127417, k = 2)
dendro.plot.127417 <- ggdendrogram(data = dendro.127417, rotate = TRUE, k = 2)
dendro.plot.127417 <- dendro.plot.127417 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127417)

melt.127417 <- melt(mat_tr127417)
head(melt.127417)
colnames(melt.127417) <- c("Targeted", "Bycatch", "TI")
melt.127417$TI_Ranges <- cut(melt.127417$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127417 <- order.dendrogram(dendro.127417)
melt.127417$Targeted <- factor(x = melt.127417$Targeted,
                               levels = melt.127417$Targeted[order.127417],
                               ordered = T)

heatmap.127417 <- ggplot(data=melt.127417, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127417)

grid.newpage()
print(heatmap.127417, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127417, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.4 2018####

tr1_27418 <- subset(tr1_274, year == 2018)

unique(tr1_27418$year)
unique(tr1_27418$foCatNat)
unique(tr1_27418$area)
tr1_27418$spp2_id <- as.numeric(as.factor(tr1_27418$spp))
spp_tr127418 <- sort(unique(tr1_27418$spp2_id))


ti_tr127418 <- data.frame() 

for(i in spp_tr127418){
  sp1 <- subset(tr1_27418, spp2_id == i) 
  for(j in spp_tr127418){
    t <- data.frame()
    sp2 <- subset(tr1_27418, spp2_id == j) 
    if(i == j){ti_tr127418[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127418[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127418)

row.names(ti_tr127418) <- sort(unique(tr1_27418$spp))
colnames(ti_tr127418) <- sort(unique(tr1_27418$spp))
View(ti_tr127418)

mat_tr127418 <- as.matrix(ti_tr127418)
dim(mat_tr127418)

##HClust with TR1 27.4 2017
dist_tr127418 <- dist(mat_tr127418, method = "euclidean") #distance matrix
hc_tr127418 <- hclust(dist_tr127418, method = "ward.D2")
plot(hc_tr127418)

coph_127418 <- cophenetic(hc_tr127418)
cor(dist_tr127418, coph_127418) #0.87, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr127418, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr127418) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr127418, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit127418 <- pvclust(t(mat_tr127418), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit127418)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit127418, alpha=.95)

dendro.127418 <- as.dendrogram(hc_tr127418, k = 2)
dendro.plot.127418 <- ggdendrogram(data = dendro.127418, rotate = TRUE, k = 2)
dendro.plot.127418 <- dendro.plot.127418 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.127418)

melt.127418 <- melt(mat_tr127418)
head(melt.127418)
colnames(melt.127418) <- c("Targeted", "Bycatch", "TI")
melt.127418$TI_Ranges <- cut(melt.127418$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.127418 <- order.dendrogram(dendro.127418)
melt.127418$Targeted <- factor(x = melt.127418$Targeted,
                               levels = melt.127418$Targeted[order.127418],
                               ordered = T)

heatmap.127418 <- ggplot(data=melt.127418, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.127418)

grid.newpage()
print(heatmap.127418, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127418, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###Compiling TR1 27.4 Heatmaps#####
pdf(file = "TR1_274_Heatmaps.pdf")

grid.newpage()
print(heatmap.1274, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1274, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127409, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127409, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127410, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127410, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127411, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127411, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127412, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127412, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127413, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127413, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127414, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127414, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127415, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127415, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127416, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127416, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127417, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127417, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.127418, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.127418, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

dev.off()
###Heatmaps for Averaging#####
ti_tr127409.2 <- data.frame() 
ti_tr127409.2[1:39, 1:39] <- -1

spp_tr127409.2 <- unique(tr1_27409$spp_id)
spp_tr127409.2 <- sort(spp_tr127409.2)

for(i in spp_tr127409.2){
  sp1 <- subset(tr1_27409, spp_id == i) 
  for(j in spp_tr127409.2){
    t <- data.frame()
    sp2 <- subset(tr1_27409, spp_id == j) 
    if(i == j){ti_tr127409.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127409.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127409.2)

row.names(ti_tr127409.2) <- sort(unique(trip382$spp))
colnames(ti_tr127409.2) <- sort(unique(trip382$spp))

mat_tr127409.2 <- as.matrix(ti_tr127409.2)
dim(mat_tr127409.2)
##2010
ti_tr127410.2 <- data.frame() 
ti_tr127410.2[1:39, 1:39] <- -1

spp_tr127410.2 <- unique(tr1_27410$spp_id)
spp_tr127410.2 <- sort(spp_tr127410.2)

for(i in spp_tr127410.2){
  sp1 <- subset(tr1_27410, spp_id == i) 
  for(j in spp_tr127410.2){
    t <- data.frame()
    sp2 <- subset(tr1_27410, spp_id == j) 
    if(i == j){ti_tr127410.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127410.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127410.2)

row.names(ti_tr127410.2) <- sort(unique(trip382$spp))
colnames(ti_tr127410.2) <- sort(unique(trip382$spp))

mat_tr127410.2 <- as.matrix(ti_tr127410.2)
dim(mat_tr127410.2)

##2011
ti_tr127411.2 <- data.frame() 
ti_tr127411.2[1:39, 1:39] <- -1

spp_tr127411.2 <- unique(tr1_27411$spp_id)
spp_tr127411.2 <- sort(spp_tr127411.2)

for(i in spp_tr127411.2){
  sp1 <- subset(tr1_27411, spp_id == i) 
  for(j in spp_tr127411.2){
    t <- data.frame()
    sp2 <- subset(tr1_27411, spp_id == j) 
    if(i == j){ti_tr127411.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127411.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127411.2)

row.names(ti_tr127411.2) <- sort(unique(trip382$spp))
colnames(ti_tr127411.2) <- sort(unique(trip382$spp))

mat_tr127411.2 <- as.matrix(ti_tr127411.2)
dim(mat_tr127411.2)

##2012
ti_tr127412.2 <- data.frame() 
ti_tr127412.2[1:39, 1:39] <- -1

spp_tr127412.2 <- unique(tr1_27412$spp_id)
spp_tr127412.2 <- sort(spp_tr127412.2)

for(i in spp_tr127412.2){
  sp1 <- subset(tr1_27412, spp_id == i) 
  for(j in spp_tr127412.2){
    t <- data.frame()
    sp2 <- subset(tr1_27412, spp_id == j) 
    if(i == j){ti_tr127412.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127412.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127410.2)

row.names(ti_tr127412.2) <- sort(unique(trip382$spp))
colnames(ti_tr127412.2) <- sort(unique(trip382$spp))

mat_tr127412.2 <- as.matrix(ti_tr127412.2)
dim(mat_tr127412.2)

##2013
ti_tr127413.2 <- data.frame() 
ti_tr127413.2[1:39, 1:39] <- -1

spp_tr127413.2 <- unique(tr1_27413$spp_id)
spp_tr127413.2 <- sort(spp_tr127413.2)

for(i in spp_tr127413.2){
  sp1 <- subset(tr1_27413, spp_id == i) 
  for(j in spp_tr127413.2){
    t <- data.frame()
    sp2 <- subset(tr1_27413, spp_id == j) 
    if(i == j){ti_tr127413.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127413.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127413.2)

row.names(ti_tr127413.2) <- sort(unique(trip382$spp))
colnames(ti_tr127413.2) <- sort(unique(trip382$spp))

mat_tr127413.2 <- as.matrix(ti_tr127413.2)
dim(mat_tr127413.2)

##2014
ti_tr127414.2 <- data.frame() 
ti_tr127414.2[1:39, 1:39] <- -1

spp_tr127414.2 <- unique(tr1_27414$spp_id)
spp_tr127414.2 <- sort(spp_tr127414.2)

for(i in spp_tr127414.2){
  sp1 <- subset(tr1_27414, spp_id == i) 
  for(j in spp_tr127414.2){
    t <- data.frame()
    sp2 <- subset(tr1_27414, spp_id == j) 
    if(i == j){ti_tr127414.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127414.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127414.2)

row.names(ti_tr127414.2) <- sort(unique(trip382$spp))
colnames(ti_tr127414.2) <- sort(unique(trip382$spp))

mat_tr127414.2 <- as.matrix(ti_tr127414.2)
dim(mat_tr127414.2)

##2015
ti_tr127415.2 <- data.frame() 
ti_tr127415.2[1:39, 1:39] <- -1

spp_tr127415.2 <- unique(tr1_27415$spp_id)
spp_tr127415.2 <- sort(spp_tr127415.2)

for(i in spp_tr127415.2){
  sp1 <- subset(tr1_27415, spp_id == i) 
  for(j in spp_tr127415.2){
    t <- data.frame()
    sp2 <- subset(tr1_27415, spp_id == j) 
    if(i == j){ti_tr127415.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127415.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127415.2)

row.names(ti_tr127415.2) <- sort(unique(trip382$spp))
colnames(ti_tr127415.2) <- sort(unique(trip382$spp))

mat_tr127415.2 <- as.matrix(ti_tr127415.2)
dim(mat_tr127410.2)

##2016
ti_tr127416.2 <- data.frame() 
ti_tr127416.2[1:39, 1:39] <- -1

spp_tr127416.2 <- unique(tr1_27416$spp_id)
spp_tr127416.2 <- sort(spp_tr127416.2)

for(i in spp_tr127416.2){
  sp1 <- subset(tr1_27416, spp_id == i) 
  for(j in spp_tr127416.2){
    t <- data.frame()
    sp2 <- subset(tr1_27416, spp_id == j) 
    if(i == j){ti_tr127416.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127416.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127416.2)

row.names(ti_tr127416.2) <- sort(unique(trip382$spp))
colnames(ti_tr127416.2) <- sort(unique(trip382$spp))

mat_tr127416.2 <- as.matrix(ti_tr127416.2)
dim(mat_tr127416.2)

##2017
ti_tr127417.2 <- data.frame() 
ti_tr127417.2[1:39, 1:39] <- -1

spp_tr127417.2 <- unique(tr1_27417$spp_id)
spp_tr127417.2 <- sort(spp_tr127417.2)

for(i in spp_tr127417.2){
  sp1 <- subset(tr1_27417, spp_id == i) 
  for(j in spp_tr127417.2){
    t <- data.frame()
    sp2 <- subset(tr1_27417, spp_id == j) 
    if(i == j){ti_tr127417.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127417.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127417.2)

row.names(ti_tr127417.2) <- sort(unique(trip382$spp))
colnames(ti_tr127417.2) <- sort(unique(trip382$spp))

mat_tr127417.2 <- as.matrix(ti_tr127417.2)
dim(mat_tr127417.2)

##2018
ti_tr127418.2 <- data.frame() 
ti_tr127418.2[1:39, 1:39] <- -1

spp_tr127418.2 <- unique(tr1_27418$spp_id)
spp_tr127418.2 <- sort(spp_tr127418.2)

for(i in spp_tr127418.2){
  sp1 <- subset(tr1_27418, spp_id == i) 
  for(j in spp_tr127418.2){
    t <- data.frame()
    sp2 <- subset(tr1_27418, spp_id == j) 
    if(i == j){ti_tr127418.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr127418.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr127418.2)

row.names(ti_tr127418.2) <- sort(unique(trip382$spp))
colnames(ti_tr127418.2) <- sort(unique(trip382$spp))

mat_tr127418.2 <- as.matrix(ti_tr127418.2)
dim(mat_tr127418.2)
##Averaging the Matrices together
matavg_tr1274 <- (mat_tr127409.2 + mat_tr127410.2 + mat_tr127411.2 + 
             mat_tr127412.2 + mat_tr127413.2 + mat_tr127414.2 + 
             mat_tr127415.2 + mat_tr127416.2 +  mat_tr127417.2 + mat_tr127418.2)/10


dist_avg1274 <- dist(matavg_tr1274, method = "euclidean") #distance matrix
hc_avg1274 <- hclust(dist_avg1274, method = "ward.D2")
plot(hc_avg1274)

coph_avg1274 <- cophenetic(hc_avg1274)
cor(dist_avg1274, coph_avg1274) #0.93, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(matavg_tr1274, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_avg1274) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_avg1274, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitavg1274 <- pvclust(t(matavg_tr1274), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fitavg1274)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitavg1274, alpha=.95)

dendro.avg1274 <- as.dendrogram(hc_avg1274, k = 2)
dendro.plot.avg1274 <- ggdendrogram(data = dendro.avg1274, rotate = TRUE, k = 2)
dendro.plot.avg1274 <- dendro.plot.avg1274 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.avg1274)

melt.avg1274 <- melt(matavg_tr1274)
head(melt.avg1274 )
colnames(melt.avg1274 ) <- c("Targeted", "Bycatch", "TI")
melt.avg1274 $TI_Ranges <- cut(melt.avg1274 $TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.avg1274 <- order.dendrogram(dendro.avg1274)
melt.avg1274$Targeted <- factor(x = melt.avg1274$Targeted,
                               levels = melt.avg1274$Targeted[order.avg1274],
                               ordered = T)

heatmap.avg1274 <- ggplot(data=melt.avg1274, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.avg1274)

grid.newpage()
print(heatmap.avg1274, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.avg1274, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

####Creating Array TR1 27.4####
a_tr1274 <- abind(mat_tr127409.2, mat_tr127410.2, mat_tr127411.2, 
              mat_tr127412.2, mat_tr127413.2, mat_tr127414.2,
              mat_tr127415.2, mat_tr127416.2, mat_tr127417.2, mat_tr127418.2, along = 3)
View(a_tr1274[,,10])

melta_tr1274 <- melt(a_tr1274)
head(melta_tr1274)
colnames(melta_tr1274) <- c("Targeted", "Bycatch", "Year", "TIS")

melta_tr1274$TI_Ranges <- cut(melta_tr1274$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
str(melta_tr1274)


melta_tr1274$Area[1:15210] <- "27.4"
melta_tr1274$Area <- as.factor(melta_tr1274$Area)
melta_tr1274$Gear[1:15210] <- "TR1"
melta_tr1274$Gear <- as.factor(melta_tr1274$Gear)

melta_tr1274$Year <- as.factor(melta_tr1274$Year)
melta_tr1274$Year <- revalue(melta_tr1274$Year, c("1"="2009", "2"="2010", "3"="2011", "4"="2012",
                             "5"="2013", "6"="2014", "7"="2015", "8"="2016",
                             "9"="2017", "10"="2018"))

melta_tr1274 <- na.omit(melta_tr1274) #removing absent species that are marked as NA in TI_Ranges

tr1274_by <- with(melta_tr1274, tapply(TIS, INDEX = list(Bycatch, Year), FUN = "mean"))
tr1274_tar <- with(melta_tr1274, tapply(TIS, INDEX = list(Targeted, Year), FUN = "mean"))
tr1274_overall <- with(melta_tr1274, tapply (TIS, INDEX = Year, FUN = "mean"))
tr1274_by <- as.data.frame(tr1274_by)
tr1274_tar <- as.data.frame(tr1274_tar)
tr1274_overall <- as.data.frame(tr1274_overall)

tr1274_overall$var <- with(melta_tr1274, tapply(TIS, INDEX = Year, FUN = "var"))
tr1274_overall$sd <- sqrt(tr1274_overall$var)
tr1274_overall$cv <- 100*((tr1274_overall$sd)/(tr1274_overall$tr1274_overall))

tr1274_by$Overall <- rowMeans(tr1274_by)
tr1274_by$var <- with(melta_tr1274, tapply(TIS, INDEX = list(Bycatch,Year), FUN = "var"))
tr1274_by$sd <- sqrt(tr1274_by$var)
tr1274_by$cv <- 100*((tr1274_by$sd)/(tr1274_by$Overall))
tr1274_by.05 <- tr1274_by[!(tr1274_by$Overall<=5),]
tr1274_by.50 <- tr1274_by[!(tr1274_by$Overall<=50),]
View(tr1274_by.05)
View(tr1274_by.50)

tr1274_tar$Overall <- rowMeans(tr1274_tar)
tr1274_tar$var <- with(melta_tr1274, tapply(TIS, INDEX = list(Targeted,Year), FUN = "var"))
tr1274_tar$sd <- sqrt(tr1274_tar$var)
tr1274_tar$cv <- 100*((tr1274_tar$sd)/(tr1274_tar$Overall))
tr1274_tar.05 <- tr1274_tar[!(tr1274_tar$Overall<=5),]
tr1274_tar.25 <- tr1274_tar[!(tr1274_tar$Overall<=25),]
View(tr1274_tar.05)
View(tr1274_tar.25)


####Centered Moving Average TR1 27.4--EACH YEAR#####
?rollmean
?select
# 
# test <- tr1274_overall %>% select(TIS_year = tr1274_overall) %>% 
#   mutate(TIS_mavg3 = rollmean(TIS_year, k = 3, fill = NA),
#          Year = 2009:2018)
# testg<- test %>% gather(metric, value, TIS_year:TIS_mavg3) %>%
#   ggplot(aes(Year, value, color = metric)) +
#   geom_line()
# print(testg)

t1274 <- data.frame() #creating dataframe that will hold centered moving average data; Centered moving average used to look at stability of each Area+Gear type combination
for(i in 2009:2018){
  selected <- c(i-1, i, i+1)
  f1 <- subset(melta_tr1274, Year %in% selected)
  if((i == 2009)| (i==2018)){
  t1274[i-2008, 1] <- NA} 
  else{t1274[i-2008,1] <- mean(f1$TIS)}
}
print(t1274)

t1274.2 <- t1274 %>% select(TIS_mavg = V1) %>% mutate(TIS_yavg = tr1274_overall$tr1274_overall,
                                                Year = 2009:2018) #Creating modified dataframe that will include numerical years for easy plotting
t1274.g <- t1274.2 %>% gather(Metric, TIS, TIS_mavg:TIS_yavg) %>%
  ggplot(aes(Year, TIS, color = Metric)) +
  geom_line() #Plotting 1274 CMA
print(t1274.g)



##########################################################################
##Area 27.4, TR2 Subsets#################
tr2_274 <- subset(trip382, gear_id == 4 & area_id == 2)
tr2_274

unique(tr2_274$year)
unique(tr2_274$foCatNat)
unique(tr2_274$area)
tr2_274$spp_id <- tr2_274$spp2_id
tr2_274$spp2_id <- as.numeric(as.factor(tr2_274$spp))
spp_tr2274 <- unique(tr2_274$spp2_id)
spp_tr2274 <- sort(spp_tr2274)

ti_tr2274 <- data.frame() 

for(i in spp_tr2274){
  sp1 <- subset(tr2_274, spp2_id == i) 
  for(j in spp_tr2274){
    t <- data.frame()
    sp2 <- subset(tr2_274, spp2_id == j) 
    if(i == j){ti_tr2274[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$year[k] %in% sp2$year)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2274[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2274)

row.names(ti_tr2274) <- sort(unique(tr2_274$spp))
colnames(ti_tr2274) <- sort(unique(tr2_274$spp))
View(ti_tr2274)

mat_tr2274 <- as.matrix(ti_tr2274)
dim(mat_tr2274)

##HClust with TR2 27.4
dist_tr2274 <- dist(mat_tr2274, method = "euclidean") #distance matrix
hc_tr2274<- hclust(dist_tr2274, method = "ward.D2")
plot(hc_tr2274)

coph_2274 <- cophenetic(hc_tr2274)
cor(dist_tr2274, coph_2274) #0.83, Greater than 0.75, therefore considered good


fviz_nbclust(mat_tr2274, FUN = hcut, method = "silhouette") #optimal number of clusters is 2
?fviz_nbclust


# Ward Hierarchical Clustering
plot(hc_tr2274) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr2274, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2 <- pvclust(t(mat_tr2274), method.hclust="ward.D2",
                method.dist="euclidean")
plot(fit2)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2, alpha=.95)

dendro.2274 <- as.dendrogram(hc_tr2274)
dendro.plot.2274 <- ggdendrogram(data = dendro.2274, rotate = TRUE, k = 2)
dendro.plot.2274 <- dendro.plot.2274 + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.2274)

melt.2274 <- melt(mat_tr2274)
head(melt.2274)
colnames(melt.2274) <- c("Targeted", "Bycatch", "TI")
melt.2274$TI_Ranges <- cut(melt.2274$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2274 <- order.dendrogram(dendro.2274)
melt.2274$Targeted <- factor(x = melt.2274$Targeted,
                             levels = melt.2274$Targeted[order.2274],
                             ordered = T)


heatmap.2274 <- ggplot(data=melt.2274, 
                       aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2274)

grid.newpage()
print(heatmap.2274, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2274, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

###
###TR2 27.4 2009########
tr2_27409 <- subset(tr2_274, year == 2009)

unique(tr2_27409$year)
unique(tr2_27409$foCatNat)
unique(tr2_27409$area)
tr2_27409$spp2_id <- as.numeric(as.factor(tr2_27409$spp))
spp_tr227409 <- unique(tr2_27409$spp2_id)
spp_tr227409 <- sort(spp_tr227409)

ti_tr227409 <- data.frame() 

for(i in spp_tr227409){
  sp1 <- subset(tr2_27409, spp2_id == i) 
  for(j in spp_tr227409){
    t <- data.frame()
    sp2 <- subset(tr2_27409, spp2_id == j) 
    if(i == j){ti_tr227409[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227409[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227409)

row.names(ti_tr227409) <- sort(unique(tr2_27409$spp))
colnames(ti_tr227409) <- sort(unique(tr2_27409$spp))
View(ti_tr227409)

mat_tr227409 <- as.matrix(ti_tr227409)
dim(mat_tr127409)
##HClust with TR2 27.4 2009
dist_tr227409 <- dist(mat_tr227409, method = "euclidean") #distance matrix
hc_tr227409<- hclust(dist_tr227409, method = "ward.D2")
plot(hc_tr227409)

coph_227409 <- cophenetic(hc_tr227409)
cor(dist_tr227409, coph_227409) #0.73, therefore not considered good


fviz_nbclust(mat_tr227409, FUN = hcut, method = "silhouette") #optimal number of clusters is 2



# Ward Hierarchical Clustering
plot(hc_tr227409) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227409, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227409 <- pvclust(t(mat_tr227409), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227409)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227409, alpha=.95)

dendro.227409 <- as.dendrogram(hc_tr227409, k = 2)
dendro.plot.227409 <- ggdendrogram(data = dendro.227409, rotate = TRUE, k = 2)
dendro.plot.227409 <- dendro.plot.227409 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227409)

melt.227409 <- melt(mat_tr227409)
head(melt.227409)
colnames(melt.227409) <- c("Targeted", "Bycatch", "TI")
melt.227409$TI_Ranges <- cut(melt.227409$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227409 <- order.dendrogram(dendro.227409)
melt.227409$Targeted <- factor(x = melt.227409$Targeted,
                               levels = melt.227409$Targeted[order.227409],
                               ordered = T)

heatmap.227409 <- ggplot(data=melt.227409, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227409)

grid.newpage()
print(heatmap.227409, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227409, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

##
###TR2 27.4 2010##########
tr2_27410 <- subset(tr2_274, year == 2010)

unique(tr2_27410$year)
unique(tr2_27410$foCatNat)
unique(tr2_27410$area)
tr2_27410$spp2_id <- as.numeric(as.factor(tr2_27410$spp))
spp_tr227410 <- unique(tr2_27410$spp2_id)
spp_tr227410 <- sort(spp_tr227410)

ti_tr227410 <- data.frame() 

for(i in spp_tr227410){
  sp1 <- subset(tr2_27410, spp2_id == i) 
  for(j in spp_tr227410){
    t <- data.frame()
    sp2 <- subset(tr2_27410, spp2_id == j) 
    if(i == j){ti_tr227410[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227410[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227410)

row.names(ti_tr227410) <- sort(unique(tr2_27410$spp))
colnames(ti_tr227410) <- sort(unique(tr2_27410$spp))
View(ti_tr227410)

mat_tr227410 <- as.matrix(ti_tr227410)
dim(mat_tr227410)


##HClust with TR2 27.4 2010
dist_tr227410 <- dist(mat_tr227410, method = "euclidean") #distance matrix
hc_tr227410 <- hclust(dist_tr227410, method = "ward.D2")
plot(hc_tr227410)

coph_227410 <- cophenetic(hc_tr227410)
cor(dist_tr227410, coph_227410) #0.82, greater than 0.75, suboptimal

fviz_nbclust(mat_tr227410, FUN = hcut, method = "silhouette") #optimal number of clusters is 7




# Ward Hierarchical Clustering
plot(hc_tr227410) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227410, k=7, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227410 <- pvclust(t(mat_tr227410), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227410)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227410, alpha=.95)

dendro.227410 <- as.dendrogram(hc_tr227410, k = 7)
dendro.plot.227410 <- ggdendrogram(data = dendro.227410, rotate = TRUE, k = 7)
dendro.plot.227410 <- dendro.plot.227410 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227410)

melt.227410 <- melt(mat_tr227410)
head(melt.227410)
colnames(melt.227410) <- c("Targeted", "Bycatch", "TI")
melt.227410$TI_Ranges <- cut(melt.227410$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227410 <- order.dendrogram(dendro.227410)
melt.227410$Targeted <- factor(x = melt.227410$Targeted,
                               levels = melt.227410$Targeted[order.227410],
                               ordered = T)

heatmap.227410 <- ggplot(data=melt.227410, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227410)

grid.newpage()
print(heatmap.227410, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227410, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))



##

###TR2 27.4 2011#######

tr2_27411 <- subset(tr2_274, year == 2011)

unique(tr2_27411$year)
unique(tr2_27411$foCatNat)
unique(tr2_27411$area)
tr2_27411$spp2_id <- as.numeric(as.factor(tr2_27411$spp))
spp_tr227411 <- sort(unique(tr2_27411$spp2_id))


ti_tr227411 <- data.frame() 

for(i in spp_tr227411){
  sp1 <- subset(tr2_27411, spp2_id == i) 
  for(j in spp_tr227411){
    t <- data.frame()
    sp2 <- subset(tr2_27411, spp2_id == j) 
    if(i == j){ti_tr227411[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227411[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227411)

row.names(ti_tr227411) <- sort(unique(tr2_27411$spp))
colnames(ti_tr227411) <- sort(unique(tr2_27411$spp))
View(ti_tr227411)

mat_tr227411 <- as.matrix(ti_tr227411)
dim(mat_tr227411)

##HClust with TR2 27.4 2011
dist_tr227411 <- dist(mat_tr227411, method = "euclidean") #distance matrix
hc_tr227411 <- hclust(dist_tr227411, method = "ward.D2")
plot(hc_tr227411)

coph_227411 <- cophenetic(hc_tr227411)
cor(dist_tr227411, coph_227411) #0.75, equal to 0.75 so is okay


#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227411, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr227411) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227411, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227411 <- pvclust(t(mat_tr227411), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227411)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227411, alpha=.95)

dendro.227411 <- as.dendrogram(hc_tr227411, k = 2)
dendro.plot.227411 <- ggdendrogram(data = dendro.227411, rotate = TRUE, k = 2)
dendro.plot.227411 <- dendro.plot.227411 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227411)

melt.227411 <- melt(mat_tr227411)
head(melt.227411)
colnames(melt.227411) <- c("Targeted", "Bycatch", "TI")
melt.227411$TI_Ranges <- cut(melt.227411$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227411 <- order.dendrogram(dendro.227411)
melt.227411$Targeted <- factor(x = melt.227411$Targeted,
                               levels = melt.227411$Targeted[order.227411],
                               ordered = T)

heatmap.227411 <- ggplot(data=melt.227411, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227411)

grid.newpage()
print(heatmap.227411, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227411, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.4 2012######

tr2_27412 <- subset(tr2_274, year == 2012)

unique(tr2_27412$year)
unique(tr2_27412$foCatNat)
unique(tr2_27412$area)
tr2_27412$spp2_id <- as.numeric(as.factor(tr2_27412$spp))
spp_tr227412 <- sort(unique(tr2_27412$spp2_id))


ti_tr227412 <- data.frame() 

for(i in spp_tr227412){
  sp1 <- subset(tr2_27412, spp2_id == i) 
  for(j in spp_tr227412){
    t <- data.frame()
    sp2 <- subset(tr2_27412, spp2_id == j) 
    if(i == j){ti_tr227412[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227412[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227412)

row.names(ti_tr227412) <- sort(unique(tr2_27412$spp))
colnames(ti_tr227412) <- sort(unique(tr2_27412$spp))
View(ti_tr227412)

mat_tr227412 <- as.matrix(ti_tr227412)
dim(mat_tr227412)

##HClust with TR2 27.4 2012
dist_tr227412 <- dist(mat_tr227412, method = "euclidean") #distance matrix
hc_tr227412 <- hclust(dist_tr227412, method = "ward.D2")
plot(hc_tr227412)

coph_227412 <- cophenetic(hc_tr227412)
cor(dist_tr227412, coph_227412) #0.78, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227412, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_tr227412) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227412, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227412 <- pvclust(t(mat_tr227412), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227412)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227412, alpha=.95)

dendro.227412 <- as.dendrogram(hc_tr227412, k = 4)
dendro.plot.227412 <- ggdendrogram(data = dendro.227412, rotate = TRUE, k = 4)
dendro.plot.227412 <- dendro.plot.227412 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227412)

melt.227412 <- melt(mat_tr227412)
head(melt.227412)
colnames(melt.227412) <- c("Targeted", "Bycatch", "TI")
melt.227412$TI_Ranges <- cut(melt.227412$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227412 <- order.dendrogram(dendro.227412)
melt.227412$Targeted <- factor(x = melt.227412$Targeted,
                               levels = melt.127412$Targeted[order.127412],
                               ordered = T)

heatmap.227412 <- ggplot(data=melt.227412, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227412)

grid.newpage()
print(heatmap.227412, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227412, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.4 2013#####

tr2_27413 <- subset(tr2_274, year == 2013)

unique(tr2_27413$year)
unique(tr2_27413$foCatNat)
unique(tr2_27413$area)
tr2_27413$spp2_id <- as.numeric(as.factor(tr2_27413$spp))
spp_tr227413 <- sort(unique(tr2_27413$spp2_id))


ti_tr227413 <- data.frame() 

for(i in spp_tr227413){
  sp1 <- subset(tr2_27413, spp2_id == i) 
  for(j in spp_tr227413){
    t <- data.frame()
    sp2 <- subset(tr2_27413, spp2_id == j) 
    if(i == j){ti_tr227413[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227413[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227413)

row.names(ti_tr227413) <- sort(unique(tr2_27413$spp))
colnames(ti_tr227413) <- sort(unique(tr2_27413$spp))
View(ti_tr227413)

mat_tr227413 <- as.matrix(ti_tr227413)
dim(mat_tr227413)

##HClust with TR2 27.4 2013
dist_tr227413 <- dist(mat_tr227413, method = "euclidean") #distance matrix
hc_tr227413 <- hclust(dist_tr227413, method = "ward.D2")
plot(hc_tr227413)

coph_227413 <- cophenetic(hc_tr227413)
cor(dist_tr227413, coph_227413) #0.74, lesser than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227413, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr227413) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227413, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227413 <- pvclust(t(mat_tr227413), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227413)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227413, alpha=.95)

dendro.227413 <- as.dendrogram(hc_tr227413, k = 2)
dendro.plot.227413 <- ggdendrogram(data = dendro.227413, rotate = TRUE, k = 2)
dendro.plot.227413 <- dendro.plot.227413 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227413)

melt.227413 <- melt(mat_tr227413)
head(melt.227413)
colnames(melt.227413) <- c("Targeted", "Bycatch", "TI")
melt.227413$TI_Ranges <- cut(melt.227413$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227413 <- order.dendrogram(dendro.227413)
melt.227413$Targeted <- factor(x = melt.227413$Targeted,
                               levels = melt.227413$Targeted[order.227413],
                               ordered = T)

heatmap.227413 <- ggplot(data=melt.227413, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227413)

grid.newpage()
print(heatmap.227413, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227413, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.4 2014#####

tr2_27414 <- subset(tr2_274, year == 2014)

unique(tr2_27414$year)
unique(tr2_27414$foCatNat)
unique(tr2_27414$area)
tr2_27414$spp2_id <- as.numeric(as.factor(tr2_27414$spp))
spp_tr227414 <- sort(unique(tr2_27414$spp2_id))


ti_tr227414 <- data.frame() 

for(i in spp_tr227414){
  sp1 <- subset(tr2_27414, spp2_id == i) 
  for(j in spp_tr227414){
    t <- data.frame()
    sp2 <- subset(tr2_27414, spp2_id == j) 
    if(i == j){ti_tr227414[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227414[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227414)

row.names(ti_tr227414) <- sort(unique(tr2_27414$spp))
colnames(ti_tr227414) <- sort(unique(tr2_27414$spp))
View(ti_tr227414)

mat_tr227414 <- as.matrix(ti_tr227414)
dim(mat_tr227414)

##HClust with TR2 27.4 2014
dist_tr227414 <- dist(mat_tr227414, method = "euclidean") #distance matrix
hc_tr227414 <- hclust(dist_tr227414, method = "ward.D2")
plot(hc_tr227414)

coph_227414 <- cophenetic(hc_tr227414)
cor(dist_tr227414, coph_227414) #0.67, less than 0.75 so less than optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227414, FUN = hcut, method = "silhouette") #optimal number of clusters is 5


# Ward Hierarchical Clustering
plot(hc_tr227414) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227414, k=5, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227414 <- pvclust(t(mat_tr227414), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227414)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227414, alpha=.95)

dendro.227414 <- as.dendrogram(hc_tr227414, k = 2)
dendro.plot.227414 <- ggdendrogram(data = dendro.227414, rotate = TRUE, k = 2)
dendro.plot.227414 <- dendro.plot.227414 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227414)

melt.227414 <- melt(mat_tr227414)
head(melt.227414)
colnames(melt.227414) <- c("Targeted", "Bycatch", "TI")
melt.227414$TI_Ranges <- cut(melt.227414$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227414 <- order.dendrogram(dendro.227414)
melt.227414$Targeted <- factor(x = melt.227414$Targeted,
                               levels = melt.227414$Targeted[order.227414],
                               ordered = T)

heatmap.227414 <- ggplot(data=melt.227414, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227414)

grid.newpage()
print(heatmap.227414, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227414, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.4 2015####

tr2_27415 <- subset(tr2_274, year == 2015)

unique(tr2_27415$year)
unique(tr2_27415$foCatNat)
unique(tr2_27415$area)
tr2_27415$spp2_id <- as.numeric(as.factor(tr2_27415$spp))
spp_tr227415 <- sort(unique(tr2_27415$spp2_id))


ti_tr227415 <- data.frame() 

for(i in spp_tr227415){
  sp1 <- subset(tr2_27415, spp2_id == i) 
  for(j in spp_tr227415){
    t <- data.frame()
    sp2 <- subset(tr2_27415, spp2_id == j) 
    if(i == j){ti_tr227415[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227415[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227415)

row.names(ti_tr227415) <- sort(unique(tr2_27415$spp))
colnames(ti_tr227415) <- sort(unique(tr2_27415$spp))
View(ti_tr227415)

mat_tr227415 <- as.matrix(ti_tr227415)
dim(mat_tr227415)

##HClust with TR2 27.4 2015
dist_tr227415 <- dist(mat_tr227415, method = "euclidean") #distance matrix
hc_tr227415 <- hclust(dist_tr227415, method = "ward.D2")
plot(hc_tr227415)

coph_227415 <- cophenetic(hc_tr227415)
cor(dist_tr227415, coph_227415) #0.69, less than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227415, FUN = hcut, method = "silhouette") #optimal number of clusters is 7


# Ward Hierarchical Clustering
plot(hc_tr227415) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227415, k=7, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227415 <- pvclust(t(mat_tr227415), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227415)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227415, alpha=.95)

dendro.227415 <- as.dendrogram(hc_tr227415, k = 7)
dendro.plot.227415 <- ggdendrogram(data = dendro.227415, rotate = TRUE, k = 7)
dendro.plot.227415 <- dendro.plot.227415 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227415)

melt.227415 <- melt(mat_tr227415)
head(melt.227415)
colnames(melt.227415) <- c("Targeted", "Bycatch", "TI")
melt.227415$TI_Ranges <- cut(melt.227415$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227415 <- order.dendrogram(dendro.227415)
melt.227415$Targeted <- factor(x = melt.227415$Targeted,
                               levels = melt.227415$Targeted[order.227415],
                               ordered = T)

heatmap.227415 <- ggplot(data=melt.227415, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227415)

grid.newpage()
print(heatmap.227415, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227415, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.4 2016####

tr2_27416 <- subset(tr2_274, year == 2016)

unique(tr2_27416$year)
unique(tr2_27416$foCatNat)
unique(tr2_27416$area)
tr2_27416$spp2_id <- as.numeric(as.factor(tr2_27416$spp))
spp_tr227416 <- sort(unique(tr2_27416$spp2_id))


ti_tr227416 <- data.frame() 

for(i in spp_tr227416){
  sp1 <- subset(tr2_27416, spp2_id == i) 
  for(j in spp_tr227416){
    t <- data.frame()
    sp2 <- subset(tr2_27416, spp2_id == j) 
    if(i == j){ti_tr227416[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227416[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227416)

row.names(ti_tr227416) <- sort(unique(tr2_27416$spp))
colnames(ti_tr227416) <- sort(unique(tr2_27416$spp))
View(ti_tr227416)

mat_tr227416 <- as.matrix(ti_tr227416)
dim(mat_tr227416)

##HClust with TR2 27.4 2016
dist_tr227416 <- dist(mat_tr227416, method = "euclidean") #distance matrix
hc_tr227416 <- hclust(dist_tr227416, method = "ward.D2")
plot(hc_tr227416)

coph_227416 <- cophenetic(hc_tr227416)
cor(dist_tr227416, coph_227416) #0.65, less than 0.75 so that's not optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227416, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_tr227416) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227416, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227416 <- pvclust(t(mat_tr227416), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227416)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227416, alpha=.95)

dendro.227416 <- as.dendrogram(hc_tr227416, k = 4)
dendro.plot.227416 <- ggdendrogram(data = dendro.227416, rotate = TRUE, k = 4)
dendro.plot.227416 <- dendro.plot.227416 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227416)

melt.227416 <- melt(mat_tr227416)
head(melt.227416)
colnames(melt.227416) <- c("Targeted", "Bycatch", "TI")
melt.227416$TI_Ranges <- cut(melt.227416$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227416 <- order.dendrogram(dendro.227416)
melt.227416$Targeted <- factor(x = melt.227416$Targeted,
                               levels = melt.227416$Targeted[order.227416],
                               ordered = T)

heatmap.227416 <- ggplot(data=melt.227416, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227416)

grid.newpage()
print(heatmap.227416, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227416, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.4 2017####

tr2_27417 <- subset(tr2_274, year == 2017)

unique(tr2_27417$year)
unique(tr2_27417$foCatNat)
unique(tr2_27417$area)
tr2_27417$spp2_id <- as.numeric(as.factor(tr2_27417$spp))
spp_tr227417 <- sort(unique(tr2_27417$spp2_id))


ti_tr227417 <- data.frame() 

for(i in spp_tr227417){
  sp1 <- subset(tr2_27417, spp2_id == i) 
  for(j in spp_tr227417){
    t <- data.frame()
    sp2 <- subset(tr2_27417, spp2_id == j) 
    if(i == j){ti_tr227417[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227417[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227417)

row.names(ti_tr227417) <- sort(unique(tr2_27417$spp))
colnames(ti_tr227417) <- sort(unique(tr2_27417$spp))
View(ti_tr227417)

mat_tr227417 <- as.matrix(ti_tr227417)
dim(mat_tr227417)

##HClust with TR1 27.4 2017
dist_tr227417 <- dist(mat_tr227417, method = "euclidean") #distance matrix
hc_tr227417 <- hclust(dist_tr227417, method = "ward.D2")
plot(hc_tr227417)

coph_227417 <- cophenetic(hc_tr227417)
cor(dist_tr227417, coph_227417) #0.77, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227417, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr227417) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227417, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227417 <- pvclust(t(mat_tr227417), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227417)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227417, alpha=.95)

dendro.227417 <- as.dendrogram(hc_tr227417, k = 2)
dendro.plot.227417 <- ggdendrogram(data = dendro.227417, rotate = TRUE, k = 2)
dendro.plot.227417 <- dendro.plot.227417 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227417)

melt.227417 <- melt(mat_tr227417)
head(melt.227417)
colnames(melt.227417) <- c("Targeted", "Bycatch", "TI")
melt.227417$TI_Ranges <- cut(melt.227417$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227417 <- order.dendrogram(dendro.227417)
melt.227417$Targeted <- factor(x = melt.227417$Targeted,
                               levels = melt.227417$Targeted[order.227417],
                               ordered = T)

heatmap.227417 <- ggplot(data=melt.227417, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227417)

grid.newpage()
print(heatmap.227417, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227417, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.4 2018####

tr2_27418 <- subset(tr2_274, year == 2018)

unique(tr2_27418$year)
unique(tr2_27418$foCatNat)
unique(tr2_27418$area)
tr2_27418$spp2_id <- as.numeric(as.factor(tr2_27418$spp))
spp_tr227418 <- sort(unique(tr2_27418$spp2_id))


ti_tr227418 <- data.frame() 

for(i in spp_tr227418){
  sp1 <- subset(tr2_27418, spp2_id == i) 
  for(j in spp_tr227418){
    t <- data.frame()
    sp2 <- subset(tr2_27418, spp2_id == j) 
    if(i == j){ti_tr227418[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227418[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227418)

row.names(ti_tr227418) <- sort(unique(tr2_27418$spp))
colnames(ti_tr227418) <- sort(unique(tr2_27418$spp))
View(ti_tr227418)

mat_tr227418 <- as.matrix(ti_tr227418)
dim(mat_tr227418)

##HClust with TR2 27.4 2017
dist_tr227418 <- dist(mat_tr227418, method = "euclidean") #distance matrix
hc_tr227418 <- hclust(dist_tr227418, method = "ward.D2")
plot(hc_tr227418)

coph_227418 <- cophenetic(hc_tr227418)
cor(dist_tr227418, coph_227418) #0.76, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr227418, FUN = hcut, method = "silhouette") #optimal number of clusters is 5


# Ward Hierarchical Clustering
plot(hc_tr227418) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr227418, k=5, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit227418 <- pvclust(t(mat_tr227418), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit227418)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit227418, alpha=.95)

dendro.227418 <- as.dendrogram(hc_tr227418, k = 2)
dendro.plot.227418 <- ggdendrogram(data = dendro.227418, rotate = TRUE, k = 2)
dendro.plot.227418 <- dendro.plot.227418 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.227418)

melt.227418 <- melt(mat_tr227418)
head(melt.227418)
colnames(melt.227418) <- c("Targeted", "Bycatch", "TI")
melt.227418$TI_Ranges <- cut(melt.227418$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.227418 <- order.dendrogram(dendro.227418)
melt.227418$Targeted <- factor(x = melt.227418$Targeted,
                               levels = melt.227418$Targeted[order.227418],
                               ordered = T)

heatmap.227418 <- ggplot(data=melt.227418, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.227418)

grid.newpage()
print(heatmap.227418, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227418, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###Compiling TR2 27.4 Heatmaps#####
pdf(file = "TR2_274_Heatmaps.pdf")

grid.newpage()
print(heatmap.2274, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2274, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227409, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227409, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227410, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227410, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227411, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227411, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227412, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227412, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227413, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227413, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227414, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227414, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227415, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227415, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227416, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227416, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227417, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227417, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.227418, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.227418, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

dev.off()

###Heatmaps for Averaging#####
ti_tr227409.2 <- data.frame() 
ti_tr227409.2[1:39, 1:39] <- -1

spp_tr227409.2 <- unique(tr2_27409$spp_id)
spp_tr227409.2 <- sort(spp_tr227409.2)

for(i in spp_tr227409.2){
  sp1 <- subset(tr2_27409, spp_id == i) 
  for(j in spp_tr227409.2){
    t <- data.frame()
    sp2 <- subset(tr2_27409, spp_id == j) 
    if(i == j){ti_tr227409.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227409.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227409.2)

row.names(ti_tr227409.2) <- sort(unique(trip382$spp))
colnames(ti_tr227409.2) <- sort(unique(trip382$spp))

mat_tr227409.2 <- as.matrix(ti_tr227409.2)
dim(mat_tr227409.2)
##2010
ti_tr227410.2 <- data.frame() 
ti_tr227410.2[1:39, 1:39] <- -1

spp_tr227410.2 <- unique(tr2_27410$spp_id)
spp_tr227410.2 <- sort(spp_tr227410.2)

for(i in spp_tr227410.2){
  sp1 <- subset(tr2_27410, spp_id == i) 
  for(j in spp_tr227410.2){
    t <- data.frame()
    sp2 <- subset(tr2_27410, spp_id == j) 
    if(i == j){ti_tr227410.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227410.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227410.2)

row.names(ti_tr227410.2) <- sort(unique(trip382$spp))
colnames(ti_tr227410.2) <- sort(unique(trip382$spp))

mat_tr227410.2 <- as.matrix(ti_tr227410.2)
dim(mat_tr227410.2)

##2011
ti_tr227411.2 <- data.frame() 
ti_tr227411.2[1:39, 1:39] <- -1

spp_tr227411.2 <- unique(tr2_27411$spp_id)
spp_tr227411.2 <- sort(spp_tr227411.2)

for(i in spp_tr227411.2){
  sp1 <- subset(tr2_27411, spp_id == i) 
  for(j in spp_tr227411.2){
    t <- data.frame()
    sp2 <- subset(tr2_27411, spp_id == j) 
    if(i == j){ti_tr227411.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227411.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227411.2)

row.names(ti_tr227411.2) <- sort(unique(trip382$spp))
colnames(ti_tr227411.2) <- sort(unique(trip382$spp))

mat_tr227411.2 <- as.matrix(ti_tr227411.2)
dim(mat_tr227411.2)

##2012
ti_tr227412.2 <- data.frame() 
ti_tr227412.2[1:39, 1:39] <- -1

spp_tr227412.2 <- unique(tr2_27412$spp_id)
spp_tr227412.2 <- sort(spp_tr227412.2)

for(i in spp_tr227412.2){
  sp1 <- subset(tr2_27412, spp_id == i) 
  for(j in spp_tr227412.2){
    t <- data.frame()
    sp2 <- subset(tr2_27412, spp_id == j) 
    if(i == j){ti_tr227412.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227412.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227410.2)

row.names(ti_tr227412.2) <- sort(unique(trip382$spp))
colnames(ti_tr227412.2) <- sort(unique(trip382$spp))

mat_tr227412.2 <- as.matrix(ti_tr227412.2)
dim(mat_tr227412.2)

##2013
ti_tr227413.2 <- data.frame() 
ti_tr227413.2[1:39, 1:39] <- -1

spp_tr227413.2 <- unique(tr2_27413$spp_id)
spp_tr227413.2 <- sort(spp_tr227413.2)

for(i in spp_tr227413.2){
  sp1 <- subset(tr2_27413, spp_id == i) 
  for(j in spp_tr227413.2){
    t <- data.frame()
    sp2 <- subset(tr2_27413, spp_id == j) 
    if(i == j){ti_tr227413.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227413.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227413.2)

row.names(ti_tr227413.2) <- sort(unique(trip382$spp))
colnames(ti_tr227413.2) <- sort(unique(trip382$spp))

mat_tr227413.2 <- as.matrix(ti_tr227413.2)
dim(mat_tr227413.2)

##2014
ti_tr227414.2 <- data.frame() 
ti_tr227414.2[1:39, 1:39] <- -1

spp_tr227414.2 <- unique(tr2_27414$spp_id)
spp_tr227414.2 <- sort(spp_tr227414.2)

for(i in spp_tr227414.2){
  sp1 <- subset(tr2_27414, spp_id == i) 
  for(j in spp_tr227414.2){
    t <- data.frame()
    sp2 <- subset(tr2_27414, spp_id == j) 
    if(i == j){ti_tr227414.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227414.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227414.2)

row.names(ti_tr227414.2) <- sort(unique(trip382$spp))
colnames(ti_tr227414.2) <- sort(unique(trip382$spp))

mat_tr227414.2 <- as.matrix(ti_tr227414.2)
dim(mat_tr227414.2)

##2015
ti_tr227415.2 <- data.frame() 
ti_tr227415.2[1:39, 1:39] <- -1

spp_tr227415.2 <- unique(tr2_27415$spp_id)
spp_tr227415.2 <- sort(spp_tr227415.2)

for(i in spp_tr227415.2){
  sp1 <- subset(tr2_27415, spp_id == i) 
  for(j in spp_tr227415.2){
    t <- data.frame()
    sp2 <- subset(tr2_27415, spp_id == j) 
    if(i == j){ti_tr227415.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227415.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227415.2)

row.names(ti_tr227415.2) <- sort(unique(trip382$spp))
colnames(ti_tr227415.2) <- sort(unique(trip382$spp))

mat_tr227415.2 <- as.matrix(ti_tr227415.2)
dim(mat_tr227410.2)

##2016
ti_tr227416.2 <- data.frame() 
ti_tr227416.2[1:39, 1:39] <- -1

spp_tr227416.2 <- unique(tr2_27416$spp_id)
spp_tr227416.2 <- sort(spp_tr227416.2)

for(i in spp_tr227416.2){
  sp1 <- subset(tr2_27416, spp_id == i) 
  for(j in spp_tr227416.2){
    t <- data.frame()
    sp2 <- subset(tr2_27416, spp_id == j) 
    if(i == j){ti_tr227416.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227416.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227416.2)

row.names(ti_tr227416.2) <- sort(unique(trip382$spp))
colnames(ti_tr227416.2) <- sort(unique(trip382$spp))

mat_tr227416.2 <- as.matrix(ti_tr227416.2)
dim(mat_tr227416.2)

##2017
ti_tr227417.2 <- data.frame() 
ti_tr227417.2[1:39, 1:39] <- -1

spp_tr227417.2 <- unique(tr2_27417$spp_id)
spp_tr227417.2 <- sort(spp_tr227417.2)

for(i in spp_tr227417.2){
  sp1 <- subset(tr2_27417, spp_id == i) 
  for(j in spp_tr227417.2){
    t <- data.frame()
    sp2 <- subset(tr2_27417, spp_id == j) 
    if(i == j){ti_tr227417.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227417.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227417.2)

row.names(ti_tr227417.2) <- sort(unique(trip382$spp))
colnames(ti_tr227417.2) <- sort(unique(trip382$spp))

mat_tr227417.2 <- as.matrix(ti_tr227417.2)
dim(mat_tr227417.2)

##2018
ti_tr227418.2 <- data.frame() 
ti_tr227418.2[1:39, 1:39] <- -1

spp_tr227418.2 <- unique(tr2_27418$spp_id)
spp_tr227418.2 <- sort(spp_tr227418.2)

for(i in spp_tr227418.2){
  sp1 <- subset(tr2_27418, spp_id == i) 
  for(j in spp_tr227418.2){
    t <- data.frame()
    sp2 <- subset(tr2_27418, spp_id == j) 
    if(i == j){ti_tr227418.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr227418.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr227418.2)

row.names(ti_tr227418.2) <- sort(unique(trip382$spp))
colnames(ti_tr227418.2) <- sort(unique(trip382$spp))

mat_tr227418.2 <- as.matrix(ti_tr227418.2)
dim(mat_tr227418.2)
##Averaging the Matrices together
matavg_tr2274 <- (mat_tr227409.2 + mat_tr227410.2 + mat_tr227411.2 + 
                    mat_tr227412.2 + mat_tr227413.2 + mat_tr227414.2 + 
                    mat_tr227415.2 + mat_tr227416.2 +  mat_tr227417.2 + mat_tr227418.2)/10


dist_avg2274 <- dist(matavg_tr2274, method = "euclidean") #distance matrix
hc_avg2274 <- hclust(dist_avg2274, method = "ward.D2")
plot(hc_avg2274)

coph_avg2274 <- cophenetic(hc_avg2274)
cor(dist_avg2274, coph_avg2274) #0.91, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(matavg_tr2274, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_avg2274) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_avg2274, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitavg2274 <- pvclust(t(matavg_tr2274), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitavg2274)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitavg2274, alpha=.95)

dendro.avg2274 <- as.dendrogram(hc_avg2274, k = 2)
dendro.plot.avg2274 <- ggdendrogram(data = dendro.avg2274, rotate = TRUE, k = 2)
dendro.plot.avg2274 <- dendro.plot.avg2274 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.avg2274)

melt.avg2274 <- melt(matavg_tr2274)
head(melt.avg2274)
colnames(melt.avg2274) <- c("Targeted", "Bycatch", "TI")
melt.avg2274$TI_Ranges <- cut(melt.avg2274$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.avg2274 <- order.dendrogram(dendro.avg2274)
melt.avg2274$Targeted <- factor(x = melt.avg2274$Targeted,
                                levels = melt.avg2274$Targeted[order.avg1274],
                                ordered = T)

heatmap.avg2274 <- ggplot(data=melt.avg2274, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.avg2274)

grid.newpage()
print(heatmap.avg2274, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.avg2274, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

####Creating Array TR2 27.4####
a_tr2274 <- abind(mat_tr227409.2, mat_tr227410.2, mat_tr227411.2, 
                  mat_tr227412.2, mat_tr227413.2, mat_tr227414.2,
                  mat_tr227415.2, mat_tr227416.2, mat_tr227417.2, mat_tr227418.2, along = 3)
View(a_tr2274[,,10])

melta_tr2274 <- melt(a_tr2274)
head(melta_tr2274)
colnames(melta_tr2274) <- c("Targeted", "Bycatch", "Year", "TIS")

melta_tr2274$TI_Ranges <- cut(melta_tr2274$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
str(melta_tr2274)


melta_tr2274$Area[1:15210] <- "27.4"
melta_tr2274$Area <- as.factor(melta_tr2274$Area)
melta_tr2274$Gear[1:15210] <- "TR2"
melta_tr2274$Gear <- as.factor(melta_tr2274$Gear)

melta_tr2274$Year <- as.factor(melta_tr2274$Year)
melta_tr2274$Year <- revalue(melta_tr2274$Year, c("1"="2009", "2"="2010", "3"="2011", "4"="2012",
                                                  "5"="2013", "6"="2014", "7"="2015", "8"="2016",
                                                  "9"="2017", "10"="2018"))

melta_tr2274 <- na.omit(melta_tr2274)

str(melta_tr2274)

tr2274_by <- with(melta_tr2274, tapply(TIS, INDEX = list(Bycatch, Year), FUN = "mean"))
tr2274_tar <- with(melta_tr2274, tapply(TIS, INDEX = list(Targeted, Year), FUN = "mean"))
tr2274_overall <- with(melta_tr2274, tapply (TIS, INDEX = Year, FUN = "mean"))
tr2274_by <- as.data.frame(tr2274_by)
tr2274_tar <- as.data.frame(tr2274_tar)
tr2274_overall <- as.data.frame(tr2274_overall)

tr2274_overall$var <- with(melta_tr2274, tapply(TIS, INDEX = Year, FUN = "var"))
tr2274_overall$sd <- sqrt(tr2274_overall$var)
tr2274_overall$cv <- 100*((tr2274_overall$sd)/(tr2274_overall$tr2274_overall))

tr2274_by$Overall <- rowMeans(tr2274_by)
tr2274_by$var <- with(melta_tr2274, tapply(TIS, INDEX = list(Bycatch,Year), FUN = "var"))
tr2274_by$sd <- sqrt(tr2274_by$var)
tr2274_by$cv <- 100*((tr2274_by$sd)/(tr1274_by$Overall))
tr2274_by.05 <- tr1274_by[!(tr2274_by$Overall<=5),]
tr2274_by.50 <- tr1274_by[!(tr2274_by$Overall<=50),]
View(tr2274_by.05)
View(tr2274_by.50)

tr2274_tar$Overall <- rowMeans(tr2274_tar)
tr2274_tar$var <- with(melta_tr2274, tapply(TIS, INDEX = list(Targeted,Year), FUN = "var"))
tr2274_tar$sd <- sqrt(tr2274_tar$var)
tr2274_tar$cv <- 100*((tr2274_tar$sd)/(tr1274_tar$Overall))
tr2274_tar.05 <- tr1274_tar[!(tr2274_tar$Overall<=5),]
tr2274_tar.25 <- tr1274_tar[!(tr2274_tar$Overall<=25),]
View(tr2274_tar.05)
View(tr2274_tar.25)


####Centered Moving Average TR2 27.4--EACH YEAR#####
?rollmean
?select
# 
# test <- tr1274_overall %>% select(TIS_year = tr1274_overall) %>% 
#   mutate(TIS_mavg3 = rollmean(TIS_year, k = 3, fill = NA),
#          Year = 2009:2018)
# testg<- test %>% gather(metric, value, TIS_year:TIS_mavg3) %>%
#   ggplot(aes(Year, value, color = metric)) +
#   geom_line()
# print(testg)

t2274 <- data.frame()
for(i in 2009:2018){
  selected <- c(i-1, i, i+1)
  f1 <- subset(melta_tr2274, Year %in% selected)
  if((i == 2009) | (i == 2018)){
    t2274[i-2008, 1] <- NA} 
  else{t2274[i-2008, 1] <- mean(f1$TIS)}
}

print(t2274)

t2274.2 <- t2274 %>% select(TIS_mavg = V1) %>% mutate(TIS_yavg = tr2274_overall$tr2274_overall,
                                                Year = 2009:2018)
t2274.g <- t2274.2 %>% gather(Metric, TIS, TIS_mavg:TIS_yavg) %>%
  ggplot(aes(Year, TIS, color = Metric)) +
  geom_line()
print(t2274.g)



##########################################################################
#####Area 27.6.a, TR1 Subsets#################

tr1_276a <- subset(trip382, gear_id == 3 & area_id == 6)
tr1_276a

unique(tr1_276a$year)
unique(tr1_276a$foCatNat)
unique(tr1_276a$area)
tr1_276a$spp_id <- tr1_276a$spp2_id
tr1_276a$spp2_id <- as.numeric(as.factor(tr1_276a$spp))
spp_tr1276a <- unique(tr1_276a$spp2_id)
spp_tr1276a <- sort(spp_tr1276a)

ti_tr1276a <- data.frame() 

for(i in spp_tr1276a){
  sp1 <- subset(tr1_276a, spp2_id == i) 
  for(j in spp_tr1276a){
    t <- data.frame()
    sp2 <- subset(tr1_276a, spp2_id == j) 
    if(i == j){ti_tr1276a[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$year[k] %in% sp2$year)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a)

row.names(ti_tr1276a) <- sort(unique(tr1_276a$spp))
colnames(ti_tr1276a) <- sort(unique(tr1_276a$spp))
View(ti_tr1276a)

mat_tr1276a <- as.matrix(ti_tr1276a)
##HClust with TR1 27.6.a
dist_tr1276a <- dist(mat_tr1276a, method = "euclidean") #distance matrix
hc_tr1276a<- hclust(dist_tr1276a, method = "ward.D2")
plot(hc_tr1276a)

coph_1276a <- cophenetic(hc_tr1276a)
cor(dist_tr1276a, coph_1276a) #0.73, Less than 0.75, therefore not good


fviz_nbclust(mat_tr1276a, FUN = hcut, method = "silhouette") #optimal number of clusters is 2
?fviz_nbclust


# Ward Hierarchical Clustering
plot(hc_tr1276a) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr1276a, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit276a <- pvclust(t(mat_tr1276a), method.hclust="ward.D2",
                method.dist="euclidean")
plot(fit276a)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit276a, alpha=.95)

dendro.1276a <- as.dendrogram(hc_tr1276a)
dendro.plot.1276a <- ggdendrogram(data = dendro.1276a, rotate = TRUE, k = 2)
dendro.plot.1276a <- dendro.plot.1276a + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a)

melt.1276a <- melt(mat_tr1276a)
head(melt.1276a)
colnames(melt.1276a) <- c("Targeted", "Bycatch", "TI")
melt.1276a$TI_Ranges <- cut(melt.1276a$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a <- order.dendrogram(dendro.1276a)
melt.1276a$Targeted <- factor(x = melt.1276a$Targeted,
                             levels = melt.1276a$Targeted[order.1276a],
                             ordered = T)


heatmap.1276a <- ggplot(data=melt.1276a, 
                       aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a)

grid.newpage()
print(heatmap.1276a, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

###TR1 27.6.a 2009########
tr1_276a09 <- subset(tr1_276a, year == 2009)

unique(tr1_276a09$year)
unique(tr1_276a09$foCatNat)
unique(tr1_276a09$area)
tr1_276a09$spp2_id <- as.numeric(as.factor(tr1_276a09$spp))
spp_tr1276a09 <- unique(tr1_276a09$spp2_id)
spp_tr1276a09 <- sort(spp_tr1276a09)

ti_tr1276a09 <- data.frame() 

for(i in spp_tr1276a09){
  sp1 <- subset(tr1_276a09, spp2_id == i) 
  for(j in spp_tr1276a09){
    t <- data.frame()
    sp2 <- subset(tr1_276a09, spp2_id == j) 
    if(i == j){ti_tr1276a09[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a09[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a09)

row.names(ti_tr1276a09) <- sort(unique(tr1_276a09$spp))
colnames(ti_tr1276a09) <- sort(unique(tr1_276a09$spp))
View(ti_tr1276a09)

mat_tr1276a09 <- as.matrix(ti_tr1276a09)
dim(mat_tr127409)
##HClust with TR1 27.6.a 2009
dist_tr1276a09 <- dist(mat_tr1276a09, method = "euclidean") #distance matrix
hc_tr1276a09<- hclust(dist_tr1276a09, method = "ward.D2")
plot(hc_tr1276a09)

coph_1276a09 <- cophenetic(hc_tr1276a09)
cor(dist_tr1276a09, coph_1276a09) #Greater than 0.75, therefore considered good


fviz_nbclust(mat_tr1276a09, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr1276a09) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276a09, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a09 <- pvclust(t(mat_tr1276a09), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a09)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a09, alpha=.95)

dendro.1276a09 <- as.dendrogram(hc_tr1276a09, k = 2)
dendro.plot.1276a09 <- ggdendrogram(data = dendro.1276a09, rotate = TRUE, k = 2)
dendro.plot.1276a09 <- dendro.plot.1276a09 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a09)

melt.1276a09 <- melt(mat_tr1276a09)
head(melt.1276a09)
colnames(melt.1276a09) <- c("Targeted", "Bycatch", "TI")
melt.1276a09$TI_Ranges <- cut(melt.1276a09$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a09 <- order.dendrogram(dendro.1276a09)
melt.1276a09$Targeted <- factor(x = melt.1276a09$Targeted,
                               levels = melt.1276a09$Targeted[order.1276a09],
                               ordered = T)

heatmap.1276a09 <- ggplot(data=melt.1276a09, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a09)

grid.newpage()
print(heatmap.1276a09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

##
###TR1 27.6.a 2010##########
tr1_276a10 <- subset(tr1_276a, year == 2010)

unique(tr1_276a10$year)
unique(tr1_276a10$foCatNat)
unique(tr1_276a10$area)
tr1_276a10$spp2_id <- as.numeric(as.factor(tr1_276a10$spp))
spp_tr1276a10 <- unique(tr1_276a10$spp2_id)
spp_tr1276a10 <- sort(spp_tr1276a10)

ti_tr1276a10 <- data.frame() 

for(i in spp_tr1276a10){
  sp1 <- subset(tr1_276a10, spp2_id == i) 
  for(j in spp_tr1276a10){
    t <- data.frame()
    sp2 <- subset(tr1_276a10, spp2_id == j) 
    if(i == j){ti_tr1276a10[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a10[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a10)

row.names(ti_tr1276a10) <- sort(unique(tr1_276a10$spp))
colnames(ti_tr1276a10) <- sort(unique(tr1_276a10$spp))
View(ti_tr1276a10)

mat_tr1276a10 <- as.matrix(ti_tr1276a10)
dim(mat_tr1276a10)


##HClust with TR1 27.6.a 2010
dist_tr1276a10 <- dist(mat_tr1276a10, method = "euclidean") #distance matrix
hc_tr1276a10 <- hclust(dist_tr1276a10, method = "ward.D2")
plot(hc_tr1276a10)

coph_1276a10 <- cophenetic(hc_tr1276a10)
cor(dist_tr1276a10, coph_1276a10) #0.87, greater than 0.75, optimal

fviz_nbclust(mat_tr1276a10, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr1276a10) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276a10, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a10 <- pvclust(t(mat_tr1276a10), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a10)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a10, alpha=.95)

dendro.1276a10 <- as.dendrogram(hc_tr1276a10, k = 2)
dendro.plot.1276a10 <- ggdendrogram(data = dendro.1276a10, rotate = TRUE, k = 2)
dendro.plot.1276a10 <- dendro.plot.1276a10 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a10)

melt.1276a10 <- melt(mat_tr1276a10)
head(melt.1276a10)
colnames(melt.1276a10) <- c("Targeted", "Bycatch", "TI")
melt.1276a10$TI_Ranges <- cut(melt.1276a10$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a10 <- order.dendrogram(dendro.1276a10)
melt.1276a10$Targeted <- factor(x = melt.1276a10$Targeted,
                               levels = melt.1276a10$Targeted[order.1276a10],
                               ordered = T)

heatmap.1276a10 <- ggplot(data=melt.1276a10, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a10)

grid.newpage()
print(heatmap.1276a10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))


##

###TR1 27.6.a 2011#######

tr1_276a11 <- subset(tr1_276a, year == 2011)

unique(tr1_276a11$year)
unique(tr1_276a11$foCatNat)
unique(tr1_276a11$area)
tr1_276a11$spp2_id <- as.numeric(as.factor(tr1_276a11$spp))
spp_tr1276a11 <- sort(unique(tr1_276a11$spp2_id))


ti_tr1276a11 <- data.frame() 

for(i in spp_tr1276a11){
  sp1 <- subset(tr1_276a11, spp2_id == i) 
  for(j in spp_tr1276a11){
    t <- data.frame()
    sp2 <- subset(tr1_276a11, spp2_id == j) 
    if(i == j){ti_tr1276a11[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a11[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a11)

row.names(ti_tr1276a11) <- sort(unique(tr1_276a11$spp))
colnames(ti_tr1276a11) <- sort(unique(tr1_276a11$spp))
View(ti_tr1276a11)

mat_tr1276a11 <- as.matrix(ti_tr1276a11)
dim(mat_tr1276a11)

##HClust with TR1 27.4 2011
dist_tr1276a11 <- dist(mat_tr1276a11, method = "euclidean") #distance matrix
hc_tr1276a11 <- hclust(dist_tr1276a11, method = "ward.D2")
plot(hc_tr1276a11)

coph_1276a11 <- cophenetic(hc_tr1276a11)
cor(dist_tr1276a11, coph_1276a11) #0.83, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a11, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_tr1276a11) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276a11, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit1276a11 <- pvclust(t(mat_tr1276a11), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a11)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a11, alpha=.95)

dendro.1276a11 <- as.dendrogram(hc_tr1276a11, k = 4)
dendro.plot.1276a11 <- ggdendrogram(data = dendro.1276a11, rotate = TRUE, k = 4)
dendro.plot.1276a11 <- dendro.plot.1276a11 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a11)

melt.1276a11 <- melt(mat_tr1276a11)
head(melt.1276a11)
colnames(melt.1276a11) <- c("Targeted", "Bycatch", "TI")
melt.1276a11$TI_Ranges <- cut(melt.1276a11$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a11 <- order.dendrogram(dendro.1276a11)
melt.1276a11$Targeted <- factor(x = melt.1276a11$Targeted,
                               levels = melt.1276a11$Targeted[order.1276a11],
                               ordered = T)

heatmap.1276a11 <- ggplot(data=melt.1276a11, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a11)

grid.newpage()
print(heatmap.1276a11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR1 27.6.a 2012######

tr1_276a12 <- subset(tr1_276a, year == 2012)

unique(tr1_276a12$year)
unique(tr1_276a12$foCatNat)
unique(tr1_276a12$area)
tr1_276a12$spp2_id <- as.numeric(as.factor(tr1_276a12$spp))
spp_tr1276a12 <- sort(unique(tr1_276a12$spp2_id))


ti_tr1276a12 <- data.frame() 

for(i in spp_tr1276a12){
  sp1 <- subset(tr1_276a12, spp2_id == i) 
  for(j in spp_tr1276a12){
    t <- data.frame()
    sp2 <- subset(tr1_276a12, spp2_id == j) 
    if(i == j){ti_tr1276a12[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a12[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a12)

row.names(ti_tr1276a12) <- sort(unique(tr1_276a12$spp))
colnames(ti_tr1276a12) <- sort(unique(tr1_276a12$spp))
View(ti_tr1276a12)

mat_tr1276a12 <- as.matrix(ti_tr1276a12)
dim(mat_tr1276a12)

##HClust with TR1 27.6.a 2012
dist_tr1276a12 <- dist(mat_tr1276a12, method = "euclidean") #distance matrix
hc_tr1276a12 <- hclust(dist_tr1276a12, method = "ward.D2")
plot(hc_tr1276a12)

coph_1276a12 <- cophenetic(hc_tr1276a12)
cor(dist_tr1276a12, coph_1276a12) #0.72, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a12, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_tr1276a12) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276a12, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a12 <- pvclust(t(mat_tr1276a12), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a12)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a12, alpha=.95)

dendro.1276a12 <- as.dendrogram(hc_tr1276a12, k = 2)
dendro.plot.1276a12 <- ggdendrogram(data = dendro.1276a12, rotate = TRUE, k = 3)
dendro.plot.1276a12 <- dendro.plot.1276a12 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a12)

melt.1276a12 <- melt(mat_tr1276a12)
head(melt.1276a12)
colnames(melt.1276a12) <- c("Targeted", "Bycatch", "TI")
melt.1276a12$TI_Ranges <- cut(melt.1276a12$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a12 <- order.dendrogram(dendro.1276a12)
melt.1276a12$Targeted <- factor(x = melt.1276a12$Targeted,
                               levels = melt.1276a12$Targeted[order.1276a12],
                               ordered = T)

heatmap.1276a12 <- ggplot(data=melt.1276a12, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a12)

grid.newpage()
print(heatmap.1276a12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR1 27.6.a 2013#####

tr1_276a13 <- subset(tr1_276a, year == 2013)

unique(tr1_276a13$year)
unique(tr1_276a13$foCatNat)
unique(tr1_276a13$area)
tr1_276a13$spp2_id <- as.numeric(as.factor(tr1_276a13$spp))
spp_tr1276a13 <- sort(unique(tr1_276a13$spp2_id))


ti_tr1276a13 <- data.frame() 

for(i in spp_tr1276a13){
  sp1 <- subset(tr1_276a13, spp2_id == i) 
  for(j in spp_tr1276a13){
    t <- data.frame()
    sp2 <- subset(tr1_276a13, spp2_id == j) 
    if(i == j){ti_tr1276a13[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a13[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a13)

row.names(ti_tr1276a13) <- sort(unique(tr1_276a13$spp))
colnames(ti_tr1276a13) <- sort(unique(tr1_276a13$spp))
View(ti_tr1276a13)

mat_tr1276a13 <- as.matrix(ti_tr1276a13)
dim(mat_tr1276a13)

##HClust with TR1 27.6.a 2013
dist_tr1276a13 <- dist(mat_tr1276a13, method = "euclidean") #distance matrix
hc_tr1276a13 <- hclust(dist_tr1276a13, method = "ward.D2")
plot(hc_tr1276a13)

coph_1276a13 <- cophenetic(hc_tr1276a13)
cor(dist_tr1276a13, coph_1276a13) #0.58, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a13, FUN = hcut, method = "silhouette") #optimal number of clusters is 10


# Ward Hierarchical Clustering
plot(hc_tr1276a13) # display dendogram
# draw dendogram with red borders around the 10 clusters 
rect.hclust(hc_tr1276a13, k=10, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a13 <- pvclust(t(mat_tr1276a13), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a13)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a13, alpha=.95)

dendro.1276a13 <- as.dendrogram(hc_tr1276a13, k = 2)
dendro.plot.1276a13 <- ggdendrogram(data = dendro.1276a13, rotate = TRUE, k = 2)
dendro.plot.1276a13 <- dendro.plot.1276a13 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a13)

melt.1276a13 <- melt(mat_tr1276a13)
head(melt.1276a13)
colnames(melt.1276a13) <- c("Targeted", "Bycatch", "TI")
melt.1276a13$TI_Ranges <- cut(melt.1276a13$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a13 <- order.dendrogram(dendro.1276a13)
melt.1276a13$Targeted <- factor(x = melt.1276a13$Targeted,
                               levels = melt.1276a13$Targeted[order.1276a13],
                               ordered = T)

heatmap.1276a13 <- ggplot(data=melt.1276a13, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a13)

grid.newpage()
print(heatmap.1276a13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR1 27.6.a 2014#####

tr1_276a14 <- subset(tr1_276a, year == 2014)

unique(tr1_276a14$year)
unique(tr1_276a14$foCatNat)
unique(tr1_276a14$area)
tr1_276a14$spp2_id <- as.numeric(as.factor(tr1_276a14$spp))
spp_tr1276a14 <- sort(unique(tr1_276a14$spp2_id))


ti_tr1276a14 <- data.frame() 

for(i in spp_tr1276a14){
  sp1 <- subset(tr1_276a14, spp2_id == i) 
  for(j in spp_tr1276a14){
    t <- data.frame()
    sp2 <- subset(tr1_276a14, spp2_id == j) 
    if(i == j){ti_tr1276a14[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a14[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a14)

row.names(ti_tr1276a14) <- sort(unique(tr1_276a14$spp))
colnames(ti_tr1276a14) <- sort(unique(tr1_276a14$spp))
View(ti_tr127414)

mat_tr1276a14 <- as.matrix(ti_tr1276a14)
dim(mat_tr1276a14)

##HClust with TR1 27.6.a 2014
dist_tr1276a14 <- dist(mat_tr1276a14, method = "euclidean") #distance matrix
hc_tr1276a14 <- hclust(dist_tr1276a14, method = "ward.D2")
plot(hc_tr1276a14)

coph_1276a14 <- cophenetic(hc_tr1276a14)
cor(dist_tr1276a14, coph_1276a14) #0.68, less than 0.75 so less than optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a14, FUN = hcut, method = "silhouette") #optimal number of clusters is 5


# Ward Hierarchical Clustering
plot(hc_tr1276a14) # display dendogram
# draw dendogram with red borders around the 5 clusters 
rect.hclust(hc_tr1276a14, k=5, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a14 <- pvclust(t(mat_tr1276a14), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a14)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a14, alpha=.95)

dendro.1276a14 <- as.dendrogram(hc_tr1276a14, k = 5)
dendro.plot.1276a14 <- ggdendrogram(data = dendro.1276a14, rotate = TRUE, k = 5)
dendro.plot.1276a14 <- dendro.plot.1276a14 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a14)

melt.1276a14 <- melt(mat_tr1276a14)
head(melt.1276a14)
colnames(melt.1276a14) <- c("Targeted", "Bycatch", "TI")
melt.1276a14$TI_Ranges <- cut(melt.1276a14$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a14 <- order.dendrogram(dendro.1276a14)
melt.1276a14$Targeted <- factor(x = melt.1276a14$Targeted,
                               levels = melt.1276a14$Targeted[order.1276a14],
                               ordered = T)

heatmap.1276a14 <- ggplot(data=melt.1276a14, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a14)

grid.newpage()
print(heatmap.1276a14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.a 2015####

tr1_276a15 <- subset(tr1_276a, year == 2015)

unique(tr1_276a15$year)
unique(tr1_276a15$foCatNat)
unique(tr1_276a15$area)
tr1_276a15$spp2_id <- as.numeric(as.factor(tr1_276a15$spp))
spp_tr1276a15 <- sort(unique(tr1_276a15$spp2_id))


ti_tr1276a15 <- data.frame() 

for(i in spp_tr1276a15){
  sp1 <- subset(tr1_276a15, spp2_id == i) 
  for(j in spp_tr1276a15){
    t <- data.frame()
    sp2 <- subset(tr1_276a15, spp2_id == j) 
    if(i == j){ti_tr1276a15[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a15[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a15)

row.names(ti_tr1276a15) <- sort(unique(tr1_276a15$spp))
colnames(ti_tr1276a15) <- sort(unique(tr1_276a15$spp))
View(ti_tr1276a15)

mat_tr1276a15 <- as.matrix(ti_tr1276a15)
dim(mat_tr1276a15)

##HClust with TR1 27.6.a 2015
dist_tr1276a15 <- dist(mat_tr1276a15, method = "euclidean") #distance matrix
hc_tr1276a15 <- hclust(dist_tr1276a15, method = "ward.D2")
plot(hc_tr1276a15)

coph_1276a15 <- cophenetic(hc_tr1276a15)
cor(dist_tr1276a15, coph_1276a15) #0.72, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a15, FUN = hcut, method = "silhouette") #optimal number of clusters is 3


# Ward Hierarchical Clustering
plot(hc_tr1276a15) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr1276a15, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a15 <- pvclust(t(mat_tr1276a15), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a15)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a15, alpha=.95)

dendro.1276a15 <- as.dendrogram(hc_tr1276a15, k = 2)
dendro.plot.1276a15 <- ggdendrogram(data = dendro.1276a15, rotate = TRUE, k = 2)
dendro.plot.1276a15 <- dendro.plot.1276a15 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a15)

melt.1276a15 <- melt(mat_tr1276a15)
head(melt.1276a15)
colnames(melt.1276a15) <- c("Targeted", "Bycatch", "TI")
melt.1276a15$TI_Ranges <- cut(melt.1276a15$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a15 <- order.dendrogram(dendro.1276a15)
melt.1276a15$Targeted <- factor(x = melt.1276a15$Targeted,
                               levels = melt.1276a15$Targeted[order.1276a15],
                               ordered = T)

heatmap.1276a15 <- ggplot(data=melt.1276a15, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a15)

grid.newpage()
print(heatmap.1276a15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.a 2016####

tr1_276a16 <- subset(tr1_276a, year == 2016)

unique(tr1_276a16$year)
unique(tr1_276a16$foCatNat)
unique(tr1_276a16$area)
tr1_276a16$spp2_id <- as.numeric(as.factor(tr1_276a16$spp))
spp_tr1276a16 <- sort(unique(tr1_276a16$spp2_id))


ti_tr1276a16 <- data.frame() 
for(i in spp_tr1276a16){
  sp1 <- subset(tr1_276a16, spp2_id == i) 
  for(j in spp_tr1276a16){
    t <- data.frame()
    sp2 <- subset(tr1_276a16, spp2_id == j) 
    if(i == j){ti_tr1276a16[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a16[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a16)

row.names(ti_tr1276a16) <- sort(unique(tr1_276a16$spp))
colnames(ti_tr1276a16) <- sort(unique(tr1_276a16$spp))
View(ti_tr1276a16)

mat_tr1276a16 <- as.matrix(ti_tr1276a16)
dim(mat_tr1276a16)

##HClust with TR1 27.6.a 2016
dist_tr1276a16 <- dist(mat_tr1276a16, method = "euclidean") #distance matrix
hc_tr1276a16 <- hclust(dist_tr1276a16, method = "ward.D2")
plot(hc_tr1276a16)

coph_1276a16 <- cophenetic(hc_tr1276a16)
cor(dist_tr1276a16, coph_1276a16) #0.62, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a16, FUN = hcut, method = "silhouette") #optimal number of clusters is 7


# Ward Hierarchical Clustering
plot(hc_tr1276a16) # display dendogram
# draw dendogram with red borders around the 7 clusters 
rect.hclust(hc_tr1276a16, k=7, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a16 <- pvclust(t(mat_tr1276a16), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a16)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a16, alpha=.95)

dendro.1276a16 <- as.dendrogram(hc_tr1276a16, k = 3)
dendro.plot.1276a16 <- ggdendrogram(data = dendro.1276a16, rotate = TRUE, k = 2)
dendro.plot.1276a16 <- dendro.plot.1276a16 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a16)

melt.1276a16 <- melt(mat_tr1276a16)
head(melt.1276a16)
colnames(melt.1276a16) <- c("Targeted", "Bycatch", "TI")
melt.1276a16$TI_Ranges <- cut(melt.1276a16$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a16 <- order.dendrogram(dendro.1276a16)
melt.1276a16$Targeted <- factor(x = melt.1276a16$Targeted,
                               levels = melt.1276a16$Targeted[order.1276a16],
                               ordered = T)

heatmap.1276a16 <- ggplot(data=melt.1276a16, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a16)

grid.newpage()
print(heatmap.1276a16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.a 2017####

tr1_276a17 <- subset(tr1_276a, year == 2017)

unique(tr1_276a17$year)
unique(tr1_276a17$foCatNat)
unique(tr1_276a17$area)
tr1_276a17$spp2_id <- as.numeric(as.factor(tr1_276a17$spp))
spp_tr1276a17 <- sort(unique(tr1_276a17$spp2_id))


ti_tr1276a17 <- data.frame() 

for(i in spp_tr1276a17){
  sp1 <- subset(tr1_276a17, spp2_id == i) 
  for(j in spp_tr1276a17){
    t <- data.frame()
    sp2 <- subset(tr1_276a17, spp2_id == j) 
    if(i == j){ti_tr1276a17[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a17[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a17)

row.names(ti_tr1276a17) <- sort(unique(tr1_276a17$spp))
colnames(ti_tr1276a17) <- sort(unique(tr1_276a17$spp))
View(ti_tr1276a17)

mat_tr1276a17 <- as.matrix(ti_tr1276a17)
dim(mat_tr1276a17)

##HClust with TR1 27.6.a 2017
dist_tr1276a17 <- dist(mat_tr1276a17, method = "euclidean") #distance matrix
hc_tr1276a17 <- hclust(dist_tr1276a17, method = "ward.D2")
plot(hc_tr1276a17)

coph_1276a17 <- cophenetic(hc_tr1276a17)
cor(dist_tr1276a17, coph_1276a17) #0.74, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a17, FUN = hcut, method = "silhouette") #optimal number of clusters is 10


# Ward Hierarchical Clustering
plot(hc_tr1276a17) # display dendogram 10 clusters 
rect.hclust(hc_tr1276a17, k=10, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a17 <- pvclust(t(mat_tr1276a17), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a17)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a17, alpha=.95)

dendro.1276a17 <- as.dendrogram(hc_tr1276a17, k = 2)
dendro.plot.1276a17 <- ggdendrogram(data = dendro.1276a17, rotate = TRUE, k = 2)
dendro.plot.1276a17 <- dendro.plot.1276a17 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a17)

melt.1276a17 <- melt(mat_tr1276a17)
head(melt.1276a17)
colnames(melt.1276a17) <- c("Targeted", "Bycatch", "TI")
melt.1276a17$TI_Ranges <- cut(melt.1276a17$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a17 <- order.dendrogram(dendro.1276a17)
melt.1276a17$Targeted <- factor(x = melt.1276a17$Targeted,
                               levels = melt.1276a17$Targeted[order.1276a17],
                               ordered = T)

heatmap.1276a17 <- ggplot(data=melt.1276a17, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a17)

grid.newpage()
print(heatmap.1276a17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.a 2018####

tr1_276a18 <- subset(tr1_276a, year == 2018)

unique(tr1_276a18$year)
unique(tr1_276a18$foCatNat)
unique(tr1_276a18$area)
tr1_276a18$spp2_id <- as.numeric(as.factor(tr1_276a18$spp))
spp_tr1276a18 <- sort(unique(tr1_276a18$spp2_id))


ti_tr1276a18 <- data.frame() 

for(i in spp_tr1276a18){
  sp1 <- subset(tr1_276a18, spp2_id == i) 
  for(j in spp_tr1276a18){
    t <- data.frame()
    sp2 <- subset(tr1_276a18, spp2_id == j) 
    if(i == j){ti_tr1276a18[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a18[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a18)

row.names(ti_tr1276a18) <- sort(unique(tr1_276a18$spp))
colnames(ti_tr1276a18) <- sort(unique(tr1_276a18$spp))
View(ti_tr1276a18)

mat_tr1276a18 <- as.matrix(ti_tr1276a18)
dim(mat_tr1276a18)

##HClust with TR1 27.6.a 2018
dist_tr1276a18 <- dist(mat_tr1276a18, method = "euclidean") #distance matrix
hc_tr1276a18 <- hclust(dist_tr1276a18, method = "ward.D2")
plot(hc_tr1276a18)

coph_1276a18 <- cophenetic(hc_tr1276a18)
cor(dist_tr1276a18, coph_1276a18) #0.88, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276a18, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr1276a18) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276a18, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276a18 <- pvclust(t(mat_tr1276a18), method.hclust="ward.D2",
                     method.dist="euclidean")
plot(fit1276a18)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276a18, alpha=.95)

dendro.1276a18 <- as.dendrogram(hc_tr1276a18, k = 2)
dendro.plot.1276a18 <- ggdendrogram(data = dendro.1276a18, rotate = TRUE, k = 2)
dendro.plot.1276a18 <- dendro.plot.1276a18 + theme(axis.text.y = element_blank(), 
                                                 axis.text.x.bottom =  element_blank())

print(dendro.plot.1276a18)

melt.1276a18 <- melt(mat_tr1276a18)
head(melt.1276a18)
colnames(melt.1276a18) <- c("Targeted", "Bycatch", "TI")
melt.1276a18$TI_Ranges <- cut(melt.1276a18$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276a18 <- order.dendrogram(dendro.1276a18)
melt.1276a18$Targeted <- factor(x = melt.1276a18$Targeted,
                               levels = melt.1276a18$Targeted[order.1276a18],
                               ordered = T)

heatmap.1276a18 <- ggplot(data=melt.1276a18, 
                         aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276a18)

grid.newpage()
print(heatmap.1276a18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###Compiling TR1 27.6.a Heatmaps#####
pdf(file = "TR1_276a_Heatmaps.pdf")

grid.newpage()
print(heatmap.1276a, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276a18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276a18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

dev.off()
###Heatmaps for Averaging#####
ti_tr1276a09.2 <- data.frame() 
ti_tr1276a09.2[1:39, 1:39] <- -1

spp_tr1276a09.2 <- unique(tr1_276a09$spp_id)
spp_tr1276a09.2 <- sort(spp_tr1276a09.2)

for(i in spp_tr1276a09.2){
  sp1 <- subset(tr1_276a09, spp_id == i) 
  for(j in spp_tr1276a09.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a09, spp_id == j) 
    if(i == j){ti_tr1276a09.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a09.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a09.2)

row.names(ti_tr1276a09.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a09.2) <- sort(unique(trip382$spp))

mat_tr1276a09.2 <- as.matrix(ti_tr1276a09.2)
dim(mat_tr1276a09.2)
##2010
ti_tr1276a10.2 <- data.frame() 
ti_tr1276a10.2[1:39, 1:39] <- -1

spp_tr1276a10.2 <- unique(tr1_276a10$spp_id)
spp_tr1276a10.2 <- sort(spp_tr1276a10.2)

for(i in spp_tr1276a10.2){
  sp1 <- subset(tr1_276a10, spp_id == i) 
  for(j in spp_tr1276a10.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a10, spp_id == j) 
    if(i == j){ti_tr1276a10.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a10.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a10.2)

row.names(ti_tr1276a10.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a10.2) <- sort(unique(trip382$spp))

mat_tr1276a10.2 <- as.matrix(ti_tr1276a10.2)
dim(mat_tr1276a10.2)

##2011
ti_tr1276a11.2 <- data.frame() 
ti_tr1276a11.2[1:39, 1:39] <- -1

spp_tr1276a11.2 <- unique(tr1_276a11$spp_id)
spp_tr1276a11.2 <- sort(spp_tr1276a11.2)

for(i in spp_tr1276a11.2){
  sp1 <- subset(tr1_276a11, spp_id == i) 
  for(j in spp_tr1276a11.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a11, spp_id == j) 
    if(i == j){ti_tr1276a11.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a11.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a11.2)

row.names(ti_tr1276a11.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a11.2) <- sort(unique(trip382$spp))

mat_tr1276a11.2 <- as.matrix(ti_tr1276a11.2)
dim(mat_tr1276a11.2)

##2012
ti_tr1276a12.2 <- data.frame() 
ti_tr1276a12.2[1:39, 1:39] <- -1

spp_tr1276a12.2 <- unique(tr1_276a12$spp_id)
spp_tr1276a12.2 <- sort(spp_tr1276a12.2)

for(i in spp_tr1276a12.2){
  sp1 <- subset(tr1_276a12, spp_id == i) 
  for(j in spp_tr1276a12.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a12, spp_id == j) 
    if(i == j){ti_tr1276a12.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a12.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a10.2)

row.names(ti_tr1276a12.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a12.2) <- sort(unique(trip382$spp))

mat_tr1276a12.2 <- as.matrix(ti_tr1276a12.2)
dim(mat_tr1276a12.2)

##2013
ti_tr1276a13.2 <- data.frame() 
ti_tr1276a13.2[1:39, 1:39] <- -1

spp_tr1276a13.2 <- unique(tr1_276a13$spp_id)
spp_tr1276a13.2 <- sort(spp_tr1276a13.2)

for(i in spp_tr1276a13.2){
  sp1 <- subset(tr1_276a13, spp_id == i) 
  for(j in spp_tr1276a13.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a13, spp_id == j) 
    if(i == j){ti_tr1276a13.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a13.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a13.2)

row.names(ti_tr1276a13.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a13.2) <- sort(unique(trip382$spp))

mat_tr1276a13.2 <- as.matrix(ti_tr1276a13.2)
dim(mat_tr1276a13.2)

##2014
ti_tr1276a14.2 <- data.frame() 
ti_tr1276a14.2[1:39, 1:39] <- -1

spp_tr1276a14.2 <- unique(tr1_276a14$spp_id)
spp_tr1276a14.2 <- sort(spp_tr1276a14.2)

for(i in spp_tr1276a14.2){
  sp1 <- subset(tr1_276a14, spp_id == i) 
  for(j in spp_tr1276a14.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a14, spp_id == j) 
    if(i == j){ti_tr1276a14.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a14.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a14.2)

row.names(ti_tr1276a14.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a14.2) <- sort(unique(trip382$spp))

mat_tr1276a14.2 <- as.matrix(ti_tr1276a14.2)
dim(mat_tr1276a14.2)

##2015
ti_tr1276a15.2 <- data.frame() 
ti_tr1276a15.2[1:39, 1:39] <- -1

spp_tr1276a15.2 <- unique(tr1_276a15$spp_id)
spp_tr1276a15.2 <- sort(spp_tr1276a15.2)

for(i in spp_tr1276a15.2){
  sp1 <- subset(tr1_276a15, spp_id == i) 
  for(j in spp_tr1276a15.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a15, spp_id == j) 
    if(i == j){ti_tr1276a15.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a15.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a15.2)

row.names(ti_tr1276a15.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a15.2) <- sort(unique(trip382$spp))

mat_tr1276a15.2 <- as.matrix(ti_tr1276a15.2)
dim(mat_tr1276a15.2)

##2016
ti_tr1276a16.2 <- data.frame() 
ti_tr1276a16.2[1:39, 1:39] <- -1

spp_tr1276a16.2 <- unique(tr1_276a16$spp_id)
spp_tr1276a16.2 <- sort(spp_tr1276a16.2)

for(i in spp_tr1276a16.2){
  sp1 <- subset(tr1_276a16, spp_id == i) 
  for(j in spp_tr1276a16.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a16, spp_id == j) 
    if(i == j){ti_tr1276a16.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a16.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a16.2)

row.names(ti_tr1276a16.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a16.2) <- sort(unique(trip382$spp))

mat_tr1276a16.2 <- as.matrix(ti_tr1276a16.2)
dim(mat_tr1276a16.2)

##2017
ti_tr1276a17.2 <- data.frame() 
ti_tr1276a17.2[1:39, 1:39] <- -1

spp_tr1276a17.2 <- unique(tr1_276a17$spp_id)
spp_tr1276a17.2 <- sort(spp_tr1276a17.2)

for(i in spp_tr1276a17.2){
  sp1 <- subset(tr1_276a17, spp_id == i) 
  for(j in spp_tr1276a17.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a17, spp_id == j) 
    if(i == j){ti_tr1276a17.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a17.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a17.2)

row.names(ti_tr1276a17.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a17.2) <- sort(unique(trip382$spp))

mat_tr1276a17.2 <- as.matrix(ti_tr1276a17.2)
dim(mat_tr1276a17.2)

##2018
ti_tr1276a18.2 <- data.frame() 
ti_tr1276a18.2[1:39, 1:39] <- -1

spp_tr1276a18.2 <- unique(tr1_276a18$spp_id)
spp_tr1276a18.2 <- sort(spp_tr1276a18.2)

for(i in spp_tr1276a18.2){
  sp1 <- subset(tr1_276a18, spp_id == i) 
  for(j in spp_tr1276a18.2){
    t <- data.frame()
    sp2 <- subset(tr1_276a18, spp_id == j) 
    if(i == j){ti_tr1276a18.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276a18.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276a18.2)

row.names(ti_tr1276a18.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276a18.2) <- sort(unique(trip382$spp))

mat_tr1276a18.2 <- as.matrix(ti_tr1276a18.2)
dim(mat_tr1276a18.2)
##Averaging the Matrices together
matavg_tr1276a <- (mat_tr1276a09.2 + mat_tr1276a10.2 + mat_tr1276a11.2 + 
                    mat_tr1276a12.2 + mat_tr1276a13.2 + mat_tr1276a14.2 + 
                    mat_tr1276a15.2 + mat_tr1276a16.2 +  mat_tr1276a17.2 + mat_tr1276a18.2)/10


dist_avg1276a <- dist(matavg_tr1276a, method = "euclidean") #distance matrix
hc_avg1276a <- hclust(dist_avg1276a, method = "ward.D2")
plot(hc_avg1276a)

coph_avg1276a <- cophenetic(hc_avg1276a)
cor(dist_avg1276a, coph_avg1276a) #0.85, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(matavg_tr1276a, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_avg1276a) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_avg1276a, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitavg1276a <- pvclust(t(matavg_tr1276a), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitavg1276a)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitavg1276a, alpha=.95)

dendro.avg1276a <- as.dendrogram(hc_avg1276a, k = 2)
dendro.plot.avg1276a <- ggdendrogram(data = dendro.avg1276a, rotate = TRUE, k = 2)
dendro.plot.avg1276a <- dendro.plot.avg1276a + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.avg1276a)

melt.avg1276a <- melt(matavg_tr1276a)
head(melt.avg1276a)
colnames(melt.avg1276a) <- c("Targeted", "Bycatch", "TI")
melt.avg1276a$TI_Ranges <- cut(melt.avg1276a$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.avg1276a <- order.dendrogram(dendro.avg1276a)
melt.avg1276a$Targeted <- factor(x = melt.avg1276a$Targeted,
                                levels = melt.avg1276a$Targeted[order.avg1276a],
                                ordered = T)

heatmap.avg1276a <- ggplot(data=melt.avg1276a, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.avg1276a)

grid.newpage()
print(heatmap.avg1276a, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.avg1276a, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

####Creating Array TR1 27.6.a####
a_tr1276a <- abind(mat_tr1276a09.2, mat_tr1276a10.2, mat_tr1276a11.2, 
                  mat_tr1276a12.2, mat_tr1276a13.2, mat_tr1276a14.2,
                  mat_tr1276a15.2, mat_tr1276a16.2, mat_tr1276a17.2, mat_tr1276a18.2, along = 3)
View(a_tr1276a[,,10])

melta_tr1276a <- melt(a_tr1276a)
head(melta_tr1276a)
colnames(melta_tr1276a) <- c("Targeted", "Bycatch", "Year", "TIS")

melta_tr1276a$TI_Ranges <- cut(melta_tr1276a$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
str(melta_tr1276a)


melta_tr1276a$Area[1:15210] <- "27.6.a"
melta_tr1276a$Area <- as.factor(melta_tr1276a$Area)
melta_tr1276a$Gear[1:15210] <- "TR1"
melta_tr1276a$Gear <- as.factor(melta_tr1276a$Gear)

melta_tr1276a$Year <- as.factor(melta_tr1276a$Year)
melta_tr1276a$Year <- revalue(melta_tr1276a$Year, c("1"="2009", "2"="2010", "3"="2011", "4"="2012",
                                                  "5"="2013", "6"="2014", "7"="2015", "8"="2016",
                                                  "9"="2017", "10"="2018"))

melta_tr1276a <- na.omit(melta_tr1276a) #removing absent species that are marked as NA in TI_Ranges

tr1276a_by <- with(melta_tr1276a, tapply(TIS, INDEX = list(Bycatch, Year), FUN = "mean"))
tr1276a_tar <- with(melta_tr1276a, tapply(TIS, INDEX = list(Targeted, Year), FUN = "mean"))
tr1276a_overall <- with(melta_tr1276a, tapply (TIS, INDEX = Year, FUN = "mean"))
tr1276a_by <- as.data.frame(tr1276a_by)
tr1276a_tar <- as.data.frame(tr1276a_tar)
tr1276a_overall <- as.data.frame(tr1276a_overall)

tr1276a_overall$var <- with(melta_tr1276a, tapply(TIS, INDEX = Year, FUN = "var"))
tr1276a_overall$sd <- sqrt(tr1276a_overall$var)
tr1276a_overall$cv <- 100*((tr1276a_overall$sd)/(tr1274_overall$tr1274_overall))

tr1276a_by$Overall <- rowMeans(tr1276a_by)
tr1276a_by$var <- with(melta_tr1276a, tapply(TIS, INDEX = list(Bycatch,Year), FUN = "var"))
tr1276a_by$sd <- sqrt(tr1276a_by$var)
tr1276a_by$cv <- 100*((tr1276a_by$sd)/(tr1276a_by$Overall))
tr1276a_by.05 <- tr1276a_by[!(tr1276a_by$Overall<=5),]
tr1276a_by.50 <- tr1276a_by[!(tr1276a_by$Overall<=50),]
View(tr1276a_by.05)
View(tr1276a_by.50)

tr1276a_tar$Overall <- rowMeans(tr1276a_tar)
tr1276a_tar$var <- with(melta_tr1276a, tapply(TIS, INDEX = list(Targeted,Year), FUN = "var"))
tr1276a_tar$sd <- sqrt(tr1276a_tar$var)
tr1276a_tar$cv <- 100*((tr1276a_tar$sd)/(tr1276a_tar$Overall))
tr1276a_tar.05 <- tr1276a_tar[!(tr1276a_tar$Overall<=5),]
tr1276a_tar.25 <- tr1276a_tar[!(tr1276a_tar$Overall<=25),]
View(tr1276a_tar.05)
View(tr1276a_tar.25)


####Centered Moving Average TR1 27.6.a--EACH YEAR#####
?rollmean
?select
# 
# test <- tr1274_overall %>% select(TIS_year = tr1274_overall) %>% 
#   mutate(TIS_mavg3 = rollmean(TIS_year, k = 3, fill = NA),
#          Year = 2009:2018)
# testg<- test %>% gather(metric, value, TIS_year:TIS_mavg3) %>%
#   ggplot(aes(Year, value, color = metric)) +
#   geom_line()
# print(testg)

t1276a <- data.frame()
for(i in 2009:2018){
  selected <- c(i-1, i, i+1)
  f1 <- subset(melta_tr1276a, Year %in% selected)
  if((i == 2009)| (i==2018)){
    t1276a[i-2008, 1] <- NA} 
  else{t1276a[i-2008,1] <- mean(f1$TIS)}
}
print(t1276a)

t1276a.2 <- t1276a %>% select(TIS_mavg = V1) %>% mutate(TIS_yavg = tr1276a_overall$tr1276a_overall,
                                                      Year = 2009:2018)
t1276a.g <- t1276a.2 %>% gather(Metric, TIS, TIS_mavg:TIS_yavg) %>%
  ggplot(aes(Year, TIS, color = Metric)) +
  geom_line()
print(t1276a.g)

##########################################################################
#####Area 27.6.a, TR2 Subsets#################

tr2_276a <- subset(trip382, gear_id == 4 & area_id == 6)
tr2_276a

unique(tr2_276a$year)
unique(tr2_276a$foCatNat)
unique(tr2_276a$area)
tr2_276a$spp_id <- tr2_276a$spp2_id
tr2_276a$spp2_id <- as.numeric(as.factor(tr2_276a$spp))
spp_tr2276a <- unique(tr2_276a$spp2_id)
spp_tr2276a <- sort(spp_tr2276a)

ti_tr2276a <- data.frame() 

for(i in spp_tr2276a){
  sp1 <- subset(tr2_276a, spp2_id == i) 
  for(j in spp_tr2276a){
    t <- data.frame()
    sp2 <- subset(tr2_276a, spp2_id == j) 
    if(i == j){ti_tr2276a[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$year[k] %in% sp2$year)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a)

row.names(ti_tr2276a) <- sort(unique(tr2_276a$spp))
colnames(ti_tr2276a) <- sort(unique(tr2_276a$spp))
View(ti_tr2276a)

mat_tr2276a <- as.matrix(ti_tr2276a)
dim(mat_tr2276a)

##HClust with TR2 27.6.a
dist_tr2276a <- dist(mat_tr2276a, method = "euclidean") #distance matrix
hc_tr2276a<- hclust(dist_tr2276a, method = "ward.D2")
plot(hc_tr2276a)

coph_2276a <- cophenetic(hc_tr2276a)
cor(dist_tr2276a, coph_2276a) #0.53, Less than 0.75, therefore not good


fviz_nbclust(mat_tr2276a, FUN = hcut, method = "silhouette") #optimal number of clusters is 6
?fviz_nbclust


# Ward Hierarchical Clustering
plot(hc_tr2276a) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr2276a, k=6, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit276a2 <- pvclust(t(mat_tr2276a), method.hclust="ward.D2",
                   method.dist="euclidean")
plot(fit276a2)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit276a2, alpha=.95)

dendro.2276a <- as.dendrogram(hc_tr2276a)
dendro.plot.2276a <- ggdendrogram(data = dendro.2276a, rotate = TRUE, k = 2)
dendro.plot.2276a <- dendro.plot.2276a + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a)

melt.2276a <- melt(mat_tr2276a)
head(melt.2276a)
colnames(melt.2276a) <- c("Targeted", "Bycatch", "TI")
melt.2276a$TI_Ranges <- cut(melt.2276a$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a <- order.dendrogram(dendro.2276a)
melt.2276a$Targeted <- factor(x = melt.2276a$Targeted,
                              levels = melt.2276a$Targeted[order.2276a],
                              ordered = T)


heatmap.2276a <- ggplot(data=melt.2276a, 
                        aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a)

grid.newpage()
print(heatmap.2276a, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

###TR2 27.6.a 2009########
tr2_276a09 <- subset(tr2_276a, year == 2009)

unique(tr2_276a09$year)
unique(tr2_276a09$foCatNat)
unique(tr2_276a09$area)
tr2_276a09$spp2_id <- as.numeric(as.factor(tr2_276a09$spp))
spp_tr2276a09 <- unique(tr2_276a09$spp2_id)
spp_tr2276a09 <- sort(spp_tr2276a09)

ti_tr2276a09 <- data.frame() 

for(i in spp_tr2276a09){
  sp1 <- subset(tr2_276a09, spp2_id == i) 
  for(j in spp_tr2276a09){
    t <- data.frame()
    sp2 <- subset(tr2_276a09, spp2_id == j) 
    if(i == j){ti_tr2276a09[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a09[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a09)

row.names(ti_tr2276a09) <- sort(unique(tr2_276a09$spp))
colnames(ti_tr2276a09) <- sort(unique(tr2_276a09$spp))
View(ti_tr2276a09)

mat_tr2276a09 <- as.matrix(ti_tr2276a09)
dim(mat_tr2276a09)

##HClust with TR2 27.6.a 2009
dist_tr2276a09 <- dist(mat_tr2276a09, method = "euclidean") #distance matrix
hc_tr2276a09<- hclust(dist_tr2276a09, method = "ward.D2")
plot(hc_tr2276a09)

coph_2276a09 <- cophenetic(hc_tr2276a09)
cor(dist_tr2276a09, coph_2276a09) #0.68, less than 0.75, therefore not good


fviz_nbclust(mat_tr2276a09, FUN = hcut, method = "silhouette") #optimal number of clusters is 6


# Ward Hierarchical Clustering
plot(hc_tr2276a09) # display dendogram
# draw dendogram with red borders around the 6 clusters 
rect.hclust(hc_tr2276a09, k=6, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a09 <- pvclust(t(mat_tr2276a09), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a09)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a09, alpha=.95)

dendro.2276a09 <- as.dendrogram(hc_tr2276a09, k = 6)
dendro.plot.2276a09 <- ggdendrogram(data = dendro.2276a09, rotate = TRUE, k = 6)
dendro.plot.2276a09 <- dendro.plot.2276a09 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a09)

melt.2276a09 <- melt(mat_tr2276a09)
head(melt.2276a09)
colnames(melt.2276a09) <- c("Targeted", "Bycatch", "TI")
melt.2276a09$TI_Ranges <- cut(melt.2276a09$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a09 <- order.dendrogram(dendro.2276a09)
melt.2276a09$Targeted <- factor(x = melt.2276a09$Targeted,
                                levels = melt.2276a09$Targeted[order.2276a09],
                                ordered = T)

heatmap.2276a09 <- ggplot(data=melt.2276a09, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a09)

grid.newpage()
print(heatmap.2276a09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

##
###TR2 27.6.a 2010##########
tr2_276a10 <- subset(tr2_276a, year == 2010)

unique(tr2_276a10$year)
unique(tr2_276a10$foCatNat)
unique(tr2_276a10$area)
tr2_276a10$spp2_id <- as.numeric(as.factor(tr2_276a10$spp))
spp_tr2276a10 <- unique(tr2_276a10$spp2_id)
spp_tr2276a10 <- sort(spp_tr2276a10)

ti_tr2276a10 <- data.frame() 

for(i in spp_tr2276a10){
  sp1 <- subset(tr2_276a10, spp2_id == i) 
  for(j in spp_tr2276a10){
    t <- data.frame()
    sp2 <- subset(tr2_276a10, spp2_id == j) 
    if(i == j){ti_tr2276a10[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a10[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a10)

row.names(ti_tr2276a10) <- sort(unique(tr2_276a10$spp))
colnames(ti_tr2276a10) <- sort(unique(tr2_276a10$spp))
View(ti_tr2276a10)

mat_tr2276a10 <- as.matrix(ti_tr2276a10)
dim(mat_tr2276a10)


##HClust with TR2 27.6.a 2010
dist_tr2276a10 <- dist(mat_tr2276a10, method = "euclidean") #distance matrix
hc_tr2276a10 <- hclust(dist_tr2276a10, method = "ward.D2")
plot(hc_tr2276a10)

coph_2276a10 <- cophenetic(hc_tr2276a10)
cor(dist_tr2276a10, coph_2276a10) #0.77, greater than 0.75, optimal

fviz_nbclust(mat_tr2276a10, FUN = hcut, method = "silhouette") #optimal number of clusters is 5


# Ward Hierarchical Clustering
plot(hc_tr2276a10) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr2276a10, k=5, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a10 <- pvclust(t(mat_tr2276a10), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a10)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a10, alpha=.95)

dendro.2276a10 <- as.dendrogram(hc_tr2276a10)
dendro.plot.2276a10 <- ggdendrogram(data = dendro.2276a10, rotate = TRUE)
dendro.plot.2276a10 <- dendro.plot.2276a10 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a10)

melt.2276a10 <- melt(mat_tr2276a10)
head(melt.2276a10)
colnames(melt.2276a10) <- c("Targeted", "Bycatch", "TI")
melt.2276a10$TI_Ranges <- cut(melt.2276a10$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a10 <- order.dendrogram(dendro.2276a10)
melt.2276a10$Targeted <- factor(x = melt.2276a10$Targeted,
                                levels = melt.2276a10$Targeted[order.2276a10],
                                ordered = T)

heatmap.2276a10 <- ggplot(data=melt.2276a10, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a10)

grid.newpage()
print(heatmap.2276a10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))


##

###TR2 27.6.a 2011#######

tr2_276a11 <- subset(tr2_276a, year == 2011)

unique(tr2_276a11$year)
unique(tr2_276a11$foCatNat)
unique(tr2_276a11$area)
tr2_276a11$spp2_id <- as.numeric(as.factor(tr2_276a11$spp))
spp_tr2276a11 <- sort(unique(tr2_276a11$spp2_id))


ti_tr2276a11 <- data.frame() 

for(i in spp_tr2276a11){
  sp1 <- subset(tr2_276a11, spp2_id == i) 
  for(j in spp_tr2276a11){
    t <- data.frame()
    sp2 <- subset(tr2_276a11, spp2_id == j) 
    if(i == j){ti_tr2276a11[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a11[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a11)

row.names(ti_tr2276a11) <- sort(unique(tr2_276a11$spp))
colnames(ti_tr2276a11) <- sort(unique(tr2_276a11$spp))
View(ti_tr2276a11)

mat_tr2276a11 <- as.matrix(ti_tr2276a11)
dim(mat_tr2276a11)

##HClust with TR2 27.6.a 2011
dist_tr2276a11 <- dist(mat_tr2276a11, method = "euclidean") #distance matrix
hc_tr2276a11 <- hclust(dist_tr2276a11, method = "ward.D2")
plot(hc_tr2276a11)

coph_2276a11 <- cophenetic(hc_tr2276a11)
cor(dist_tr2276a11, coph_2276a11) #0.65, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a11, FUN = hcut, method = "silhouette") #optimal number of clusters is 5


# Ward Hierarchical Clustering
plot(hc_tr2276a11) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr2276a11, k=5, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit2276a11 <- pvclust(t(mat_tr2276a11), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a11)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a11, alpha=.95)

dendro.2276a11 <- as.dendrogram(hc_tr2276a11)
dendro.plot.2276a11 <- ggdendrogram(data = dendro.2276a11, rotate = TRUE)
dendro.plot.2276a11 <- dendro.plot.2276a11 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a11)

melt.2276a11 <- melt(mat_tr2276a11)
head(melt.2276a11)
colnames(melt.2276a11) <- c("Targeted", "Bycatch", "TI")
melt.2276a11$TI_Ranges <- cut(melt.2276a11$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a11 <- order.dendrogram(dendro.2276a11)
melt.2276a11$Targeted <- factor(x = melt.2276a11$Targeted,
                                levels = melt.2276a11$Targeted[order.2276a11],
                                ordered = T)

heatmap.2276a11 <- ggplot(data=melt.2276a11, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a11)

grid.newpage()
print(heatmap.2276a11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.6.a 2012######

tr2_276a12 <- subset(tr2_276a, year == 2012)

unique(tr2_276a12$year)
unique(tr2_276a12$foCatNat)
unique(tr2_276a12$area)
tr2_276a12$spp2_id <- as.numeric(as.factor(tr2_276a12$spp))
spp_tr2276a12 <- sort(unique(tr2_276a12$spp2_id))


ti_tr2276a12 <- data.frame() 

for(i in spp_tr2276a12){
  sp1 <- subset(tr2_276a12, spp2_id == i) 
  for(j in spp_tr2276a12){
    t <- data.frame()
    sp2 <- subset(tr2_276a12, spp2_id == j) 
    if(i == j){ti_tr2276a12[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a12[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a12)

row.names(ti_tr2276a12) <- sort(unique(tr2_276a12$spp))
colnames(ti_tr2276a12) <- sort(unique(tr2_276a12$spp))
View(ti_tr2276a12)

mat_tr2276a12 <- as.matrix(ti_tr2276a12)
dim(mat_tr2276a12)

##HClust with TR2 27.6.a 2012
dist_tr2276a12 <- dist(mat_tr2276a12, method = "euclidean") #distance matrix
hc_tr2276a12 <- hclust(dist_tr2276a12, method = "ward.D2")
plot(hc_tr2276a12)

coph_2276a12 <- cophenetic(hc_tr2276a12)
cor(dist_tr2276a12, coph_2276a12) #0.76, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a12, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr2276a12) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr2276a12, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a12 <- pvclust(t(mat_tr2276a12), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a12)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a12, alpha=.95)

dendro.2276a12 <- as.dendrogram(hc_tr2276a12)
dendro.plot.2276a12 <- ggdendrogram(data = dendro.2276a12, rotate = TRUE)
dendro.plot.2276a12 <- dendro.plot.2276a12 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a12)

melt.2276a12 <- melt(mat_tr2276a12)
head(melt.2276a12)
colnames(melt.2276a12) <- c("Targeted", "Bycatch", "TI")
melt.2276a12$TI_Ranges <- cut(melt.2276a12$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a12 <- order.dendrogram(dendro.2276a12)
melt.2276a12$Targeted <- factor(x = melt.2276a12$Targeted,
                                levels = melt.2276a12$Targeted[order.2276a12],
                                ordered = T)

heatmap.2276a12 <- ggplot(data=melt.2276a12, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a12)

grid.newpage()
print(heatmap.2276a12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.6.a 2013#####

tr2_276a13 <- subset(tr2_276a, year == 2013)

unique(tr2_276a13$year)
unique(tr2_276a13$foCatNat)
unique(tr2_276a13$area)
tr2_276a13$spp2_id <- as.numeric(as.factor(tr2_276a13$spp))
spp_tr2276a13 <- sort(unique(tr2_276a13$spp2_id))


ti_tr2276a13 <- data.frame() 

for(i in spp_tr2276a13){
  sp1 <- subset(tr2_276a13, spp2_id == i) 
  for(j in spp_tr2276a13){
    t <- data.frame()
    sp2 <- subset(tr2_276a13, spp2_id == j) 
    if(i == j){ti_tr2276a13[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a13[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a13)

row.names(ti_tr2276a13) <- sort(unique(tr2_276a13$spp))
colnames(ti_tr2276a13) <- sort(unique(tr2_276a13$spp))
View(ti_tr2276a13)

mat_tr2276a13 <- as.matrix(ti_tr2276a13)
dim(mat_tr2276a13)

##HClust with TR2 27.6.a 2013
dist_tr2276a13 <- dist(mat_tr2276a13, method = "euclidean") #distance matrix
hc_tr2276a13 <- hclust(dist_tr2276a13, method = "ward.D2")
plot(hc_tr2276a13)

coph_2276a13 <- cophenetic(hc_tr2276a13)
cor(dist_tr2276a13, coph_2276a13) #0.76, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a13, FUN = hcut, method = "silhouette") #optimal number of clusters is 8


# Ward Hierarchical Clustering
plot(hc_tr2276a13) # display dendogram
# draw dendogram with red borders around the clusters 
rect.hclust(hc_tr2276a13, k=8, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a13 <- pvclust(t(mat_tr2276a13), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a13)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a13, alpha=.95)

dendro.2276a13 <- as.dendrogram(hc_tr2276a13)
dendro.plot.2276a13 <- ggdendrogram(data = dendro.2276a13, rotate = TRUE)
dendro.plot.2276a13 <- dendro.plot.2276a13 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a13)

melt.2276a13 <- melt(mat_tr2276a13)
head(melt.2276a13)
colnames(melt.2276a13) <- c("Targeted", "Bycatch", "TI")
melt.2276a13$TI_Ranges <- cut(melt.2276a13$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a13 <- order.dendrogram(dendro.2276a13)
melt.2276a13$Targeted <- factor(x = melt.2276a13$Targeted,
                                levels = melt.2276a13$Targeted[order.2276a13],
                                ordered = T)

heatmap.2276a13 <- ggplot(data=melt.2276a13, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a13)

grid.newpage()
print(heatmap.2276a13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.6.a 2014#####

tr2_276a14 <- subset(tr2_276a, year == 2014)

unique(tr2_276a14$year)
unique(tr2_276a14$foCatNat)
unique(tr2_276a14$area)
tr2_276a14$spp2_id <- as.numeric(as.factor(tr2_276a14$spp))
spp_tr2276a14 <- sort(unique(tr2_276a14$spp2_id))


ti_tr2276a14 <- data.frame() 

for(i in spp_tr2276a14){
  sp1 <- subset(tr2_276a14, spp2_id == i) 
  for(j in spp_tr2276a14){
    t <- data.frame()
    sp2 <- subset(tr2_276a14, spp2_id == j) 
    if(i == j){ti_tr2276a14[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a14[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a14)

row.names(ti_tr2276a14) <- sort(unique(tr2_276a14$spp))
colnames(ti_tr2276a14) <- sort(unique(tr2_276a14$spp))
View(ti_tr227414)

mat_tr2276a14 <- as.matrix(ti_tr2276a14)
dim(mat_tr2276a14)

##HClust with TR2 27.6.a 2014
dist_tr2276a14 <- dist(mat_tr2276a14, method = "euclidean") #distance matrix
hc_tr2276a14 <- hclust(dist_tr2276a14, method = "ward.D2")
plot(hc_tr2276a14)

coph_2276a14 <- cophenetic(hc_tr2276a14)
cor(dist_tr2276a14, coph_2276a14) #0.74, less than 0.75 so less than optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a14, FUN = hcut, method = "silhouette") #optimal number of clusters is 3


# Ward Hierarchical Clustering
plot(hc_tr2276a14) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr2276a14, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a14 <- pvclust(t(mat_tr2276a14), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a14)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a14, alpha=.95)

dendro.2276a14 <- as.dendrogram(hc_tr2276a14)
dendro.plot.2276a14 <- ggdendrogram(data = dendro.2276a14, rotate = TRUE)
dendro.plot.2276a14 <- dendro.plot.2276a14 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a14)

melt.2276a14 <- melt(mat_tr2276a14)
head(melt.2276a14)
colnames(melt.2276a14) <- c("Targeted", "Bycatch", "TI")
melt.2276a14$TI_Ranges <- cut(melt.2276a14$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a14 <- order.dendrogram(dendro.2276a14)
melt.2276a14$Targeted <- factor(x = melt.2276a14$Targeted,
                                levels = melt.2276a14$Targeted[order.2276a14],
                                ordered = T)

heatmap.2276a14 <- ggplot(data=melt.2276a14, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a14)

grid.newpage()
print(heatmap.2276a14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.6.a 2015####

tr2_276a15 <- subset(tr2_276a, year == 2015)

unique(tr2_276a15$year)
unique(tr2_276a15$foCatNat)
unique(tr2_276a15$area)
tr2_276a15$spp2_id <- as.numeric(as.factor(tr2_276a15$spp))
spp_tr2276a15 <- sort(unique(tr2_276a15$spp2_id))


ti_tr2276a15 <- data.frame() 

for(i in spp_tr2276a15){
  sp1 <- subset(tr2_276a15, spp2_id == i) 
  for(j in spp_tr2276a15){
    t <- data.frame()
    sp2 <- subset(tr2_276a15, spp2_id == j) 
    if(i == j){ti_tr2276a15[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a15[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a15)

row.names(ti_tr2276a15) <- sort(unique(tr2_276a15$spp))
colnames(ti_tr2276a15) <- sort(unique(tr2_276a15$spp))
View(ti_tr2276a15)

mat_tr2276a15 <- as.matrix(ti_tr2276a15)
dim(mat_tr2276a15)

##HClust with TR2 27.6.a 2015
dist_tr2276a15 <- dist(mat_tr2276a15, method = "euclidean") #distance matrix
hc_tr2276a15 <- hclust(dist_tr2276a15, method = "ward.D2")
plot(hc_tr2276a15)

coph_2276a15 <- cophenetic(hc_tr2276a15)
cor(dist_tr2276a15, coph_2276a15) #0.55, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a15, FUN = hcut, method = "silhouette") #optimal number of clusters is 6


# Ward Hierarchical Clustering
plot(hc_tr2276a15) # display dendogram
# draw dendogram with red borders around the 6 clusters 
rect.hclust(hc_tr2276a15, k=6, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a15 <- pvclust(t(mat_tr2276a15), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a15)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a15, alpha=.95)

dendro.2276a15 <- as.dendrogram(hc_tr2276a15, k = 2)
dendro.plot.2276a15 <- ggdendrogram(data = dendro.2276a15, rotate = TRUE, k = 2)
dendro.plot.2276a15 <- dendro.plot.2276a15 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a15)

melt.2276a15 <- melt(mat_tr2276a15)
head(melt.2276a15)
colnames(melt.2276a15) <- c("Targeted", "Bycatch", "TI")
melt.2276a15$TI_Ranges <- cut(melt.2276a15$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a15 <- order.dendrogram(dendro.2276a15)
melt.2276a15$Targeted <- factor(x = melt.2276a15$Targeted,
                                levels = melt.2276a15$Targeted[order.2276a15],
                                ordered = T)

heatmap.2276a15 <- ggplot(data=melt.2276a15, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a15)

grid.newpage()
print(heatmap.2276a15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.6.a 2016####

tr2_276a16 <- subset(tr2_276a, year == 2016)

unique(tr2_276a16$year)
unique(tr2_276a16$foCatNat)
unique(tr2_276a16$area)
tr2_276a16$spp2_id <- as.numeric(as.factor(tr2_276a16$spp))
spp_tr2276a16 <- sort(unique(tr2_276a16$spp2_id))


ti_tr2276a16 <- data.frame() 
for(i in spp_tr2276a16){
  sp1 <- subset(tr2_276a16, spp2_id == i) 
  for(j in spp_tr2276a16){
    t <- data.frame()
    sp2 <- subset(tr2_276a16, spp2_id == j) 
    if(i == j){ti_tr2276a16[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a16[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a16)

row.names(ti_tr2276a16) <- sort(unique(tr2_276a16$spp))
colnames(ti_tr2276a16) <- sort(unique(tr2_276a16$spp))
View(ti_tr2276a16)

mat_tr2276a16 <- as.matrix(ti_tr2276a16)
dim(mat_tr2276a16)

##HClust with TR1 27.6.a 2016
dist_tr2276a16 <- dist(mat_tr2276a16, method = "euclidean") #distance matrix
hc_tr2276a16 <- hclust(dist_tr2276a16, method = "ward.D2")
plot(hc_tr2276a16)

coph_2276a16 <- cophenetic(hc_tr2276a16)
cor(dist_tr2276a16, coph_2276a16) #0.51, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a16, FUN = hcut, method = "silhouette") #optimal number of clusters is 3


# Ward Hierarchical Clustering
plot(hc_tr2276a16) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr2276a16, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a16 <- pvclust(t(mat_tr2276a16), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a16)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a16, alpha=.95)

dendro.2276a16 <- as.dendrogram(hc_tr2276a16, k = 3)
dendro.plot.2276a16 <- ggdendrogram(data = dendro.2276a16, rotate = TRUE, k = 2)
dendro.plot.2276a16 <- dendro.plot.2276a16 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a16)

melt.2276a16 <- melt(mat_tr2276a16)
head(melt.2276a16)
colnames(melt.2276a16) <- c("Targeted", "Bycatch", "TI")
melt.2276a16$TI_Ranges <- cut(melt.2276a16$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a16 <- order.dendrogram(dendro.2276a16)
melt.2276a16$Targeted <- factor(x = melt.2276a16$Targeted,
                                levels = melt.2276a16$Targeted[order.2276a16],
                                ordered = T)

heatmap.2276a16 <- ggplot(data=melt.2276a16, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a16)

grid.newpage()
print(heatmap.2276a16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR2 27.6.a 2017####

tr2_276a17 <- subset(tr2_276a, year == 2017)

unique(tr2_276a17$year)
unique(tr2_276a17$foCatNat)
unique(tr2_276a17$area)
tr2_276a17$spp2_id <- as.numeric(as.factor(tr2_276a17$spp))
spp_tr2276a17 <- sort(unique(tr2_276a17$spp2_id))


ti_tr2276a17 <- data.frame() 

for(i in spp_tr2276a17){
  sp1 <- subset(tr2_276a17, spp2_id == i) 
  for(j in spp_tr2276a17){
    t <- data.frame()
    sp2 <- subset(tr2_276a17, spp2_id == j) 
    if(i == j){ti_tr2276a17[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a17[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a17)

row.names(ti_tr2276a17) <- sort(unique(tr2_276a17$spp))
colnames(ti_tr2276a17) <- sort(unique(tr2_276a17$spp))
View(ti_tr2276a17)

mat_tr2276a17 <- as.matrix(ti_tr2276a17)
dim(mat_tr2276a17)

##HClust with TR2 27.6.a 2017
dist_tr2276a17 <- dist(mat_tr2276a17, method = "euclidean") #distance matrix
hc_tr2276a17 <- hclust(dist_tr2276a17, method = "ward.D2")
plot(hc_tr2276a17)

coph_2276a17 <- cophenetic(hc_tr2276a17)
cor(dist_tr2276a17, coph_2276a17) #0.65, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a17, FUN = hcut, method = "silhouette") #optimal number of clusters is 9

# Ward Hierarchical Clustering
plot(hc_tr2276a17) # display dendogram 9 clusters 
rect.hclust(hc_tr2276a17, k=9, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a17 <- pvclust(t(mat_tr2276a17), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a17)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a17, alpha=.95)

dendro.2276a17 <- as.dendrogram(hc_tr2276a17, k = 2)
dendro.plot.2276a17 <- ggdendrogram(data = dendro.2276a17, rotate = TRUE, k = 2)
dendro.plot.2276a17 <- dendro.plot.2276a17 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a17)

melt.2276a17 <- melt(mat_tr2276a17)
head(melt.2276a17)
colnames(melt.2276a17) <- c("Targeted", "Bycatch", "TI")
melt.2276a17$TI_Ranges <- cut(melt.2276a17$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a17 <- order.dendrogram(dendro.2276a17)
melt.2276a17$Targeted <- factor(x = melt.2276a17$Targeted,
                                levels = melt.2276a17$Targeted[order.2276a17],
                                ordered = T)

heatmap.2276a17 <- ggplot(data=melt.2276a17, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a17)

grid.newpage()
print(heatmap.2276a17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR2 27.6.a 2018####

tr2_276a18 <- subset(tr2_276a, year == 2018)

unique(tr2_276a18$year)
unique(tr2_276a18$foCatNat)
unique(tr2_276a18$area)
tr2_276a18$spp2_id <- as.numeric(as.factor(tr2_276a18$spp))
spp_tr2276a18 <- sort(unique(tr2_276a18$spp2_id))


ti_tr2276a18 <- data.frame() 

for(i in spp_tr2276a18){
  sp1 <- subset(tr2_276a18, spp2_id == i) 
  for(j in spp_tr2276a18){
    t <- data.frame()
    sp2 <- subset(tr2_276a18, spp2_id == j) 
    if(i == j){ti_tr2276a18[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a18[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a18)

row.names(ti_tr2276a18) <- sort(unique(tr2_276a18$spp))
colnames(ti_tr2276a18) <- sort(unique(tr2_276a18$spp))
View(ti_tr2276a18)

mat_tr2276a18 <- as.matrix(ti_tr2276a18)
dim(mat_tr2276a18)

##HClust with TR2 27.6.a 2018
dist_tr2276a18 <- dist(mat_tr2276a18, method = "euclidean") #distance matrix
hc_tr2276a18 <- hclust(dist_tr2276a18, method = "ward.D2")
plot(hc_tr2276a18)

coph_2276a18 <- cophenetic(hc_tr2276a18)
cor(dist_tr2276a18, coph_2276a18) #0.88, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr2276a18, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr2276a18) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr2276a18, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit2276a18 <- pvclust(t(mat_tr2276a18), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit2276a18)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2276a18, alpha=.95)

dendro.2276a18 <- as.dendrogram(hc_tr2276a18, k = 2)
dendro.plot.2276a18 <- ggdendrogram(data = dendro.2276a18, rotate = TRUE, k = 2)
dendro.plot.2276a18 <- dendro.plot.2276a18 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.2276a18)

melt.2276a18 <- melt(mat_tr2276a18)
head(melt.2276a18)
colnames(melt.2276a18) <- c("Targeted", "Bycatch", "TI")
melt.2276a18$TI_Ranges <- cut(melt.2276a18$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.2276a18 <- order.dendrogram(dendro.2276a18)
melt.2276a18$Targeted <- factor(x = melt.2276a18$Targeted,
                                levels = melt.2276a18$Targeted[order.2276a18],
                                ordered = T)

heatmap.2276a18 <- ggplot(data=melt.2276a18, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.2276a18)

grid.newpage()
print(heatmap.2276a18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###Compiling TR2 27.6.a Heatmaps#####
pdf(file = "TR2_276a_Heatmaps.pdf")

grid.newpage()
print(heatmap.2276a, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.2276a18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.2276a18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

dev.off()
###Heatmaps for Averaging#####
ti_tr2276a09.2 <- data.frame() 
ti_tr2276a09.2[1:39, 1:39] <- -1

spp_tr2276a09.2 <- unique(tr2_276a09$spp_id)
spp_tr2276a09.2 <- sort(spp_tr2276a09.2)

for(i in spp_tr2276a09.2){
  sp1 <- subset(tr2_276a09, spp_id == i) 
  for(j in spp_tr2276a09.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a09, spp_id == j) 
    if(i == j){ti_tr2276a09.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a09.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a09.2)

row.names(ti_tr2276a09.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a09.2) <- sort(unique(trip382$spp))

mat_tr2276a09.2 <- as.matrix(ti_tr2276a09.2)
dim(mat_tr2276a09.2)
##2010
ti_tr2276a10.2 <- data.frame() 
ti_tr2276a10.2[1:39, 1:39] <- -1

spp_tr2276a10.2 <- unique(tr2_276a10$spp_id)
spp_tr2276a10.2 <- sort(spp_tr2276a10.2)

for(i in spp_tr2276a10.2){
  sp1 <- subset(tr2_276a10, spp_id == i) 
  for(j in spp_tr2276a10.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a10, spp_id == j) 
    if(i == j){ti_tr2276a10.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a10.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a10.2)

row.names(ti_tr2276a10.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a10.2) <- sort(unique(trip382$spp))

mat_tr2276a10.2 <- as.matrix(ti_tr2276a10.2)
dim(mat_tr2276a10.2)

##2011
ti_tr2276a11.2 <- data.frame() 
ti_tr2276a11.2[1:39, 1:39] <- -1

spp_tr2276a11.2 <- unique(tr2_276a11$spp_id)
spp_tr2276a11.2 <- sort(spp_tr2276a11.2)

for(i in spp_tr2276a11.2){
  sp1 <- subset(tr2_276a11, spp_id == i) 
  for(j in spp_tr2276a11.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a11, spp_id == j) 
    if(i == j){ti_tr2276a11.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a11.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a11.2)

row.names(ti_tr2276a11.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a11.2) <- sort(unique(trip382$spp))

mat_tr2276a11.2 <- as.matrix(ti_tr2276a11.2)
dim(mat_tr2276a11.2)

##2012
ti_tr2276a12.2 <- data.frame() 
ti_tr2276a12.2[1:39, 1:39] <- -1

spp_tr2276a12.2 <- unique(tr2_276a12$spp_id)
spp_tr2276a12.2 <- sort(spp_tr2276a12.2)

for(i in spp_tr2276a12.2){
  sp1 <- subset(tr2_276a12, spp_id == i) 
  for(j in spp_tr2276a12.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a12, spp_id == j) 
    if(i == j){ti_tr2276a12.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a12.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a10.2)

row.names(ti_tr2276a12.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a12.2) <- sort(unique(trip382$spp))

mat_tr2276a12.2 <- as.matrix(ti_tr2276a12.2)
dim(mat_tr2276a12.2)

##2013
ti_tr2276a13.2 <- data.frame() 
ti_tr2276a13.2[1:39, 1:39] <- -1

spp_tr2276a13.2 <- unique(tr2_276a13$spp_id)
spp_tr2276a13.2 <- sort(spp_tr2276a13.2)

for(i in spp_tr2276a13.2){
  sp1 <- subset(tr2_276a13, spp_id == i) 
  for(j in spp_tr2276a13.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a13, spp_id == j) 
    if(i == j){ti_tr2276a13.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a13.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a13.2)

row.names(ti_tr2276a13.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a13.2) <- sort(unique(trip382$spp))

mat_tr2276a13.2 <- as.matrix(ti_tr2276a13.2)
dim(mat_tr2276a13.2)

##2014
ti_tr2276a14.2 <- data.frame() 
ti_tr2276a14.2[1:39, 1:39] <- -1

spp_tr2276a14.2 <- unique(tr2_276a14$spp_id)
spp_tr2276a14.2 <- sort(spp_tr2276a14.2)

for(i in spp_tr2276a14.2){
  sp1 <- subset(tr2_276a14, spp_id == i) 
  for(j in spp_tr2276a14.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a14, spp_id == j) 
    if(i == j){ti_tr2276a14.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a14.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a14.2)

row.names(ti_tr2276a14.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a14.2) <- sort(unique(trip382$spp))

mat_tr2276a14.2 <- as.matrix(ti_tr2276a14.2)
dim(mat_tr2276a14.2)

##2015
ti_tr2276a15.2 <- data.frame() 
ti_tr2276a15.2[1:39, 1:39] <- -1

spp_tr2276a15.2 <- unique(tr2_276a15$spp_id)
spp_tr2276a15.2 <- sort(spp_tr2276a15.2)

for(i in spp_tr2276a15.2){
  sp1 <- subset(tr2_276a15, spp_id == i) 
  for(j in spp_tr2276a15.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a15, spp_id == j) 
    if(i == j){ti_tr2276a15.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a15.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a15.2)

row.names(ti_tr2276a15.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a15.2) <- sort(unique(trip382$spp))

mat_tr2276a15.2 <- as.matrix(ti_tr2276a15.2)
dim(mat_tr2276a15.2)

##2016
ti_tr2276a16.2 <- data.frame() 
ti_tr2276a16.2[1:39, 1:39] <- -1

spp_tr2276a16.2 <- unique(tr2_276a16$spp_id)
spp_tr2276a16.2 <- sort(spp_tr2276a16.2)

for(i in spp_tr2276a16.2){
  sp1 <- subset(tr2_276a16, spp_id == i) 
  for(j in spp_tr2276a16.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a16, spp_id == j) 
    if(i == j){ti_tr2276a16.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a16.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a16.2)

row.names(ti_tr2276a16.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a16.2) <- sort(unique(trip382$spp))

mat_tr2276a16.2 <- as.matrix(ti_tr2276a16.2)
dim(mat_tr2276a16.2)

##2017
ti_tr2276a17.2 <- data.frame() 
ti_tr2276a17.2[1:39, 1:39] <- -1

spp_tr2276a17.2 <- unique(tr2_276a17$spp_id)
spp_tr2276a17.2 <- sort(spp_tr2276a17.2)

for(i in spp_tr2276a17.2){
  sp1 <- subset(tr2_276a17, spp_id == i) 
  for(j in spp_tr2276a17.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a17, spp_id == j) 
    if(i == j){ti_tr2276a17.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a17.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a17.2)

row.names(ti_tr2276a17.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a17.2) <- sort(unique(trip382$spp))

mat_tr2276a17.2 <- as.matrix(ti_tr2276a17.2)
dim(mat_tr2276a17.2)

##2018
ti_tr2276a18.2 <- data.frame() 
ti_tr2276a18.2[1:39, 1:39] <- -1

spp_tr2276a18.2 <- unique(tr2_276a18$spp_id)
spp_tr2276a18.2 <- sort(spp_tr2276a18.2)

for(i in spp_tr2276a18.2){
  sp1 <- subset(tr2_276a18, spp_id == i) 
  for(j in spp_tr2276a18.2){
    t <- data.frame()
    sp2 <- subset(tr2_276a18, spp_id == j) 
    if(i == j){ti_tr2276a18.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr2276a18.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr2276a18.2)

row.names(ti_tr2276a18.2) <- sort(unique(trip382$spp))
colnames(ti_tr2276a18.2) <- sort(unique(trip382$spp))

mat_tr2276a18.2 <- as.matrix(ti_tr2276a18.2)
dim(mat_tr2276a18.2)
##Averaging the Matrices together
matavg_tr2276a <- (mat_tr2276a09.2 + mat_tr2276a10.2 + mat_tr2276a11.2 + 
                     mat_tr2276a12.2 + mat_tr2276a13.2 + mat_tr2276a14.2 + 
                     mat_tr2276a15.2 + mat_tr2276a16.2 +  mat_tr2276a17.2 + mat_tr2276a18.2)/10


dist_avg2276a <- dist(matavg_tr2276a, method = "euclidean") #distance matrix
hc_avg2276a <- hclust(dist_avg2276a, method = "ward.D2")
plot(hc_avg2276a)

coph_avg2276a <- cophenetic(hc_avg2276a)
cor(dist_avg2276a, coph_avg2276a) #0.93, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(matavg_tr2276a, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_avg2276a) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_avg2276a, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitavg2276a <- pvclust(t(matavg_tr2276a), method.hclust="ward.D2",
                       method.dist="euclidean")
plot(fitavg2276a)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitavg2276a, alpha=.95)

dendro.avg2276a <- as.dendrogram(hc_avg2276a, k = 2)
dendro.plot.avg2276a <- ggdendrogram(data = dendro.avg2276a, rotate = TRUE, k = 2)
dendro.plot.avg2276a <- dendro.plot.avg2276a + theme(axis.text.y = element_blank(), 
                                                     axis.text.x.bottom =  element_blank())

print(dendro.plot.avg2276a)

melt.avg2276a <- melt(matavg_tr2276a)
head(melt.avg2276a)
colnames(melt.avg2276a) <- c("Targeted", "Bycatch", "TI")
melt.avg2276a$TI_Ranges <- cut(melt.avg2276a$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.avg2276a <- order.dendrogram(dendro.avg2276a)
melt.avg2276a$Targeted <- factor(x = melt.avg2276a$Targeted,
                                 levels = melt.avg2276a$Targeted[order.avg2276a],
                                 ordered = T)

heatmap.avg2276a <- ggplot(data=melt.avg2276a, 
                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.avg2276a)

grid.newpage()
print(heatmap.avg2276a, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.avg2276a, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

####Creating Array TR2 27.6.a####
a_tr2276a <- abind(mat_tr2276a09.2, mat_tr2276a10.2, mat_tr2276a11.2, 
                   mat_tr2276a12.2, mat_tr2276a13.2, mat_tr2276a14.2,
                   mat_tr2276a15.2, mat_tr2276a16.2, mat_tr2276a17.2, mat_tr2276a18.2, along = 3)
View(a_tr2276a[,,10])

melta_tr2276a <- melt(a_tr2276a)
head(melta_tr2276a)
colnames(melta_tr2276a) <- c("Targeted", "Bycatch", "Year", "TIS")

melta_tr2276a$TI_Ranges <- cut(melta_tr2276a$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
str(melta_tr2276a)


melta_tr2276a$Area[1:15210] <- "27.6.a"
melta_tr2276a$Area <- as.factor(melta_tr2276a$Area)
melta_tr2276a$Gear[1:15210] <- "TR2"
melta_tr2276a$Gear <- as.factor(melta_tr2276a$Gear)

melta_tr2276a$Year <- as.factor(melta_tr2276a$Year)
melta_tr2276a$Year <- revalue(melta_tr2276a$Year, c("1"="2009", "2"="2010", "3"="2011", "4"="2012",
                                                    "5"="2013", "6"="2014", "7"="2015", "8"="2016",
                                                    "9"="2017", "10"="2018"))

melta_tr2276a <- na.omit(melta_tr2276a) #removing absent species that are marked as NA in TI_Ranges

tr2276a_by <- with(melta_tr2276a, tapply(TIS, INDEX = list(Bycatch, Year), FUN = "mean"))
tr2276a_tar <- with(melta_tr2276a, tapply(TIS, INDEX = list(Targeted, Year), FUN = "mean"))
tr2276a_overall <- with(melta_tr2276a, tapply (TIS, INDEX = Year, FUN = "mean"))
tr2276a_by <- as.data.frame(tr2276a_by)
tr2276a_tar <- as.data.frame(tr2276a_tar)
tr2276a_overall <- as.data.frame(tr2276a_overall)

tr2276a_overall$var <- with(melta_tr2276a, tapply(TIS, INDEX = Year, FUN = "var"))
tr2276a_overall$sd <- sqrt(tr2276a_overall$var)
tr2276a_overall$cv <- 100*((tr2276a_overall$sd)/(tr1274_overall$tr1274_overall))

tr2276a_by$Overall <- rowMeans(tr2276a_by)
tr2276a_by$var <- with(melta_tr2276a, tapply(TIS, INDEX = list(Bycatch,Year), FUN = "var"))
tr2276a_by$sd <- sqrt(tr2276a_by$var)
tr2276a_by$cv <- 100*((tr2276a_by$sd)/(tr2276a_by$Overall))
tr2276a_by.05 <- tr2276a_by[!(tr2276a_by$Overall<=5),]
tr2276a_by.50 <- tr2276a_by[!(tr2276a_by$Overall<=50),]
View(tr2276a_by.05)
View(tr2276a_by.50)

tr2276a_tar$Overall <- rowMeans(tr2276a_tar)
tr2276a_tar$var <- with(melta_tr2276a, tapply(TIS, INDEX = list(Targeted,Year), FUN = "var"))
tr2276a_tar$sd <- sqrt(tr2276a_tar$var)
tr2276a_tar$cv <- 100*((tr2276a_tar$sd)/(tr2276a_tar$Overall))
tr2276a_tar.05 <- tr2276a_tar[!(tr2276a_tar$Overall<=5),]
tr2276a_tar.25 <- tr2276a_tar[!(tr2276a_tar$Overall<=25),]
View(tr2276a_tar.05)
View(tr2276a_tar.25)


####Centered Moving Average TR2 27.6.a--EACH YEAR#####
?rollmean
?select
# 
# test <- tr1274_overall %>% select(TIS_year = tr1274_overall) %>% 
#   mutate(TIS_mavg3 = rollmean(TIS_year, k = 3, fill = NA),
#          Year = 2009:2018)
# testg<- test %>% gather(metric, value, TIS_year:TIS_mavg3) %>%
#   ggplot(aes(Year, value, color = metric)) +
#   geom_line()
# print(testg)

t2276a <- data.frame()
for(i in 2009:2018){
  selected <- c(i-1, i, i+1)
  f1 <- subset(melta_tr2276a, Year %in% selected)
  if((i == 2009)| (i==2018)){
    t2276a[i-2008, 1] <- NA} 
  else{t2276a[i-2008,1] <- mean(f1$TIS)}
}
print(t2276a)

t2276a.2 <- t2276a %>% select(TIS_mavg = V1) %>% mutate(TIS_yavg = tr2276a_overall$tr2276a_overall,
                                                        Year = 2009:2018)
t2276a.g <- t2276a.2 %>% gather(Metric, TIS, TIS_mavg:TIS_yavg) %>%
  ggplot(aes(Year, TIS, color = Metric)) +
  geom_line()
print(t2276a.g)

##########################################################################
#####Area 27.6.b, TR1 Subsets#################

tr1_276b <- subset(trip382, gear_id == 3 & area_id == 7)
tr1_276b

unique(tr1_276b$year)
unique(tr1_276b$foCatNat)
unique(tr1_276b$area)
tr1_276b$spp_id <- tr1_276b$spp2_id
tr1_276b$spp2_id <- as.numeric(as.factor(tr1_276b$spp))
spp_tr1276b <- unique(tr1_276b$spp2_id)
spp_tr1276b <- sort(spp_tr1276b)

ti_tr1276b <- data.frame() 

for(i in spp_tr1276b){
  sp1 <- subset(tr1_276b, spp2_id == i) 
  for(j in spp_tr1276b){
    t <- data.frame()
    sp2 <- subset(tr1_276b, spp2_id == j) 
    if(i == j){ti_tr1276b[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$year[k] %in% sp2$year)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b)

row.names(ti_tr1276b) <- sort(unique(tr1_276b$spp))
colnames(ti_tr1276b) <- sort(unique(tr1_276b$spp))
View(ti_tr1276b)

mat_tr1276b <- as.matrix(ti_tr1276b)
##HClust with TR1 27.6.b
dist_tr1276b <- dist(mat_tr1276b, method = "euclidean") #distance matrix
hc_tr1276b<- hclust(dist_tr1276b, method = "ward.D2")
plot(hc_tr1276b)

coph_1276b <- cophenetic(hc_tr1276b)
cor(dist_tr1276b, coph_1276b) #0.80, greater than 0.75, therefore good

fviz_nbclust(mat_tr1276b, FUN = hcut, method = "silhouette") #optimal number of clusters is 3

# Ward Hierarchical Clustering
plot(hc_tr1276b) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_tr1276b, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit276b <- pvclust(t(mat_tr1276b), method.hclust="ward.D2",
                   method.dist="euclidean")
plot(fit276b)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit276b, alpha=.95)

dendro.1276b <- as.dendrogram(hc_tr1276b)
dendro.plot.1276b <- ggdendrogram(data = dendro.1276b, rotate = TRUE)
dendro.plot.1276b <- dendro.plot.1276b + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b)

melt.1276b <- melt(mat_tr1276b)
head(melt.1276b)
colnames(melt.1276b) <- c("Targeted", "Bycatch", "TI")
melt.1276b$TI_Ranges <- cut(melt.1276b$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b <- order.dendrogram(dendro.1276b)
melt.1276b$Targeted <- factor(x = melt.1276b$Targeted,
                              levels = melt.1276b$Targeted[order.1276b],
                              ordered = T)


heatmap.1276b <- ggplot(data=melt.1276b, 
                        aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b)

grid.newpage()
print(heatmap.1276b, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

###TR1 27.6.b 2009########
tr1_276b09 <- subset(tr1_276b, year == 2009)

unique(tr1_276b09$year)
unique(tr1_276b09$foCatNat)
unique(tr1_276b09$area)
tr1_276b09$spp2_id <- as.numeric(as.factor(tr1_276b09$spp))
spp_tr1276b09 <- unique(tr1_276b09$spp2_id)
spp_tr1276b09 <- sort(spp_tr1276b09)

ti_tr1276b09 <- data.frame() 

for(i in spp_tr1276b09){
  sp1 <- subset(tr1_276b09, spp2_id == i) 
  for(j in spp_tr1276b09){
    t <- data.frame()
    sp2 <- subset(tr1_276b09, spp2_id == j) 
    if(i == j){ti_tr1276b09[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b09[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b09)

row.names(ti_tr1276b09) <- sort(unique(tr1_276b09$spp))
colnames(ti_tr1276b09) <- sort(unique(tr1_276b09$spp))
View(ti_tr1276b09)

mat_tr1276b09 <- as.matrix(ti_tr1276b09)
dim(mat_tr1276b09)
##HClust with TR1 27.6.a 2009
dist_tr1276b09 <- dist(mat_tr1276b09, method = "euclidean") #distance matrix
hc_tr1276b09<- hclust(dist_tr1276b09, method = "ward.D2")
plot(hc_tr1276b09)

coph_1276b09 <- cophenetic(hc_tr1276b09)
cor(dist_tr1276b09, coph_1276b09) #Greater than 0.75, therefore considered good


fviz_nbclust(mat_tr1276b09, FUN = hcut, method = "silhouette") #optimal number of clusters is 6


# Ward Hierarchical Clustering
plot(hc_tr1276b09) # display dendogram
# draw dendogram with red borders around the 6 clusters 
rect.hclust(hc_tr1276b09, k=6, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b09 <- pvclust(t(mat_tr1276b09), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b09)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b09, alpha=.95)

dendro.1276b09 <- as.dendrogram(hc_tr1276b09)
dendro.plot.1276b09 <- ggdendrogram(data = dendro.1276b09, rotate = TRUE)
dendro.plot.1276b09 <- dendro.plot.1276b09 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b09)

melt.1276b09 <- melt(mat_tr1276b09)
head(melt.1276b09)
colnames(melt.1276b09) <- c("Targeted", "Bycatch", "TI")
melt.1276b09$TI_Ranges <- cut(melt.1276b09$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b09 <- order.dendrogram(dendro.1276b09)
melt.1276b09$Targeted <- factor(x = melt.1276b09$Targeted,
                                levels = melt.1276b09$Targeted[order.1276b09],
                                ordered = T)

heatmap.1276b09 <- ggplot(data=melt.1276b09, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b09)

grid.newpage()
print(heatmap.1276b09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

##
###NO DATA FOR 27.6.b 2010##########
# tr1_276a10 <- subset(tr1_276a, year == 2010)
# 
# unique(tr1_276a10$year)
# unique(tr1_276a10$foCatNat)
# unique(tr1_276a10$area)
# tr1_276a10$spp2_id <- as.numeric(as.factor(tr1_276a10$spp))
# spp_tr1276a10 <- unique(tr1_276a10$spp2_id)
# spp_tr1276a10 <- sort(spp_tr1276a10)
# 
# ti_tr1276a10 <- data.frame() 
# 
# for(i in spp_tr1276a10){
#   sp1 <- subset(tr1_276a10, spp2_id == i) 
#   for(j in spp_tr1276a10){
#     t <- data.frame()
#     sp2 <- subset(tr1_276a10, spp2_id == j) 
#     if(i == j){ti_tr1276a10[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_tr1276a10[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_tr1276a10)
# 
# row.names(ti_tr1276a10) <- sort(unique(tr1_276a10$spp))
# colnames(ti_tr1276a10) <- sort(unique(tr1_276a10$spp))
# View(ti_tr1276a10)
# 
# mat_tr1276a10 <- as.matrix(ti_tr1276a10)
# dim(mat_tr1276a10)
# 
# 
# ##HClust with TR1 27.6.a 2010
# dist_tr1276a10 <- dist(mat_tr1276a10, method = "euclidean") #distance matrix
# hc_tr1276a10 <- hclust(dist_tr1276a10, method = "ward.D2")
# plot(hc_tr1276a10)
# 
# coph_1276a10 <- cophenetic(hc_tr1276a10)
# cor(dist_tr1276a10, coph_1276a10) #0.87, greater than 0.75, optimal
# 
# fviz_nbclust(mat_tr1276a10, FUN = hcut, method = "silhouette") #optimal number of clusters is 2
# 
# 
# # Ward Hierarchical Clustering
# plot(hc_tr1276a10) # display dendogram
# # draw dendogram with red borders around the 2 clusters 
# rect.hclust(hc_tr1276a10, k=2, border="red") 
# 
# # Ward Hierarchical Clustering with Bootstrapped p values
# fit1276a10 <- pvclust(t(mat_tr1276a10), method.hclust="ward.D2",
#                       method.dist="euclidean")
# plot(fit1276a10)# dendogram with p values
# # add rectangles around groups highly supported by the data
# pvrect(fit1276a10, alpha=.95)
# 
# dendro.1276a10 <- as.dendrogram(hc_tr1276a10, k = 2)
# dendro.plot.1276a10 <- ggdendrogram(data = dendro.1276a10, rotate = TRUE, k = 2)
# dendro.plot.1276a10 <- dendro.plot.1276a10 + theme(axis.text.y = element_blank(), 
#                                                    axis.text.x.bottom =  element_blank())
# 
# print(dendro.plot.1276a10)
# 
# melt.1276a10 <- melt(mat_tr1276a10)
# head(melt.1276a10)
# colnames(melt.1276a10) <- c("Targeted", "Bycatch", "TI")
# melt.1276a10$TI_Ranges <- cut(melt.1276a10$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# order.1276a10 <- order.dendrogram(dendro.1276a10)
# melt.1276a10$Targeted <- factor(x = melt.1276a10$Targeted,
#                                 levels = melt.1276a10$Targeted[order.1276a10],
#                                 ordered = T)
# 
# heatmap.1276a10 <- ggplot(data=melt.1276a10, 
#                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)")
# print(heatmap.1276a10)
# 
# grid.newpage()
# print(heatmap.1276a10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot.1276a10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))


##

###TR1 27.6.b 2011#######

tr1_276b11 <- subset(tr1_276b, year == 2011)

unique(tr1_276b11$year)
unique(tr1_276b11$foCatNat)
unique(tr1_276b11$area)
tr1_276b11$spp2_id <- as.numeric(as.factor(tr1_276b11$spp))
spp_tr1276b11 <- sort(unique(tr1_276b11$spp2_id))


ti_tr1276b11 <- data.frame() 

for(i in spp_tr1276b11){
  sp1 <- subset(tr1_276b11, spp2_id == i) 
  for(j in spp_tr1276b11){
    t <- data.frame()
    sp2 <- subset(tr1_276b11, spp2_id == j) 
    if(i == j){ti_tr1276b11[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b11[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b11)

row.names(ti_tr1276b11) <- sort(unique(tr1_276b11$spp))
colnames(ti_tr1276b11) <- sort(unique(tr1_276b11$spp))
View(ti_tr1276b11)

mat_tr1276b11 <- as.matrix(ti_tr1276b11)
dim(mat_tr1276b11)

##HClust with TR1 27.4 2011
dist_tr1276b11 <- dist(mat_tr1276b11, method = "euclidean") #distance matrix
hc_tr1276b11 <- hclust(dist_tr1276b11, method = "ward.D2")
plot(hc_tr1276b11)

coph_1276b11 <- cophenetic(hc_tr1276b11)
cor(dist_tr1276b11, coph_1276b11) #0.99, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b11, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr1276b11) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276b11, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b11 <- pvclust(t(mat_tr1276b11), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b11)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b11, alpha=.95)

dendro.1276b11 <- as.dendrogram(hc_tr1276b11)
dendro.plot.1276b11 <- ggdendrogram(data = dendro.1276b11, rotate = TRUE)
dendro.plot.1276b11 <- dendro.plot.1276b11 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b11)

melt.1276b11 <- melt(mat_tr1276b11)
head(melt.1276b11)
colnames(melt.1276b11) <- c("Targeted", "Bycatch", "TI")
melt.1276b11$TI_Ranges <- cut(melt.1276b11$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b11 <- order.dendrogram(dendro.1276b11)
melt.1276b11$Targeted <- factor(x = melt.1276b11$Targeted,
                                levels = melt.1276b11$Targeted[order.1276b11],
                                ordered = T)

heatmap.1276b11 <- ggplot(data=melt.1276b11, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b11)

grid.newpage()
print(heatmap.1276b11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR1 27.6.b 2012######

tr1_276b12 <- subset(tr1_276b, year == 2012)

unique(tr1_276b12$year)
unique(tr1_276b12$foCatNat)
unique(tr1_276b12$area)
tr1_276b12$spp2_id <- as.numeric(as.factor(tr1_276b12$spp))
spp_tr1276b12 <- sort(unique(tr1_276b12$spp2_id))


ti_tr1276b12 <- data.frame() 

for(i in spp_tr1276b12){
  sp1 <- subset(tr1_276b12, spp2_id == i) 
  for(j in spp_tr1276b12){
    t <- data.frame()
    sp2 <- subset(tr1_276b12, spp2_id == j) 
    if(i == j){ti_tr1276b12[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b12[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b12)

row.names(ti_tr1276b12) <- sort(unique(tr1_276b12$spp))
colnames(ti_tr1276b12) <- sort(unique(tr1_276b12$spp))
View(ti_tr1276b12)

mat_tr1276b12 <- as.matrix(ti_tr1276b12)
dim(mat_tr1276b12)

##HClust with TR1 27.6.b 2012
dist_tr1276b12 <- dist(mat_tr1276b12, method = "euclidean") #distance matrix
hc_tr1276b12 <- hclust(dist_tr1276b12, method = "ward.D2")
plot(hc_tr1276b12)

coph_1276b12 <- cophenetic(hc_tr1276b12)
cor(dist_tr1276b12, coph_1276b12) #0.76, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b12, FUN = hcut, method = "silhouette") #optimal number of clusters is 9


# Ward Hierarchical Clustering
plot(hc_tr1276b12) # display dendogram
# draw dendogram with red borders around the 9 clusters 
rect.hclust(hc_tr1276b12, k=9, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b12 <- pvclust(t(mat_tr1276b12), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b12)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b12, alpha=.95)

dendro.1276b12 <- as.dendrogram(hc_tr1276b12, k = 2)
dendro.plot.1276b12 <- ggdendrogram(data = dendro.1276b12, rotate = TRUE, k = 3)
dendro.plot.1276b12 <- dendro.plot.1276b12 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b12)

melt.1276b12 <- melt(mat_tr1276b12)
head(melt.1276b12)
colnames(melt.1276b12) <- c("Targeted", "Bycatch", "TI")
melt.1276b12$TI_Ranges <- cut(melt.1276b12$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b12 <- order.dendrogram(dendro.1276b12)
melt.1276b12$Targeted <- factor(x = melt.1276b12$Targeted,
                                levels = melt.1276b12$Targeted[order.1276b12],
                                ordered = T)

heatmap.1276b12 <- ggplot(data=melt.1276b12, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b12)

grid.newpage()
print(heatmap.1276b12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###ALL 100% TR1 27.6.b 2013#####

tr1_276b13 <- subset(tr1_276b, year == 2013)

unique(tr1_276b13$year)
unique(tr1_276b13$foCatNat)
unique(tr1_276b13$area)
tr1_276b13$spp2_id <- as.numeric(as.factor(tr1_276b13$spp))
spp_tr1276b13 <- sort(unique(tr1_276b13$spp2_id))


ti_tr1276b13 <- data.frame() 

for(i in spp_tr1276b13){
  sp1 <- subset(tr1_276b13, spp2_id == i) 
  for(j in spp_tr1276b13){
    t <- data.frame()
    sp2 <- subset(tr1_276b13, spp2_id == j) 
    if(i == j){ti_tr1276b13[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b13[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b13)

row.names(ti_tr1276b13) <- sort(unique(tr1_276b13$spp))
colnames(ti_tr1276b13) <- sort(unique(tr1_276b13$spp))
View(ti_tr1276b13)

mat_tr1276b13 <- as.matrix(ti_tr1276b13)
dim(mat_tr1276b13)

##HClust with TR1 27.6.b 2013
dist_tr1276b13 <- dist(mat_tr1276b13, method = "euclidean") #distance matrix
hc_tr1276b13 <- hclust(dist_tr1276b13, method = "ward.D2")
plot(hc_tr1276b13)

coph_1276b13 <- cophenetic(hc_tr1276b13)
cor(dist_tr1276b13, coph_1276b13) #0.58, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b13, FUN = hcut, method = "silhouette") #optimal number of clusters is 10


# Ward Hierarchical Clustering
plot(hc_tr1276b13) # display dendogram
# draw dendogram with red borders around the 10 clusters 
rect.hclust(hc_tr1276b13, k=10, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b13 <- pvclust(t(mat_tr1276b13), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b13)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b13, alpha=.95)

dendro.1276b13 <- as.dendrogram(hc_tr1276b13, k = 2)
dendro.plot.1276b13 <- ggdendrogram(data = dendro.1276b13, rotate = TRUE, k = 2)
dendro.plot.1276b13 <- dendro.plot.1276b13 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b13)

melt.1276b13 <- melt(mat_tr1276b13)
head(melt.1276b13)
colnames(melt.1276b13) <- c("Targeted", "Bycatch", "TI")
melt.1276b13$TI_Ranges <- cut(melt.1276b13$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b13 <- order.dendrogram(dendro.1276b13)
melt.1276b13$Targeted <- factor(x = melt.1276b13$Targeted,
                                levels = melt.1276b13$Targeted[order.1276b13],
                                ordered = T)

heatmap.1276b13 <- ggplot(data=melt.1276b13, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b13)

grid.newpage()
print(heatmap.1276b13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###TR1 27.6.b 2014#####

tr1_276b14 <- subset(tr1_276b, year == 2014)

unique(tr1_276b14$year)
unique(tr1_276b14$foCatNat)
unique(tr1_276b14$area)
tr1_276b14$spp2_id <- as.numeric(as.factor(tr1_276b14$spp))
spp_tr1276b14 <- sort(unique(tr1_276b14$spp2_id))


ti_tr1276b14 <- data.frame() 

for(i in spp_tr1276b14){
  sp1 <- subset(tr1_276b14, spp2_id == i) 
  for(j in spp_tr1276b14){
    t <- data.frame()
    sp2 <- subset(tr1_276b14, spp2_id == j) 
    if(i == j){ti_tr1276b14[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b14[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b14)

row.names(ti_tr1276b14) <- sort(unique(tr1_276b14$spp))
colnames(ti_tr1276b14) <- sort(unique(tr1_276b14$spp))
View(ti_tr127414)

mat_tr1276b14 <- as.matrix(ti_tr1276b14)
dim(mat_tr1276b14)

##HClust with TR1 27.6.a 2014
dist_tr1276b14 <- dist(mat_tr1276b14, method = "euclidean") #distance matrix
hc_tr1276b14 <- hclust(dist_tr1276b14, method = "ward.D2")
plot(hc_tr1276b14)

coph_1276b14 <- cophenetic(hc_tr1276b14)
cor(dist_tr1276b14, coph_1276b14) #0.95, greater than 0.75, optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b14, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr1276b14) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_tr1276b14, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b14 <- pvclust(t(mat_tr1276b14), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b14)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b14, alpha=.95)

dendro.1276b14 <- as.dendrogram(hc_tr1276b14)
dendro.plot.1276b14 <- ggdendrogram(data = dendro.1276b14, rotate = TRUE)
dendro.plot.1276b14 <- dendro.plot.1276b14 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b14)

melt.1276b14 <- melt(mat_tr1276b14)
head(melt.1276b14)
colnames(melt.1276b14) <- c("Targeted", "Bycatch", "TI")
melt.1276b14$TI_Ranges <- cut(melt.1276b14$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b14 <- order.dendrogram(dendro.1276b14)
melt.1276b14$Targeted <- factor(x = melt.1276b14$Targeted,
                                levels = melt.1276b14$Targeted[order.1276b14],
                                ordered = T)

heatmap.1276b14 <- ggplot(data=melt.1276b14, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b14)

grid.newpage()
print(heatmap.1276b14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.b 2015####

tr1_276b15 <- subset(tr1_276b, year == 2015)

unique(tr1_276b15$year)
unique(tr1_276b15$foCatNat)
unique(tr1_276b15$area)
tr1_276b15$spp2_id <- as.numeric(as.factor(tr1_276b15$spp))
spp_tr1276b15 <- sort(unique(tr1_276b15$spp2_id))


ti_tr1276b15 <- data.frame() 

for(i in spp_tr1276b15){
  sp1 <- subset(tr1_276b15, spp2_id == i) 
  for(j in spp_tr1276b15){
    t <- data.frame()
    sp2 <- subset(tr1_276b15, spp2_id == j) 
    if(i == j){ti_tr1276b15[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b15[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b15)

row.names(ti_tr1276b15) <- sort(unique(tr1_276b15$spp))
colnames(ti_tr1276b15) <- sort(unique(tr1_276b15$spp))
View(ti_tr1276b15)

mat_tr1276b15 <- as.matrix(ti_tr1276b15)
dim(mat_tr1276b15)

##HClust with TR1 27.6.b 2015
dist_tr1276b15 <- dist(mat_tr1276b15, method = "euclidean") #distance matrix
hc_tr1276b15 <- hclust(dist_tr1276b15, method = "ward.D2")
plot(hc_tr1276b15)

coph_1276b15 <- cophenetic(hc_tr1276b15)
cor(dist_tr1276b15, coph_1276b15) #0.85, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b15, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_tr1276b15) # display dendogram
# draw dendogram with red borders around the 4 clusters 
rect.hclust(hc_tr1276b15, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b15 <- pvclust(t(mat_tr1276b15), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b15)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b15, alpha=.95)

dendro.1276b15 <- as.dendrogram(hc_tr1276b15)
dendro.plot.1276b15 <- ggdendrogram(data = dendro.1276b15, rotate = TRUE)
dendro.plot.1276b15 <- dendro.plot.1276b15 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b15)

melt.1276b15 <- melt(mat_tr1276b15)
head(melt.1276b15)
colnames(melt.1276b15) <- c("Targeted", "Bycatch", "TI")
melt.1276b15$TI_Ranges <- cut(melt.1276b15$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b15 <- order.dendrogram(dendro.1276b15)
melt.1276b15$Targeted <- factor(x = melt.1276b15$Targeted,
                                levels = melt.1276b15$Targeted[order.1276b15],
                                ordered = T)

heatmap.1276b15 <- ggplot(data=melt.1276b15, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b15)

grid.newpage()
print(heatmap.1276b15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.b 2016####

tr1_276b16 <- subset(tr1_276b, year == 2016)

unique(tr1_276b16$year)
unique(tr1_276b16$foCatNat)
unique(tr1_276b16$area)
tr1_276b16$spp2_id <- as.numeric(as.factor(tr1_276b16$spp))
spp_tr1276b16 <- sort(unique(tr1_276b16$spp2_id))


ti_tr1276b16 <- data.frame() 
for(i in spp_tr1276b16){
  sp1 <- subset(tr1_276b16, spp2_id == i) 
  for(j in spp_tr1276b16){
    t <- data.frame()
    sp2 <- subset(tr1_276b16, spp2_id == j) 
    if(i == j){ti_tr1276b16[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b16[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b16)

row.names(ti_tr1276b16) <- sort(unique(tr1_276b16$spp))
colnames(ti_tr1276b16) <- sort(unique(tr1_276b16$spp))
View(ti_tr1276b16)

mat_tr1276b16 <- as.matrix(ti_tr1276b16)
dim(mat_tr1276b16)

##HClust with TR1 27.6.a 2016
dist_tr1276b16 <- dist(mat_tr1276b16, method = "euclidean") #distance matrix
hc_tr1276b16 <- hclust(dist_tr1276b16, method = "ward.D2")
plot(hc_tr1276b16)

coph_1276b16 <- cophenetic(hc_tr1276b16)
cor(dist_tr1276b16, coph_1276b16) #0.77, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b16, FUN = hcut, method = "silhouette") #optimal number of clusters is 5


# Ward Hierarchical Clustering
plot(hc_tr1276b16) # display dendogram
# draw dendogram with red borders around the 5 clusters 
rect.hclust(hc_tr1276b16, k=5, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b16 <- pvclust(t(mat_tr1276b16), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b16)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b16, alpha=.95)

dendro.1276b16 <- as.dendrogram(hc_tr1276b16, k = 3)
dendro.plot.1276b16 <- ggdendrogram(data = dendro.1276b16, rotate = TRUE, k = 2)
dendro.plot.1276b16 <- dendro.plot.1276b16 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b16)

melt.1276b16 <- melt(mat_tr1276b16)
head(melt.1276b16)
colnames(melt.1276b16) <- c("Targeted", "Bycatch", "TI")
melt.1276b16$TI_Ranges <- cut(melt.1276b16$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b16 <- order.dendrogram(dendro.1276b16)
melt.1276b16$Targeted <- factor(x = melt.1276b16$Targeted,
                                levels = melt.1276b16$Targeted[order.1276b16],
                                ordered = T)

heatmap.1276b16 <- ggplot(data=melt.1276b16, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b16)

grid.newpage()
print(heatmap.1276b16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.b 2017####

tr1_276b17 <- subset(tr1_276b, year == 2017)

unique(tr1_276b17$year)
unique(tr1_276b17$foCatNat)
unique(tr1_276b17$area)
tr1_276b17$spp2_id <- as.numeric(as.factor(tr1_276b17$spp))
spp_tr1276b17 <- sort(unique(tr1_276b17$spp2_id))


ti_tr1276b17 <- data.frame() 

for(i in spp_tr1276b17){
  sp1 <- subset(tr1_276b17, spp2_id == i) 
  for(j in spp_tr1276b17){
    t <- data.frame()
    sp2 <- subset(tr1_276b17, spp2_id == j) 
    if(i == j){ti_tr1276b17[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b17[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b17)

row.names(ti_tr1276b17) <- sort(unique(tr1_276b17$spp))
colnames(ti_tr1276b17) <- sort(unique(tr1_276b17$spp))
View(ti_tr1276b17)

mat_tr1276b17 <- as.matrix(ti_tr1276b17)
dim(mat_tr1276b17)

##HClust with TR1 27.6.a 2017
dist_tr1276b17 <- dist(mat_tr1276b17, method = "euclidean") #distance matrix
hc_tr1276b17 <- hclust(dist_tr1276b17, method = "ward.D2")
plot(hc_tr1276b17)

coph_1276b17 <- cophenetic(hc_tr1276b17)
cor(dist_tr1276b17, coph_1276b17) #0.90, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b17, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_tr1276b17) # display dendogram 2 clusters 
rect.hclust(hc_tr1276b17, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b17 <- pvclust(t(mat_tr1276b17), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b17)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b17, alpha=.95)

dendro.1276b17 <- as.dendrogram(hc_tr1276b17, k = 2)
dendro.plot.1276b17 <- ggdendrogram(data = dendro.1276b17, rotate = TRUE, k = 2)
dendro.plot.1276b17 <- dendro.plot.1276b17 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b17)

melt.1276b17 <- melt(mat_tr1276b17)
head(melt.1276b17)
colnames(melt.1276b17) <- c("Targeted", "Bycatch", "TI")
melt.1276b17$TI_Ranges <- cut(melt.1276b17$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b17 <- order.dendrogram(dendro.1276b17)
melt.1276b17$Targeted <- factor(x = melt.1276b17$Targeted,
                                levels = melt.1276b17$Targeted[order.1276b17],
                                ordered = T)

heatmap.1276b17 <- ggplot(data=melt.1276b17, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b17)

grid.newpage()
print(heatmap.1276b17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###TR1 27.6.b 2018####

tr1_276b18 <- subset(tr1_276b, year == 2018)

unique(tr1_276b18$year)
unique(tr1_276b18$foCatNat)
unique(tr1_276b18$area)
tr1_276b18$spp2_id <- as.numeric(as.factor(tr1_276b18$spp))
spp_tr1276b18 <- sort(unique(tr1_276b18$spp2_id))


ti_tr1276b18 <- data.frame() 

for(i in spp_tr1276b18){
  sp1 <- subset(tr1_276b18, spp2_id == i) 
  for(j in spp_tr1276b18){
    t <- data.frame()
    sp2 <- subset(tr1_276b18, spp2_id == j) 
    if(i == j){ti_tr1276b18[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b18[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b18)

row.names(ti_tr1276b18) <- sort(unique(tr1_276b18$spp))
colnames(ti_tr1276b18) <- sort(unique(tr1_276b18$spp))
View(ti_tr1276b18)

mat_tr1276b18 <- as.matrix(ti_tr1276b18)
dim(mat_tr1276b18)

##HClust with TR1 27.6.a 2018
dist_tr1276b18 <- dist(mat_tr1276b18, method = "euclidean") #distance matrix
hc_tr1276b18 <- hclust(dist_tr1276b18, method = "ward.D2")
plot(hc_tr1276b18)

coph_1276b18 <- cophenetic(hc_tr1276b18)
cor(dist_tr1276b18, coph_1276b18) #0.72, less than 0.75 so that's not good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_tr1276b18, FUN = hcut, method = "silhouette") #optimal number of clusters is 7


# Ward Hierarchical Clustering
plot(hc_tr1276b18) # display dendogram
# draw dendogram with red borders around the 7 clusters 
rect.hclust(hc_tr1276b18, k=7, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit1276b18 <- pvclust(t(mat_tr1276b18), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fit1276b18)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit1276b18, alpha=.95)

dendro.1276b18 <- as.dendrogram(hc_tr1276b18)
dendro.plot.1276b18 <- ggdendrogram(data = dendro.1276b18, rotate = TRUE)
dendro.plot.1276b18 <- dendro.plot.1276b18 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.1276b18)

melt.1276b18 <- melt(mat_tr1276b18)
head(melt.1276b18)
colnames(melt.1276b18) <- c("Targeted", "Bycatch", "TI")
melt.1276b18$TI_Ranges <- cut(melt.1276b18$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.1276b18 <- order.dendrogram(dendro.1276b18)
melt.1276b18$Targeted <- factor(x = melt.1276b18$Targeted,
                                levels = melt.1276b18$Targeted[order.1276b18],
                                ordered = T)

heatmap.1276b18 <- ggplot(data=melt.1276b18, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.1276b18)

grid.newpage()
print(heatmap.1276b18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###Compiling TR1 27.6.b Heatmaps#####
pdf(file = "TR1_276b_Heatmaps.pdf")

grid.newpage()
print(heatmap.1276b, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

dev.off()
###Heatmaps for Averaging#####
ti_tr1276b09.2 <- data.frame() 
ti_tr1276b09.2[1:39, 1:39] <- -1

spp_tr1276b09.2 <- unique(tr1_276b09$spp_id)
spp_tr1276b09.2 <- sort(spp_tr1276b09.2)

for(i in spp_tr1276b09.2){
  sp1 <- subset(tr1_276b09, spp_id == i) 
  for(j in spp_tr1276b09.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b09, spp_id == j) 
    if(i == j){ti_tr1276b09.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b09.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b09.2)

row.names(ti_tr1276b09.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b09.2) <- sort(unique(trip382$spp))

mat_tr1276b09.2 <- as.matrix(ti_tr1276b09.2)
dim(mat_tr1276b09.2)

##2010
ti_tr1276b10.2 <- data.frame()
ti_tr1276b10.2[1:39, 1:39] <- -1
# 
# spp_tr1276a10.2 <- unique(tr1_276a10$spp_id)
# spp_tr1276a10.2 <- sort(spp_tr1276a10.2)
# 
# for(i in spp_tr1276a10.2){
#   sp1 <- subset(tr1_276a10, spp_id == i) 
#   for(j in spp_tr1276a10.2){
#     t <- data.frame()
#     sp2 <- subset(tr1_276a10, spp_id == j) 
#     if(i == j){ti_tr1276a10.2[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_tr1276a10.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_tr1276a10.2)
# 
row.names(ti_tr1276b10.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b10.2) <- sort(unique(trip382$spp))

mat_tr1276b10.2 <- as.matrix(ti_tr1276b10.2)
dim(mat_tr1276b10.2)

##2011
ti_tr1276b11.2 <- data.frame() 
ti_tr1276b11.2[1:39, 1:39] <- -1

spp_tr1276b11.2 <- unique(tr1_276b11$spp_id)
spp_tr1276b11.2 <- sort(spp_tr1276b11.2)

for(i in spp_tr1276b11.2){
  sp1 <- subset(tr1_276b11, spp_id == i) 
  for(j in spp_tr1276b11.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b11, spp_id == j) 
    if(i == j){ti_tr1276b11.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b11.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b11.2)

row.names(ti_tr1276b11.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b11.2) <- sort(unique(trip382$spp))

mat_tr1276b11.2 <- as.matrix(ti_tr1276b11.2)
dim(mat_tr1276b11.2)

##2012
ti_tr1276b12.2 <- data.frame() 
ti_tr1276b12.2[1:39, 1:39] <- -1

spp_tr1276b12.2 <- unique(tr1_276b12$spp_id)
spp_tr1276b12.2 <- sort(spp_tr1276b12.2)

for(i in spp_tr1276b12.2){
  sp1 <- subset(tr1_276b12, spp_id == i) 
  for(j in spp_tr1276b12.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b12, spp_id == j) 
    if(i == j){ti_tr1276b12.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b12.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b12.2)

row.names(ti_tr1276b12.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b12.2) <- sort(unique(trip382$spp))

mat_tr1276b12.2 <- as.matrix(ti_tr1276b12.2)
dim(mat_tr1276b12.2)

##2013
ti_tr1276b13.2 <- data.frame() 
ti_tr1276b13.2[1:39, 1:39] <- -1

spp_tr1276b13.2 <- unique(tr1_276b13$spp_id)
spp_tr1276b13.2 <- sort(spp_tr1276b13.2)

for(i in spp_tr1276b13.2){
  sp1 <- subset(tr1_276b13, spp_id == i) 
  for(j in spp_tr1276b13.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b13, spp_id == j) 
    if(i == j){ti_tr1276b13.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b13.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b13.2)

row.names(ti_tr1276b13.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b13.2) <- sort(unique(trip382$spp))

mat_tr1276b13.2 <- as.matrix(ti_tr1276b13.2)
dim(mat_tr1276b13.2)

##2014
ti_tr1276b14.2 <- data.frame() 
ti_tr1276b14.2[1:39, 1:39] <- -1

spp_tr1276b14.2 <- unique(tr1_276b14$spp_id)
spp_tr1276b14.2 <- sort(spp_tr1276b14.2)

for(i in spp_tr1276b14.2){
  sp1 <- subset(tr1_276b14, spp_id == i) 
  for(j in spp_tr1276b14.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b14, spp_id == j) 
    if(i == j){ti_tr1276b14.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b14.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b14.2)

row.names(ti_tr1276b14.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b14.2) <- sort(unique(trip382$spp))

mat_tr1276b14.2 <- as.matrix(ti_tr1276b14.2)
dim(mat_tr1276b14.2)

##2015
ti_tr1276b15.2 <- data.frame() 
ti_tr1276b15.2[1:39, 1:39] <- -1

spp_tr1276b15.2 <- unique(tr1_276b15$spp_id)
spp_tr1276b15.2 <- sort(spp_tr1276b15.2)

for(i in spp_tr1276b15.2){
  sp1 <- subset(tr1_276b15, spp_id == i) 
  for(j in spp_tr1276b15.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b15, spp_id == j) 
    if(i == j){ti_tr1276b15.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b15.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b15.2)

row.names(ti_tr1276b15.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b15.2) <- sort(unique(trip382$spp))

mat_tr1276b15.2 <- as.matrix(ti_tr1276b15.2)
dim(mat_tr1276b15.2)

##2016
ti_tr1276b16.2 <- data.frame() 
ti_tr1276b16.2[1:39, 1:39] <- -1

spp_tr1276b16.2 <- unique(tr1_276b16$spp_id)
spp_tr1276b16.2 <- sort(spp_tr1276b16.2)

for(i in spp_tr1276b16.2){
  sp1 <- subset(tr1_276b16, spp_id == i) 
  for(j in spp_tr1276b16.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b16, spp_id == j) 
    if(i == j){ti_tr1276b16.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b16.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b16.2)

row.names(ti_tr1276b16.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b16.2) <- sort(unique(trip382$spp))

mat_tr1276b16.2 <- as.matrix(ti_tr1276b16.2)
dim(mat_tr1276b16.2)

##2017
ti_tr1276b17.2 <- data.frame() 
ti_tr1276b17.2[1:39, 1:39] <- -1

spp_tr1276b17.2 <- unique(tr1_276b17$spp_id)
spp_tr1276b17.2 <- sort(spp_tr1276b17.2)

for(i in spp_tr1276b17.2){
  sp1 <- subset(tr1_276b17, spp_id == i) 
  for(j in spp_tr1276b17.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b17, spp_id == j) 
    if(i == j){ti_tr1276b17.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b17.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b17.2)

row.names(ti_tr1276b17.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b17.2) <- sort(unique(trip382$spp))

mat_tr1276b17.2 <- as.matrix(ti_tr1276b17.2)
dim(mat_tr1276b17.2)

##2018
ti_tr1276b18.2 <- data.frame() 
ti_tr1276b18.2[1:39, 1:39] <- -1

spp_tr1276b18.2 <- unique(tr1_276b18$spp_id)
spp_tr1276b18.2 <- sort(spp_tr1276b18.2)

for(i in spp_tr1276b18.2){
  sp1 <- subset(tr1_276b18, spp_id == i) 
  for(j in spp_tr1276b18.2){
    t <- data.frame()
    sp2 <- subset(tr1_276b18, spp_id == j) 
    if(i == j){ti_tr1276b18.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_tr1276b18.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_tr1276b18.2)

row.names(ti_tr1276b18.2) <- sort(unique(trip382$spp))
colnames(ti_tr1276b18.2) <- sort(unique(trip382$spp))

mat_tr1276b18.2 <- as.matrix(ti_tr1276b18.2)
dim(mat_tr1276b18.2)
##Averaging the Matrices together
matavg_tr1276b <- (mat_tr1276b09.2 + mat_tr1276b11.2 + 
                     mat_tr1276b12.2 + mat_tr1276b13.2 + mat_tr1276b14.2 + 
                     mat_tr1276b15.2 + mat_tr1276b16.2 +  mat_tr1276b17.2 + mat_tr1276b18.2)/9


dist_avg1276b <- dist(matavg_tr1276b, method = "euclidean") #distance matrix
hc_avg1276b <- hclust(dist_avg1276b, method = "ward.D2")
plot(hc_avg1276b)

coph_avg1276b <- cophenetic(hc_avg1276b)
cor(dist_avg1276b, coph_avg1276b) #0.85, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(matavg_tr1276b, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_avg1276b) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_avg1276b, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitavg1276b <- pvclust(t(matavg_tr1276b), method.hclust="ward.D2",
                       method.dist="euclidean")
plot(fitavg1276b)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitavg1276b, alpha=.95)

dendro.avg1276b <- as.dendrogram(hc_avg1276b, k = 2)
dendro.plot.avg1276b <- ggdendrogram(data = dendro.avg1276b, rotate = TRUE, k = 2)
dendro.plot.avg1276b <- dendro.plot.avg1276b + theme(axis.text.y = element_blank(), 
                                                     axis.text.x.bottom =  element_blank())

print(dendro.plot.avg1276b)

melt.avg1276b <- melt(matavg_tr1276b)
head(melt.avg1276b)
colnames(melt.avg1276b) <- c("Targeted", "Bycatch", "TI")
melt.avg1276b$TI_Ranges <- cut(melt.avg1276b$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.avg1276b <- order.dendrogram(dendro.avg1276b)
melt.avg1276b$Targeted <- factor(x = melt.avg1276b$Targeted,
                                 levels = melt.avg1276b$Targeted[order.avg1276b],
                                 ordered = T)

heatmap.avg1276b <- ggplot(data=melt.avg1276b, 
                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.avg1276b)

grid.newpage()
print(heatmap.avg1276b, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.avg1276b, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

####Creating Array TR1 27.6.b####
a_tr1276b <- abind(mat_tr1276b09.2, mat_tr1276b10.2, mat_tr1276b11.2, 
                   mat_tr1276b12.2, mat_tr1276b13.2, mat_tr1276b14.2,
                   mat_tr1276b15.2, mat_tr1276b16.2, mat_tr1276b17.2, mat_tr1276b18.2, along = 3)
View(a_tr1276b[,,9])

melta_tr1276b <- melt(a_tr1276b)
head(melta_tr1276b)
colnames(melta_tr1276b) <- c("Targeted", "Bycatch", "Year", "TIS")

melta_tr1276b$TI_Ranges <- cut(melta_tr1276b$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
str(melta_tr1276b)


melta_tr1276b$Area[1:15210] <- "27.6.b"
melta_tr1276b$Area <- as.factor(melta_tr1276b$Area)
melta_tr1276b$Gear[1:15210] <- "TR1"
melta_tr1276b$Gear <- as.factor(melta_tr1276b$Gear)

melta_tr1276b$Year <- as.factor(melta_tr1276b$Year)
melta_tr1276b$Year <- revalue(melta_tr1276b$Year, c("1"="2009", "2"="2010", "3"="2011", "4"="2012",
                                                    "5"="2013", "6"="2014", "7"="2015", "8"="2016",
                                                    "9"="2017", "10"="2018"))

melta_tr1276b <- na.omit(melta_tr1276b) #removing absent species that are marked as NA in TI_Ranges

tr1276b_by <- with(melta_tr1276b, tapply(TIS, INDEX = list(Bycatch, Year), FUN = "mean"))
tr1276b_tar <- with(melta_tr1276b, tapply(TIS, INDEX = list(Targeted, Year), FUN = "mean"))
tr1276b_overall <- with(melta_tr1276b, tapply (TIS, INDEX = Year, FUN = "mean"))
tr1276b_by <- as.data.frame(tr1276b_by)
tr1276b_tar <- as.data.frame(tr1276b_tar)
tr1276b_overall <- as.data.frame(tr1276b_overall)

tr1276b_overall$var <- with(melta_tr1276b, tapply(TIS, INDEX = Year, FUN = "var"))
tr1276b_overall$sd <- sqrt(tr1276b_overall$var)
tr1276b_overall$cv <- 100*((tr1276b_overall$sd)/(tr1276b_overall$tr1276b_overall))

tr1276b_by$Overall <- rowMeans(tr1276b_by)
tr1276b_by$var <- with(melta_tr1276b, tapply(TIS, INDEX = list(Bycatch,Year), FUN = "var"))
tr1276b_by$sd <- sqrt(tr1276b_by$var)
tr1276b_by$cv <- 100*((tr1276b_by$sd)/(tr1276b_by$Overall))
tr1276b_by.05 <- tr1276b_by[!(tr1276b_by$Overall<=5),]
tr1276b_by.50 <- tr1276b_by[!(tr1276b_by$Overall<=50),]
View(tr1276b_by.05)
View(tr1276b_by.50)

tr1276b_tar$Overall <- rowMeans(tr1276b_tar)
tr1276b_tar$var <- with(melta_tr1276b, tapply(TIS, INDEX = list(Targeted,Year), FUN = "var"))
tr1276b_tar$sd <- sqrt(tr1276b_tar$var)
tr1276b_tar$cv <- 100*((tr1276b_tar$sd)/(tr1276b_tar$Overall))
tr1276b_tar.05 <- tr1276b_tar[!(tr1276b_tar$Overall<=5),]
tr1276b_tar.25 <- tr1276b_tar[!(tr1276b_tar$Overall<=25),]
View(tr1276b_tar.05)
View(tr1276b_tar.25)


####Centered Moving Average TR1 27.6.b--EACH YEAR#####
?rollmean
?select
# 
# test <- tr1274_overall %>% select(TIS_year = tr1274_overall) %>% 
#   mutate(TIS_mavg3 = rollmean(TIS_year, k = 3, fill = NA),
#          Year = 2009:2018)
# testg<- test %>% gather(metric, value, TIS_year:TIS_mavg3) %>%
#   ggplot(aes(Year, value, color = metric)) +
#   geom_line()
# print(testg)

t1276b <- data.frame()
for(i in 2009:2018){
  selected <- c(i-1, i, i+1)
  f1 <- subset(melta_tr1276b, Year %in% selected)
  if((i == 2009)| (i==2018)){
    t1276b[i-2008, 1] <- NA} 
  else{t1276b[i-2008,1] <- mean(f1$TIS)}
}
print(t1276b) #The 2010 Mavg is just the avg of 2009 and 2011

print(tr1276b_overall)
tr1276b_overall[2,] <- c(NA,0,0,0)


t1276b.2 <- t1276b %>% select(TIS_mavg = V1) %>% mutate(TIS_yavg = tr1276b_overall$tr1276b_overall,
                                                        Year = 2009:2018)
t1276b.g <- t1276b.2 %>% gather(Metric, TIS, TIS_mavg:TIS_yavg) %>%
  ggplot(aes(Year, TIS, color = Metric)) +
  geom_line()
print(t1276b.g)

##########################################################################
#####Area 27.6.b, TRO Subsets#################

trO_276b <- subset(trip382, gear_id == 5 & area_id == 7)
trO_276b

unique(trO_276b$year)
unique(trO_276b$foCatNat)
unique(trO_276b$area)
trO_276b$spp_id <- trO_276b$spp2_id
trO_276b$spp2_id <- as.numeric(as.factor(trO_276b$spp))
spp_trO276b <- unique(trO_276b$spp2_id)
spp_trO276b <- sort(spp_trO276b)

ti_trO276b <- data.frame() 

for(i in spp_trO276b){
  sp1 <- subset(trO_276b, spp2_id == i) 
  for(j in spp_trO276b){
    t <- data.frame()
    sp2 <- subset(trO_276b, spp2_id == j) 
    if(i == j){ti_trO276b[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$year[k] %in% sp2$year)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b)

row.names(ti_trO276b) <- sort(unique(trO_276b$spp))
colnames(ti_trO276b) <- sort(unique(trO_276b$spp))
View(ti_trO276b)

mat_trO276b <- as.matrix(ti_trO276b)
##HClust with trO 27.6.b
dist_trO276b <- dist(mat_trO276b, method = "euclidean") #distance matrix
hc_trO276b<- hclust(dist_trO276b, method = "ward.D2")
plot(hc_trO276b)

coph_O276b <- cophenetic(hc_trO276b)
cor(dist_trO276b, coph_O276b) #0.94, greater than 0.75, therefore good

fviz_nbclust(mat_trO276b, FUN = hcut, method = "silhouette") #optimal number of clusters is 4

# Ward Hierarchical Clustering
plot(hc_trO276b) # display dendogram
# draw dendogram with red borders around the 4 clusters 
rect.hclust(hc_trO276b, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fit276b <- pvclust(t(mat_trO276b), method.hclust="ward.D2",
                   method.dist="euclidean")
plot(fit276b)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit276b, alpha=.95)

dendro.O276b <- as.dendrogram(hc_trO276b)
dendro.plot.O276b <- ggdendrogram(data = dendro.O276b, rotate = TRUE)
dendro.plot.O276b <- dendro.plot.O276b + theme(axis.text.y = element_blank(), axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b)

melt.O276b <- melt(mat_trO276b)
head(melt.O276b)
colnames(melt.O276b) <- c("Targeted", "Bycatch", "TI")
melt.O276b$TI_Ranges <- cut(melt.O276b$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b <- order.dendrogram(dendro.O276b)
melt.O276b$Targeted <- factor(x = melt.O276b$Targeted,
                              levels = melt.O276b$Targeted[order.O276b],
                              ordered = T)


heatmap.O276b <- ggplot(data=melt.O276b, 
                        aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b)

grid.newpage()
print(heatmap.O276b, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

###trO 27.6.b 2009########
trO_276b09 <- subset(trO_276b, year == 2009)

unique(trO_276b09$year)
unique(trO_276b09$foCatNat)
unique(trO_276b09$area)
trO_276b09$spp2_id <- as.numeric(as.factor(trO_276b09$spp))
spp_trO276b09 <- unique(trO_276b09$spp2_id)
spp_trO276b09 <- sort(spp_trO276b09)

ti_trO276b09 <- data.frame() 

for(i in spp_trO276b09){
  sp1 <- subset(trO_276b09, spp2_id == i) 
  for(j in spp_trO276b09){
    t <- data.frame()
    sp2 <- subset(trO_276b09, spp2_id == j) 
    if(i == j){ti_trO276b09[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b09[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b09)

row.names(ti_trO276b09) <- sort(unique(trO_276b09$spp))
colnames(ti_trO276b09) <- sort(unique(trO_276b09$spp))
View(ti_trO276b09)

mat_trO276b09 <- as.matrix(ti_trO276b09)
dim(mat_trO276b09)
##HClust with trO 27.6.a 2009
dist_trO276b09 <- dist(mat_trO276b09, method = "euclidean") #distance matrix
hc_trO276b09<- hclust(dist_trO276b09, method = "ward.D2")
plot(hc_trO276b09)

coph_O276b09 <- cophenetic(hc_trO276b09)
cor(dist_trO276b09, coph_O276b09) #Greater than 0.99, therefore considered good


fviz_nbclust(mat_trO276b09, FUN = hcut, method = "silhouette") #optimal number of clusters is 6


# Ward Hierarchical Clustering
plot(hc_trO276b09) # display dendogram
# draw dendogram with red borders around the 6 clusters 
rect.hclust(hc_trO276b09, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitO276b09 <- pvclust(t(mat_trO276b09), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitO276b09)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitO276b09, alpha=.95)

dendro.O276b09 <- as.dendrogram(hc_trO276b09)
dendro.plot.O276b09 <- ggdendrogram(data = dendro.O276b09, rotate = TRUE)
dendro.plot.O276b09 <- dendro.plot.O276b09 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b09)

melt.O276b09 <- melt(mat_trO276b09)
head(melt.O276b09)
colnames(melt.O276b09) <- c("Targeted", "Bycatch", "TI")
melt.O276b09$TI_Ranges <- cut(melt.O276b09$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b09 <- order.dendrogram(dendro.O276b09)
melt.O276b09$Targeted <- factor(x = melt.O276b09$Targeted,
                                levels = melt.O276b09$Targeted[order.O276b09],
                                ordered = T)

heatmap.O276b09 <- ggplot(data=melt.O276b09, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b09)

grid.newpage()
print(heatmap.O276b09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

##
###NO DATA FOR 27.6.b 2010##########
# trO_276a10 <- subset(trO_276a, year == 2010)
# 
# unique(trO_276a10$year)
# unique(trO_276a10$foCatNat)
# unique(trO_276a10$area)
# trO_276a10$spp2_id <- as.numeric(as.factor(trO_276a10$spp))
# spp_trO276a10 <- unique(trO_276a10$spp2_id)
# spp_trO276a10 <- sort(spp_trO276a10)
# 
# ti_trO276a10 <- data.frame() 
# 
# for(i in spp_trO276a10){
#   sp1 <- subset(trO_276a10, spp2_id == i) 
#   for(j in spp_trO276a10){
#     t <- data.frame()
#     sp2 <- subset(trO_276a10, spp2_id == j) 
#     if(i == j){ti_trO276a10[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_trO276a10[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_trO276a10)
# 
# row.names(ti_trO276a10) <- sort(unique(trO_276a10$spp))
# colnames(ti_trO276a10) <- sort(unique(trO_276a10$spp))
# View(ti_trO276a10)
# 
# mat_trO276a10 <- as.matrix(ti_trO276a10)
# dim(mat_trO276a10)
# 
# 
# ##HClust with trO 27.6.a 2010
# dist_trO276a10 <- dist(mat_trO276a10, method = "euclidean") #distance matrix
# hc_trO276a10 <- hclust(dist_trO276a10, method = "ward.D2")
# plot(hc_trO276a10)
# 
# coph_1276a10 <- cophenetic(hc_trO276a10)
# cor(dist_trO276a10, coph_1276a10) #0.87, greater than 0.75, optimal
# 
# fviz_nbclust(mat_trO276a10, FUN = hcut, method = "silhouette") #optimal number of clusters is 2
# 
# 
# # Ward Hierarchical Clustering
# plot(hc_trO276a10) # display dendogram
# # draw dendogram with red borders around the 2 clusters 
# rect.hclust(hc_trO276a10, k=2, border="red") 
# 
# # Ward Hierarchical Clustering with Bootstrapped p values
# fit1276a10 <- pvclust(t(mat_trO276a10), method.hclust="ward.D2",
#                       method.dist="euclidean")
# plot(fit1276a10)# dendogram with p values
# # add rectangles around groups highly supported by the data
# pvrect(fit1276a10, alpha=.95)
# 
# dendro.1276a10 <- as.dendrogram(hc_trO276a10, k = 2)
# dendro.plot.1276a10 <- ggdendrogram(data = dendro.1276a10, rotate = TRUE, k = 2)
# dendro.plot.1276a10 <- dendro.plot.1276a10 + theme(axis.text.y = element_blank(), 
#                                                    axis.text.x.bottom =  element_blank())
# 
# print(dendro.plot.1276a10)
# 
# melt.1276a10 <- melt(mat_trO276a10)
# head(melt.1276a10)
# colnames(melt.1276a10) <- c("Targeted", "Bycatch", "TI")
# melt.1276a10$TI_Ranges <- cut(melt.1276a10$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# order.1276a10 <- order.dendrogram(dendro.1276a10)
# melt.1276a10$Targeted <- factor(x = melt.1276a10$Targeted,
#                                 levels = melt.1276a10$Targeted[order.1276a10],
#                                 ordered = T)
# 
# heatmap.1276a10 <- ggplot(data=melt.1276a10, 
#                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)")
# print(heatmap.1276a10)
# 
# grid.newpage()
# print(heatmap.1276a10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot.1276a10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))


##

###trO 27.6.b 2011#######

trO_276b11 <- subset(trO_276b, year == 2011)

unique(trO_276b11$year)
unique(trO_276b11$foCatNat)
unique(trO_276b11$area)
trO_276b11$spp2_id <- as.numeric(as.factor(trO_276b11$spp))
spp_trO276b11 <- sort(unique(trO_276b11$spp2_id))


ti_trO276b11 <- data.frame() 

for(i in spp_trO276b11){
  sp1 <- subset(trO_276b11, spp2_id == i) 
  for(j in spp_trO276b11){
    t <- data.frame()
    sp2 <- subset(trO_276b11, spp2_id == j) 
    if(i == j){ti_trO276b11[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b11[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b11)

row.names(ti_trO276b11) <- sort(unique(trO_276b11$spp))
colnames(ti_trO276b11) <- sort(unique(trO_276b11$spp))
View(ti_trO276b11)

mat_trO276b11 <- as.matrix(ti_trO276b11)
dim(mat_trO276b11)

##HClust with trO 27.4 2011
dist_trO276b11 <- dist(mat_trO276b11, method = "euclidean") #distance matrix
hc_trO276b11 <- hclust(dist_trO276b11, method = "ward.D2")
plot(hc_trO276b11)

coph_O276b11 <- cophenetic(hc_trO276b11)
cor(dist_trO276b11, coph_O276b11) #0.87, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_trO276b11, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_trO276b11) # display dendogram
# draw dendogram with red borders around the 4 clusters 
rect.hclust(hc_trO276b11, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitO276b11 <- pvclust(t(mat_trO276b11), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitO276b11)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitO276b11, alpha=.95)

dendro.O276b11 <- as.dendrogram(hc_trO276b11)
dendro.plot.O276b11 <- ggdendrogram(data = dendro.O276b11, rotate = TRUE)
dendro.plot.O276b11 <- dendro.plot.O276b11 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b11)

melt.O276b11 <- melt(mat_trO276b11)
head(melt.O276b11)
colnames(melt.O276b11) <- c("Targeted", "Bycatch", "TI")
melt.O276b11$TI_Ranges <- cut(melt.O276b11$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b11 <- order.dendrogram(dendro.O276b11)
melt.O276b11$Targeted <- factor(x = melt.O276b11$Targeted,
                                levels = melt.O276b11$Targeted[order.O276b11],
                                ordered = T)

heatmap.O276b11 <- ggplot(data=melt.O276b11, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b11)

grid.newpage()
print(heatmap.O276b11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###trO 27.6.b 2012######

trO_276b12 <- subset(trO_276b, year == 2012)

unique(trO_276b12$year)
unique(trO_276b12$foCatNat)
unique(trO_276b12$area)
trO_276b12$spp2_id <- as.numeric(as.factor(trO_276b12$spp))
spp_trO276b12 <- sort(unique(trO_276b12$spp2_id))


ti_trO276b12 <- data.frame() 

for(i in spp_trO276b12){
  sp1 <- subset(trO_276b12, spp2_id == i) 
  for(j in spp_trO276b12){
    t <- data.frame()
    sp2 <- subset(trO_276b12, spp2_id == j) 
    if(i == j){ti_trO276b12[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b12[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b12)

row.names(ti_trO276b12) <- sort(unique(trO_276b12$spp))
colnames(ti_trO276b12) <- sort(unique(trO_276b12$spp))
View(ti_trO276b12)

mat_trO276b12 <- as.matrix(ti_trO276b12)
dim(mat_trO276b12)

##HClust with trO 27.6.b 2012
dist_trO276b12 <- dist(mat_trO276b12, method = "euclidean") #distance matrix
hc_trO276b12 <- hclust(dist_trO276b12, method = "ward.D2")
plot(hc_trO276b12)

coph_O276b12 <- cophenetic(hc_trO276b12)
cor(dist_trO276b12, coph_O276b12) #0.79, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_trO276b12, FUN = hcut, method = "silhouette") #optimal number of clusters is 9


# Ward Hierarchical Clustering
plot(hc_trO276b12) # display dendogram
# draw dendogram with red borders around the 9 clusters 
rect.hclust(hc_trO276b12, k=9, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitO276b12 <- pvclust(t(mat_trO276b12), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitO276b12)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitO276b12, alpha=.95)

dendro.O276b12 <- as.dendrogram(hc_trO276b12, k = 2)
dendro.plot.O276b12 <- ggdendrogram(data = dendro.O276b12, rotate = TRUE, k = 3)
dendro.plot.O276b12 <- dendro.plot.O276b12 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b12)

melt.O276b12 <- melt(mat_trO276b12)
head(melt.O276b12)
colnames(melt.O276b12) <- c("Targeted", "Bycatch", "TI")
melt.O276b12$TI_Ranges <- cut(melt.O276b12$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b12 <- order.dendrogram(dendro.O276b12)
melt.O276b12$Targeted <- factor(x = melt.O276b12$Targeted,
                                levels = melt.O276b12$Targeted[order.O276b12],
                                ordered = T)

heatmap.O276b12 <- ggplot(data=melt.O276b12, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b12)

grid.newpage()
print(heatmap.O276b12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###trO 27.6.b 2013#####

trO_276b13 <- subset(trO_276b, year == 2013)

unique(trO_276b13$year)
unique(trO_276b13$foCatNat)
unique(trO_276b13$area)
trO_276b13$spp2_id <- as.numeric(as.factor(trO_276b13$spp))
spp_trO276b13 <- sort(unique(trO_276b13$spp2_id))


ti_trO276b13 <- data.frame() 

for(i in spp_trO276b13){
  sp1 <- subset(trO_276b13, spp2_id == i) 
  for(j in spp_trO276b13){
    t <- data.frame()
    sp2 <- subset(trO_276b13, spp2_id == j) 
    if(i == j){ti_trO276b13[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b13[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b13)

row.names(ti_trO276b13) <- sort(unique(trO_276b13$spp))
colnames(ti_trO276b13) <- sort(unique(trO_276b13$spp))
View(ti_trO276b13)

mat_trO276b13 <- as.matrix(ti_trO276b13)
dim(mat_trO276b13)

##HClust with trO 27.6.b 2013
dist_trO276b13 <- dist(mat_trO276b13, method = "euclidean") #distance matrix
hc_trO276b13 <- hclust(dist_trO276b13, method = "ward.D2")
plot(hc_trO276b13)

coph_O276b13 <- cophenetic(hc_trO276b13)
cor(dist_trO276b13, coph_O276b13) #0.77, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_trO276b13, FUN = hcut, method = "silhouette") #optimal number of clusters is 3


# Ward Hierarchical Clustering
plot(hc_trO276b13) # display dendogram
# draw dendogram with red borders around the 3 clusters 
rect.hclust(hc_trO276b13, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitO276b13 <- pvclust(t(mat_trO276b13), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitO276b13)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitO276b13, alpha=.95)

dendro.O276b13 <- as.dendrogram(hc_trO276b13, k = 2)
dendro.plot.O276b13 <- ggdendrogram(data = dendro.O276b13, rotate = TRUE, k = 2)
dendro.plot.O276b13 <- dendro.plot.O276b13 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b13)

melt.O276b13 <- melt(mat_trO276b13)
head(melt.O276b13)
colnames(melt.O276b13) <- c("Targeted", "Bycatch", "TI")
melt.O276b13$TI_Ranges <- cut(melt.O276b13$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b13 <- order.dendrogram(dendro.O276b13)
melt.O276b13$Targeted <- factor(x = melt.O276b13$Targeted,
                                levels = melt.O276b13$Targeted[order.O276b13],
                                ordered = T)

heatmap.O276b13 <- ggplot(data=melt.O276b13, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b13)

grid.newpage()
print(heatmap.O276b13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

###ALL 100% trO 27.6.b 2014#####

trO_276b14 <- subset(trO_276b, year == 2014)

unique(trO_276b14$year)
unique(trO_276b14$foCatNat)
unique(trO_276b14$area)
trO_276b14$spp2_id <- as.numeric(as.factor(trO_276b14$spp))
spp_trO276b14 <- sort(unique(trO_276b14$spp2_id))


ti_trO276b14 <- data.frame() 

for(i in spp_trO276b14){
  sp1 <- subset(trO_276b14, spp2_id == i) 
  for(j in spp_trO276b14){
    t <- data.frame()
    sp2 <- subset(trO_276b14, spp2_id == j) 
    if(i == j){ti_trO276b14[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b14[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b14)

row.names(ti_trO276b14) <- sort(unique(trO_276b14$spp))
colnames(ti_trO276b14) <- sort(unique(trO_276b14$spp))
View(ti_trO27414)

mat_trO276b14 <- as.matrix(ti_trO276b14)
dim(mat_trO276b14)

##HClust with trO 27.6.a 2014
dist_trO276b14 <- dist(mat_trO276b14, method = "euclidean") #distance matrix
hc_trO276b14 <- hclust(dist_trO276b14, method = "ward.D2")
plot(hc_trO276b14)

coph_O276b14 <- cophenetic(hc_trO276b14)
cor(dist_trO276b14, coph_O276b14) #0.95, greater than 0.75, optimal

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_trO276b14, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_trO276b14) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_trO276b14, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitO276b14 <- pvclust(t(mat_trO276b14), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitO276b14)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitO276b14, alpha=.95)

dendro.O276b14 <- as.dendrogram(hc_trO276b14)
dendro.plot.O276b14 <- ggdendrogram(data = dendro.O276b14, rotate = TRUE)
dendro.plot.O276b14 <- dendro.plot.O276b14 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b14)

melt.O276b14 <- melt(mat_trO276b14)
head(melt.O276b14)
colnames(melt.O276b14) <- c("Targeted", "Bycatch", "TI")
melt.O276b14$TI_Ranges <- cut(melt.O276b14$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b14 <- order.dendrogram(dendro.O276b14)
melt.O276b14$Targeted <- factor(x = melt.O276b14$Targeted,
                                levels = melt.O276b14$Targeted[order.O276b14],
                                ordered = T)

heatmap.O276b14 <- ggplot(data=melt.O276b14, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b14)

grid.newpage()
print(heatmap.O276b14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###ALL 100% trO 27.6.b 2015####

trO_276b15 <- subset(trO_276b, year == 2015)

unique(trO_276b15$year)
unique(trO_276b15$foCatNat)
unique(trO_276b15$area)
trO_276b15$spp2_id <- as.numeric(as.factor(trO_276b15$spp))
spp_trO276b15 <- sort(unique(trO_276b15$spp2_id))


ti_trO276b15 <- data.frame() 

for(i in spp_trO276b15){
  sp1 <- subset(trO_276b15, spp2_id == i) 
  for(j in spp_trO276b15){
    t <- data.frame()
    sp2 <- subset(trO_276b15, spp2_id == j) 
    if(i == j){ti_trO276b15[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b15[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b15)

row.names(ti_trO276b15) <- sort(unique(trO_276b15$spp))
colnames(ti_trO276b15) <- sort(unique(trO_276b15$spp))
View(ti_trO276b15)

mat_trO276b15 <- as.matrix(ti_trO276b15)
dim(mat_trO276b15)

##HClust with trO 27.6.b 2015
dist_trO276b15 <- dist(mat_trO276b15, method = "euclidean") #distance matrix
hc_trO276b15 <- hclust(dist_trO276b15, method = "ward.D2")
plot(hc_trO276b15)

coph_O276b15 <- cophenetic(hc_trO276b15)
cor(dist_trO276b15, coph_O276b15) #0.85, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(mat_trO276b15, FUN = hcut, method = "silhouette") #optimal number of clusters is 4


# Ward Hierarchical Clustering
plot(hc_trO276b15) # display dendogram
# draw dendogram with red borders around the 4 clusters 
rect.hclust(hc_trO276b15, k=4, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitO276b15 <- pvclust(t(mat_trO276b15), method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitO276b15)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitO276b15, alpha=.95)

dendro.O276b15 <- as.dendrogram(hc_trO276b15)
dendro.plot.O276b15 <- ggdendrogram(data = dendro.O276b15, rotate = TRUE)
dendro.plot.O276b15 <- dendro.plot.O276b15 + theme(axis.text.y = element_blank(), 
                                                   axis.text.x.bottom =  element_blank())

print(dendro.plot.O276b15)

melt.O276b15 <- melt(mat_trO276b15)
head(melt.O276b15)
colnames(melt.O276b15) <- c("Targeted", "Bycatch", "TI")
melt.O276b15$TI_Ranges <- cut(melt.O276b15$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.O276b15 <- order.dendrogram(dendro.O276b15)
melt.O276b15$Targeted <- factor(x = melt.O276b15$Targeted,
                                levels = melt.O276b15$Targeted[order.O276b15],
                                ordered = T)

heatmap.O276b15 <- ggplot(data=melt.O276b15, 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.O276b15)

grid.newpage()
print(heatmap.O276b15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.O276b15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###NO DATA trO 27.6.b 2016####

# trO_276b16 <- subset(tr1_276b, year == 2016)
# 
# unique(tr1_276b16$year)
# unique(tr1_276b16$foCatNat)
# unique(tr1_276b16$area)
# tr1_276b16$spp2_id <- as.numeric(as.factor(tr1_276b16$spp))
# spp_tr1276b16 <- sort(unique(tr1_276b16$spp2_id))
# 
# 
# ti_tr1276b16 <- data.frame() 
# for(i in spp_tr1276b16){
#   sp1 <- subset(tr1_276b16, spp2_id == i) 
#   for(j in spp_tr1276b16){
#     t <- data.frame()
#     sp2 <- subset(tr1_276b16, spp2_id == j) 
#     if(i == j){ti_tr1276b16[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_tr1276b16[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_tr1276b16)
# 
# row.names(ti_tr1276b16) <- sort(unique(tr1_276b16$spp))
# colnames(ti_tr1276b16) <- sort(unique(tr1_276b16$spp))
# View(ti_tr1276b16)
# 
# mat_tr1276b16 <- as.matrix(ti_tr1276b16)
# dim(mat_tr1276b16)
# 
# ##HClust with TR1 27.6.a 2016
# dist_tr1276b16 <- dist(mat_tr1276b16, method = "euclidean") #distance matrix
# hc_tr1276b16 <- hclust(dist_tr1276b16, method = "ward.D2")
# plot(hc_tr1276b16)
# 
# coph_1276b16 <- cophenetic(hc_tr1276b16)
# cor(dist_tr1276b16, coph_1276b16) #0.77, greater than 0.75 so that's good
# 
# #Using Silhouette method to determine optimal numer of clusters
# fviz_nbclust(mat_tr1276b16, FUN = hcut, method = "silhouette") #optimal number of clusters is 5
# 
# 
# # Ward Hierarchical Clustering
# plot(hc_tr1276b16) # display dendogram
# # draw dendogram with red borders around the 5 clusters 
# rect.hclust(hc_tr1276b16, k=5, border="red") 
# 
# # Ward Hierarchical Clustering with Bootstrapped p values
# fit1276b16 <- pvclust(t(mat_tr1276b16), method.hclust="ward.D2",
#                       method.dist="euclidean")
# plot(fit1276b16)# dendogram with p values
# # add rectangles around groups highly supported by the data
# pvrect(fit1276b16, alpha=.95)
# 
# dendro.1276b16 <- as.dendrogram(hc_tr1276b16, k = 3)
# dendro.plot.1276b16 <- ggdendrogram(data = dendro.1276b16, rotate = TRUE, k = 2)
# dendro.plot.1276b16 <- dendro.plot.1276b16 + theme(axis.text.y = element_blank(), 
#                                                    axis.text.x.bottom =  element_blank())
# 
# print(dendro.plot.1276b16)
# 
# melt.1276b16 <- melt(mat_tr1276b16)
# head(melt.1276b16)
# colnames(melt.1276b16) <- c("Targeted", "Bycatch", "TI")
# melt.1276b16$TI_Ranges <- cut(melt.1276b16$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# order.1276b16 <- order.dendrogram(dendro.1276b16)
# melt.1276b16$Targeted <- factor(x = melt.1276b16$Targeted,
#                                 levels = melt.1276b16$Targeted[order.1276b16],
#                                 ordered = T)
# 
# heatmap.1276b16 <- ggplot(data=melt.1276b16, 
#                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)")
# print(heatmap.1276b16)
# 
# grid.newpage()
# print(heatmap.1276b16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot.1276b16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
####NO DATA TRO27.6.b 2017####
# 
# tr1_276b17 <- subset(tr1_276b, year == 2017)
# 
# unique(tr1_276b17$year)
# unique(tr1_276b17$foCatNat)
# unique(tr1_276b17$area)
# tr1_276b17$spp2_id <- as.numeric(as.factor(tr1_276b17$spp))
# spp_tr1276b17 <- sort(unique(tr1_276b17$spp2_id))
# 
# 
# ti_tr1276b17 <- data.frame() 
# 
# for(i in spp_tr1276b17){
#   sp1 <- subset(tr1_276b17, spp2_id == i) 
#   for(j in spp_tr1276b17){
#     t <- data.frame()
#     sp2 <- subset(tr1_276b17, spp2_id == j) 
#     if(i == j){ti_tr1276b17[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_tr1276b17[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_tr1276b17)
# 
# row.names(ti_tr1276b17) <- sort(unique(tr1_276b17$spp))
# colnames(ti_tr1276b17) <- sort(unique(tr1_276b17$spp))
# View(ti_tr1276b17)
# 
# mat_tr1276b17 <- as.matrix(ti_tr1276b17)
# dim(mat_tr1276b17)
# 
# ##HClust with TR1 27.6.a 2017
# dist_tr1276b17 <- dist(mat_tr1276b17, method = "euclidean") #distance matrix
# hc_tr1276b17 <- hclust(dist_tr1276b17, method = "ward.D2")
# plot(hc_tr1276b17)
# 
# coph_1276b17 <- cophenetic(hc_tr1276b17)
# cor(dist_tr1276b17, coph_1276b17) #0.90, greater than 0.75 so that's good
# 
# #Using Silhouette method to determine optimal numer of clusters
# fviz_nbclust(mat_tr1276b17, FUN = hcut, method = "silhouette") #optimal number of clusters is 2
# 
# 
# # Ward Hierarchical Clustering
# plot(hc_tr1276b17) # display dendogram 2 clusters 
# rect.hclust(hc_tr1276b17, k=2, border="red") 
# 
# # Ward Hierarchical Clustering with Bootstrapped p values
# fit1276b17 <- pvclust(t(mat_tr1276b17), method.hclust="ward.D2",
#                       method.dist="euclidean")
# plot(fit1276b17)# dendogram with p values
# # add rectangles around groups highly supported by the data
# pvrect(fit1276b17, alpha=.95)
# 
# dendro.1276b17 <- as.dendrogram(hc_tr1276b17, k = 2)
# dendro.plot.1276b17 <- ggdendrogram(data = dendro.1276b17, rotate = TRUE, k = 2)
# dendro.plot.1276b17 <- dendro.plot.1276b17 + theme(axis.text.y = element_blank(), 
#                                                    axis.text.x.bottom =  element_blank())
# 
# print(dendro.plot.1276b17)
# 
# melt.1276b17 <- melt(mat_tr1276b17)
# head(melt.1276b17)
# colnames(melt.1276b17) <- c("Targeted", "Bycatch", "TI")
# melt.1276b17$TI_Ranges <- cut(melt.1276b17$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# order.1276b17 <- order.dendrogram(dendro.1276b17)
# melt.1276b17$Targeted <- factor(x = melt.1276b17$Targeted,
#                                 levels = melt.1276b17$Targeted[order.1276b17],
#                                 ordered = T)
# 
# heatmap.1276b17 <- ggplot(data=melt.1276b17, 
#                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)")
# print(heatmap.1276b17)
# 
# grid.newpage()
# print(heatmap.1276b17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot.1276b17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
####NO DATA TRO 27.6.b 2018####
# 
# tr1_276b18 <- subset(tr1_276b, year == 2018)
# 
# unique(tr1_276b18$year)
# unique(tr1_276b18$foCatNat)
# unique(tr1_276b18$area)
# tr1_276b18$spp2_id <- as.numeric(as.factor(tr1_276b18$spp))
# spp_tr1276b18 <- sort(unique(tr1_276b18$spp2_id))
# 
# 
# ti_tr1276b18 <- data.frame() 
# 
# for(i in spp_tr1276b18){
#   sp1 <- subset(tr1_276b18, spp2_id == i) 
#   for(j in spp_tr1276b18){
#     t <- data.frame()
#     sp2 <- subset(tr1_276b18, spp2_id == j) 
#     if(i == j){ti_tr1276b18[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_tr1276b18[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_tr1276b18)
# 
# row.names(ti_tr1276b18) <- sort(unique(tr1_276b18$spp))
# colnames(ti_tr1276b18) <- sort(unique(tr1_276b18$spp))
# View(ti_tr1276b18)
# 
# mat_tr1276b18 <- as.matrix(ti_tr1276b18)
# dim(mat_tr1276b18)
# 
# ##HClust with TR1 27.6.a 2018
# dist_tr1276b18 <- dist(mat_tr1276b18, method = "euclidean") #distance matrix
# hc_tr1276b18 <- hclust(dist_tr1276b18, method = "ward.D2")
# plot(hc_tr1276b18)
# 
# coph_1276b18 <- cophenetic(hc_tr1276b18)
# cor(dist_tr1276b18, coph_1276b18) #0.72, less than 0.75 so that's not good
# 
# #Using Silhouette method to determine optimal numer of clusters
# fviz_nbclust(mat_tr1276b18, FUN = hcut, method = "silhouette") #optimal number of clusters is 7
# 
# 
# # Ward Hierarchical Clustering
# plot(hc_tr1276b18) # display dendogram
# # draw dendogram with red borders around the 7 clusters 
# rect.hclust(hc_tr1276b18, k=7, border="red") 
# 
# # Ward Hierarchical Clustering with Bootstrapped p values
# fit1276b18 <- pvclust(t(mat_tr1276b18), method.hclust="ward.D2",
#                       method.dist="euclidean")
# plot(fit1276b18)# dendogram with p values
# # add rectangles around groups highly supported by the data
# pvrect(fit1276b18, alpha=.95)
# 
# dendro.1276b18 <- as.dendrogram(hc_tr1276b18)
# dendro.plot.1276b18 <- ggdendrogram(data = dendro.1276b18, rotate = TRUE)
# dendro.plot.1276b18 <- dendro.plot.1276b18 + theme(axis.text.y = element_blank(), 
#                                                    axis.text.x.bottom =  element_blank())
# 
# print(dendro.plot.1276b18)
# 
# melt.1276b18 <- melt(mat_tr1276b18)
# head(melt.1276b18)
# colnames(melt.1276b18) <- c("Targeted", "Bycatch", "TI")
# melt.1276b18$TI_Ranges <- cut(melt.1276b18$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# order.1276b18 <- order.dendrogram(dendro.1276b18)
# melt.1276b18$Targeted <- factor(x = melt.1276b18$Targeted,
#                                 levels = melt.1276b18$Targeted[order.1276b18],
#                                 ordered = T)
# 
# heatmap.1276b18 <- ggplot(data=melt.1276b18, 
#                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)")
# print(heatmap.1276b18)
# 
# grid.newpage()
# print(heatmap.1276b18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot.1276b18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))
###Compiling TRO 27.6.b Heatmaps#####
pdf(file = "TRO_276b_Heatmaps.pdf")

grid.newpage()
print(heatmap.1276b, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b09, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b09, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b10, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b10, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b11, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b12, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b12, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b13, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b13, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b14, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b14, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b15, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b15, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b16, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b16, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b17, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b17, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

grid.newpage()
print(heatmap.1276b18, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.1276b18, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.738))

dev.off()
###Heatmaps for Averaging#####
ti_trO276b09.2 <- data.frame() 
ti_trO276b09.2[1:39, 1:39] <- -1

spp_trO276b09.2 <- unique(trO_276b09$spp_id)
spp_trO276b09.2 <- sort(spp_trO276b09.2)

for(i in spp_trO276b09.2){
  sp1 <- subset(trO_276b09, spp_id == i) 
  for(j in spp_trO276b09.2){
    t <- data.frame()
    sp2 <- subset(trO_276b09, spp_id == j) 
    if(i == j){ti_trO276b09.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b09.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b09.2)

row.names(ti_trO276b09.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b09.2) <- sort(unique(trip382$spp))

mat_trO276b09.2 <- as.matrix(ti_trO276b09.2)
dim(mat_trO276b09.2)

##2010
ti_trO276b10.2 <- data.frame()
ti_trO276b10.2[1:39, 1:39] <- -1
# 
# spp_trO276a10.2 <- unique(trO_276a10$spp_id)
# spp_trO276a10.2 <- sort(spp_trO276a10.2)
# 
# for(i in spp_trO276a10.2){
#   sp1 <- subset(trO_276a10, spp_id == i) 
#   for(j in spp_trO276a10.2){
#     t <- data.frame()
#     sp2 <- subset(trO_276a10, spp_id == j) 
#     if(i == j){ti_trO276a10.2[i,j] <- 100} 
#     else{for(k in 1:nrow(sp1)){ 
#       if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
#         t[k,1] <- sp1$wt[k]
#       } else{t[k,1] <- 0}}
#       ti_trO276a10.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
# }                                              
# 
# anyNA(ti_trO276a10.2)
# 
row.names(ti_trO276b10.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b10.2) <- sort(unique(trip382$spp))

mat_trO276b10.2 <- as.matrix(ti_trO276b10.2)
dim(mat_trO276b10.2)

##2011
ti_trO276b11.2 <- data.frame() 
ti_trO276b11.2[1:39, 1:39] <- -1

spp_trO276b11.2 <- unique(trO_276b11$spp_id)
spp_trO276b11.2 <- sort(spp_trO276b11.2)

for(i in spp_trO276b11.2){
  sp1 <- subset(trO_276b11, spp_id == i) 
  for(j in spp_trO276b11.2){
    t <- data.frame()
    sp2 <- subset(trO_276b11, spp_id == j) 
    if(i == j){ti_trO276b11.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b11.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b11.2)

row.names(ti_trO276b11.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b11.2) <- sort(unique(trip382$spp))

mat_trO276b11.2 <- as.matrix(ti_trO276b11.2)
dim(mat_trO276b11.2)

##2012
ti_trO276b12.2 <- data.frame() 
ti_trO276b12.2[1:39, 1:39] <- -1

spp_trO276b12.2 <- unique(trO_276b12$spp_id)
spp_trO276b12.2 <- sort(spp_trO276b12.2)

for(i in spp_trO276b12.2){
  sp1 <- subset(trO_276b12, spp_id == i) 
  for(j in spp_trO276b12.2){
    t <- data.frame()
    sp2 <- subset(trO_276b12, spp_id == j) 
    if(i == j){ti_trO276b12.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b12.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b12.2)

row.names(ti_trO276b12.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b12.2) <- sort(unique(trip382$spp))

mat_trO276b12.2 <- as.matrix(ti_trO276b12.2)
dim(mat_trO276b12.2)

##2013
ti_trO276b13.2 <- data.frame() 
ti_trO276b13.2[1:39, 1:39] <- -1

spp_trO276b13.2 <- unique(trO_276b13$spp_id)
spp_trO276b13.2 <- sort(spp_trO276b13.2)

for(i in spp_trO276b13.2){
  sp1 <- subset(trO_276b13, spp_id == i) 
  for(j in spp_trO276b13.2){
    t <- data.frame()
    sp2 <- subset(trO_276b13, spp_id == j) 
    if(i == j){ti_trO276b13.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b13.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b13.2)

row.names(ti_trO276b13.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b13.2) <- sort(unique(trip382$spp))

mat_trO276b13.2 <- as.matrix(ti_trO276b13.2)
dim(mat_trO276b13.2)

##2014
ti_trO276b14.2 <- data.frame() 
ti_trO276b14.2[1:39, 1:39] <- -1

spp_trO276b14.2 <- unique(trO_276b14$spp_id)
spp_trO276b14.2 <- sort(spp_trO276b14.2)

for(i in spp_trO276b14.2){
  sp1 <- subset(trO_276b14, spp_id == i) 
  for(j in spp_trO276b14.2){
    t <- data.frame()
    sp2 <- subset(trO_276b14, spp_id == j) 
    if(i == j){ti_trO276b14.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b14.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b14.2)

row.names(ti_trO276b14.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b14.2) <- sort(unique(trip382$spp))

mat_trO276b14.2 <- as.matrix(ti_trO276b14.2)
dim(mat_trO276b14.2)

##2015
ti_trO276b15.2 <- data.frame() 
ti_trO276b15.2[1:39, 1:39] <- -1

spp_trO276b15.2 <- unique(trO_276b15$spp_id)
spp_trO276b15.2 <- sort(spp_trO276b15.2)

for(i in spp_trO276b15.2){
  sp1 <- subset(trO_276b15, spp_id == i) 
  for(j in spp_trO276b15.2){
    t <- data.frame()
    sp2 <- subset(trO_276b15, spp_id == j) 
    if(i == j){ti_trO276b15.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b15.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b15.2)

row.names(ti_trO276b15.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b15.2) <- sort(unique(trip382$spp))

mat_trO276b15.2 <- as.matrix(ti_trO276b15.2)
dim(mat_trO276b15.2)

##2016
ti_trO276b16.2 <- data.frame() 
ti_trO276b16.2[1:39, 1:39] <- -1

spp_trO276b16.2 <- unique(trO_276b16$spp_id)
spp_trO276b16.2 <- sort(spp_trO276b16.2)

for(i in spp_trO276b16.2){
  sp1 <- subset(trO_276b16, spp_id == i) 
  for(j in spp_trO276b16.2){
    t <- data.frame()
    sp2 <- subset(trO_276b16, spp_id == j) 
    if(i == j){ti_trO276b16.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b16.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b16.2)

row.names(ti_trO276b16.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b16.2) <- sort(unique(trip382$spp))

mat_trO276b16.2 <- as.matrix(ti_trO276b16.2)
dim(mat_trO276b16.2)

##2017
ti_trO276b17.2 <- data.frame() 
ti_trO276b17.2[1:39, 1:39] <- -1

spp_trO276b17.2 <- unique(trO_276b17$spp_id)
spp_trO276b17.2 <- sort(spp_trO276b17.2)

for(i in spp_trO276b17.2){
  sp1 <- subset(trO_276b17, spp_id == i) 
  for(j in spp_trO276b17.2){
    t <- data.frame()
    sp2 <- subset(trO_276b17, spp_id == j) 
    if(i == j){ti_trO276b17.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b17.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b17.2)

row.names(ti_trO276b17.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b17.2) <- sort(unique(trip382$spp))

mat_trO276b17.2 <- as.matrix(ti_trO276b17.2)
dim(mat_trO276b17.2)

##2018
ti_trO276b18.2 <- data.frame() 
ti_trO276b18.2[1:39, 1:39] <- -1

spp_trO276b18.2 <- unique(trO_276b18$spp_id)
spp_trO276b18.2 <- sort(spp_trO276b18.2)

for(i in spp_trO276b18.2){
  sp1 <- subset(trO_276b18, spp_id == i) 
  for(j in spp_trO276b18.2){
    t <- data.frame()
    sp2 <- subset(trO_276b18, spp_id == j) 
    if(i == j){ti_trO276b18.2[i,j] <- 100} 
    else{for(k in 1:nrow(sp1)){ 
      if((sp1$quarter[k] %in% sp2$quarter)&(sp1$trpCode[k] %in% sp2$trpCode)){
        t[k,1] <- sp1$wt[k]
      } else{t[k,1] <- 0}}
      ti_trO276b18.2[i,j] <- ((sum(t)/sum(sp1$wt))*100)}} 
}                                              

anyNA(ti_trO276b18.2)

row.names(ti_trO276b18.2) <- sort(unique(trip382$spp))
colnames(ti_trO276b18.2) <- sort(unique(trip382$spp))

mat_trO276b18.2 <- as.matrix(ti_trO276b18.2)
dim(mat_trO276b18.2)
##Averaging the Matrices together
matavg_trO276b <- (mat_trO276b09.2 + mat_trO276b11.2 + 
                     mat_trO276b12.2 + mat_trO276b13.2 + mat_trO276b14.2 + 
                     mat_trO276b15.2)/6


dist_avgO276b <- dist(matavg_trO276b, method = "euclidean") #distance matrix
hc_avgO276b <- hclust(dist_avgO276b, method = "ward.D2")
plot(hc_avgO276b)

coph_avgO276b <- cophenetic(hc_avgO276b)
cor(dist_avgO276b, coph_avgO276b) #0.95, greater than 0.75 so that's good

#Using Silhouette method to determine optimal numer of clusters
fviz_nbclust(matavg_trO276b, FUN = hcut, method = "silhouette") #optimal number of clusters is 2


# Ward Hierarchical Clustering
plot(hc_avgO276b) # display dendogram
# draw dendogram with red borders around the 2 clusters 
rect.hclust(hc_avgO276b, k=2, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
fitavgO276b <- pvclust(t(matavg_trO276b), method.hclust="ward.D2",
                       method.dist="euclidean")
plot(fitavgO276b)# dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fitavgO276b, alpha=.95)

dendro.avgO276b <- as.dendrogram(hc_avgO276b, k = 2)
dendro.plot.avgO276b <- ggdendrogram(data = dendro.avgO276b, rotate = TRUE, k = 2)
dendro.plot.avgO276b <- dendro.plot.avgO276b + theme(axis.text.y = element_blank(), 
                                                     axis.text.x.bottom =  element_blank())

print(dendro.plot.avgO276b)

melt.avgO276b <- melt(matavg_trO276b)
head(melt.avgO276b)
colnames(melt.avgO276b) <- c("Targeted", "Bycatch", "TI")
melt.avgO276b$TI_Ranges <- cut(melt.avgO276b$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
order.avgO276b <- order.dendrogram(dendro.avgO276b)
melt.avgO276b$Targeted <- factor(x = melt.avgO276b$Targeted,
                                 levels = melt.avgO276b$Targeted[order.avgO276b],
                                 ordered = T)

heatmap.avgO276b <- ggplot(data=melt.avgO276b, 
                           aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)")
print(heatmap.avgO276b)

grid.newpage()
print(heatmap.avgO276b, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot.avgO276b, vp = viewport(x = 0.89, y = 0.58, width = 0.2, height = 0.739))

####Creating Array trO 27.6.b####
a_trO276b <- abind(mat_trO276b09.2, mat_trO276b10.2, mat_trO276b11.2, 
                   mat_trO276b12.2, mat_trO276b13.2, mat_trO276b14.2,
                   mat_trO276b15.2, mat_trO276b16.2, mat_trO276b17.2, mat_trO276b18.2, along = 3)
View(a_trO276b[,,9])

melta_trO276b <- melt(a_trO276b)
head(melta_trO276b)
colnames(melta_trO276b) <- c("Targeted", "Bycatch", "Year", "TIS")

melta_trO276b$TI_Ranges <- cut(melta_trO276b$TI, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
str(melta_trO276b)


melta_trO276b$Area[1:15210] <- "27.6.b"
melta_trO276b$Area <- as.factor(melta_trO276b$Area)
melta_trO276b$Gear[1:15210] <- "TRO"
melta_trO276b$Gear <- as.factor(melta_trO276b$Gear)

melta_trO276b$Year <- as.factor(melta_trO276b$Year)
melta_trO276b$Year <- revalue(melta_trO276b$Year, c("1"="2009", "2"="2010", "3"="2011", "4"="2012",
                                                    "5"="2013", "6"="2014", "7"="2015", "8"="2016",
                                                    "9"="2017", "10"="2018"))

melta_trO276b <- na.omit(melta_trO276b) #removing absent species that are marked as NA in TI_Ranges

trO276b_by <- with(melta_trO276b, tapply(TIS, INDEX = list(Bycatch, Year), FUN = "mean"))
trO276b_tar <- with(melta_trO276b, tapply(TIS, INDEX = list(Targeted, Year), FUN = "mean"))
trO276b_overall <- with(melta_trO276b, tapply (TIS, INDEX = Year, FUN = "mean"))
trO276b_by <- as.data.frame(trO276b_by)
trO276b_tar <- as.data.frame(trO276b_tar)
trO276b_overall <- as.data.frame(trO276b_overall)

trO276b_overall$var <- with(melta_trO276b, tapply(TIS, INDEX = Year, FUN = "var"))
trO276b_overall$sd <- sqrt(trO276b_overall$var)
trO276b_overall$cv <- 100*((trO276b_overall$sd)/(trO276b_overall$trO276b_overall))

trO276b_by$Overall <- rowMeans(trO276b_by)
trO276b_by$var <- with(melta_trO276b, tapply(TIS, INDEX = list(Bycatch,Year), FUN = "var"))
trO276b_by$sd <- sqrt(trO276b_by$var)
trO276b_by$cv <- 100*((trO276b_by$sd)/(trO276b_by$Overall))
trO276b_by.05 <- trO276b_by[!(trO276b_by$Overall<=5),]
trO276b_by.50 <- trO276b_by[!(trO276b_by$Overall<=50),]
View(trO276b_by.05)
View(trO276b_by.50)

trO276b_tar$Overall <- rowMeans(trO276b_tar)
trO276b_tar$var <- with(melta_trO276b, tapply(TIS, INDEX = list(Targeted,Year), FUN = "var"))
trO276b_tar$sd <- sqrt(trO276b_tar$var)
trO276b_tar$cv <- 100*((trO276b_tar$sd)/(trO276b_tar$Overall))
trO276b_tar.05 <- trO276b_tar[!(trO276b_tar$Overall<=5),]
trO276b_tar.25 <- trO276b_tar[!(trO276b_tar$Overall<=25),]
View(trO276b_tar.05)
View(trO276b_tar.25)


####Centered Moving Average trO 27.6.b--EACH YEAR#####
?rollmean
?select
# 
# test <- trO274_overall %>% select(TIS_year = trO274_overall) %>% 
#   mutate(TIS_mavg3 = rollmean(TIS_year, k = 3, fill = NA),
#          Year = 2009:2018)
# testg<- test %>% gather(metric, value, TIS_year:TIS_mavg3) %>%
#   ggplot(aes(Year, value, color = metric)) +
#   geom_line()
# print(testg)

tO276b <- data.frame()
for(i in 2009:2018){
  selected <- c(i-1, i, i+1)
  f1 <- subset(melta_trO276b, Year %in% selected)
  if((i == 2009)| (i==2018)){
    tO276b[i-2008, 1] <- NA} 
  else{tO276b[i-2008,1] <- mean(f1$TIS)}
}
print(tO276b) #The 2010 Mavg is just the avg of 2009 and 2011

print(trO276b_overall)
trO276b_overall[2,] <- c(NA,0,0,0)
trO276b_overall[8,] <- c(NA,0,0,0)
trO276b_overall[9,] <- c(NA,0,0,0)
trO276b_overall[10,] <- c(NA,0,0,0)

tO276b.2 <- tO276b %>% select(TIS_mavg = V1) %>% mutate(TIS_yavg = trO276b_overall$trO276b_overall,
                                                        Year = 2009:2018)
tO276b.g <- tO276b.2 %>% gather(Metric, TIS, TIS_mavg:TIS_yavg) %>%
  ggplot(aes(Year, TIS, color = Metric)) +
  geom_line()
print(tO276b.g)

########################################################################################
#####################CREATING MASTER DATAFRAME WITH ALL GEAR+YEAR+AREA COMBOS###########
master_melt <- rbind(melta_tr1274, melta_tr2274, melta_tr1276a, melta_tr2276a, melta_tr1276b, melta_trO276b)#Combining all dataframes that were results of SPP38 dataframes
head(master_melt)
area_gear <- tapply(master_melt$TIS, INDEX = list(master_melt$Gear, master_melt$Area), FUN = "mean") #finding average for each combo
area_yearly <- tapply(master_melt$TIS, INDEX = list(master_melt$Area, master_melt$Year), FUN = "mean") #finding average for each area
gear_yearly <- tapply(master_melt$TIS, INDEX = list(master_melt$Gear, master_melt$Year), FUN = "mean") #finding average for each gear

hist(master_melt$TIS, freq = F) #bimodal distrubtion, 0 and 100
boxplot(master_melt$TIS ~ master_melt$Area)
boxplot(master_melt$TIS ~ master_melt$Gear)
boxplot(master_melt$TIS ~ master_melt$Year) #Distribution and Averages different year to year


##Number of Interactions per Area+Gear per Year
inter_count <- tapply(master_melt$TIS, INDEX = list(master_melt$Gear, master_melt$Area, master_melt$Year), FUN = "length") #counting number of interactions for each combo per year
ic2 <- melt(inter_count)
colnames(ic2) <- c("Gear", "Area", "Year", "Interactions")
ic2 <- na.omit(ic2)#removing NA's as these disrupt further calculations and are not relevant
ic2$Gear <- as.factor(ic2$Gear)
ic2$Area <- as.factor(ic2$Area)


ic2$Combos <- paste(ic2$Area, "-", ic2$Gear)
ic2$Combos <- as.factor(ic2$Combos)
str(ic2)
ggplot(ic2)+geom_line(aes(y=Interactions, x=Year, color = Combos))

ic2_summ <- summarySE(data = ic2, measurevar = "Interactions", groupvars = "Combos")#finding average interactions for each combo

mic.g <- ggplot(ic2_summ, aes(x = Combos, y = Interactions)) + geom_col() + 
  geom_errorbar(aes(ymin=Interactions-ci, ymax=Interactions+ci), width = 0.5, size = 1) + 
  scale_y_continuous(limits = c(0,1000), expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = c(0,0)) + 
  theme(text = element_text(size = 13, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##Average TIS per Area+Gear 
TIS_avg <- summarySE(data = master_melt, measurevar = "TIS", groupvars = c("Area", "Gear"))
TIS_avg

TIS_avg$Combos <- paste(TIS_avg$Area, "-", TIS_avg$Gear)
TIS_avg$Combos <- as.factor(TIS_avg$Combos)

mtis.g <- ggplot(TIS_avg, aes(x = Combos, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), width = 0.5, size = 1) + 
  scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  scale_x_discrete(expand = c(0, 0)) + 
  theme(text = element_text(size = 13, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 



pdf(file = "Area and Gear Type Plots.pdf")
dev.off()

##Looking at TIS vs number of interactions
TIS_count <- tapply(master_melt$TIS, INDEX = list(master_melt$Gear, master_melt$Area, master_melt$Year), FUN = "mean")
tc2 <- melt(TIS_count)
tc2

colnames(tc2) <- c("Gear", "Area", "Year", "TIS")
tc2 <- na.omit(tc2)
tc2$Gear <- as.factor(tc2$Gear)
tc2$Area <- as.factor(ic2$Area)


tc2$Combos <- paste(tc2$Area, "-", tc2$Gear)
tc2$Combos <- as.factor(tc2$Combos)

nrow(tc2)
tic <- data.frame()
tic[1:55,1] <- 0
tic$TIS <- tc2$TIS
tic$inter <- ic2$Interactions
tic$combos  <- ic2$Combos

tic <- tic[,-1]

ggplot(tic, (aes(x= inter, y = TIS, color = combos))) + 
  geom_point() + geom_abline(intercept = 87.767565, slope = -0.038118, color = "red") + 
  labs(x = "Number of Interactions per Year", y = "Yearly Average TIS (%)") +
  guides(color=guide_legend(title="Area-Gear"))
pdf(file = "TIS vs Interactions.pdf")
dev.off()

resiplot(tic.lm)
summary(tic.lm)

# Call:
#   lm(formula = TIS ~ inter, data = tic)
# 
# Residuals:
#  Min      1Q      Median  3Q     Max 
# -17.316  -5.394   1.007   4.877  20.809 
# 
# Coefficients:
#                Estimate   Std. Error t value   Pr(>|t|)    
# (Intercept)    87.767565  2.814081   31.189  < 2e-16 ***
#   inter       -0.038118   0.004884  -7.804     2.3e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 8.905 on 53 degrees of freedom
# Multiple R-squared:  0.5347,	Adjusted R-squared:  0.5259 
# F-statistic:  60.9 on 1 and 53 DF,  p-value: 2.303e-10

anova(tic.lm)
# Analysis of Variance Table
# 
# Response: TIS
#            Df Sum Sq  Mean Sq  F value    Pr(>F)    
# inter      1  4829.4  4829.4   60.902     2.303e-10 ***
# Residuals 53  4202.8    79.3                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##Number of Trips per Area+Gear per Year
trips6_year <- tapply(triparea6$trpCode, 
                      list(triparea6$year, triparea6$foCatNat), length)#Area 27.6.a
trips6_year <- trips6_year[,3:4]
t6y_melt <- melt(trips6_year)
t6y_melt$Area <- "27.6.a"
t6y_melt$combos <- paste(t6y_melt$Area, "-", t6y_melt$Var2)
colnames(t6y_melt) <- c("Year", "Gear", "Trips", "Area", "Combos")

trips2_year <- tapply(triparea2$trpCode, 
                      list(triparea2$year, triparea2$foCatNat), length)#Area 27.4
trips2_year <- trips2_year[,2:3]
t2y_melt <- melt(trips2_year)
t2y_melt$Area <- "27.4"
t2y_melt$combos <- paste(t2y_melt$Area, "-", t2y_melt$Var2)
colnames(t2y_melt) <- c("Year", "Gear", "Trips", "Area", "Combos")

trips7_year <- tapply(triparea7$trpCode, 
                      list(triparea7$year, triparea7$foCatNat), length)#Area 27.6.b 
trips7_year#No trips in 2010
t7y_melt <- melt(trips7_year)
t7y_melt$Var2 <- revalue(t7y_melt$Var2, c("TRother" = "TRO", "TR1" = "TR1"))
t7y_melt$Area <- "27.6.b"
t7y_melt$combos <- paste(t7y_melt$Area, "-", t7y_melt$Var2)
colnames(t7y_melt) <- c("Year", "Gear", "Trips", "Area", "Combos")

alltrps <- rbind(t2y_melt, t6y_melt, t7y_melt)
alltrps
nrow(alltrps)
alltrps <- na.omit(alltrps)
nrow(alltrps)

merg1 <- merge(tc2, alltrps, by = c("Year", "Combos"))
merg1 <- merg1[,-c(6,8)]
colnames(merg1) <- c("Year", "Combos", "Gear", "Area", "TIS", "Trips")

merg2 <- merge(merg1, ic2, by = c("Year", "Combos"))
merg2 <- merg2[,-c(7:8)]
colnames(merg2) <- c("Year", "Combos", "Gear", "Area", "TIS", "Trips", "Interactions")
merg2$fYear <- as.factor(merg2$Year)


stis.g <- ggplot(merg2, (aes(x= log(Trips), y = TIS, color = Combos))) + 
  geom_point() + 
  labs(x = "Natural Log of Trips per Year", y = "Yearly Average TIS (%)") +
  guides(color=guide_legend(title="Area-Gear")) + 
  theme(legend.position = "none") + scale_color_brewer(palette = "Dark2") + 
  theme(text = element_text(size = 13, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


sint.g <- ggplot(merg2, (aes(x= log(Trips), y = Interactions, color = Combos))) + 
  geom_point() + labs(x = "Natural Log of Trips per Year", y = "Interactions per Year") +
  guides(color=guide_legend(title="Area-Gear"))+ theme(legend.position = "right") + 
  scale_color_brewer(palette = "Dark2") + 
  theme(text = element_text(size = 13, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(merg2, aes(x = Combos, y = Trips, fill = fYear)) + geom_col(color = "black") + 
  labs(y="Total Trips", x = "Area-Gear", fill = "Year") + 
  scale_fill_brewer(palette="RdGy")

ggplot(merg2, aes(x = Year, y = Trips, fill = Combos)) +
  geom_col(color = "black") + labs(y = "Total Trips", x = "Year", fill = "Area-Gear") + 
  scale_fill_brewer(palette="RdGy")

pdf(file = "Trip Plots per AG.pdf")
dev.off()

blankplot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
slegend <- get_legend(sint.g)
sint.g2 <- sint.g + theme(legend.position = "none")
grid.arrange((stis.g+labs(tag = "A")), (sint.g2+labs(tag = "B")), slegend, nrow = 1, widths = c(2.1, 2.0, 0.8))

#####################LOOKING AT SPECIFIC SPECIES INTERACTIONS###########################
table(master_melt$Targeted)

##These are the species that either bring in tens of millions of pounds 
###or are marked as Critically Endangered TARGETED & BYCATCH
#10: Dipturus batis
#11: Gadus morhua
#21: Melanogrammus aeglefinnus
#28: Nephrops norvegicus

#34: Scomber scombrus


mm2 <- master_melt
mm2$sppcode <- as.numeric(mm2$Targeted)
mm2$sppcodeb <- as.numeric(mm2$Bycatch)
mm2_guide <- unique(mm2$sppcode)
mm2_guide <- as.data.frame(mm2_guide)
mm2_guide$Targeted <- unique(mm2$Targeted)
View(mm2_guide)

mm3 <- subset(mm2, sppcode %in% c(10, 11, 21, 28, 34) | sppcodeb %in% c(10, 11, 21, 28, 34))
table(mm3$Bycatch)

head(mm3)
nrow(mm3)

mm3.1274 <- subset(mm3, Area %in% "27.4" & Gear %in% "TR1")
mm3.2274 <- subset(mm3, Area %in% "27.4" & Gear %in% "TR2")
mm3.1276a <- subset(mm3, Area %in% "27.6.a" & Gear %in% "TR1")
mm3.2276a <- subset(mm3, Area %in% "27.6.a" & Gear %in% "TR2")
mm3.1276b <- subset(mm3, Area %in% "27.6.b" & Gear %in% "TR1")
mm3.O276b <- subset(mm3, Area %in% "27.6.b" & Gear %in% "TRO")

####1274 ssp6####

unique(mm3.1274$Area)
unique(mm3.1274$Gear)
e
heatmap.mm3127409 <- ggplot(data=subset(mm3.1274, Year == 2009), 
                          aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127409)

heatmap.mm3127410 <- ggplot(data=subset(mm3.1274, Year == 2010), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127410)

heatmap.mm3127411 <- ggplot(data=subset(mm3.1274, Year == 2011), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127411)

heatmap.mm3127412 <- ggplot(data=subset(mm3.1274, Year == 2012), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127412)

heatmap.mm3127413 <- ggplot(data=subset(mm3.1274, Year == 2013), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127413)

heatmap.mm3127414 <- ggplot(data=subset(mm3.1274, Year == 2014), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127414)

heatmap.mm3127415 <- ggplot(data=subset(mm3.1274, Year == 2015), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127415)

heatmap.mm3127416 <- ggplot(data=subset(mm3.1274, Year == 2016), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127416)

heatmap.mm3127417 <- ggplot(data=subset(mm3.1274, Year == 2017), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127417)

heatmap.mm3127418 <- ggplot(data=subset(mm3.1274, Year == 2018), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3127418)

####2274 ssp6####
unique(mm3.2274$Area)
unique(mm3.2274$Gear)

heatmap.mm3227409 <- ggplot(data=subset(mm3.2274, Year == 2009), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227409)

heatmap.mm3227410 <- ggplot(data=subset(mm3.2274, Year == 2010), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227410)

heatmap.mm3227411 <- ggplot(data=subset(mm3.2274, Year == 2011), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227411)

heatmap.mm3227412 <- ggplot(data=subset(mm3.2274, Year == 2012), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227412)

heatmap.mm3227413 <- ggplot(data=subset(mm3.2274, Year == 2013), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227413)

heatmap.mm3227414 <- ggplot(data=subset(mm3.2274, Year == 2014), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227414)

heatmap.mm3227415 <- ggplot(data=subset(mm3.2274, Year == 2015), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227415)

heatmap.mm3227416 <- ggplot(data=subset(mm3.2274, Year == 2016), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227416)

heatmap.mm3227417 <- ggplot(data=subset(mm3.2274, Year == 2017), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227417)

heatmap.mm3227418 <- ggplot(data=subset(mm3.2274, Year == 2018), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3227418)


####1276a ssp6####

unique(mm3.1276a$Area)
unique(mm3.1276a$Gear)

heatmap.mm31276a09 <- ggplot(data=subset(mm3.1276a, Year == 2009), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a09)

heatmap.mm31276a10 <- ggplot(data=subset(mm3.1276a, Year == 2010), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a10)

heatmap.mm31276a11 <- ggplot(data=subset(mm3.1276a, Year == 2011), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a11)

heatmap.mm31276a12 <- ggplot(data=subset(mm3.1276a, Year == 2012), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a12)

heatmap.mm31276a13 <- ggplot(data=subset(mm3.1276a, Year == 2013), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a13)

heatmap.mm31276a14 <- ggplot(data=subset(mm3.1276a, Year == 2014), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a14)

heatmap.mm31276a15 <- ggplot(data=subset(mm3.1276a, Year == 2015), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a15)

heatmap.mm31276a16 <- ggplot(data=subset(mm3.1276a, Year == 2016), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a16)

heatmap.mm31276a17 <- ggplot(data=subset(mm3.1276a, Year == 2017), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a17)

heatmap.mm31276a18 <- ggplot(data=subset(mm3.1276a, Year == 2018), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276a18)

####2276a spp6####

unique(mm3.2276a$Area)
unique(mm3.2276a$Gear)

heatmap.mm32276a09 <- ggplot(data=subset(mm3.2276a, Year == 2009), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a09)

heatmap.mm32276a10 <- ggplot(data=subset(mm3.2276a, Year == 2010), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a10)

heatmap.mm32276a11 <- ggplot(data=subset(mm3.2276a, Year == 2011), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a11)

heatmap.mm32276a12 <- ggplot(data=subset(mm3.2276a, Year == 2012), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a12)

heatmap.mm32276a13 <- ggplot(data=subset(mm3.2276a, Year == 2013), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a13)

heatmap.mm32276a14 <- ggplot(data=subset(mm3.2276a, Year == 2014), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a14)

heatmap.mm32276a15 <- ggplot(data=subset(mm3.2276a, Year == 2015), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a15)

heatmap.mm32276a16 <- ggplot(data=subset(mm3.2276a, Year == 2016), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a16)

heatmap.mm32276a17 <- ggplot(data=subset(mm3.2276a, Year == 2017), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a17)

heatmap.mm32276a18 <- ggplot(data=subset(mm3.2276a, Year == 2018), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm32276a18)
####1276b spp6####

unique(mm3.1276b$Area)
unique(mm3.1276b$Gear)

heatmap.mm31276b09 <- ggplot(data=subset(mm3.1276b, Year == 2009), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b09)

heatmap.mm31276b11 <- ggplot(data=subset(mm3.1276b, Year == 2011), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b11)

heatmap.mm31276b12 <- ggplot(data=subset(mm3.1276b, Year == 2012), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b12)

heatmap.mm31276b13 <- ggplot(data=subset(mm3.1276b, Year == 2013), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b13)

heatmap.mm31276b14 <- ggplot(data=subset(mm3.1276b, Year == 2014), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b14)

heatmap.mm31276b15 <- ggplot(data=subset(mm3.1276b, Year == 2015), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b15)

heatmap.mm31276b16 <- ggplot(data=subset(mm3.1276b, Year == 2016), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b16)

heatmap.mm31276b17 <- ggplot(data=subset(mm3.1276b, Year == 2017), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b17)

heatmap.mm31276b18 <- ggplot(data=subset(mm3.1276b, Year == 2018), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm31276b18)
####O276b spp6####

unique(mm3.O276b$Area)
unique(mm3.O276b$Gear)

heatmap.mm3O276b09 <- ggplot(data=subset(mm3.O276b, Year == 2009), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3O276b09)

heatmap.mm3O276b11 <- ggplot(data=subset(mm3.O276b, Year == 2011), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3O276b11)

heatmap.mm3O276b12 <- ggplot(data=subset(mm3.O276b, Year == 2012), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3O276b12)

heatmap.mm3O276b13 <- ggplot(data=subset(mm3.O276b, Year == 2013), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3O276b13)

heatmap.mm3O276b14 <- ggplot(data=subset(mm3.O276b, Year == 2014), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3O276b14)

heatmap.mm3O276b15 <- ggplot(data=subset(mm3.O276b, Year == 2015), 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.mm3O276b15)
####Compiling heatmaps into 6 pdfs####
##1274
pdf(file = "Spp6 1274.pdf")
print(heatmap.mm3127409)
print(heatmap.mm3127410)
print(heatmap.mm3127411)
print(heatmap.mm3127412)
print(heatmap.mm3127413)
print(heatmap.mm3127414)
print(heatmap.mm3127415)
print(heatmap.mm3127416)
print(heatmap.mm3127417)
print(heatmap.mm3127418)
dev.off()

##2274
pdf(file = "Spp6 2274.pdf")
print(heatmap.mm3227409)
print(heatmap.mm3227410)
print(heatmap.mm3227411)
print(heatmap.mm3227412)
print(heatmap.mm3227413)
print(heatmap.mm3227414)
print(heatmap.mm3227415)
print(heatmap.mm3227416)
print(heatmap.mm3227417)
print(heatmap.mm3227418)
dev.off()

##1276a
pdf(file = "Spp6 1276a.pdf")
print(heatmap.mm31276a09)
print(heatmap.mm31276a10)
print(heatmap.mm31276a11)
print(heatmap.mm31276a12)
print(heatmap.mm31276a13)
print(heatmap.mm31276a14)
print(heatmap.mm31276a15)
print(heatmap.mm31276a16)
print(heatmap.mm31276a17)
print(heatmap.mm31276a18)
dev.off()

##2276a
pdf(file = "Spp6 2276a.pdf")
print(heatmap.mm32276a09)
print(heatmap.mm32276a10)
print(heatmap.mm32276a11)
print(heatmap.mm32276a12)
print(heatmap.mm32276a13)
print(heatmap.mm32276a14)
print(heatmap.mm32276a15)
print(heatmap.mm32276a16)
print(heatmap.mm32276a17)
print(heatmap.mm32276a18)
dev.off()

##1276b
pdf(file = "Spp6 1276b.pdf")
print(heatmap.mm31276b09)
#No 2010
print(heatmap.mm31276b11)
print(heatmap.mm31276b12)
print(heatmap.mm31276b13)
print(heatmap.mm31276b14)
print(heatmap.mm31276b15)
print(heatmap.mm31276b16)
print(heatmap.mm31276b17)
print(heatmap.mm31276b18)
dev.off()

##O276b
pdf(file = "Spp6 O276b.pdf")
print(heatmap.mm3O276b09)
#No 2010
print(heatmap.mm3O276b11)
print(heatmap.mm3O276b12)
print(heatmap.mm3O276b13)
print(heatmap.mm3O276b14)
print(heatmap.mm3O276b15)
#No 2016, 2017, 2018
dev.off()
####Averages per Spp6#####
#1274
avg1274yb <- summarySE(data = mm3.1274, measurevar = "TIS", groupvars = c("Year","Bycatch"))
avg1274yb$sppcode <- as.numeric(as.factor(avg1274yb$Bycatch))
avg1274yb6 <- subset(avg1274yb, sppcode %in% c(10,11,21,28,34))

avg1274b <- summarySE(data = mm3.1274, measurevar = "TIS", groupvars = "Bycatch")
avg1274b$sppcode <- as.numeric(as.factor(avg1274b$Bycatch))
avg1274b6 <- subset(avg1274b, sppcode %in% c(10,11,21,28,34))

avg1274yt <- summarySE(data = mm3.1274, measurevar = "TIS", groupvars = c("Year","Targeted"))
avg1274yt$sppcode <- as.numeric(as.factor(avg1274yt$Targeted))
avg1274yt6 <- subset(avg1274yt, sppcode %in% c(10,11,21,28,34))

avg1274t <- summarySE(data = mm3.1274, measurevar = "TIS", groupvars = "Targeted")
avg1274t$sppcode <- as.numeric(as.factor(avg1274t$Targeted))
avg1274t6 <- subset(avg1274t, sppcode %in% c(10,11,21,28,34))

#2274
avg2274yb <- summarySE(data = mm3.2274, measurevar = "TIS", groupvars = c("Year","Bycatch"))
avg2274yb$sppcode <- as.numeric(as.factor(avg2274yb$Bycatch))
avg2274yb6 <- subset(avg2274yb, sppcode %in% c(10,11,21,28,34))

avg2274b <- summarySE(data = mm3.2274, measurevar = "TIS", groupvars = "Bycatch")
avg2274b$sppcode <- as.numeric(as.factor(avg2274b$Bycatch))
avg2274b6 <- subset(avg2274b, sppcode %in% c(10,11,21,28,34))

avg2274yt <- summarySE(data = mm3.2274, measurevar = "TIS", groupvars = c("Year","Targeted"))
avg2274yt$sppcode <- as.numeric(as.factor(avg2274yt$Targeted))
avg2274yt6 <- subset(avg2274yt, sppcode %in% c(10,11,21,28,34))
unique(avg2274yt6$Targeted)

avg2274t <- summarySE(data = mm3.2274, measurevar = "TIS", groupvars = "Targeted")
avg2274t$sppcode <- as.numeric(as.factor(avg2274t$Targeted))
avg2274t6 <- subset(avg2274t, sppcode %in% c(10,11,21,28,34))
unique(avg2274t6$Targeted)

#1276a
avg1276ayb <- summarySE(data = mm3.1276a, measurevar = "TIS", groupvars = c("Year","Bycatch"))
avg1276ayb$sppcode <- as.numeric(as.factor(avg1276ayb$Bycatch))
avg1276ayb6 <- subset(avg1276ayb, sppcode %in% c(10,11,21,28,34))

avg1276ab <- summarySE(data = mm3.1276a, measurevar = "TIS", groupvars = "Bycatch")
avg1276ab$sppcode <- as.numeric(as.factor(avg1276ab$Bycatch))
avg1276ab6 <- subset(avg1276ab, sppcode %in% c(10,11,21,28,34))

avg1276ayt <- summarySE(data = mm3.1276a, measurevar = "TIS", groupvars = c("Year","Targeted"))
avg1276ayt$sppcode <- as.numeric(as.factor(avg1276ayt$Targeted))
avg1276ayt6 <- subset(avg1276ayt, sppcode %in% c(10,11,21,28,34))

avg1276at <- summarySE(data = mm3.1276a, measurevar = "TIS", groupvars = "Targeted")
avg1276at$sppcode <- as.numeric(as.factor(avg1276at$Targeted))
avg1276at6 <- subset(avg1276at, sppcode %in% c(10,11,21,28,34))

#2276a
avg2276ayb <- summarySE(data = mm3.2276a, measurevar = "TIS", groupvars = c("Year","Bycatch"))
avg2276ayb$sppcode <- as.numeric(as.factor(avg2276ayb$Bycatch))
avg2276ayb6 <- subset(avg2276ayb, sppcode %in% c(10,11,21,28,34))

avg2276ab <- summarySE(data = mm3.2276a, measurevar = "TIS", groupvars = "Bycatch")
avg2276ab$sppcode <- as.numeric(as.factor(avg2276ab$Bycatch))
avg2276ab6 <- subset(avg2276ab, sppcode %in% c(10,11,21,28,34))

avg2276ayt <- summarySE(data = mm3.2276a, measurevar = "TIS", groupvars = c("Year","Targeted"))
avg2276ayt$sppcode <- as.numeric(as.factor(avg2276ayt$Targeted))
avg2276ayt6 <- subset(avg2276ayt, sppcode %in% c(10,11,21,28,34))

avg2276at <- summarySE(data = mm3.2276a, measurevar = "TIS", groupvars = "Targeted")
avg2276at$sppcode <- as.numeric(as.factor(avg2276at$Targeted))
avg2276at6 <- subset(avg2276at, sppcode %in% c(10,11,21,28,34))

#1276b
avg1276byb <- summarySE(data = mm3.1276b, measurevar = "TIS", groupvars = c("Year","Bycatch"))
avg1276byb$sppcode <- as.numeric(as.factor(avg1276byb$Bycatch))
avg1276byb6 <- subset(avg1276byb, sppcode %in% c(10,11,21,28,34))

avg1276bb <- summarySE(data = mm3.1276b, measurevar = "TIS", groupvars = "Bycatch")
avg1276bb$sppcode <- as.numeric(as.factor(avg1276bb$Bycatch))
avg1276bb6 <- subset(avg1276bb, sppcode %in% c(10,11,21,28,34))

avg1276byt <- summarySE(data = mm3.1276b, measurevar = "TIS", groupvars = c("Year","Targeted"))
avg1276byt$sppcode <- as.numeric(as.factor(avg1276byt$Targeted))
avg1276byt6 <- subset(avg1276byt, sppcode %in% c(10,11,21,28,34))

avg1276bt <- summarySE(data = mm3.1276b, measurevar = "TIS", groupvars = "Targeted")
avg1276bt$sppcode <- as.numeric(as.factor(avg1276bt$Targeted))
avg1276bt6 <- subset(avg1276bt, sppcode %in% c(10,11,21,28,34))

#O276b
avgO276byb <- summarySE(data = mm3.O276b, measurevar = "TIS", groupvars = c("Year","Bycatch"))
avgO276byb$sppcode <- as.numeric(as.factor(avgO276byb$Bycatch))
avgO276byb6 <- subset(avgO276byb, sppcode %in% c(10,11,21,28,34))

avgO276bb <- summarySE(data = mm3.O276b, measurevar = "TIS", groupvars = "Bycatch")
avgO276bb$sppcode <- as.numeric(as.factor(avgO276bb$Bycatch))
avgO276bb6 <- subset(avgO276bb, sppcode %in% c(10,11,21,28,34))

avgO276byt <- summarySE(data = mm3.O276b, measurevar = "TIS", groupvars = c("Year","Targeted"))
avgO276byt$sppcode <- as.numeric(as.factor(avgO276byt$Targeted))
avgO276byt6 <- subset(avgO276byt, sppcode %in% c(10,11,21,28,34))

avgO276bt <- summarySE(data = mm3.O276b, measurevar = "TIS", groupvars = "Targeted")
avgO276bt$sppcode <- as.numeric(as.factor(avgO276bt$Targeted))
avgO276bt6 <- subset(avgO276bt, sppcode %in% c(10,11,21,28,34))


####Graphing Averages#####
#List of Average Dataframes
##1274
avg1274b6 #10yr avg for each species in 27.4, TR1 BYCATCH
avg1274yb6 #Yearly avg for each species in 27.4, TR1 BYCATCH
avg1274t6  #10yr avg for each species in 27.4, TR1 TARGETED
avg1274yt6 #Yearly avg for each species in 27.4, TR1 TARGETED 


avg1274yb6$nYear <- as.numeric(levels(avg1274yb6$Year))[avg1274yb6$Year]
ggplot(avg1274yb6, aes(x = nYear, y = TIS, color = Bycatch))+geom_line() + ylim(0,100)
str(avg1274yb6)

avg1274yt6$nYear <- as.numeric(levels(avg1274yt6$Year))[avg1274yt6$Year]
ggplot(avg1274yt6, aes(x = nYear, y = TIS, color = Targeted))+geom_line()+ ylim(0,100)
str(avg1274yb6)

ggplot(avg1274b6, aes(x = Bycatch, y = TIS)) + geom_col() + 
    geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


ggplot(avg1274t6, aes(x = Targeted, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


pdf(file = "Spp 6 1274 Averages.pdf")
dev.off()

##2274
avg2274b6 
avg2274yb6
avg2274t6
avg2274yt6

avg2274yb6$nYear <- as.numeric(levels(avg2274yb6$Year))[avg2274yb6$Year]
ggplot(avg2274yb6, aes(x = nYear, y = TIS, color = Bycatch))+geom_line()+ ylim(0,100)

avg2274yt6$nYear <- as.numeric(levels(avg2274yt6$Year))[avg2274yt6$Year]
ggplot(avg2274yt6, aes(x = nYear, y = TIS, color = Targeted))+geom_line()+ ylim(0,100)

ggplot(avg2274b6, aes(x = Bycatch, y = TIS)) + geom_col() + ylim(0,100) + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci))  + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


ggplot(avg2274t6, aes(x = Targeted, y = TIS)) + geom_col() + ylim(0,100) + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci))  + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


pdf(file = "Spp 6 2274 Averages.pdf")
dev.off()
##1276a
avg1276ab6
avg1276ayb6
avg1276at6
avg1276ayt6

avg1276ayb6$nYear <- as.numeric(levels(avg1276ayb6$Year))[avg1276ayb6$Year]
ggplot(avg1276ayb6, aes(x = nYear, y = TIS, color = Bycatch))+geom_line() + ylim(0,100)

avg1276ayt6$nYear <- as.numeric(levels(avg1276ayt6$Year))[avg1276ayt6$Year]
ggplot(avg1276ayt6, aes(x = nYear, y = TIS, color = Targeted))+geom_line() + ylim(0,100)

ggplot(avg1276ab6, aes(x = Bycatch, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


ggplot(avg1276at6, aes(x = Targeted, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


pdf(file = "Spp 6 1276a Averages.pdf")
dev.off()

##2276b
avg2276ab6
avg2276ayb6
avg2276at6
avg2276ayt6

avg2276ayb6$nYear <- as.numeric(levels(avg2276ayb6$Year))[avg2276ayb6$Year]
ggplot(avg2276ayb6, aes(x = nYear, y = TIS, color = Bycatch))+geom_line() + ylim(0,100)

avg2276ayt6$nYear <- as.numeric(levels(avg2276ayt6$Year))[avg2276ayt6$Year]
ggplot(avg2276ayt6, aes(x = nYear, y = TIS, color = Targeted))+geom_line() + ylim(0,100)

ggplot(avg2276ab6, aes(x = Bycatch, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


ggplot(avg2276at6, aes(x = Targeted, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))



pdf(file = "Spp 6 2276a Averages.pdf")
dev.off()
##1276b
avg1276bb6
avg1276byb6
avg1276bt6
avg1276byt6

avg1276byb6$nYear <- as.numeric(levels(avg1276byb6$Year))[avg1276byb6$Year]
ggplot(avg1276byb6, aes(x = nYear, y = TIS, color = Bycatch))+geom_line() + ylim(0,100)

avg1276byt6$nYear <- as.numeric(levels(avg1276byt6$Year))[avg1276byt6$Year]
ggplot(avg1276byt6, aes(x = nYear, y = TIS, color = Targeted))+geom_line() + ylim(0,100)

ggplot(avg1276bb6, aes(x = Bycatch, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


ggplot(avg1276bt6, aes(x = Targeted, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,102) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


pdf(file = "Spp 6 1276b Averages.pdf")
dev.off()

##O276b
avgO276bb6
avgO276byb6
avgO276bt6
avgO276byt6

avgO276byb6$nYear <- as.numeric(levels(avgO276byb6$Year))[avgO276byb6$Year]
ggplot(avgO276byb6, aes(x = nYear, y = TIS, color = Bycatch))+geom_line() + ylim(0,100)

avgO276byt6$nYear <- as.numeric(levels(avgO276byt6$Year))[avgO276byt6$Year]
ggplot(avgO276byt6, aes(x = nYear, y = TIS, color = Targeted))+geom_line() + ylim(0,100)

ggplot(avgO276bb6, aes(x = Bycatch, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,100) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


ggplot(avgO276bt6, aes(x = Targeted, y = TIS)) + geom_col() + 
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci)) + ylim(0,104) + scale_x_discrete(labels=c("Dipturus batis" = "D. batis", "Gadus morhua" = "G. morhua", "Melanogrammus aeglefinus" = "M. aeglefinus", "Nephrops norvegicus" = "N. norvegicus", "Scomber scombrus" = "S. scombrus"))


pdf(file = "Spp 6 O276b Averages.pdf")
dev.off()

####27.4 v 27.6.a########
##Bycatch
mm3.274 <- rbind(mm3.1274, mm3.2274)
head(mm3.274)

mm3.276a <- rbind(mm3.1276a, mm3.2276a)
head(mm3.276a)

mm3.274by <- subset(mm3.274, sppcodeb %in% c(10,11,21,28,34)) #Taking Rostroraja alba out as it rarely appears
unique(mm3.274by$Bycatch)
mm3.276aby <- subset(mm3.276a, sppcodeb %in% c(10,11,21,28,34)) #Taking rostroraja alba out as it rarely appears
unique(mm3.276aby$Bycatch)

mm3.4_6aby <- rbind(mm3.274by, mm3.276aby)
head(mm3.4_6aby)
avg.4_6aby <- summarySE(data = mm3.4_6aby, measurevar = "TIS", groupvars = c("Area", "Bycatch"))

ggplot(avg.4_6aby, aes(y = TIS, x = Bycatch, fill = Area)) + 
  geom_col(position = position_dodge()) + geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge())

##Targeted
mm3.274ta <- subset(mm3.274, sppcode %in% c(10,11,21,28,34)) #Taking R. alba out
unique(mm3.274ta$Targeted)
mm3.276ata <- subset(mm3.276a, sppcode %in% c(10,11,21,28,34)) #Taking R. alba out
unique(mm3.276ata$Targeted)

mm3.4_6ata <- rbind(mm3.274ta, mm3.276ata)
head(mm3.4_6ata)
avg.4_6ata <- summarySE(data = mm3.4_6ata, measurevar = "TIS", groupvars = c("Area", "Targeted"))

ggplot(avg.4_6ata, aes(y = TIS, x = Targeted, fill = Area)) + 
  geom_col(position = position_dodge()) + geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge())


####TR1 v TR2#######
head(mm3.4_6aby)
avg.tr1_2by <- summarySE(data = mm3.4_6aby, measurevar = "TIS", groupvars = c("Gear", "Bycatch"))

ggplot(avg.tr1_2by, aes(y = TIS, x = Bycatch, fill = Gear)) + 
  geom_col(position = position_dodge()) + geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge())

##Targeted
head(mm3.4_6ata)
avg.tr1_2ta <- summarySE(data = mm3.4_6ata, measurevar = "TIS", groupvars = c("Gear", "Targeted"))

ggplot(avg.tr1_2ta, aes(y = TIS, x = Targeted, fill = Gear)) + 
  geom_col(position = position_dodge()) + geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge())

####1274 v 2274 v 1276a v 2276a##########
mm3.4_6aby$Combos <- paste(mm3.4_6aby$Area, "-", mm3.4_6aby$Gear)
head(mm3.4_6aby)
avg.4_6abyco <- summarySE(data = mm3.4_6aby, measurevar = "TIS", groupvars = c("Combos", "Bycatch"))

ggplot(avg.4_6abyco, aes(y = TIS, x = Bycatch, fill = Combos)) + 
  geom_col(position = position_dodge()) + geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge())

##Targeted
mm3.4_6ata$Combos <- paste(mm3.4_6ata$Area, "-", mm3.4_6ata$Gear)
head(mm3.4_6ata)
avg.4_6ataco <- summarySE(data = mm3.4_6ata, measurevar = "TIS", groupvars = c("Combos", "Targeted"))

ggplot(avg.4_6ataco, aes(y = TIS, x = Targeted, fill = Combos)) + 
  geom_col(position = position_dodge()) + geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge())

####Compiling Comparative Graphs####
pdf("Spp 6 Comparative Bar Graphs.pdf")

qby.g <- ggplot(avg.4_6aby, aes(y = TIS, x = Bycatch, fill = Area)) + 
  geom_col(position = "dodge") + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge(0.8), width = 0.5, size = 1)+ 
  scale_x_discrete(labels=c("Dipturus batis" = "*D. batis*", 
                            "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
                            "Nephrops norvegicus" = "*N. norvegicus*", 
                            "Scomber scombrus" = "*S. scombrus*"), expand = c(0, 0)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(text = element_text(size = 13, face = "bold")) + scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


qtar.g <- ggplot(avg.4_6ata, aes(y = TIS, x = Targeted, fill = Area)) + 
  geom_col(position = "dodge") + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge(0.8), width = 0.5, size = 1)+ 
  scale_x_discrete(labels=c("Dipturus batis" = "*D. batis*", 
                            "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
                            "Nephrops norvegicus" = "*N. norvegicus*", 
                            "Scomber scombrus" = "*S. scombrus*"), expand = c(0, 0)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(text = element_text(size = 13, face = "bold")) + scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
##For Report
qlegend <- get_legend(qby.g)
qby.g2 <- qby.g + theme(legend.position = "none")
qtar.g <- qtar.g + theme(legend.position = "none")
grid.arrange((qtar.g+labs(tag="A")+ylab("TIS (%)")), blankplot, (qby.g2+labs(tag="B")+ylab("TIS (%)")), qlegend, ncol = 2, widths = c(2,0.5))

#
kby.g <- ggplot(avg.tr1_2by, aes(y = TIS, x = Bycatch, fill = Gear)) + 
  geom_col(position = "dodge") + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge(0.8), width = 0.5, size = 1)+ 
  scale_x_discrete(labels=c("Dipturus batis" = "*D. batis*", 
                            "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
                            "Nephrops norvegicus" = "*N. norvegicus*", 
                            "Scomber scombrus" = "*S. scombrus*"), expand = c(0, 0)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(text = element_text(size = 13, face = "bold")) + scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ktar.g <- ggplot(avg.tr1_2ta, aes(y = TIS, x = Targeted, fill = Gear)) + 
  geom_col(position = "dodge") + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge(0.8), width = 0.5, size = 1)+ 
  scale_x_discrete(labels=c("Dipturus batis" = "*D. batis*", 
                            "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
                            "Nephrops norvegicus" = "*N. norvegicus*", 
                            "Scomber scombrus" = "*S. scombrus*"), expand = c(0, 0)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(text = element_text(size = 13, face = "bold")) + scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##For Report
klegend <- get_legend(kby.g)
kby.g2 <- kby.g + theme(legend.position = "none")
ktar.g <- ktar.g + theme(legend.position = "none")
grid.arrange((ktar.g+labs(tag="A")+ylab("TIS (%)")), blankplot, (kby.g2+labs(tag="B")+ylab("TIS (%)")), klegend, ncol = 2, widths = c(2,0.5))


#
dby.g <- ggplot(avg.4_6abyco, aes(y = TIS, x = Bycatch, fill = Combos)) + 
  geom_col(position = position_dodge()) + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge(0.9), width = 0.5, size = 1)+ 
  scale_x_discrete(labels=c("Dipturus batis" = "*D. batis*", 
                            "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
                            "Nephrops norvegicus" = "*N. norvegicus*", 
                            "Scomber scombrus" = "*S. scombrus*"), expand = c(0, 0)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(text = element_text(size = 13, face = "bold")) + scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


dtar.g <- ggplot(avg.4_6ataco, aes(y = TIS, x = Targeted, fill = Combos)) + 
  geom_col(position = position_dodge()) + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, .1))) + ylab("TIS (%)") +
  geom_errorbar(aes(ymin=TIS-ci, ymax=TIS+ci), position = position_dodge(0.9), width = 0.5, size = 1)+ 
  scale_x_discrete(labels=c("Dipturus batis" = "*D. batis*", 
                            "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
                            "Nephrops norvegicus" = "*N. norvegicus*", 
                            "Scomber scombrus" = "*S. scombrus*"), expand = c(0, 0)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(text = element_text(size = 13, face = "bold")) + scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##For Report
dlegend <- get_legend((dby.g+labs(fill = "Area-Gear")))
dby.g2 <- dby.g + theme(legend.position = "none")
dtar.g <- dtar.g + theme(legend.position = "none")
grid.arrange((dtar.g+labs(tag="A")+ylab("TIS (%)")), blankplot, (dby.g2+labs(tag="B")+ylab("TIS (%)")), dlegend, ncol = 2, widths = c(2, 0.7))






####Moving Averages Spp6#####

mavg.6b <- function(data){
  t <- data.frame()
  t[1:10, 1:39] <- 0
for(j in c(10,11,21,28,34)){
  for(i in 2009:2018){
    selected <- c(i-1, i, i+1)
    f1 <- subset(data, Year %in% selected & sppcodeb %in% j)
    if((i == 2009 | (i == 2018))){
      t[i-2008, j] <- NA}
    else{t[i-2008, j] <- mean(f1$TIS)}
  }
}
t <- t[,-c(1:9,12:20,22:27,29:33,35:39)]
colnames(t) <- c("Dipturus batis", "Gadus morhua", "Melanogrammus aeglefinus", "Nephrops norvegicus", "Scomber scombrus")
row.names(t) <- c(2009:2018)
mt <- melt(t, id.vars = NULL)
mt[,3] <- c(2009:2018)
colnames(mt) <- c("Bycatch", "TIS", "Year")
mt$Bycatch <- as.factor(mt$Bycatch)
print(mt)
}

mavg.6t <- function(data){
  t <- data.frame()
  t[1:10, 1:39] <- 0
  for(j in c(10,11,21,28,34)){
    for(i in 2009:2018){
      selected <- c(i-1, i, i+1)
      f1 <- subset(data, Year %in% selected & sppcode %in% j)
      if((i == 2009 | (i == 2018))){
        t[i-2008, j] <- NA}
      else{t[i-2008, j] <- mean(f1$TIS)}
    }
  }
  t <- t[,-c(1:9,12:20,22:27,29:33,35:39)]
  colnames(t) <- c("Dipturus batis", "Gadus morhua", "Melanogrammus aeglefinus", "Nephrops norvegicus", "Scomber scombrus")
  row.names(t) <- c(2009:2018)
  mt <- melt(t,id.vars = NULL)
  mt[,3] <- c(2009:2018)
  colnames(mt) <- c("Targeted", "TIS", "Year")
  mt$Targeted <- as.factor(mt$Targeted)
  print(mt)
}

##1274
mavg1274by6 <- mavg.6b(data = mm3.1274)
mavg1274ta6 <- mavg.6t(data = mm3.1274)

pdf("12746 Mavg.pdf")
zb1 <- ggplot(mavg1274by6, aes(x = Year, y = TIS, color = Bycatch)) +
  geom_line()+ylim(0,100) + ylab("TIS (%)") + theme(legend.position = "none") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


zt1 <- ggplot(mavg1274ta6, aes(x = Year, y = TIS, color = Targeted))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##2274
mavg2274by6 <- mavg.6b(data = mm3.2274)
mavg2274ta6 <- mavg.6t(data = mm3.2274)

pdf("22746 Mavg.pdf")
zb2 <- ggplot(mavg2274by6, aes(x = Year, y = TIS, color = Bycatch))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


zt2 <- ggplot(mavg2274ta6, aes(x = Year, y = TIS, color = Targeted))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##1276a
mavg1276aby6 <- mavg.6b(data = mm3.1276a)
mavg1276ata6 <- mavg.6t(data = mm3.1276a)

pdf("1276a6 Mavg.pdf")
zb3 <- ggplot(mavg1276aby6, aes(x = Year, y = TIS, color = Bycatch))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

zt3 <- ggplot(mavg1276ata6, aes(x = Year, y = TIS, color = Targeted))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
##2276a
mavg2276aby6 <- mavg.6b(data = mm3.2276a)
mavg2276ata6 <- mavg.6t(data = mm3.2276a)

pdf("22746a6 Mavg.pdf")
zb4 <- ggplot(mavg2276aby6, aes(x = Year, y = TIS, color = Bycatch))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

zt4 <- ggplot(mavg2276ata6, aes(x = Year, y = TIS, color = Targeted))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##1276b
mavg1276bby6 <- mavg.6b(data = mm3.1276b)
mavg1276bta6 <- mavg.6t(data = mm3.1276b)

pdf("1276b6 Mavg.pdf")
zb5 <- ggplot(mavg1276bby6, aes(x = Year, y = TIS, color = Bycatch))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

zt5 <- ggplot(mavg1276bta6, aes(x = Year, y = TIS, color = Targeted))+geom_line()+ylim(0,100) +
  theme(legend.position = "none") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##O276b
mavgO276bby6 <- mavg.6b(data = mm3.O276b)
mavgO276bta6 <- mavg.6t(data = mm3.O276b)

pdf("O276b6 Mavg.pdf")
zb6 <- ggplot(mavgO276bby6, aes(x = Year, y = TIS, color = Bycatch))+geom_line()+ylim(0,100) +
  theme(legend.position = "right") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2", 
      labels=c("Dipturus batis" = "*D. batis*", 
        "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
        "Nephrops norvegicus" = "*N. norvegicus*", 
        "Scomber scombrus" = "*S. scombrus*")) + 
  theme(legend.text = ggtext::element_markdown()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

zt6 <- ggplot(mavgO276bta6, aes(x = Year, y = TIS, color = Targeted))+geom_line()+ylim(0,100) +
  theme(legend.position = "right") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2", 
      labels=c("Dipturus batis" = "*D. batis*", 
        "Gadus morhua" = "*G. morhua*", "Melanogrammus aeglefinus" = "*M. aeglefinus*", 
        "Nephrops norvegicus" = "*N. norvegicus*", 
        "Scomber scombrus" = "*S. scombrus*")) + 
  theme(legend.text = ggtext::element_markdown()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

###Multipanel plot
##Bycatch Spp6
#Code found on http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


blankplot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
zblegend <- get_legend(zb6)
zb6.2 <- zb6 + theme(legend.position = "none")
grid.arrange((zb1+labs(tag = "A", title = "27.4-TR1")+theme(axis.text.x = element_text(angle = 90))), (zb2+labs(tag = "B", title = "27.4-TR2")+theme(axis.text.x = element_text(angle = 90))), blankplot, (zb3+labs(tag = "C", title = "27.6.a-TR1")+theme(axis.text.x = element_text(angle = 90))),
             (zb4+labs(tag = "D", title = "27.6.a-TR2")+theme(axis.text.x = element_text(angle = 90))), blankplot, 
             (zb5+labs(tag = "E", title = "27.6.b-TR1")+theme(axis.text.x = element_text(angle = 90))), (zb6.2+labs(tag = "F", title = "27.6.b-TROther")+theme(axis.text.x = element_text(angle = 90))), 
             zblegend, ncol = 3)


##Targeted Spp6
blankplot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
ztlegend <- get_legend(zt6)
zt6.2 <- zt6 + theme(legend.position = "none")
grid.arrange((zt1+labs(tag = "A", title = "27.4-TR1")+theme(axis.text.x = element_text(angle = 90))), (zt2+labs(tag = "B", title = "27.4-TR2")+theme(axis.text.x = element_text(angle = 90))), blankplot, (zt3+labs(tag = "C", title = "27.6.a-TR1")+theme(axis.text.x = element_text(angle = 90))),
             (zt4+labs(tag = "D", title = "27.6.a-TR2")+theme(axis.text.x = element_text(angle = 90))), blankplot, 
             (zt5+labs(tag = "E", title = "27.6.b-TR1")+theme(axis.text.x = element_text(angle = 90))), (zt6.2+labs(tag = "F", title = "27.6.b-TROther")+theme(axis.text.x = element_text(angle = 90))), 
             ztlegend, ncol = 3)




####Weighted Averages for Overall AG Combos######
##Will be finding weighted averages using the following method:
  ######(Avg_2009*NumTrips_2009)+(Avg_2010*NumTrips_2010)...+(Avg_2018*NumTrips_2018)/Total Trips

head(merg2)# this has the avg for each year for each combination as well as the number of trips in that year

wavg_10yr <- data.frame()
wavg_10yr[1,1:6] <- 0

colnames(wavg_10yr) <- unique(merg2$Combos)
row.names(wavg_10yr) <- "Weighted_Average"


##Code found from http://flovv.github.io/Weighted_Mean_Standard_Errors/
weighted.summarySE <- function(data=NULL, measurevar,  groupvars=NULL, weights, na.rm=FALSE,
                               conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  #weighted - SD function!
  w.sd <- function(x, w,na.rm=TRUE )  ( (sum(w*x*x, na.rm=na.rm)/sum(w, na.rm=na.rm)) - weighted.mean(x,w, na.rm=na.rm)^2 )^.5
  
  
  # This does the summary. For each group's data frame, return a vector with
  datac <- ddply(data, groupvars,
                 .fun = function(xx, col, weights) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = weighted.mean(xx[[col]], xx[[weights]], na.rm=na.rm),
                     sd   = w.sd(xx[[col]], xx[[weights]], na.rm=na.rm)
                   )
                 },
                 measurevar, weights
  )
  

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


test <- weighted.summarySE(data = merg2, measurevar = "TIS", groupvars = "Combos", weights = "Trips")

##comparing non-weighted averages 
test
TIS_avg
##Basically the same


####Moving Averages AG Combos Compiled/Plotted#######
t1274
t2274
t1276a
t2276a
t1276b
tO276b

toverall_mavg <- cbind(t1274, t2274, t1276a, t2276a, t1276b, tO276b)
print(toverall_mavg)
colnames(toverall_mavg) <- c("27.4-TR1", "27.4-TR2", "27.6.a-TR1", "27.6.a-TR2", "27.6.b-TR1", "27.6.b-TROther")
toverall_mavg

?gather

toverall_mavg.2 <- toverall_mavg %>%  mutate(Year = 2009:2018)
toverall_mavg.g <- toverall_mavg.2 %>% gather(Area_Gear, TIS, "27.4-TR1":"27.6.b-TROther") %>%
  ggplot(aes(Year, TIS, color = Area_Gear)) +
  geom_line() + xlim(2009,2018) + labs(tag = "A", colour = "Area-Gear Combinations")

toverall_mavg.g <- toverall_mavg.g + geom_line() +
  theme(legend.position = "top") + ylab("TIS (%)") +
  theme(text = element_text(size = 13, face = "bold")) + scale_color_brewer(palette = "Dark2") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks = seq(40, 100, by = 5))


print(toverall_mavg.g)

pdf("Area-Gear MAvg Compiled.pdf")
dev.off()

grid.arrange(blankplot, (mtis.g+labs(tag = "B")+xlab("Area-Gear Combinations")), (toverall_mavg.g+labs(tag="A", colour = "Area-Gear") + 
                                                                                    theme(legend.position = "top") + ylab("TIS (%)")),
             (mic.g+labs(tag="C")+xlab("Area-Gear Combinations")), ncol = 2)

      
####Averaged Heatmaps for 27.4 and 27.6.a######
head(mm3.1274) 
head(mm3.2274)
head(mm3.1276a)
head(mm3.2276a)
head(mm3.1276b)
head(mm3.O276b)

##
##27.4 TR1
avghm1274 <- tapply(mm3.1274$TIS, list(mm3.1274$Targeted, mm3.1274$Bycatch), mean)
View(avghm1274)

avghm1274[is.na(avghm1274)] <- -1

##HClust 
dist_avghm1274b <- dist(t(avghm1274), method = "euclidean")#distance matrix
hc_avghm1274 <- hclust(dist_avghm1274, method = "ward.D2")


dendro.avghm1274 <- as.dendrogram(hc_avghm1274)


mmavghm1274 <- melt(avghm1274)
colnames(mmavghm1274) <- c("Targeted", "Bycatch", "TIS")
head(mmavghm1274)

mmavghm1274$TI_Ranges <- cut(mmavghm1274$TIS, breaks = c(0,5,50,75,100), right = F, include.lowest = T)

order.avghm1274 <- order.dendrogram(dendro.avghm1274)
mmavghm1274$Bycatch <- factor(x = mmavghm1274$Bycatch,
                               levels = unique(mmavghm1274$Bycatch)[order.avghm1274],
                               ordered = T)

mmavghm1274$Targeted <- factor(x = mmavghm1274$Targeted,
                               levels = mmavghm1274$Targeted[order.avghm1274],
                               ordered = T)

heatmap.avghm1274 <- ggplot(data=mmavghm1274, 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "top")+labs(fill = "TIS Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.avghm1274)

##
##27.4 TR2
avghm2274 <- tapply(mm3.2274$TIS, list(mm3.2274$Targeted, mm3.2274$Bycatch), mean)
View(avghm2274)

avghm2274[is.na(avghm2274)] <- -1

##HClust 
dist_avghm2274 <- dist(t(avghm2274), method = "euclidean") #distance matrix
hc_avghm2274 <- hclust(dist_avghm2274, method = "ward.D2")

dendro.avghm2274 <- as.dendrogram(hc_avghm2274)

mmavghm2274 <- melt(avghm2274)
colnames(mmavghm2274) <- c("Targeted", "Bycatch", "TIS")
head(mmavghm2274)

mmavghm2274$TI_Ranges <- cut(mmavghm2274$TIS, breaks = c(0,5,50,75,100), right = F, include.lowest = T)

order.avghm2274 <- order.dendrogram(dendro.avghm2274)
mmavghm2274$Bycatch <- factor(x = mmavghm2274$Bycatch,
                              levels = unique(mmavghm2274$Bycatch)[order.avghm2274],
                              ordered = T)

mmavghm2274$Targeted <- factor(x = mmavghm2274$Targeted,
                               levels = mmavghm2274$Targeted[order.avghm2274],
                               ordered = T)

heatmap.avghm2274 <- ggplot(data=mmavghm2274, 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "top")+labs(fill = "TIS Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.avghm2274)

##
##27.6.a TR1
avghm1276a <- tapply(mm3.1276a$TIS, list(mm3.1276a$Targeted, mm3.1276a$Bycatch), mean)
View(avghm1276a)

avghm1276a[is.na(avghm1276a)] <- -1

##HClust 
dist_avghm1276a <- dist(avghm1276a, method = "euclidean") #distance matrix
hc_avghm1276a <- hclust(dist_avghm1276a, method = "ward.D2")

dendro.avghm1276a <- as.dendrogram(hc_avghm1276a)

mmavghm1276a <- melt(avghm1276a)
colnames(mmavghm1276a) <- c("Targeted", "Bycatch", "TIS")
head(mmavghm1276a)

mmavghm1276a$TI_Ranges <- cut(mmavghm1276a$TIS, breaks = c(0,5,50,75,100), right = F, include.lowest = T)

order.avghm1276a <- order.dendrogram(dendro.avghm1276a)
mmavghm1276a$Bycatch <- factor(x = mmavghm1276a$Bycatch,
                              levels = unique(mmavghm1276a$Bycatch)[order.avghm1276a],
                              ordered = T)

mmavghm1276a$Targeted <- factor(x = mmavghm1276a$Targeted,
                               levels = mmavghm1276a$Targeted[order.avghm1276a],
                               ordered = T)

heatmap.avghm1276a <- ggplot(data=mmavghm1276a, 
                            aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "top")+labs(fill = "TIS Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.avghm1276a)

##
##27.6.a TR2
avghm2276a <- tapply(mm3.2276a$TIS, list(mm3.2276a$Targeted, mm3.2276a$Bycatch), mean)
View(avghm2276a)

avghm2276a[is.na(avghm2276a)] <- -1

##HClust 
dist_avghm2276a <- dist(t(avghm2276a), method = "euclidean") #distance matrix
hc_avghm2276a <- hclust(dist_avghm2276a, method = "ward.D2")

dendro.avghm2276a <- as.dendrogram(hc_avghm2276a)

mmavghm2276a <- melt(avghm2276a)
colnames(mmavghm2276a) <- c("Targeted", "Bycatch", "TIS")
head(mmavghm2276a)

mmavghm2276a$TI_Ranges <- cut(mmavghm2276a$TIS, breaks = c(0,5,50,75,100), right = F, include.lowest = T)

order.avghm2276a <- order.dendrogram(dendro.avghm2276a)
mmavghm2276a$Bycatch <- factor(x = mmavghm2276a$Bycatch,
                               levels = unique(mmavghm2276a$Bycatch)[order.avghm2276a],
                               ordered = T)

mmavghm2276a$Targeted <- factor(x = mmavghm2276a$Targeted,
                                levels = mmavghm2276a$Targeted[order.avghm2276a],
                                ordered = T)

heatmap.avghm2276a <- ggplot(data=mmavghm2276a, 
                             aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
  scale_fill_brewer(palette = "YlOrRd" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "top")+labs(fill = "TIS Ranges (%)") + theme(panel.background = element_blank())
print(heatmap.avghm2276a)

# ##
# ##27.6.b TR1
# avghm1276b <- tapply(mm3.1276b$TIS, list(mm3.1276b$Targeted, mm3.1276b$Bycatch), mean)
# View(avghm1276b)
# mmavghm1276b <- melt(avghm1276b)
# colnames(mmavghm1276b) <- c("Targeted", "Bycatch", "TIS")
# head(mmavghm1276b)
# 
# mmavghm1276b$TI_Ranges <- cut(mmavghm1276b$TIS, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# 
# heatmap.avg1276b6 <- ggplot(data=mmavghm1276b, 
#                             aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
# print(heatmap.avg1276b6)
# 
# ##
# ##27.6.b TRO
# avghmO276b <- tapply(mm3.O276b$TIS, list(mm3.O276b$Targeted, mm3.O276b$Bycatch), mean)
# View(avghmO276b)
# mmavghmO276b <- melt(avghmO276b)
# colnames(mmavghmO276b) <- c("Targeted", "Bycatch", "TIS")
# head(mmavghmO276b)
# 
# mmavghmO276b$TI_Ranges <- cut(mmavghmO276b$TIS, breaks = c(0,5,50,75,100), right = F, include.lowest = T)
# 
# heatmap.avgO276b6 <- ggplot(data=mmavghmO276b, 
#                             aes(y=Targeted, x=Bycatch, fill = TI_Ranges)) + geom_tile(color = "grey") + 
#   scale_fill_brewer(palette = "YlOrRd" ) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#         legend.position = "top")+labs(fill = "TI Ranges (%)") + theme(panel.background = element_blank())
# print(heatmap.avgO276b6)
# 
# unique(mmavghm1274$Targeted)
##
##Making the Figure
l1 <- heatmap.avghm1274 + theme(legend.position = "none") + labs(tag = "A") + 
  theme(text = element_text(face = "bold"), axis.text.x = element_text(face = c("italic")), 
        axis.text.y = element_text(face = c("italic"))) + 
  scale_y_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata")) +
  scale_x_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata"))
l2 <- heatmap.avghm2274 + labs(tag = "B") + theme(text = element_text(face = "bold")) + 
  theme(text = element_text(face = "bold"), axis.text.x = element_text(face = c("italic")), 
        axis.text.y = element_text(face = c("italic")))  + 
  scale_y_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata")) +
  scale_x_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata")) 
l3 <-heatmap.avghm1276a + theme(legend.position = "none") + labs(tag = "C") + theme(text = element_text(face = "bold")) + 
  theme(text = element_text(face = "bold"), axis.text.x = element_text(face = c("italic")), 
        axis.text.y = element_text(face = c("italic")))  + 
  scale_y_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata")) +
  scale_x_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata"))
l4 <-heatmap.avghm2276a + theme(legend.position = "none") + labs(tag = "D") + theme(text = element_text(face = "bold")) + 
  theme(text = element_text(face = "bold"), axis.text.x = element_text(face = c("italic")), 
        axis.text.y = element_text(face = c("italic")))  + 
  scale_y_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata")) +
  scale_x_discrete(labels = c("Trachurus trachurus" = "T. trachurus", "Squalus acanthias" = "S. acanthias", "Sebastes norvegicus" = "S. norvegicus",
                              "Sebastes mentella" = "S. mentella", "Scophthalmus maximus" = "S. maximus", "Scomber scombrus" = "S. scombrus", 
                              "Rostroraja alba " = "R. alba", "Pollachius virens" = "P. virens", "Pollachius pollachius" = "P. pollachius", 
                              "Pleuronectes platessa" = "P. platessa", "Other" = "Other", "Nephrops norvegicus" = "N. norvegicus", "Molva molva" = "M. movla",
                              "Molva dypterygia" = "M. dypterygia", "Microstomus kitt" = "M. kitt", "Micromesistius poutassou" = "M. poutassou", 
                              "Merluccius merluccius" = "M. merluccius", "Merlangius merlangus" = "M. merlangus", "Melanogrammus aeglefinus" = "M. aeglefinus",
                              "Lophius piscatorius" = "L. piscatorius", "Lophius budegassa" = "L. budegassa", "Leucoraja fullonica" = "L. fullonica", 
                              "Leucoraja circularis" = "L. circularis", "Lepidorhombus whiffiagonis" = "L. whiffiagonis", "Hoplostethus atlanticus" = "H. atlanticus",
                              "Hippoglossus hippoglossus" = "H. hippoglossus", "Glyptocephalus cynoglossus" = "G. cynoglossus", "Galeorhinus galeus" = "G. galeus", 
                              "Gadus morhua" = "G. morhua", "Dipturus batis" = "D. batis", "Deania calcea" = "D. calcea", "Dasyatis pastinaca" = "D. pastinaca",
                              "Dalatias licha" = "D. licha", "Coryphaenoides rupestris" = "C. rupestris", "Clupea harengus" = "C. harengus", "Cetorhinus maximus" = "C. maximus",
                              "Centroscymnus coelolepis" = "C. coelolepis", "Centrophorus squamosus" = "C. squamosus", "Amblyraja radiata" = "A. radiata"))
#l5 <-heatmap.avg1276b6 + theme(legend.position = "none")
#l6 <-heatmap.avgO276b6 + theme(legend.position = "none")

print(l1)
print(l2)
print(l3)
print(l4)

llegend <- get_legend(l2)
l2.2 <- l2+theme(legend.position = "none")
grid.arrange(blankplot, llegend, l1, l2.2, nrow = 2, heights = c(0.1,1))
grid.arrange(l3, l4, nrow = 1)
###
###