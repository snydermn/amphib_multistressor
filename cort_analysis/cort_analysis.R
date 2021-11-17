
#This document combines the corticosterone data from the individual EIAs into one large dataset. That 
#dataset is used to investigate the effect of the treatments on the corticosterone level of individual frogs. 
#Corticosterone was measured in the water at the beginning of the 48 hour exposure and at the end of the 
#exposure. A standard curve was created for each EIA plate individually as abundance units vs. log (concentration) values. 


library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(AICcmodavg)

#Import files used:

# setwd("~/git/amphib_multistressor/cort_analysis")
#import sample identification key
#key_022916 <- read.csv("~/Dropbox/amphib_metabolomics/PIP/data/CORT/final/key_022916.csv")
#import cort values for plate1 02/22/16
cort_plate5_022216 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate5_022216.csv")
plate1<-cort_plate5_022216
dim(plate1)
#import cort values for plate2 02/23/16
cort_plate6_022316 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate6_022316.csv")
plate2<-cort_plate6_022316
#import cort values for plate3 02/26/16
cort_plate14_022616 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate14_022616.csv")
plate3<-cort_plate14_022616
#import cort values for plate4 02/29/16
cort_plate7_022916 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate7_022916.csv")
plate4<-cort_plate7_022916
#import cort values for plate5 03/02/16
cort_plate5_030216 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate5_030216.csv")
plate5<-cort_plate5_030216
#import key for plate 1
EIA_Cort_022216_412_KEY <- read.csv("c://git/amphib_multistressor/cort_analysis/EIA_Cort_022216_412_KEY.csv")
key1<-EIA_Cort_022216_412_KEY
#import key for plate 2
EIA_Cort_022316_412_key <- read.csv("c://git/amphib_multistressor/cort_analysis/EIA_Cort_022316_412_key.csv")
key2<-EIA_Cort_022316_412_key
#import key for plate 3
EIA_Cort_022616_key <- read.csv("c://git/amphib_multistressor/cort_analysis/EIA_Cort_022616_key.csv")
key3<-EIA_Cort_022616_key
#import key for plate 4
key_022916 <- read.csv("c://git/amphib_multistressor/cort_analysis/key_022916.csv")
key4<-key_022916
#import key for plate 5
key_030216 <- read.csv("c://git/amphib_multistressor/cort_analysis/key_030216.csv")
key5<-key_030216


# Average duplicate values for samples and merge
#### plate 1 ####
#convert to lower case
plate1$cell<-tolower(plate1$cell)
plate1a<-merge(plate1, key1, by.x='cell', by.y='plate_key')
#View(plate1a)
#get mean of sample duplicate values
plate1.mean<-ddply(plate1a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate1b<-merge(key1, plate1.mean, by.x='sample_name', by.y='sample_name')
#View(plate1b)
#remove duplicate rows
plate1c<-plate1b[!duplicated(plate1b$sample_name),]
plate1c$plate<-'plate1'
#View(plate1c)

#### plate 2 ####
#convert to lower case
plate2$cell<-tolower(plate2$cell)
plate2a<-merge(plate2, key2, by.x='cell', by.y='plate_key')
#View(plate2a)
#get mean of sample duplicate values
plate2.mean<-ddply(plate2a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate2b<-merge(key2, plate2.mean, by.x='sample_name', by.y='sample_name')
#View(plate2b)
#remove duplicate rows
plate2c<-plate2b[!duplicated(plate2b$sample_name),]
plate2c$plate<-'plate2'
#View(plate2c)

#### plate 3 ####
#convert to lower case
plate3$cell<-tolower(plate3$cell)
plate3a<-merge(plate3, key3, by.x='cell', by.y='plate_key')
#View(plate3a)
#get mean of sample duplicate values
plate3.mean<-ddply(plate3a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate3b<-merge(key3, plate3.mean, by.x='sample_name', by.y='sample_name')
#View(plate3b)
#remove duplicate rows
plate3c<-plate3b[!duplicated(plate3b$sample_name),]
plate3c$plate<-'plate3'
#View(plate3c)

#### plate 4 ####
#convert to lower case
plate4$cell<-tolower(plate4$cell)
plate4a<-merge(plate4, key4, by.x='cell', by.y='plate_key')
#View(plate4a)
#get mean of sample duplicate values
plate4.mean<-ddply(plate4a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate4b<-merge(key4, plate4.mean, by.x='sample_name', by.y='sample_name')
#View(plate4b)
#remove duplicate rows
plate4c<-plate4b[!duplicated(plate4b$sample_name),]
plate4c$plate<-'plate4'
#View(plate4c)

#### plate 5 ####
#convert to lower case
plate5$cell<-tolower(plate5$cell)
plate5a<-merge(plate5, key5, by.x='cell', by.y='plate_key')
#View(plate5a)
#get mean of sample duplicate values
plate5.mean<-ddply(plate5a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate5b<-merge(key5, plate5.mean, by.x='sample_name', by.y='sample_name')
#View(plate5b)
#remove duplicate rows
plate5c<-plate5b[!duplicated(plate5b$sample_name),]
plate5c$plate<-'plate5'
#View(plate5c)

#### combine dataframes from each plate into 1 df #### 
# convert numeric to factor for experiment column
plate2c$experiment<-as.factor(as.character(plate2c$experiment))
plate3c$experiment<-as.factor(as.character(plate3c$experiment))
plate4c$experiment<-as.factor(as.character(plate4c$experiment))
plate5c$experiment<-as.factor(as.character(plate5c$experiment))
#combine all plates into 1 df
cort_all<-rbind(plate1c, plate2c, plate3c, plate4c, plate5c)
#View(cort_all)
# sort by sample number
cort_all_sort<-cort_all[order(cort_all$sample_no),]
#sum samples with multiple cort cells
c49<-(1.284533e+03+2.092193e+01) 
c54<-(5.855812e+01+4.860801e+00)
s26<-(2.235718e+02+6.634500e+00)
s15<-(30.35141825+1.357286e+02)
#add1<-c(as.factor(as.character("exp4_C49base")),as.factor(as.character("a11g10")),as.factor(as.character("4")), as.factor(as.character("1")), "C49", as.numeric(c49))
cort_all_sort$sample_name<-as.factor(as.character(cort_all_sort$sample_name))
#cort_all_sort1<-rbind(cort_all_sort, add1)
write.csv(cort_all_sort, "cort_out.csv")


# Anova for experimental day effects
#reading court_out3.csv--I assume some data corrections here
cort_out2 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_out3.csv")
#View(cort_out2)
cort1<-cort_out2
cort1$experiment<-as.factor(as.character(cort1$experiment))
summary(cort1)
anova1<-lm(V1~experiment, data=subset(cort1, time=="1"))
summary(anova1)
anova2<-lm(V1~experiment, data=subset(cort1, time=="2"))
summary(anova2)

boxplot(V1~experiment, data=subset(cort1, time=="2"), main="Cort End")
boxplot(V1~experiment, data=subset(cort1, time=="1"), main= "Cort Beginning")

#create dataframe separate for initial and final time periods
cort1.base<-cort1[cort1$time=="1",]
#View(cort1.base)
cort1.end<-cort1[cort1$time=="2",]
#View(cort1.end)

#merge based on experiment number and sample name
cort2<-merge(cort1.base, cort1.end, by=c("sample_no", "experiment"))
summary(cort2)
#View(cort2)
dim(cort2)
#subtract final measurement from initial cort measurement
cort2$cort.diff<-cort2$V1.y-cort2$V1.x

anova3<-lm(cort.diff~experiment, data=cort2)
summary(anova3)
boxplot(cort.diff~experiment, data=cort2, main="Cort difference")


#Treatment 1= control
#Treatment 2= pesticide
#Treatment 3= bullfrog
#Treatment 4= bullfrog+pesticide

#Anova and boxplots for treatment effects:

#adds metadata
sample_list <- read.csv("c://git/amphib_multistressor/cort_analysis/sample_list.csv")
cort3<-merge(cort1, sample_list, by=c("sample_no", "experiment"))
#View(cort3)
cort3$treatment<-as.factor(as.character(cort3$treatment))
cort3$time<-as.factor(as.character(cort3$time))
#anovas for effect of treatment on cort level at beginning and end
m1<-lm(V1~treatment, data=subset(cort3, time=="1"))
summary(m1)
m2<-lm(V1~treatment, data=subset(cort3, time=="2"))
summary(m2)

boxplot(V1~treatment, data=subset(cort3, time=="2"), xlab="treatment", names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="End measurement")
boxplot(V1~treatment, data=subset(cort3, time=="1"), xlab="treatment", names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="Base measurement")

#create dataframe separate for initial and final time periods
cort3.base<-cort3[cort3$time=="1",]
#View(cort3.base)
cort3.end<-cort3[cort3$time=="2",]
#View(cort3.end)
cort4<-merge(cort3.base, cort3.end, by=c("sample_no", "experiment"))
#cort4<-merge(cort3, sample_list, by=c("sample_no", "experiment"))
cort4$cort.diff<-cort4$V1.y-cort4$V1.x
#View(cort4)

#anova for treatment effect on cort difference
m3<-lm(cort.diff~treatment.x, data=cort4)
summary(m3)

boxplot(cort.diff~treatment.x, data=cort4, xlab="treatment.x", names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="difference measurement")
#anova for EIA plate effect on cort difference
model.plate.beg<-lm(V1~plate, data=subset(cort3, time=="1"))
summary(model.plate.beg)
model.plate.end<-lm(V1~plate, data=subset(cort3, time=="2"))
summary(model.plate.end)


# Linear models with more than one predictor factor
m3<-lm(V1.y~treatment.x+plate.x, data=cort4)
summary(m3)

model.plate.end2<-lm(V1~treatment+plate, data=subset(cort3, time=="2"))
summary(model.plate.end2)


#When comparing the experiment time period, the significant differences are in the comparison of 
#cort levels to experimental time period is that the beginning base line cort level between 
#experiment 3 and 4. For end time period intercept and experimental time 2 is significantly 
#different. When comparing treatment levels to corticosterone levels, Treatment 2 and 4 is significantly 
#different in end time period Cort levels. In difference cort levels or beginning cort levels there is no 
#signficant treatment effects. Treatment 2 is carbaryl only. Treatment 4 is carbaryl and bullfrog predation 
#stress. When looking for effects of plate, plate4 is significantly different than the other plates. When 
#we went back and examined the standard curve on plate 4 the lower abundance values which correspond to the 
#higher cort levels are quite different than for the other plates. We compared slope and intercepts for the 
#standard curve for all the plates and found that plate4 is the most different so we applied an average of 
#the slope and intercept from all the other plates to the plate 4 data. We also ran the models 
#without plate 4 data. 

#ANOVAS without plate 4 for cort beginning and end time periods:

#create dataset without plate4
cort5<-cort3[cort3$plate!="plate4",]
summary(cort5)
#anova without plate4 for end time period
model.end.woplate4<-lm(V1~treatment+plate, data=subset(cort5, time=="2"))
summary(model.end.woplate4)
model.end.woplate5<-lm(V1~Pesticide+predator+plate, data=subset(cort5, time=="2"))
summary(model.end.woplate5)

boxplot(V1~treatment, data=subset(cort5, time=="2"))
model.end.woplate4.2<-lm(V1~treatment, data=subset(cort5, time=="2"))
summary(model.end.woplate4.2)
model.end.woplate4.3<-lm(V1~plate, data=subset(cort5, time=="2"))
summary(model.end.woplate4.3)

#anova without plate4 for base time period
model.end.woplate4<-lm(V1~treatment+plate, data=subset(cort5, time=="1"))
summary(model.end.woplate4)
boxplot(V1~treatment, data=subset(cort5, time=="1"))
model.end.woplate4.2<-lm(V1~treatment, data=subset(cort5, time=="1"))
summary(model.end.woplate4.2)


#Exploratory plots:
plot(cort4$V1.x, cort4$cort.diff)
abline(0, 0)
plot(cort4$V1.y, cort4$cort.diff)
abline(0,0)
summary(cort4$cort.diff)
summary(cort4$V1.x)
summary(cort4$V1.y)
#make id for experiment sample number to compare base and end measurements for same frog 
cort4$id <- with(cort4, interaction(experiment,  sample_no))
plot(cort4$V1.x, cort4$V1.y)
#plot base and end cort measurements for same frog on same plot
plot(cort4$V1.y~cort4$id, pch=1, col="red" )
points(cort4$V1.x~cort4$id, col="red")



################################################################
### Starting over?
#Re-calibrated plate4 standard curve with average slope and intercept from other plates. Reimporting and re-analyzing data.
###################################################################
#Import files used:

#import sample identification key
#key_022916 <- read.csv("~/Dropbox/amphib_metabolomics/PIP/data/CORT/final/key_022916.csv")
#import cort values for plate1 02/22/16
cort_plate5_022216 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate5_022216.csv")
plate1<-cort_plate5_022216
#import cort values for plate2 02/23/16
cort_plate6_022316 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate6_022316.csv")
plate2<-cort_plate6_022316
#import cort values for plate3 02/26/16
cort_plate14_022616 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate14_022616.csv")
plate3<-cort_plate14_022616
#import cort values for plate4 02/29/16
cort_plate7_022916 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate7_022916_2.csv")
plate4<-cort_plate7_022916
#import cort values for plate5 03/02/16
cort_plate5_030216 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_plate5_030216.csv")
plate5<-cort_plate5_030216
#import key for plate 1
EIA_Cort_022216_412_KEY <- read.csv("c://git/amphib_multistressor/cort_analysis/EIA_Cort_022216_412_KEY.csv")
key1<-EIA_Cort_022216_412_KEY
#import key for plate 2
EIA_Cort_022316_412_key <- read.csv("c://git/amphib_multistressor/cort_analysis/EIA_Cort_022316_412_key.csv")
key2<-EIA_Cort_022316_412_key
#import key for plate 3
EIA_Cort_022616_key <- read.csv("c://git/amphib_multistressor/cort_analysis/EIA_Cort_022616_key.csv")
key3<-EIA_Cort_022616_key
#import key for plate 4
key_022916 <- read.csv("c://git/amphib_multistressor/cort_analysis/key_022916.csv")
key4<-key_022916
#import key for plate 5
key_030216 <- read.csv("c://git/amphib_multistressor/cort_analysis/key_030216.csv")
key5<-key_030216


#Average duplicate values for samples and merge


#### plate 1 ####
#convert to lower case
plate1$cell<-tolower(plate1$cell)
plate1a<-merge(plate1, key1, by.x='cell', by.y='plate_key')
#View(plate1a)
#get mean of sample duplicate values
plate1.mean<-ddply(plate1a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate1b<-merge(key1, plate1.mean, by.x='sample_name', by.y='sample_name')
#View(plate1b)
#remove duplicate rows
plate1c<-plate1b[!duplicated(plate1b$sample_name),]
plate1c$plate<-'plate1'
#View(plate1c)

#### plate 2 ####
#convert to lower case
plate2$cell<-tolower(plate2$cell)
plate2a<-merge(plate2, key2, by.x='cell', by.y='plate_key')
#View(plate2a)
#get mean of sample duplicate values
plate2.mean<-ddply(plate2a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate2b<-merge(key2, plate2.mean, by.x='sample_name', by.y='sample_name')
#View(plate2b)
#remove duplicate rows
plate2c<-plate2b[!duplicated(plate2b$sample_name),]
plate2c$plate<-'plate2'
#View(plate2c)

#### plate 3 ####
#convert to lower case
plate3$cell<-tolower(plate3$cell)
plate3a<-merge(plate3, key3, by.x='cell', by.y='plate_key')
#View(plate3a)
#get mean of sample duplicate values
plate3.mean<-ddply(plate3a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate3b<-merge(key3, plate3.mean, by.x='sample_name', by.y='sample_name')
#View(plate3b)
#remove duplicate rows
plate3c<-plate3b[!duplicated(plate3b$sample_name),]
plate3c$plate<-'plate3'
#View(plate3c)

#### plate 4 ####
#convert to lower case
plate4$cell<-tolower(plate4$cell)
plate4a<-merge(plate4, key4, by.x='cell', by.y='plate_key')
#View(plate4a)
#get mean of sample duplicate values
plate4.mean<-ddply(plate4a, .(sample_name), function(x) mean(x$pg.ml))
#merge df with means with original plate key
plate4b<-merge(key4, plate4.mean, by.x='sample_name', by.y='sample_name')
#View(plate4b)
#remove duplicate rows
plate4c<-plate4b[!duplicated(plate4b$sample_name),]
plate4c$plate<-'plate4'
#View(plate4c)

#### plate 5 ####
#convert to lower case
plate5$cell<-tolower(plate5$cell)
plate5a<-merge(plate5, key5, by.x='cell', by.y='plate_key')
#View(plate5a)
#get mean of sample duplicate values
plate5.mean<-ddply(plate5a, .(sample_name), function(x) mean(x$cort.pg.ml))
#merge df with means with original plate key
plate5b<-merge(key5, plate5.mean, by.x='sample_name', by.y='sample_name')
#View(plate5b)
#remove duplicate rows
plate5c<-plate5b[!duplicated(plate5b$sample_name),]
plate5c$plate<-'plate5'
#View(plate5c)

#### combine dataframes from each plate into 1 df #### 
# convert numeric to factor for experiment column
plate2c$experiment<-as.factor(as.character(plate2c$experiment))
plate3c$experiment<-as.factor(as.character(plate3c$experiment))
plate4c$experiment<-as.factor(as.character(plate4c$experiment))
plate5c$experiment<-as.factor(as.character(plate5c$experiment))
#combine all plates into 1 df
cort_all<-rbind(plate1c, plate2c, plate3c, plate4c, plate5c)
#View(cort_all)
# sort by sample number
cort_all_sort<-cort_all[order(cort_all$sample_no),]
#sum samples with multiple cort cells
c49<-(1.284533e+03+2.092193e+01) 
c54<-(5.855812e+01+4.860801e+00)
s26<-(2.235718e+02+6.634500e+00)
s15<-(30.35141825+1.357286e+02)
#add1<-c(as.factor(as.character("exp4_C49base")),as.factor(as.character("a11g10")),as.factor(as.character("4")), as.factor(as.character("1")), "C49", as.numeric(c49))
cort_all_sort$sample_name<-as.factor(as.character(cort_all_sort$sample_name))
#cort_all_sort1<-rbind(cort_all_sort, add1)
write.csv(cort_all_sort, "cort_out.csv")



#Anova for experimental day effects

cort_out4 <- read.csv("c://git/amphib_multistressor/cort_analysis/cort_out4.csv")
#View(cort_out4)
scort1<-cort_out4
scort1$experiment<-as.factor(as.character(scort1$experiment))
summary(scort1)
anova1<-lm(V1~experiment, data=subset(scort1, time=="1"))
summary(anova1)
anova2<-lm(V1~experiment, data=subset(scort1, time=="2"))
summary(anova2)

boxplot(V1~experiment, data=subset(scort1, time=="2"), main="Cort End")
boxplot(V1~experiment, data=subset(scort1, time=="1"), main= "Cort Beginning")

#create dataframe separate for initial and final time periods
scort1.base<-scort1[scort1$time=="1",]
#View(scort1.base)
scort1.end<-scort1[scort1$time=="2",]
#View(scort1.end)

#merge based on experiment number and sample name
scort2<-merge(scort1.base, scort1.end, by=c("sample_no", "experiment"))
summary(scort2)
#View(scort2)
dim(scort2)
#subtract final measurement from initial cort measurement
scort2$cort.diff<-scort2$V1.y-scort2$V1.x

anova3<-lm(cort.diff~experiment, data=scort2)
summary(anova3)
boxplot(cort.diff~experiment, data=scort2, main="Cort difference")

sample_list <- read.csv("c://git/amphib_multistressor/cort_analysis/sample_list.csv")
scort3<-merge(scort1, sample_list, by=c("sample_no", "experiment"))
#View(scort3)
scort3$treatment<-as.factor(as.character(scort3$treatment))
scort3$time<-as.factor(as.character(scort3$time))
#anovas for effect of treatment on cort level at beginning and end
m1<-lm(V1~treatment, data=subset(scort3, time=="1"))
summary(m1)
m2<-lm(V1~treatment, data=subset(scort3, time=="2"))
summary(m2)

# this is the one MNS originally had in the draft manuscript
boxplot(V1~treatment, data=subset(scort3, time=="2"), xlab="treatment", names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="End measurement")
boxplot(V1~treatment, data=subset(scort3, time=="1"), xlab="treatment", names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="Base measurement")

#Merge dfs together to create cort difference column

#create dataframe separate for initial and final time periods
scort3.base<-scort3[scort3$time=="1",]
#View(scort3.base)
scort3.end<-scort3[scort3$time=="2",]
#View(scort3.end)
scort4<-merge(scort3.base, scort3.end, by=c("sample_no", "experiment"))
#cort4<-merge(cort3, sample_list, by=c("sample_no", "experiment"))

# transform difference cocentrations to hourly rates
# by dividing before concentration by 24 and after concentration by 48
scort4$cort.diff<-(scort4$V1.y/48)-(scort4$V1.x/24)


# this is (not) the anova for the manuscript
#View(scort4)
#anova for treatment effect on cort difference
m3<-lm(cort.diff~treatment.x, data=scort4)
summary(m3)

#original figure but without absolute value for differences
boxplot(cort.diff~treatment.x, data=scort4, xlab="treatment.x", 
        names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="")


#cort begin time
m4<-lm(V1.x~treatment.x+experiment, data=scort4)
summary(m4)

#cort end time models
m4<-lm(V1.y~treatment.x+experiment, data=scort4)
summary(m4)
m4<-lm(V1.y~treatment.x+plate.x+experiment, data=scort4)
summary(m4)
m5<-lm(V1.y~Pesticide.x+predator.x+experiment+plate.x, data=scort4)
summary(m5)
m5<-lm(V1.y~experiment, data=scort4)
summary(m5)
m7<-lm(V1.y~Pesticide.x*predator.x+experiment, data=scort4)
summary(m7)
m8<-lm(V1.y~Pesticide.x*predator.x, data=scort4)
summary(m8)
#boxplot cort endtime
boxplot(V1.y~Pesticide.x+predator.x, data=scort4, main="Cort final")

#cort difference 
m9<-lm(cort.diff~Pesticide.x+predator.x+experiment, data=scort4)
summary(m9)
m10<-lm(cort.diff~treatment.x+experiment, data=scort4)
summary(m10)
m11<-lm(cort.diff~Pesticide.x+predator.x, data=scort4)
summary(m11)
m12<-lm(cort.diff~predator.x+experiment, data=scort4)
summary(m12)
plot(m12)


boxplot(cort.diff~Pesticide.x+predator.x, data=scort4, main="Cort diff")



# The boxplots indicate some skew of the data. To help offset some of the skew I will run models with natural log of cort values

#histograms of cort values
#cort diff
hist(scort4$cort.diff)
hist(log(scort4$cort.diff+1000))
hist(log10(scort4$cort.diff+1000))
#cort end
hist(scort4$V1.y)
hist(log(scort4$V1.y))
hist(log10(scort4$V1.y))

#models of cort diff ln
m9<-lm(log(cort.diff+1000)~Pesticide.x+predator.x+experiment, data=scort4)
summary(m9) #sig pesticide and experiment
boxplot(log(cort.diff+1000)~Pesticide.x+predator.x, data=scort4, main="Cort diff")
m10<-lm(log(cort.diff+1000)~treatment.x+experiment, data=scort4)
summary(m10)
m11<-lm(log(cort.diff+1000)~Pesticide.x+predator.x, data=scort4)
summary(m11)
m12<-lm(log(cort.diff+1000)~predator.x+experiment, data=scort4)
summary(m12)

# model of end cort ln
m4<-lm(log(V1.y)~treatment.x+experiment, data=scort4)
summary(m4)
m4<-lm(log(V1.y)~treatment.x+plate.x+experiment, data=scort4)
summary(m4)
m5<-lm(log(V1.y)~Pesticide.x+predator.x+experiment+plate.x, data=scort4)
summary(m5)
m5<-lm(log(V1.y)~experiment, data=scort4)
summary(m5)
m7<-lm(log(V1.y)~Pesticide.x*predator.x+experiment, data=scort4)
summary(m7)
m8<-lm(log(V1.y)~Pesticide.x*predator.x, data=scort4)
summary(m8)
m9<-lm(log(V1.y)~Pesticide.x*predator.x+plate.x, data=scort4)
summary(m9)
boxplot(log(V1.y)~Pesticide.x+predator.x, data=scort4, main="Cort final")

######################################################
#### FINAL FIGURE AND ANOVAS
View(scort4)
dim(scort3) # long format
summary(scort3) 
dim(scort4) # tall format
summary(scort4)

#final figure
#boxplot of begin and end cort by treatment
bp_contrast <- ggplot(scort3, aes(x=treatment, y=V1, fill=time)) + 
  geom_boxplot() +
  labs(title="",x="", y = "Cort (pg/mL)") +
  scale_fill_discrete(labels = c("Start", "End")) +
  theme_classic() +
  theme(legend.position="top")
bp_contrast

bp_diff <- ggplot(scort4, aes(x=treatment.x, y=cort.diff)) + 
  geom_boxplot(fill='plum4') +
  labs(title="",x="Treatment", y = "Cort Difference (pg/mL/hr)") +
  scale_x_discrete(labels=c("Control","Carbaryl","Predator","Carbaryl+Predator")) +
  theme_classic()
bp_diff

bp_final_figure <- ggarrange(bp_contrast, bp_diff, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
bp_final_figure


cort_barplot <- paste("c://git/amphib_multistressor/cort_analysis/snyder_cort_figure.jpg",sep="")
jpeg(cort_barplot, width = 5, height = 7, units = "in",res=600)
  bp_final_figure
dev.off()

#boxplot(cort.diff~treatment.x, data=scort4, xlab="treatment.x", 
#        names=c("control","pesticide","predator","pesticide+predator"), ylab="Cort pg/mL", main="")


# final aov wo and with interaction, plus a blocking term for plate
one_way_predator <- aov(cort.diff ~ predator.x, data = scort4)
summary(one_way_predator)

one_way_pesticide <- aov(cort.diff ~ Pesticide.y, data = scort4)
summary(one_way_pesticide)

one_way_plate <- aov(cort.diff ~ plate.y, data = scort4)
summary(one_way_plate)

two_way_cort <- aov(cort.diff ~ predator.x + Pesticide.y, data = scort4)
summary(two_way_cort)

two_way_cort_block <- aov(cort.diff ~ predator.x + Pesticide.y + plate.y, data = scort4)
summary(two_way_cort_block)

two_way_cort_interaction <- aov(cort.diff ~ predator.x * Pesticide.y, data = scort4)
summary(two_way_cort_interaction)

two_way_cort_interaction_block <- aov(cort.diff ~ predator.x * Pesticide.y + plate.y, data = scort4)
summary(two_way_cort_interaction_block)

# AIC R code to determine best fit model
library(AICcmodavg)

model.set <- list(one_way_predator, one_way_pesticide, one_way_plate,
                  two_way_cort, two_way_cort_block, 
                  two_way_cort_interaction, two_way_cort_interaction_block)
model.names <- c("one_way_predator", "one_way_pesticide", "one_way_plate",
                 "two_way_cort", "two_way_cort_block",
                "two_way_cort_interaction", "two_way_cort_interaction_block")

aictab(model.set, modnames = model.names)
#Model selection based on AICc:
#  
#                              K   AICc Delta_AICc AICcWt Cum.Wt      LL
#one_way_pesticide              3 513.22       0.00   0.42   0.42 -253.41
#two_way_cort                   4 513.88       0.66   0.30   0.73 -252.59
#two_way_cort_block             7 515.41       2.19   0.14   0.87 -249.67
#two_way_cort_interaction       5 516.21       2.99   0.09   0.96 -252.57
#two_way_cort_interaction_block 8 518.04       4.82   0.04   1.00 -249.66
#one_way_predator               3 527.57      14.34   0.00   1.00 -260.58
#one_way_plate                  5 531.94      18.72   0.00   1.00 -260.43

summary(one_way_pesticide)
#Df Sum Sq Mean Sq F value  Pr(>F)    
#Pesticide.y  1   3944    3944   18.37 6.7e-05 ***
#  Residuals   60  12882     215                    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

TukeyHSD(one_way_pesticide)
#Tukey multiple comparisons of means
#95% family-wise confidence level
#
#Fit: aov(formula = cort.diff ~ Pesticide.y, data = scort4)
#
#$Pesticide.y
#diff      lwr       upr   p adj
#control-carbaryl -15.98387 -23.4442 -8.523542 6.7e-05

summary(two_way_cort)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#predator.x   1    590     590   2.775 0.101037    
#Pesticide.y  1   3688    3688  17.340 0.000103 ***
#  Residuals   59  12548     213                     
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

TukeyHSD(two_way_cort)
#Tukey multiple comparisons of means
#95% family-wise confidence level
#
#Fit: aov(formula = cort.diff ~ predator.x + Pesticide.y, data = scort4)
#
#$predator.x
#diff       lwr      upr     p adj
#yes-no 6.17078 -1.241291 13.58285 0.1010365
#
#$Pesticide.y
#diff       lwr       upr     p adj
#control-carbaryl -15.3842 -22.81175 -7.956658 0.0001103

