############################################################################################
# script used to subset peaks based on magnitude of difference betweeen controls and treatments
# MNSnyder
# Oct. 22, 2019
##############################################################################################

# libraries
library(tidyverse)

#### set up directory structure ####
# to run script on different computers specify sys.info()[4] and
# specify path to model results folder and
# folder containing look up tables

#### class and treatment key ####
# 1 = control
# 2 = pesticide
# 3 = bullfrog
# 4 = pesticide + bullfrog


#### Tom's laptop ####
if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  datadir<-path.expand("c:/git/PIP/GCMS_data_XCMS/") # path to git repo
}

# Tom's linux
if(Sys.info()[4]=="d2626ut7920d.rtpnc.epa.gov"){
  datadir<-path.expand("~/git/PIP/GCMS_data_XCMS/") # path to git repo
}

# Marcia's epa computer -> WORKS
if(Sys.info()[4]=="LZ2626XMSNYDE02"){
  datadir<-path.expand("D:/git_repo/PIP/GCMS_data_XCMS/") # path to git repo
}

setwd(datadir)


#### read in auto-scaled xcms aligned data #### 
# output from Analysis_101119.RMD #
# set up file names
liver1_class_filename <- paste0(datadir,"WorkDir_liver1/Preprocessing_Data_a/class.csv")
liver1_processedtable_filename <- paste0(datadir,"WorkDir_liver1/Preprocessing_Data_a/ProcessedTable.csv")
# read in data using filenames 
liver_a<-read.csv(liver1_processedtable_filename)
liver_a_class<-read.csv(liver1_class_filename)
# combine class with data
liver<-cbind(liver_a,liver_a_class[2] )
# sort by class
liver_o<-liver[order(liver$V1),]

#### data prep for comparing values ####
# create DF for each treatment of scaled peaks 
liver_t1<-dplyr::filter(liver_o, V1 == 1)
liver_t2<-dplyr::filter(liver_o, V1 == 2)
liver_t3<-dplyr::filter(liver_o, V1 == 3)
liver_t4<-dplyr::filter(liver_o, V1 == 4)

# drop non-numeric columns
liver_t1<-dplyr::select(liver_t1, -X, -V1)
liver_t2<-dplyr::select(liver_t2, -X, -V1)
liver_t3<-dplyr::select(liver_t3, -X, -V1)
liver_t4<-dplyr::select(liver_t4, -X, -V1)

#### get control means for peaks and compare to each treatment ####

# get mean value of each peak across all control treatments
liver_t1_mean<-colMeans(liver_t1)

# each treatment minus the control means

# control - treatment 2
# create empty DF from column names from original
diff1<-liver_t2[FALSE,]
for (i in 1:nrow(liver_t2)){
  temp<-liver_t2[i,] - liver_t1_mean
  diff1[i,]<-temp
}

# control - treatment 3
# create empty DF from column names from original
diff2<-liver_t3[FALSE,]
for (i in 1:nrow(liver_t3)){
  temp<-liver_t3[i,] - liver_t1_mean
  diff2[i,]<-temp
}

# control - treatment 4
# create empty DF from column names from original
diff3<-liver_t4[FALSE,]
for (i in 1:nrow(liver_t4)){
  temp<-liver_t4[i,] - liver_t1_mean
  diff3[i,]<-temp
}

# plot difference histograms
par(mfrow=c(3,1))
par(mar = c(4, 4, 1, 4))
hist(colMeans(abs(diff1)), xlab="Difference (Pesticide - Control)", main="", col="red")
hist(colMeans(abs(diff2)), xlab="Difference (Bullfrog - Control)", main="", col="blue")
hist(colMeans(abs(diff3)), xlab="Difference (Pesticide & Bullfrog - Control)", main="", col="purple")

#### subset columns based on column means abs(difference) > X ####

###
## pesticide only
# take absolute value of difference first
test1<-as.data.frame(colMeans(abs(diff1)))
names(test1)<-"peak_value"
#test1$V2<-abs(test1$peak_value)
# create column for peak no
peak_no<-seq(1,1584, 1)
test1_no<-cbind(test1, peak_no)
# drop peaks
threshold<-0.90
test1_f<-filter(test1_no, test1 > threshold)
# check to see ~20% dropped
dim(test1_f)[1] / dim(test1)[1]

###
## bullfrog only
test2<-as.data.frame(colMeans(abs(diff2)))
names(test2)<-"peak_value"
#test1$V2<-abs(test1$peak_value)
# create column for peak no
peak_no<-seq(1,1584, 1)
test2_no<-cbind(test2, peak_no)
# drop peaks
threshold<-0.56
test2_f<-filter(test2_no, test2 > threshold)
# check to see ~20% dropped
dim(test2_f)[1] / dim(test2)[1]

###
## pesticide and bullfrog
test3<-as.data.frame(colMeans(abs(diff2)))
names(test3)<-"peak_value"
#test1$V2<-abs(test1$peak_value)
# create column for peak no
peak_no<-seq(1,1584, 1)
test3_no<-cbind(test3, peak_no)
# drop peaks
threshold<-0.56
test3_f<-filter(test3_no, test3 > threshold)
# check to see ~20% dropped
dim(test3_f)[1] / dim(test3)[1]


###################################
#### import in SVM-RFE ranked peaks ####
rfe_pest<-read.csv(paste0(datadir, "output/RFE_ctrlvspest.csv"))
rfe_frog<-read.csv(paste0(datadir,"output/RFE_ctrlvbullfrog.csv"))
rfe_combined<-read.csv(paste0(datadir,"output/RFE_ctrlvpest_and_bullfrog_allbins.csv"))

#### import ranked and identified list for each treatment ####
rfeid_pest<-read.csv(paste0(datadir, "output/pest_metabolite_id_all_treatments_combined_wmh_final_110119.csv"))
rfeid_combined<-read.csv(paste0(datadir,"output/combined_metabolite_id_all_treatments_combined_wmh_final_110119.csv"))
rfeid_frog<-read.csv(paste0(datadir,"output/bullfrog_metabolite_id_all_treatments_combined_wmh_final_110119.csv"))

########################################
#### format SVM-RFE ranked peaks ####

## PEST
# drop first 2 rows
rfe_pest<-rfe_pest[-1:-2,]
# drop all columns except the peak time and the rank
rfe_pest<-rfe_pest[,1:2]
# add peak_no column for join
peak_no<-seq(1,1585, 1)
rfe_pest2<-cbind(rfe_pest, peak_no)

## FROG
# drop first 2 rows
rfe_frog<-rfe_frog[-1:-2,]
# drop all columns except the peak time and the rank
rfe_frog<-rfe_frog[,1:2]
# add peak_no column for join
peak_no<-seq(1,1585, 1)
rfe_frog2<-cbind(rfe_frog, peak_no)

## COMBINED
# drop first 2 rows
rfe_combined<-rfe_combined[-1:-2,]
# drop all columns except the peak time and the rank
rfe_combined<-rfe_combined[,1:2]
# add peak_no column for join
peak_no<-seq(1,1585, 1)
rfe_combined2<-cbind(rfe_combined, peak_no)

#################################################
#### merge with mean filtered peak list ####
# pest #
pest.m<-merge(test1_f, rfe_pest2, by.x="peak_no", by.y="peak_no")
pest.f<-filter(pest.m, featureRankedList4join2 < 100)
dim(pest.f) # 73 of 100 ranked peaks remaining
# frog #
frog.m<-merge(test2_f, rfe_frog2, by.x="peak_no", by.y="peak_no")
frog.f<-filter(frog.m, featureRankedList4join2 < 100)
dim(frog.f) # 74 of 100 ranked peaks remaining
# combined #
combined.m<-merge(test3_f, rfe_combined2, by.x="peak_no", by.y="peak_no")
combined.f<-filter(combined.m, featureRankedList4join2 < 100)
dim(combined.f) # 75 of 100 ranked peaks remaining

#### merge with identified peak list from Matthew ####
# pest #
rfeid_pest_sub<-merge(rfeid_pest, pest.f,by.x="RT",by.y="X")
View(rfeid_pest_sub)
# frog # 
rfeid_frog_sub<-merge(rfeid_frog, frog.f,by.x="RT",by.y="X")
View(rfeid_frog_sub)
# combined #
rfeid_combined_sub<-merge(rfeid_combined, combined.f,by.x="RT",by.y="X")
View(rfeid_combined_sub)

# write out combined CSV
write.csv(rfeid_pest_sub, paste0(datadir, "output/rfeid_pest_sub.csv"))
write.csv(rfeid_frog_sub, paste0(datadir, "output/rfeid_frog_sub.csv"))
write.csv(rfeid_combined_sub, paste0(datadir, "output/rfeid_combined_sub.csv"))


##########################################################################
######### get sum of differences #############################

## control minus treatment 2

array1 <- array(data=NA, dim=c(10,1584,10))

for (i in 1:nrow(liver_t1)){
  for (j in 1:nrow(liver_t2)){
  temp<-liver_t1[i,] - liver_t2[j,]
  a<- abs(temp)
  array1[i,,j]<-as.matrix(a)
}}

#z <- array(1:20, dim = c(3,4,5))
#dim(z)
#apply(z, MARGIN = c(3), sum)

# add up 1st and 3rd dimensions
sum_out1<-apply(array1, MARGIN=c(2), sum)

## control minus treatment 3
array2 <- array(data=NA, dim=c(10,1584,10))

for (i in 1:nrow(liver_t1)){
  for (j in 1:nrow(liver_t3)){
    temp<-liver_t1[i,] - liver_t3[j,]
    a<- abs(temp)
    array2[i,,j]<-as.matrix(a)
  }}
# add up 1st and 3rd dimensions
sum_out2<-apply(array2, MARGIN=c(2), sum)

## control minus treatment 4
array3 <- array(data=NA, dim=c(10,1584,10))

for (i in 1:nrow(liver_t1)){
  for (j in 1:nrow(liver_t4)){
    temp<-liver_t1[i,] - liver_t4[j,]
    a<- abs(temp)
    array3[i,,j]<-as.matrix(a)
  }}
# add up 1st and 3rd dimensions
sum_out3<-apply(array3, MARGIN=c(2), sum)

# combine treatment comparison sums
sum_out_all<-cbind.data.frame(sum_out1, sum_out2, sum_out3)

# plot histograms of total sums for each trtmt comparison
par(mfrow=c(3,1))
par(mar = c(4, 4, 1, 4))
hist(sum_out_all$sum_out1, xlab="Sum Difference (Control - Pesticide)", main="",col="red", xlim=c(0,150))
hist(sum_out_all$sum_out2, xlab="Sum Difference (Control - Bullfrog)", main="",col="blue", xlim=c(0,150))
hist(sum_out_all$sum_out3, xlab="Sum Difference (Control - Pesticide & Bullfrog)", main="", col="purple", xlim=c(0,150))

##########################################################################
#### subset based on threshold for each trtmt control comparison ####

## control vs pest
# create column for peak no
peak_no<-seq(1,1584, 1)
sumout1_no<-cbind.data.frame(sum_out1, peak_no)
# drop peaks
threshold<-113
sumout1_f<-filter(sumout1_no, sum_out1 > threshold)
# check to see ~20% dropped
dim(sumout1_f)[1] / dim(sumout1_no)[1]

## control vs frog
# create column for peak no
peak_no<-seq(1,1584, 1)
sumout2_no<-cbind.data.frame(sum_out2, peak_no)
# drop peaks
threshold<-82
sumout2_f<-filter(sumout2_no, sum_out2 > threshold)
# check to see ~20% dropped
dim(sumout2_f)[1] / dim(sumout2_no)[1]

## control vs pest plus frog
# create column for peak no
peak_no<-seq(1,1584, 1)
sumout3_no<-cbind.data.frame(sum_out3, peak_no)
# drop peaks
threshold<-76
sumout3_f<-filter(sumout3_no, sum_out3 > threshold)
# check to see ~20% dropped
dim(sumout3_f)[1] / dim(sumout3_no)[1]

#################################################
#### merge sum filtered peak list with top 100 ranked peaks from RFE- SVM ####

# pest #
pest.m2<-merge(sumout1_f, rfe_pest2, by.x="peak_no", by.y="peak_no")
pest.f2<-filter(pest.m2, featureRankedList4join2 < 100)
dim(pest.f2) # 73 of 100 ranked peaks remaining
# frog #
frog.m2<-merge(sumout1_f, rfe_frog2, by.x="peak_no", by.y="peak_no")
frog.f2<-filter(frog.m2, featureRankedList4join2 < 100)
dim(frog.f2) # 74 of 100 ranked peaks remaining
# combined #
combined.m2<-merge(sumout1_f, rfe_combined2, by.x="peak_no", by.y="peak_no")
combined.f2<-filter(combined.m2, featureRankedList4join2 < 100)
dim(combined.f2) # 75 of 100 ranked peaks remaining

#### merge with identified peak list from Matthew ####
# pest #
rfeid_pest_sub2<-merge(rfeid_pest, pest.f2,by.x="RT",by.y="X")
View(rfeid_pest_sub2)
# frog # 
rfeid_frog_sub2<-merge(rfeid_frog, frog.f2,by.x="RT",by.y="X")
View(rfeid_frog_sub2)
# combined #
rfeid_combined_sub2<-merge(rfeid_combined, combined.f2,by.x="RT",by.y="X")
View(rfeid_combined_sub2)

# write out combined CSV
write.csv(rfeid_pest_sub2, paste0(datadir, "output/rfeid_pest_sub2.csv"))
write.csv(rfeid_frog_sub2, paste0(datadir, "output/rfeid_frog_sub2.csv"))
write.csv(rfeid_combined_sub2, paste0(datadir, "output/rfeid_combined_sub2.csv"))


#############  NOT DONE  ##################################################3



