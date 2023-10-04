#######################################################
# Analyse impact of alternative deferral strategy
#######################################################
# initialize R
rm(list=ls())
while (!is.null(dev.list())) dev.off()
cat("\014")  

# set the working path
library(this.path)
setwd(this.dir()) # set the active working directory to the directory of this file
getwd()

# Set code date 
maincodedatestamp<-"20231004"

#############################
# Define user input
#############################

# set name of the datafile to use
FileToUse<-"testdata.RDS" # to be set by the USER 

# "FileToUSe" should contain the following data columns (data type in brackets):
# KeyID   : Unique identifier for each donor (integer)
# Sex     : indicator for Male (M) or Female (F) donor (factor)
# DonDate : date of donation (Date)
# Hb      : donor Hb at donation date (numeric)

# All relevant input and output information is stored in the "tosave" variable
# This information is stored in a file called "SavedDeferralData_DD-MM-YYYY.RDS"

# parameter to determine whether the plots go to a pdf file or to screen.
plot_to_pdf<-F # to be set by the USER 

# Set minimum acceptable Hb levels for males and females
dtm<-8.4  # to be set by the USER 
dtf<-7.8  # to be set by the USER 

# Are Hb levels given in g/L (T) or in mmol/L (F)
Hb_in_gpl<-F # to be set by the USER 

# European threshold (g/L)
135*0.06206 # 135g/L for males = 8.3781 mmol/L
round(135*0.06206,1) # 8.4 mmol/L
125*0.06206 # 125g/L for females = 7.7575 mmol/L
round(125*0.06206,1) # 7.8 mmol/L

# Is a change in anonymous donor IDs required?
changeIDs<-T # to be set by the USER 

# Set deferral percentile level
cutoffperc<-0.99 # to be set by the USER 

# set minimum size of the groups for aggregated data
mingroupsize<-20 # to be set by the USER 

# Note that the analyses may take a substantial amount of time. Therefore an analysis file 
# ("donations_analysis_data.RDS") will be created, which allows speeding up the analyses
# by reading in this file whenever the analyses are repeated (e.g. for making additional graphs)
# Do not forget to remove this file in case the original data file is updated!

#############################
# Initialize R
#############################
# load various R libraries required for the analyses
library("zoo")       # required for rollmean function
library("lubridate") # required for extracting year info from a date variable
library("Hmisc")     # to enable calculating splitpoints
library("fitdistrplus") # to fit and plot normal distribution fits to the data
library("stringr")   # to manipulate strings

# restore old color palette
palette("R3")

# load support functions from the second code file
source("General_functions.R")

#############################
# read data
#############################
# but first write some of the input data to the tosave variable
# the first item is the name of the datafile used for the analyses
tosave<-list(FileToUse=FileToUse)

# save code datestamps
tosave<-append(tosave, list(maincodedatestamp=maincodedatestamp))
tosave<-append(tosave, list(generalfunctionscodedatestamp=generalfunctionscodedatestamp))

# save cutoff percentage applied
tosave<-append(tosave, list(cutoffperc=cutoffperc))

# now read in the data
data<-readRDS(FileToUse)
(classes<-sapply(data,class))
#     KeyID         Sex     DonDate          Hb 
# "integer" "character"      "Date"   "numeric" 

# when was the first donation
daterange<-min(data$DonDate)
(daterange<-c(daterange,max(data$DonDate)))
nrrecs<-nrow(data)

# save info
tosave<-append(tosave, list(daterange=daterange))

# remove all missing records
sum(is.na(data$Hb))
data<-data[!is.na(data$Hb),]
nrrecs<-c(nrrecs,nrow(data))

data<-data[!is.na(data$KeyID),]
nrrecs<-c(nrrecs,nrow(data))
data<-data[!is.na(data$DonDate),]
nrrecs<-c(nrrecs,nrow(data))
nrrecs<-c(nrrecs,length(unique(data$KeyID)))

# save changes in dataset
tosave<-append(tosave, list(nrrecs=nrrecs))

# Sort by date per donor
data<-data[order(data$KeyID,data$DonDate),] 
# Create index for nr of donations
data$numdons <- sequence(rle(data$KeyID)$lengths)


# nr of donors
length(unique(data$KeyID)) 
# nr of donations
length((data$KeyID))       

# change KeyIDs if requested by the user
if (changeIDs) {
  donIDs<-unique(data$KeyID)
  donIDs<-cbind(donIDs, 1:length(donIDs))
  #donIDs<-as.data.frame(cbind(donIDs, 1:length(donIDs)))
  colnames(donIDs)<-c("KeyID", "NewID")
  data<-merge(data,donIDs, by="KeyID")
  data$KeyID<-NULL
  colnames(data)[which(colnames(data)=="NewID")]<-"KeyID"
}

# calculate distribution of nr of donations per donor
table(data$numdons, data$Sex) 

# convert Hb levels if so required
if (!Hb_in_gpl){
  dtm<-dtm/0.06206 -1e-6 # subtract a small margin to compensate for rounding errors
  dtf<-dtf/0.06206 -1e-6
  data$Hb<-data$Hb/0.06206
}

# calculate sex, mean Hb, Sd and nr of donations per donor
meanHb<-aggregate(data$Hb, by=list(data$KeyID), mean)
SdHb<-aggregate(data$Hb, by=list(data$KeyID), sd)
nHb<-aggregate(data$Hb, by=list(data$KeyID, data$Sex), length)
# aggregate these results
adata<-merge(meanHb,SdHb, by="Group.1")
adata<-merge(adata,nHb, by="Group.1")
adata$Group.1<-NULL
colnames(adata)<-c("Hb", "sd", "Sex", "Nrdon")

# calculate distributions for various subsets of nr of donations

if(plot_to_pdf) pdf(file="Hb_distribution_male.pdf")
# nrofquantiles are to be set by the USER 
malefits  <-fitHbdistributions(adata[adata$Sex=="M",],nrofquantiles=20)
if(plot_to_pdf) dev.off()

if(plot_to_pdf) pdf(file="Hb_distribution_female.pdf")
# nrofquantiles are to be set by the USER 
femalefits<-fitHbdistributions(adata[adata$Sex=="F",],nrofquantiles=20)
if(plot_to_pdf) dev.off()

# stop execution of groupsize is larger than required by the user
if(malefits$minsubset<mingroupsize | femalefits$minsubset<mingroupsize) {
  print("Code stopped because aggregated group size is smaller than specified by the user")
  print("Please decrease the nrofquantiles parameter in the fitHbdistributions functions (line 167/172)")
  print("or increase the mingroupsize (line 56)")
  stop("Change analysis code")
}

# set parameter for maximum follow-up
maxDons<-max(nHb$x)

# save distribution fits and maxDons
tosave<-append(tosave, list(malefits=malefits))
tosave<-append(tosave, list(femalefits=femalefits))
tosave<-append(tosave, list(maxDons=maxDons))

# Set indicator for deferral
data$def<-ifelse(data$Hb<dtf,1,0)
data$def[data$Sex=="M"]<-ifelse(data$Hb[data$Sex=="M"]<dtm,1,0)

#create year variable
data$year<-year(data$DonDate)

# table deferrals per year
with(data,table(year, Sex))
with(data,table(Sex,year,def))
deff<-with(data[data$Sex=="F",],table(year,def))
proportions(deff, margin=1)
defm<-with(data[data$Sex=="M",],table(year,def))
proportions(defm, margin=1)
tosave<-append(tosave, list(defm=defm))
tosave<-append(tosave, list(deff=deff))

# plot deferrals per year
maxdef<-max(c(proportions(defm, margin=1)[,2],proportions(deff, margin=1)[,2]))
if(plot_to_pdf) pdf(file="Deferrals_per_year.pdf")
plot(rownames(defm), proportions(defm, margin=1)[,2], ylim=c(0,maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(deff), proportions(deff, margin=1)[,2], col="red")
points(rownames(defm), proportions(defm, margin=1)[,2], pch=1, col="blue")
points(rownames(deff), proportions(deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

if(plot_to_pdf) pdf(file="Distribution_of_donation_counts.pdf")
# plot distribution of number of donations per sex
with(data[data$Sex=="F",], plot(as.numeric(table(numdons)), type="l", col="red", xlim=c(1,maxDons), log='y',
                                xlab="Number of donations", ylab="number of donors"))
with(data[data$Sex=="M",], lines(as.numeric(table(numdons)), type="l", col="blue"))
with(data[data$Sex=="M",], points(as.numeric(table(numdons)), pch=1, col="blue"))
with(data[data$Sex=="F",], points(as.numeric(table(numdons)), pch=2, col="red"))
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

##################################
# Estimate measurement variation
##################################

# create an index for pointing to the previous donation
idx<-1:nrow(data)
precursor<-idx-1
precursor[1]<-nrow(data)

# calculate time since last donation
data$dt<-NA
data$dt<-as.numeric(data$DonDate-data$DonDate[precursor])
data$dt[data$KeyID != data$KeyID[precursor]]<-NA
# calculate log10 of the time difference
data$ldt<-log10(data$dt)

# calculate change in Hb from previous donation
data$dHb<-NA
data$dHb<-data$Hb-data$Hb[precursor]
data$dHb[data$KeyID != data$KeyID[precursor]]<-NA

# set plot parameters # to be set by the USER 
intervaltoshow<-c(100, 730) # time interval to show
ylim<-c(-40, 40)           # interval in dHb to show
nrtoprint<-15000           # nr of observations to select for printing
rollmeanWidth<-1000        # width of the rolling window

# plot association between time and Hb change for females
set.seed(1)
self<-which(data$Sex=="F" & data$dt>=intervaltoshow[1] & !is.na(data$ldt))
sel<-sample(self, nrtoprint, replace=F)
sel<-sel[order(data$dt[sel])]
if(plot_to_pdf) pdf(file="Change_in_Hb_level_by_interval_female.pdf")
with(data[sel,], plot(dt, jitter(dHb, amount=.5), col="red", log="x", xlim=intervaltoshow, ylim=ylim,
     xlab="Days between donations", ylab="Change in Hb level [g/L]"))
abline(h=0,col=8)
# add rolling mean
with(data[sel,], lines(dt, rollmean(dHb, rollmeanWidth, fill = list(NA, NULL, NA)),
  col = 3, lwd = 3 ))
linfitf<-lm(dHb~ldt, data=data[self,])
summary(linfitf)
cx<-as.data.frame(log10(intervaltoshow))
colnames(cx)<-"ldt"
lines(intervaltoshow, predict(linfitf,cx), lwd=2)
if(plot_to_pdf) dev.off()
mean(data$Hb[self])   
sd(data$Hb[self])     
sd(data$dHb[self])    
sd(linfitf$residuals) 

# save info
tosave<-append(tosave, list(coeff=linfitf$coefficients))
tosave<-append(tosave, list(sdf=c(mean(data$Hb[self]), sd(data$Hb[self]), 
                            sd(data$dHb[self]), sd(linfitf$residuals))))

# plot association between time and Hb change for males
intervaltoshow<-c(50, 730) # time interval to show
set.seed(1)
selm<-which(data$Sex=="M" & data$dt>=intervaltoshow[1] & !is.na(data$ldt))
sel<-sample(selm, nrtoprint, replace=F)
sel<-sel[order(data$dt[sel])]
if(plot_to_pdf) pdf(file="Change_in_Hb_level_by_interval_male.pdf")
with(data[sel,], plot(dt, jitter(dHb, amount=.5), col="blue", log="x", xlim=intervaltoshow, ylim=ylim,
                      xlab="Days between donations", ylab="Change in Hb level [g/L]"))
abline(h=0,col=8)
# add rolling mean
with(data[sel,], lines(dt, rollmean(dHb, rollmeanWidth, fill = list(NA, NULL, NA)),
                       col = 3, lwd = 3 ))
linfitm<-lm(dHb~ldt, data=data[selm,])
summary(linfitm)
lines(intervaltoshow, predict(linfitm,cx), lwd=2)

if(plot_to_pdf) dev.off()
mean(data$Hb[selm])   
sd(data$Hb[selm])     
sd(data$dHb[selm])    
sd(linfitm$residuals) 


# save info
tosave<-append(tosave, list(coefm=linfitm$coefficients))
tosave<-append(tosave, list(sdm=c(mean(data$Hb[selm]), sd(data$Hb[selm]), 
                                  sd(data$dHb[selm]), sd(linfitm$residuals))))

###########################################################
# Create and save analysis file if it doesn't exist yet
###########################################################
# the analysis file contains per donor
# 1) KeyID
# 2) Sex
# 3) per donor per visit the 
#     Hb level (Hbi), 
#     Mean Hb level (MeanHbi) of all past donations, 
#     Nr of visits (nHbi), and 
#     Nr of measurements that were above the threshold value (HbOki)
# the new dataset is called datt

if(!file.exists("donations_analysis_data.RDS")){ 

  # tabulate correct and incorrect donations
  data$HbOk<-1-data$def
  table(data$HbOk, data$Sex )
  data$Hbres<-1
  
  # create dataset with first donations
  dats<-data[data$numdons<=1,]
  dat1<-dats[ ,c("KeyID", "Sex", "HbOk", "Hb", "Hbres")]
  colnames(dat1)<-c("KeyID", "Sex", "HbOk1", "MeanHb1", "nHb1")
  dat1<-cbind(dat1,dat1$MeanHb1)
  colnames(dat1)<-c("KeyID", "Sex", "HbOk1", "MeanHb1", "nHb1", "Hb1")
  dat1<-dat1[,c("KeyID", "Sex", "Hb1", "MeanHb1", "nHb1", "HbOk1")]
  
  # create dataset with information on second donations
  dats<-data[data$numdons<=2,]
  dat2<-aggregate(dats$HbOk, by=list(dats$KeyID), sum)
  colnames(dat2)<-c("KeyID", "HbOk2")
  datm<-aggregate(dats$Hb, by=list(dats$KeyID), mean, na.rm=T)
  colnames(datm)<-c("KeyID", "MeanHb2")
  datn<-aggregate(dats$Hbres, by=list(dats$KeyID), sum)
  colnames(datn)<-c("KeyID", "nHb2")
  dat2<-merge(dat2, datm, by="KeyID")
  dat2<-merge(dat2, datn, by="KeyID")
  dats<-dats[dats$numdons==2, c("KeyID", "Hb")]
  colnames(dats)<-c("KeyID", "Hb2")
  dat2<-merge(dat2, dats, by="KeyID", all=T)
  dat2<-dat2[,c("KeyID", "Hb2", "MeanHb2", "nHb2", "HbOk2")]

  # merge these files
  datt<-merge(dat1,dat2, by="KeyID")
  rm(dat1,dat2)
  
  # Now repeat this last step for all subsequent donations
  for (i in 3:maxDons) {
    print(paste("Add data for donation nr",i))
    #dats<-data[data$numdons<=2,]
    eval(parse(text=paste0("dats<-data[data$numdons<=",i,",]")))
    #dat2<-aggregate(dats$HbOk, by=list(dats$KeyID), sum)
    dat<-aggregate(dats$HbOk, by=list(dats$KeyID), sum)
    #colnames(dat2)<-c("KeyID", "HbOk2")
    eval(parse(text=paste0("colnames(dat)<-c(\"KeyID\", \"HbOk",i,"\")")))
    #datm<-aggregate(dats$Hb, by=list(dats$KeyID), mean)
    eval(parse(text=paste0("datm<-aggregate(dats$Hb, by=list(dats$KeyID), mean)")))
    #colnames(datm)<-c("KeyID", "MeanHb2")
    eval(parse(text=paste0("colnames(datm)<-c(\"KeyID\", \"MeanHb",i,"\")")))
    #datn<-aggregate(dats$Hbres, by=list(dats$KeyID), sum)
    eval(parse(text=paste0("datn<-aggregate(dats$Hbres, by=list(dats$KeyID), sum)")))
    #colnames(datn)<-c("KeyID", "nHb2")
    eval(parse(text=paste0("colnames(datn)<-c(\"KeyID\", \"nHb",i,"\")")))
    #dat2<-merge(dat2, datm, by="KeyID")
    dat<-merge(dat, datm, by="KeyID")
    #dat2<-merge(dat2, datn, by="KeyID")
    dat<-merge(dat, datn, by="KeyID")
    #dats<-dats[dats$numdons==2, c("KeyID", "Hb")]
    eval(parse(text=paste0("dats<-dats[dats$numdons==",i,", c(\"KeyID\", \"Hb\")]")))
    #colnames(dats)<-c("KeyID", "Hb2")
    eval(parse(text=paste0("colnames(dats)<-c(\"KeyID\", \"Hb",i,"\")")))
    #dat2<-merge(dat2, dats, by="KeyID", all=T)
    dat<-merge(dat, dats, by="KeyID", all=T)
    #dat2<-dat2[,c("KeyID", "Hb2", "MeanHb2", "nHb2", "HbOk2")]
    eval(parse(text=paste0("dat<-dat[,c(\"KeyID\", \"Hb",i,"\", \"MeanHb",i,"\", \"nHb",i,"\", \"HbOk",i,"\")]")))
    
    # Now merge the data for this number of donations to the total dataset
    datt<-merge(datt,dat, by="KeyID")
  }

  # remove intermediate datasets
  rm(list=c("dat", "datm", "datn", "dats"))
  
  # save file if doesn't exist (wich was already checked)
  if (!file.exists("donations_analysis_data.RDS")) saveRDS(datt, file="donations_analysis_data.RDS")

} else {
  # if the analysis file already exists, open it
  datt<-readRDS("donations_analysis_data.RDS")
}

############################
# Analyse deferrals
############################

# standard deviation of individual measurements
(malesd<-sd(linfitm$residuals)/sqrt(2))
(femalesd<-sd(linfitf$residuals)/sqrt(2))

# Calculate (un)acceptable deviation per sex
datt$d<-qnorm(cutoffperc)*femalesd
datt$d[datt$Sex=="M"]<-qnorm(cutoffperc)*malesd
# set thresholds per gender
datt$th<-dtf
datt$th[datt$Sex=="M"]<-dtm
table(datt$th-datt$d, datt$Sex)

###################################################################
# now analyse what the new donor deferral policy would achieve 
# in terms of number of donors now allowed to donate
###################################################################

# detach potentially attached objects
while ("datt" %in% search()) detach(datt)
while ("dattf" %in% search()) detach(dattf)
while ("dattm" %in% search()) detach(dattm)

# attach the datt file to enable analysis of policy impact
attach(datt)
analysisresults<-AnalysePolicyImpact()
outputsummarytable<-analysisresults$outputsummarytable
tosave<-append(tosave, list(analysisresults=analysisresults))
while ("datt" %in% search()) detach(datt)

dattf<-datt[datt$Sex=="F",]
attach(dattf)
analysisresults_f<-AnalysePolicyImpact()
tosave<-append(tosave, list(analysisresults_f=analysisresults_f))
while ("dattf" %in% search()) detach(dattf)

dattm<-datt[datt$Sex=="M",]
attach(dattm)
analysisresults_m<-AnalysePolicyImpact()
tosave<-append(tosave, list(analysisresults_m=analysisresults_m))
while ("dattm" %in% search()) detach(dattm)

###################################################################
# Write tosave data to datafile
###################################################################
saveRDS(tosave, paste0("SavedDeferralData_",Sys.Date(),".RDS"))

###################################################################
# now calculate various statistics of the updated policy 
###################################################################

# calculate sums per column
sms<-colSums(outputsummarytable)
(totn<-sms[1]+sms[2]) # number of measurements included

sms2<-colSums(outputsummarytable[2:nrow(outputsummarytable),])
(totn2<-sms2[1]+sms2[2]) # number of measurements included
sms2[1]/totn2   # proportion deferred 
sms2[3]/totn2   # proportion that should not have been deferred 
(sms2[1]-sms2[3])/totn2   # proportion with deviant Hb values 
sms2[3]/sms2[1] # proportion unnecessary deferrals within deviation of mean
sms2[7]/totn2   # proportion that Should not have been deferred 
(sms2[1]-sms2[7])/totn2   # proportion that should have been deferred 
sms2[7]/sms2[1] # proportion unnecessary deferrals
sms2[5]         # number of donations 
sms2[5]/totn2   # proportion that should not have donated 
sms2[6]         # number missed by deferred donors
sms2[6]/totn2   # proportion missed by deferred donors 
sms2[8]/totn2   # Reviewed for low relative Hb 

# however, if you only want to start this rule after two donations
sms3<-colSums(outputsummarytable[3:nrow(outputsummarytable),])
(totn3<-sms3[1]+sms3[2])  # number of measurements included
sms3[1]/totn3   # proportion deferred 
sms3[3]/totn3   # proportion that should not have been deferred 
(sms3[1]-sms3[3])/totn3   # proportion with deviant Hb values 
sms3[3]/sms3[1] # proportion unnecessary deferrals within deviation of mean
sms3[7]/totn3   # proportion that should not have been deferred 
(sms3[1]-sms3[7])/totn3   # proportion that should have been deferred 
sms3[7]/sms3[1] # proportion unnecessary deferrals
sms3[5]         # number of donations 
sms3[5]/totn3   # proportion that should not have donated
sms3[6]         # number missed by deferred donors
sms3[6]/totn3   # proportion missed by deferred donors
sms3[8]/totn3   # Reviewed for low relative Hb 



#######################################################
# plot some individual donor profiles
#######################################################
# FOR INTERNAL USE ONLY

maxplots<-3  # USER: Set the maximum number of graphs to plot in row/column of a matrix
# plotdonorprofile(KeyID[250], leg=T, ylim=c(0,190)) # nice illustration with a range of 10 unnecessary deferrals 
# plotdonorprofile(KeyID[250], ylim=c(75,210)) # nice illustration with a range of 10 unnecessary deferrals 

# for plotting datt dataset needs to be attached
attach(datt)

###########################
# Create a selection of donors that should have been deferred at donation 'def' but did donate 
# at least n times
def<-5 # to be set by the USER 
n<-7   # to be set by the USER 
# Nr of donors that fit the criterion
eval(parse(text=paste0("sum(MeanHb",def,"+d/sqrt(",def,")<th & Hb",def,">th & !is.na(Hb",n,"))")))
if (eval(parse(text=paste0("sum(MeanHb",def,"+d/sqrt(",def,")<th & Hb",def,">th & !is.na(Hb",n,"))")))>0) {
  # set selection of donors
  eval(parse(text=paste0("selID<-KeyID[MeanHb",def,"+d/sqrt(",def,")<th & Hb",def,">th & !is.na(Hb",n,")]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at donation ",def," minimum of ",n," donations.pdf")))
  plotmatrix(selID,maxplots=maxplots, ylim=c(70,160),)
  if(plot_to_pdf) dev.off()
}
###########################
# select donors with at least ndef deferrals at donation n and an average Hb level above the deferral threshold
n<-29     # to be set by the USER 
ndef<-8  # to be set by the USER 
# Nr of donors selected
eval(parse(text=paste0("sum( nHb",n,"-HbOk",n,">",ndef-1," & MeanHb",n,">th & !is.na(Hb",n,"))")))
if(eval(parse(text=paste0("sum( nHb",n,"-HbOk",n,">",ndef-1," & MeanHb",n,">th & !is.na(Hb",n,"))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-KeyID[nHb",n,"-HbOk",n,">",ndef-1," & MeanHb",n,">th& !is.na(Hb",n,")]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at ",n,", but with an average Hb level above the threshold.pdf")))
  plotmatrix(selID,maxplots=maxplots)
  if(plot_to_pdf) dev.off()
}
###########################
# Create a selection of donors that were deferred at donation n but have 
# an average Hb level well above (a value "delta") the deferral threshold
n<-30    # to be set by the USER 
delta<-5 # to be set by the USER 
# Nr of donors selected
eval(parse(text=paste0("sum(MeanHb",n,">th+delta & Hb",n,"<th & !is.na(Hb",n,"))")))
if (eval(parse(text=paste0("sum(MeanHb",n,">th+delta & Hb",n,"<th & !is.na(Hb",n,"))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-KeyID[MeanHb",n,">th+delta & Hb",n,"<th & !is.na(Hb",n,")]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at ",n," but with an average of ",delta," above the threshold.pdf")))
  plotmatrix(selID, maxplots=maxplots)
  if(plot_to_pdf) dev.off()
  # plot another subset
  # plotmatrix(selID, 2, seedvalue=2)
}
