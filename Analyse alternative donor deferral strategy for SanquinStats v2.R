#######################################################
# Analyse impact of alternative deferral strategy
#######################################################
# initialize R
rm(list=ls())
while (!is.null(dev.list())) dev.off()
cat("\014")  

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
changeIDs<-F # to be set by the USER 

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
data<-data[data$Hb>5,]
nrrecs<-c(nrrecs,nrow(data))
data<-data[data$Hb<15,]
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
malefits  <-fitHbdistributions(adata[adata$Sex=="M",],nrofquantiles=20)
femalefits<-fitHbdistributions(adata[adata$Sex=="F",],nrofquantiles=20)

# stop execution of groupsize is larger than required by the user
if(malefits$minsubset<mingroupsize | malefits$minsubset<mingroupsize) {
  print("Code stopped because aggregated group size is smaller than specified by the user")
  print("Please decrease the nrofquantiles parameter in the fitHbdistributions functions (line 144/145)")
  print("or increase the mingroupsize (line 44)")
  stop("Change analysis code")
}

# set parameter for maximum follow-up
maxDons<-max(nHb$x)

# save distribution fits and maxDons
tosave<-append(tosave, list(malefits=malefits))
tosave<-append(tosave, list(femalefits=femalefits))
tosave<-append(tosave, list(maxDons=maxDons))

# convert Hb levels if so required
if (!Hb_in_gpl){
  dtm<-dtm/0.06206 -1e-6 # subtract a small margin to compensate for rounding errors
  dtf<-dtf/0.06206 -1e-6
  data$Hb<-data$Hb/0.06206
}

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
plot(rownames(defm), proportions(defm, margin=1)[,2], ylim=c(0,maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(deff), proportions(deff, margin=1)[,2], col="red")
points(rownames(defm), proportions(defm, margin=1)[,2], pch=1, col="blue")
points(rownames(deff), proportions(deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))

# plot distribution of number of donations per sex
with(data[data$Sex=="F",], plot(as.numeric(table(numdons)), type="l", col="red", xlim=c(1,maxDons),# log='y',
                                xlab="Number of donations", ylab="number of donors"))
with(data[data$Sex=="M",], lines(as.numeric(table(numdons)), type="l", col="blue"))
legend("topright", c("Males", "Females"), lty=c(1,1), col=c("blue","red"))

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

# set plot parameters
intervaltoshow<-c(100, 730) # time interval to show
ylim<-c(-40, 40)           # interval in dHb to show
nrtoprint<-15000           # nr of observations to select for printing
rollmeanWidth<-1000        # width of the rolling window

# plot association between time and Hb change for females
set.seed(1)
self<-which(data$Sex=="F" & data$dt>=intervaltoshow[1] & !is.na(data$ldt))
sel<-sample(self, nrtoprint, replace=F)
sel<-sel[order(data$dt[sel])]
with(data[sel,], plot(dt, dHb, col="red", log="x", xlim=intervaltoshow, ylim=ylim,
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
with(data[sel,], plot(dt, dHb, col="blue", log="x", xlim=intervaltoshow, ylim=ylim,
                      xlab="Days between donations", ylab="Change in Hb level [g/L]"))
abline(h=0,col=8)
# add rolling mean
with(data[sel,], lines(dt, rollmean(dHb, rollmeanWidth, fill = list(NA, NULL, NA)),
                       col = 3, lwd = 3 ))
linfitm<-lm(dHb~ldt, data=data[selm,])
summary(linfitm)
lines(intervaltoshow, predict(linfitm,cx), lwd=2)
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
(malesd<-sd(linfitm$residuals))   
(femalesd<-sd(linfitf$residuals)) 

# Calculate (un)acceptable deviation per sex
datt$d<-qnorm(cutoffperc)*femalesd
datt$d[datt$Sex=="M"]<-qnorm(cutoffperc)*malesd
# set thresholds per gender
datt$th<-dtf
datt$th[datt$Sex=="M"]<-dtm
table(datt$th-datt$d, datt$Sex)
# 95% Thresholds: 109.619737077851 (F) 119.111266716166 (M)
# 99% Thresholds: 103.247400641088 (F) 112.528261305199 (M)

# attach the datt file
while ("datt" %in% search()) detach(datt)
attach(datt)

###################################################################
# now analyse what the new donor deferral policy would achieve 
# in terms of number of donors now allowed to donate
###################################################################
# define an output array with 8 columns, one row per subsequent donation
# the items stored in each row is explained below
outputsummarytable<-as.data.frame(matrix(0,maxDons,8))
colnames(outputsummarytable)<-c("Deferred", "Non-deferred", "ShouldNotDeferred", "ShouldNotDonate", 
        "Should_not_have_donated","Missed_by_stopped_donor", "ShouldNotDeferred2", "RequiresReview")

stopped<-rep(F, length(Hb1)) # indicator for whether a donors has stopped or not
                             # is set when the mean Hb level was demonstrably below 
                             # the eligibility threshold at previous donation 
stopped2<-rep(F, length(Hb1)) # indicator for whether a donors has stopped or not
                              # is set when the mean Hb level is below the 
                              # eligibility threshold at previous donation
stopafter<-3 # stop donating after significant evidence only after stopafter donations have been made

for (i in 1:maxDons ){
  
  # 1 - Deferred donors
  # The number of donors deferred at step i are those with a Hb value that is  
  # below the deferral threshold
  eval(parse(text=paste0("outputsummarytable[",i,", 1]<-sum(!is.na(Hb",i,") & Hb",i,"<th)")))
  
  # 2 - non-Deferred donors
  # The number of donors not deferred at step i are those with a Hb value that   
  # is equal or larger than the deferral threshold
  eval(parse(text=paste0("outputsummarytable[",i,", 2]<-sum(!is.na(Hb",i,") & Hb",i,">=th)")))
  
  # 3 - Donors that should not have been deferred as the Hb deviation relative to
  #     their mean Hb value does not provide sufficient evidence against donation
  # only count events from second donation onwards
  if(i>1) eval(parse(text=paste0("outputsummarytable[",i,", 3]<-sum(!is.na(Hb",i,") & Hb",i,"<th & Hb",i,">=MeanHb",i-1,"-d)")))

  # 4 - Donors that should not donate as their Hb is demonstrably below the eligibility threshold
  eval(parse(text=paste0("outputsummarytable[",i,",4]<-sum(!is.na(Hb",i,") & MeanHb",i,"<th-d/sqrt(",i,"))")))
  
  # 5 - Donors that should not have donated
  if(i>1) eval(parse(text=paste0("outputsummarytable[",i,",5]<-sum(!is.na(Hb",i,") & MeanHb",i-1,"<th-d/sqrt(",i-1,") & Hb",i,">=th)")))
  
  # set new index for (previously) stopped donors
  # index stopped indicates that the mean Hb level was demonstrably below the eligibility threshold at previous donation 
  if (i>stopafter) eval(parse(text=paste0("stopped <-stopped  | (!is.na(Hb",i,") & MeanHb",i-1,"<th-d/sqrt(",i-1,"))")))
  # index stopped2 indicates that the mean Hb level is below the eligibility threshold at previous donation
  if (i>stopafter) eval(parse(text=paste0("stopped2<-stopped2 | (!is.na(Hb",i,") & MeanHb",i-1,"<th)")))

  # 6 - donations missed as a result of new deferral rule
  eval(parse(text=paste0("outputsummarytable[",i,",6]<-sum(!is.na(Hb",i,") & Hb",i,">=th & stopped)")))
  
  # 7 - Donors that should not have been deferred as the Hb deviation relative to 
  #     the absolute Hb threshold is insufficient (see also evaluation 3 above)
  # this basically presumes that anyone with a Hb level over th-d may donate
  eval(parse(text=        paste0("outputsummarytable[",i,", 7]<-sum(!is.na(Hb",i,") & Hb",i,"<th & Hb",i,">=th-d)")))
  
  # 8 - Identified as outlier, but not deferred
  if(i>1) eval(parse(text=paste0("outputsummarytable[",i,",8]<-sum(!is.na(Hb",i,") & Hb",i,"< MeanHb",i-1,"-d & Hb",i,">=th)")))
}
sum(stopped) 
sum(stopped2)

sum(stopped)/length(stopped) # proportion of stopped donors
sum(stopped2)/length(stopped) # proportion of stopped2 donors

outputsummarytable$defprop<-outputsummarytable$Deferred/(outputsummarytable$Deferred+outputsummarytable$`Non-deferred`)
outputsummarytable$nondefprop<-outputsummarytable$ShouldNotDeferred/outputsummarytable$Deferred
outputsummarytable
tosave<-append(tosave, list(outputsummarytable=outputsummarytable))

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
# for internal use only
# note that the number of 

maxplots<-4  # USER: Set the maximum number of graphs to plot in row/column of a matrix
# plotdonorprofile(KeyID[250], leg=T, ylim=c(0,190)) # nice illustration with a range of 10 unnecessary deferrals 
# plotdonorprofile(KeyID[250], ylim=c(75,210)) # nice illustration with a range of 10 unnecessary deferrals 

# parameter to determine whether the plots go to a pdf file or to screen.
plot_to_pdf<-F # to be set by the USER 

###########################
# Create a selection of donors that should have been deferred at donation i but did donate 
# at least i times
def<-5 # to be set by the USER 
i<-7   # to be set by the USER 
# Nr of donors that fit the criterion
eval(parse(text=paste0("sum(MeanHb",def,"+d/sqrt(",def,")<th & Hb",def,">th & !is.na(Hb",i,"))")))
if (eval(parse(text=paste0("sum(MeanHb",def,"+d/sqrt(",def,")<th & Hb",def,">th & !is.na(Hb",i,"))")))>0) {
  # set selection of donors
  eval(parse(text=paste0("selID<-KeyID[MeanHb",def,"+d/sqrt(",def,")<th & Hb",def,">th & !is.na(Hb",i,")]")))
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=paste0("Deferral at donation ",def,", minimum of ",i," donations.pdf"))
  plotmatrix(selID, 2, ylim=c(70,160))
  if(plot_to_pdf) dev.off()
}
###########################
# select donors with ndef deferrals at donation i and an average Hb level above the threshold
i<-30     # to be set by the USER 
ndef<-10  # to be set by the USER 
# Nr of donors selected
eval(parse(text=paste0("sum( nHb",i,"-HbOk",i,">",ndef-1," & MeanHb",i,">th )")))
if(eval(parse(text=paste0("sum( nHb",i,"-HbOk",i,">",ndef-1," & MeanHb",i,">th )")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-KeyID[nHb",i,"-HbOk",i,">",ndef-1," & MeanHb",i,">th]")))
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=paste0("Deferral at ",i,", but with an average Hb level above the threshold.pdf"))
  plotmatrix(selID, 2)
  if(plot_to_pdf) dev.off()
}
###########################
# Create a selection of donors that were deferred at donation i but have 
# an average Hb level well above (a value "delta") the deferral threshold
i<-30    # to be set by the USER 
delta<-5 # to be set by the USER 
# Nr of donors selected
eval(parse(text=paste0("sum(MeanHb",i,">th+delta & Hb",i,"<th & !is.na(Hb",i,"))")))
if (eval(parse(text=paste0("sum(MeanHb",i,">th+delta & Hb",i,"<th & !is.na(Hb",i,"))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-KeyID[MeanHb",i,">th+delta & Hb",i,"<th & !is.na(Hb",i,")]")))
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=paste0("Deferral at ",i,", but with an average of ",delta," above the threshold.pdf"))
  plotmatrix(selID, 2)
  if(plot_to_pdf) dev.off()
  # plot another subset
  plotmatrix(selID, 2, seedvalue=2)
}
