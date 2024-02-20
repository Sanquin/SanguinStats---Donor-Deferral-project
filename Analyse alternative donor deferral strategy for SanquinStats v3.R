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
maincodedatestamp<-"20240209"

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
library("dplyr")     # data wrangling

# restore old color palette
palette("R3")

# load support functions from the second code file
source("General_functions.R")

#############################
# read data
#############################
# but first write some of the input data to the tosave variable
# the first item is the name of the datafile used for the analyses
folder <- paste0("cutoff_", cutoffperc, "/")
dir.create(folder)
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
numdons<-table(data$numdons, data$Sex) 
numdons

#convert Hb levels if so required
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

# calculate average Hb levels
adata_mean_m<-with(adata[adata$Sex=="M",],c(sum(Nrdon), sum(Hb*Nrdon)/sum(Nrdon), mean(Hb)))
adata_mean_m<-c(adata_mean_m, sum(data$Sex=="M" & !is.na(data$Hb)), mean(data$Hb[data$Sex=="M" & !is.na(data$Hb)]))
adata_mean_f<-with(adata[adata$Sex=="F",],c(sum(Nrdon), sum(Hb*Nrdon)/sum(Nrdon), mean(Hb)))
adata_mean_f<-c(adata_mean_f, sum(data$Sex=="F" & !is.na(data$Hb)), mean(data$Hb[data$Sex=="F" & !is.na(data$Hb)]))

# calculate and plot distributions for various subsets of nr of donations
if(plot_to_pdf) pdf(file=paste0(folder, "Hb_distribution_male.pdf"), paper = "a4", height = 20, width =20)
malefits_hb<-FitDistributions(adata[adata$Sex=="M",],"Hb")
if(plot_to_pdf) dev.off()

if(plot_to_pdf) pdf(file=paste0(folder, "Hb_distribution_female.pdf"),paper = "a4", height = 20, width =20)
femalefits_hb<-FitDistributions(adata[adata$Sex=="F",],"Hb")
if(plot_to_pdf) dev.off()

# set parameter for maximum follow-up
maxDons<-max(nHb$x)

# save distribution fits and maxDons
tosave<-append(tosave, list(numdons=numdons))
tosave<-append(tosave, list(malefits_hb=malefits_hb))
tosave<-append(tosave, list(malemeanHb=adata_mean_m))
tosave<-append(tosave, list(femalefits_hb=femalefits_hb))
tosave<-append(tosave, list(femalemeanHb=adata_mean_f))
tosave<-append(tosave, list(maxDons=maxDons))

# calculate and plot distributions for various subsets of donation intervals

#create donation interval variable
bdata <- data %>% group_by(KeyID)%>% mutate(interval = as.numeric(DonDate - lag(DonDate)), meanInt = mean(interval, na.rm=T), sdInt = sd(interval, na.rm=T), Nrdon = max(numdons)) %>% ungroup() %>% select(Sex, meanInt, sdInt, Nrdon, KeyID) %>% distinct(KeyID, .keep_all = T) %>% select(-KeyID)
colnames(bdata) <- c("Sex", "interval", "sd", "Nrdon")
bdata <- bdata[!is.na(bdata$interval),]

# calculate average donation levels
bdata_mean_m<-with(bdata[bdata$Sex=="M",],c(sum(Nrdon), sum(interval*Nrdon)/sum(Nrdon), mean(interval)))
bdata_mean_m<-c(bdata_mean_m, sum(data$Sex=="M"), sum(data$Sex=="M" & data$numdons>1))
bdata_mean_f<-with(bdata[bdata$Sex=="F",],c(sum(Nrdon), sum(interval*Nrdon)/sum(Nrdon), mean(interval)))
bdata_mean_f<-c(bdata_mean_f, sum(data$Sex=="F"), sum(data$Sex=="F" & data$numdons>1))

# plot distributions
if(plot_to_pdf) pdf(file=paste0(folder, "DonInterval_distribution_male.pdf"),paper = "a4", height = 20, width =20)
malefits_int <-FitDistributions(bdata[bdata$Sex=="M",],"interval")
if(plot_to_pdf) dev.off()
if(plot_to_pdf) pdf(file=paste0(folder,"DonInterval_distribution_female.pdf"),paper = "a4", height = 20, width =20)
femalefits_int<-FitDistributions(bdata[bdata$Sex=="F",],"interval")
if(plot_to_pdf) dev.off()

# save distribution fits and maxDons
tosave<-append(tosave, list(malefits_int=malefits_int))
tosave<-append(tosave, list(malemeanint=bdata_mean_m))
tosave<-append(tosave, list(femalefits_int=femalefits_int))
tosave<-append(tosave, list(femalemeanint=bdata_mean_f))

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
if(plot_to_pdf) pdf(file=paste0(folder,"Deferrals_per_year.pdf"))
plot(rownames(defm), proportions(defm, margin=1)[,2], ylim=c(0,maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(deff), proportions(deff, margin=1)[,2], col="red")
points(rownames(defm), proportions(defm, margin=1)[,2], pch=1, col="blue")
points(rownames(deff), proportions(deff, margin=1)[,2], pch=2, col="red")
legend("topright", c("Males", "Females"), pch=c(1,2), lty=c(1,1), col=c("blue","red"))
if(plot_to_pdf) dev.off()

if(plot_to_pdf) pdf(file=paste0(folder,"Distribution_of_donation_counts.pdf"))

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
if(plot_to_pdf) pdf(file=paste0(folder,"Change_in_Hb_level_by_interval_female.pdf"))
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
if(plot_to_pdf) pdf(file=paste0(folder,"Change_in_Hb_level_by_interval_male.pdf"))
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

data <- data %>%
  group_by(KeyID) %>%
  mutate(cum_Hb= cumsum(Hb)) %>% ungroup() %>% mutate(meanHb = cum_Hb/numdons)


#make a column for the previous Mean Hb
idx<-1:nrow(data)
precursor<-idx-1
precursor[1]<-nrow(data)

data$prevMeanHb<-NA
data$prevMeanHb<-data$meanHb[precursor]
data$prevMeanHb[data$KeyID != data$KeyID[precursor]]<-NA

data$prevMeanHb[is.na(data$prevMeanHb)] <- data$meanHb[is.na(data$prevMeanHb)]

############################
# Analyse deferrals
############################

# standard deviation of individual measurements
(malesd<-sd(linfitm$residuals)/sqrt(2))
(femalesd<-sd(linfitf$residuals)/sqrt(2))

# Calculate (un)acceptable deviation per sex
data$d<-qnorm(cutoffperc)*femalesd
data$d[data$Sex=="M"]<-qnorm(cutoffperc)*malesd
# set thresholds per gender
data$th<-dtf
data$th[data$Sex=="M"]<-dtm
table(data$th-data$d, data$Sex)

###################################################################
# now analyse what the new donor deferral policy would achieve 
# in terms of number of donors now allowed to donate
###################################################################

analysisresults<-AnalysePolicyImpact(dataframe=data)
outputsummarytable<-analysisresults$outputsummarytable
tosave<-append(tosave, list(analysisresults=analysisresults))


dataf<-data[data$Sex=="F",]
analysisresults_f<-AnalysePolicyImpact(dataframe=dataf)
tosave<-append(tosave, list(analysisresults_f=analysisresults_f))


datam<-data[data$Sex=="M",]
analysisresults_m<-AnalysePolicyImpact(dataframe=datam)
tosave<-append(tosave, list(analysisresults_m=analysisresults_m))

###################################################################
# Write tosave data to datafile
###################################################################
saveRDS(tosave, paste0(folder, "SavedDeferralData_",cutoffperc,"_",Sys.Date(),".RDS"))

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
data <- data %>% group_by(KeyID) %>% mutate(totaldon = n())

maxplots<-3  # USER: Set the maximum number of graphs to plot in row/column of a matrix
# plotdonorprofile(KeyID[250], leg=T, ylim=c(0,190)) # nice illustration with a range of 10 unnecessary deferrals 
# plotdonorprofile(KeyID[250], ylim=c(75,210)) # nice illustration with a range of 10 unnecessary deferrals 

###########################
# Create a selection of donors that should have been deferred at donation 'def' but did donate 
# at least n times
def<-6 # to be set by the USER 
n<-9  # to be set by the USER 
# Nr of donors that fit the criterion
(totnr<-eval(parse(text=paste0("sum(data$numdons==",def, "& data$meanHb + data$d/sqrt(",def,")<data$th & data$Hb > data$th & !is.na(data$Hb) & data$totaldon >",n,", na.rm=T)"))))
if (eval(parse(text=paste0("sum(data$numdons==",def, "& data$meanHb + data$d/sqrt(",def,")<data$th & data$Hb > data$th & !is.na(data$Hb)  & data$totaldon >",n,", na.rm=T)")))>0) {
  # set selection of donors
  eval(parse(text=paste0("selID<-data$KeyID[data$numdons==",def, "& data$meanHb + data$d/sqrt(",def,")<data$th & data$Hb > data$th & !is.na(data$Hb)  & data$totaldon >",n,"]")))
  print(selID[!is.na(selID)])
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at donation ",def," minimum of ",n," donations; ",totnr," in total.pdf")))
  plotmatrix(selID,maxplots=maxplots, ylim=c(70,160),)
  if(plot_to_pdf) dev.off()
}
###########################
# select donors with at least ndef deferrals at donation n and an average Hb level above the deferral threshold
data <- data %>% group_by(KeyID) %>% mutate(def_count = cumsum(def))

n<-29     # to be set by the USER 
ndef<-8  # to be set by the USER 
# Nr of donors selected
(totnr<-eval(parse(text=paste0("sum(data$numdons ==",n," & data$def_count >= ndef & data$meanHb > data$th & !is.na(data$Hb))"))))
if(eval(parse(text=paste0("sum(data$numdons ==",n," & data$def_count >= ndef & data$meanHb > data$th & !is.na(data$Hb))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-data$KeyID[data$numdons ==",n," & data$def_count >= ndef & data$meanHb > data$th & !is.na(data$Hb)]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at ",n," but with an average Hb level above the threshold; ",totnr," in total.pdf")))
  plotmatrix(selID,maxplots=maxplots)
  if(plot_to_pdf) dev.off()
}
###########################
# Create a selection of donors that were deferred at donation n but have 
# an average Hb level well above (a value "delta") the deferral threshold
n<-30    # to be set by the USER 
delta<-5 # to be set by the USER 
# Nr of donors selected
(totnr<-eval(parse(text=paste0("sum(data$numdons ==",n,"& data$meanHb>data$th+delta & data$Hb<data$th & !is.na(data$Hb))"))))
if (eval(parse(text=paste0("sum(data$numdons ==",n,"& data$meanHb>data$th+delta & data$Hb<data$th & !is.na(data$Hb))")))>0){
  # set selection of donors
  eval(parse(text=paste0("selID<-data$KeyID[data$numdons ==",n,"& data$meanHb>data$th+delta & data$Hb<data$th & !is.na(data$Hb)]")))
  print(selID)
  # plot the donor profile in a matrix
  if(plot_to_pdf) pdf(file=gsub(" ","_",paste0("Deferral at ",n," but with an average of ",delta," above the threshold; ",totnr," in total.pdf")))
  plotmatrix(selID, maxplots=maxplots)
  if(plot_to_pdf) dev.off()
  # plot another subset
  # plotmatrix(selID, 2, seedvalue=2)
}
