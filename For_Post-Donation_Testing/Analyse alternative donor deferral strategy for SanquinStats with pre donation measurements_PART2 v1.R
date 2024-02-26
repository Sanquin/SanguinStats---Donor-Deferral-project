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
maincodedatestamp<-"20240226"

#############################
# Define user input
#############################

# set name of the datafile to use
FileToUse<-"testdata.RDS" # to be set by the USER
PostDonationResultsFile <- "cutoff_0.99/SavedDeferralData_0.99_2024-02-25.RDS" #to be set by the USER after running the post donation analysis

# "FileToUSe" should contain the following data columns (data type in brackets):
# KeyID   : Unique identifier for each donor (integer)
# Sex     : indicator for Male (M) or Female (F) donor (factor)
# DonDate : date of donation (Date)
# Hb      : donor Hb at donation date (numeric)
# pre_Hb  : donor Hb at donation date (numeric) measured pre-donation

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
library("car")      # regression

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

data<-readRDS(FileToUse)
(classes<-sapply(data,class))
#     KeyID         Sex     DonDate          Hb    pre_Hb
# "integer" "character"      "Date"   "numeric"   "numeric"

PostDonationResults <- readRDS(PostDonationResultsFile)

tosave<-list(FileToUse=FileToUse)
tosave <- append(tosave, list(PostDonationResultsFile=PostDonationResultsFile))
# save code datestamps
tosave<-append(tosave, list(maincodedatestamp=maincodedatestamp))
tosave<-append(tosave, list(generalfunctionscodedatestamp=generalfunctionscodedatestamp))

# save cutoff percentage applied
tosave<-append(tosave, list(cutoffperc=cutoffperc))


# when was the first donation
daterange<-min(data$DonDate)
(daterange<-c(daterange,max(data$DonDate)))
nrrecs<-nrow(data)

# save info
tosave<-append(tosave, list(daterange=daterange))

# remove all missing records
data<-data[!is.na(data$KeyID),]
nrrecs<-c(nrrecs,nrow(data))
data<-data[!is.na(data$DonDate),]
nrrecs<-c(nrrecs,nrow(data))
nrrecs<-c(nrrecs,length(unique(data$KeyID)))

sum(is.na(data$Hb))
data_complete <- data #Amber: hier is/was een probleem, namelijk donaties met een pre-donatiescreening onder de treshold hebben een NA Hb waarde (omdat er na de predonatie screening geen donatie was), maar we hebben in de verdere analyse wel beide nodig
data_prescreening <- data[!is.na(data$pre_Hb),] #Amber: make a separate dataset with all records where a pre-screening was conducted
data<-data[!is.na(data$Hb),] 
nrrecs<-c(nrrecs,nrow(data), nrow(data_complete), nrow(data_prescreening), sum(unique(data_prescreening$KeyID) %in% data_complete$KeyID))

# save changes in dataset
tosave<-append(tosave, list(nrrecs=nrrecs))

# Sort by date per donor
data<-data[order(data$KeyID,data$DonDate),] 
data_complete<-data_complete[order(data_complete$KeyID,data_complete$DonDate),] 
# Create index for nr of donations
data$numdons <- sequence(rle(data$KeyID)$lengths)
data_complete$numdons <- sequence(rle(data_complete$KeyID)$lengths)

# nr of donors
length(unique(data$KeyID)) 
# nr of donations
length((data$KeyID))       


#convert Hb levels if so required
if (!Hb_in_gpl){
  dtm<-dtm/0.06206 -1e-6 # subtract a small margin to compensate for rounding errors
  dtf<-dtf/0.06206 -1e-6
  data$Hb<-data$Hb/0.06206
  data$pre_Hb<-data$pre_Hb/0.06206
  data_complete$Hb <- data_complete$Hb/0.06206
  data_complete$pre_Hb <- data_complete$pre_Hb/0.06206
  data_prescreening$Hb <- data_prescreening$Hb/0.06206
  data_prescreening$pre_Hb <- data_prescreening$pre_Hb/0.06206
} 

# Analysis of predonation screenings 
# Set indicator for deferral
data_complete$predef<-ifelse(data_complete$pre_Hb<dtf,1,0)
data_complete$predef[data_complete$Sex=="M"]<-ifelse(data_complete$pre_Hb[data_complete$Sex=="M"]<dtm,1,0)
data_complete$predef[is.na(data_complete$predef)] <- 0 

data_prescreening$def<-ifelse(data_prescreening$pre_Hb<dtf,1,0)
data_prescreening$def[data_prescreening$Sex=="M"]<-ifelse(data_prescreening$pre_Hb[data_prescreening$Sex=="M"]<dtm,1,0)

#create year variable
data_prescreening$year<-year(data_prescreening$DonDate)

# table deferrals per year 
with(data_prescreening,table(year, Sex))
with(data_prescreening,table(Sex,year,def))
pre_deff<-with(data_prescreening[data_prescreening$Sex=="F",],table(year,def))
proportions(pre_deff, margin=1)
pre_defm<-with(data_prescreening[data_prescreening$Sex=="M",],table(year,def))
proportions(pre_defm, margin=1)
tosave<-append(tosave, list(pre_defm=pre_defm))
tosave<-append(tosave, list(pre_deff=pre_deff))

# plot deferrals per year
pre_maxdef<-max(c(proportions(pre_defm, margin=1)[,2],proportions(pre_deff, margin=1)[,2]))
if(plot_to_pdf) pdf(file=paste0(folder,"Predonationscreening_Deferrals_per_year.pdf"))
plot(rownames(pre_defm), proportions(pre_defm, margin=1)[,2], ylim=c(0,pre_maxdef*1.2), col="blue", type="l",
     ylab="Proportion deferred", xlab="Year")
lines(rownames(pre_deff), proportions(pre_deff, margin=1)[,2], col="red")
points(rownames(pre_defm), proportions(pre_defm, margin=1)[,2], pch=1, col="blue")
points(rownames(pre_deff), proportions(pre_deff, margin=1)[,2], pch=2, col="red")
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

# also make a column for the previous Hb value
data$prevHb<-NA
data$prevHb<-data$Hb[precursor]
data$prevHb[data$KeyID != data$KeyID[precursor]]<-NA

###########################################################
# Hb recovery extra analysis
###########################################################

# set plot parameters # to be set by the USER 
intervaltoshow<-c(100, 730) # time interval to show
ylim<-c(-40, 40)           # interval in dHb to show
nrtoprint<-15000           # nr of observations to select for printing
rollmeanWidth<-1000        # width of the rolling window

#index for the previous donation
idx<-1:nrow(data_complete)
precursor<-idx-1
precursor[1]<-nrow(data_complete)

# check whether the previous observation was indeed a donation
idx_noprev_donation<-idx[data_complete$KeyID == data_complete$KeyID[precursor] & is.na(data_complete$Hb[precursor])]
precursor2<-precursor[idx_noprev_donation]-1
if(F){
  ### change some records to test 
  idx2<-which(data_complete$KeyID[idx_noprev_donation]==data_complete$KeyID[precursor2])
  data_complete[precursor2[idx2[1:10]],]<-data_complete[(precursor[idx_noprev_donation])[idx2[1:10]],]
  View(data_complete[data_complete$KeyID %in% data_complete$KeyID[idx_noprev_donation[idx2][1:10]],])
}
precisvalid<-data_complete$KeyID[idx_noprev_donation] == data_complete$KeyID[precursor2]
precisHb<-!is.na(data_complete$Hb[precursor2])
sum(precisvalid) # previous record was from the same donor
sum(precisHb) # previous record had a valid Hb
precisset<-precisvalid & precisHb # reference is set if both conditions are satisfied
sum(precisset)
i<-0
while ( sum(!precisvalid | precisset)!=length(precursor2)){ # loop until all previous records found are either invalid or contain an Hb value 
  i<-i+1
  precursor2[!precisset & precisvalid]<-precursor2[!precisset & precisvalid]-1
  precisvalid[!precisset & precisvalid]<-(data_complete$KeyID[idx_noprev_donation] == data_complete$KeyID[precursor2])[!precisset & precisvalid]
  precisHb[!precisset & precisvalid]<-!is.na(data_complete$Hb[precursor2])[!precisset & precisvalid]
  precisset[!precisset & precisvalid]<-(precisvalid & precisHb)[!precisset & precisvalid]
}
paste("Nr of loops to find a valid Hb:",i)
precursor[idx_noprev_donation]<-precursor2
precursor[idx_noprev_donation][!precisset]<-NA
sum(is.na(precursor))

# also make a column for the previous Hb value
data_complete$prevHb<-NA
data_complete$prevHb<-data_complete$Hb[precursor]
data_complete$prevHb[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

data_complete$prevHb_prescr<-NA
data_complete$prevHb_prescr<-data_complete$pre_Hb[precursor]
data_complete$prevHb_prescr[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

# calculate time since last donation
data_complete$dt<-NA
data_complete$dt<-as.numeric(data_complete$DonDate-data_complete$DonDate[precursor])
data_complete$dt[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

# calculate log10 of the time difference
data_complete$ldt<-log10(data_complete$dt)

# calculate change in Hb from previous donation
data_complete$dHb<-NA
data_complete$dHb<-data_complete$Hb-data_complete$prevHb

data_complete$dHb_prescr <- NA
data_complete$dHb_prescr <- data_complete$pre_Hb-data_complete$prevHb

#association for females:
#select all females with donations above the cut off value
self<-which(data_complete$Sex=="F" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb>dtf)
linfitf_notdeferred<-lm(dHb~ldt, data=data_complete[self,])
summary(linfitf_notdeferred)
tosave<-append(tosave, list(coeff_notdeferred=c(linfitf_notdeferred$coefficients, length(self), sd(linfitf_notdeferred$residuals))))

#select all females with donations below the cut off value and the capillary follow up measurement
self<-which(data_complete$Sex=="F" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf & !is.na(data_complete$pre_Hb))
linfitf_deferred_capillary<-lm(dHb_prescr~ldt, data=data_complete[self,])
summary(linfitf_deferred_capillary)
tosave<-append(tosave, list(coeff_deferred_capillary=c(linfitf_deferred_capillary$coefficients, length(self), sd(linfitf_deferred_capillary$residuals))))

#select all females with donations below the cut off value and the venous follow up measurement
self<-which(data_complete$Sex=="F" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf)
linfitf_deferred_venous<-lm(dHb~ldt, data=data_complete[self,])
summary(linfitf_deferred_venous)
tosave<-append(tosave, list(coeff_deferred_venous=c(linfitf_deferred_venous$coefficients, length(self), sd(linfitf_deferred_venous$residuals))))

#associations for males
#select all males with donations above the cut off value
selm<-which(data_complete$Sex=="M" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb>dtf)
linfitm_notdeferred<-lm(dHb~ldt, data=data_complete[selm,])
summary(linfitm_notdeferred)
tosave<-append(tosave, list(coefm_notdeferred=c(linfitm_notdeferred$coefficients, length(selm), sd(linfitm_notdeferred$residuals))))

#select all males with donations below the cut off value and the capillary follow up measurement
selm<-which(data_complete$Sex=="M" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf & !is.na(data_complete$pre_Hb))
linfitm_deferred_capillary<-lm(dHb_prescr~ldt, data=data_complete[selm,])
summary(linfitm_deferred_capillary)
tosave<-append(tosave, list(coefm_deferred_capillary=c(linfitm_deferred_capillary$coefficients, length(selm), sd(linfitm_deferred_capillary$residuals))))

#select all males with donations below the cut off value and the venous follow up measurement
selm<-which(data_complete$Sex=="M" & data_complete$dt>=intervaltoshow[1] & data_complete$prevHb<dtf)
linfitm_deferred_venous<-lm(dHb~ldt, data=data_complete[selm,])
summary(linfitm_deferred_venous)
tosave<-append(tosave, list(coefm_deferred_venous=c(linfitm_deferred_venous$coefficients, length(selm), sd(linfitm_deferred_venous$residuals))))

#save the dataset with the prescreening
data_prescreening <- data_complete[!is.na(data_complete$pre_Hb),]

data_complete$Hbfilled<-data_complete$Hb
data_complete$Hbfilled[is.na(data_complete$Hbfilled)]<-data_complete$pre_Hb[is.na(data_complete$Hbfilled)]
data_complete <- data_complete %>%
  group_by(KeyID) %>%
  mutate(cum_Hb= cumsum(Hbfilled)) %>% ungroup() %>% mutate(meanHb = cum_Hb/numdons)

idx<-1:nrow(data_complete)
precursor<-idx-1
precursor[1]<-nrow(data_complete)

data_complete$prevMeanHb<-NA
data_complete$prevMeanHb<-data_complete$meanHb[precursor]
data_complete$prevMeanHb[data_complete$KeyID != data_complete$KeyID[precursor]]<-NA

data_complete$prevMeanHb[is.na(data_complete$prevMeanHb)] <- data_complete$meanHb[is.na(data_complete$prevMeanHb)]

###########################################################
# Analyse capillary variation 
###########################################################

data_cap<-data_complete[!is.na(data_complete$pre_Hb),]
data_cap$cap_dHb<-data_cap$pre_Hb-data_cap$prevMeanHb
(sd_cap<-sd(data_cap$cap_dHb))
(sd_capf<-sd(data_cap$cap_dHb[data_cap$Sex=="F"]))
(sd_capm<-sd(data_cap$cap_dHb[data_cap$Sex=="M"]))
(mean_cap<-mean(data_cap$cap_dHb))
(mean_capf<-mean(data_cap$cap_dHb[data_cap$Sex=="F"]))
(mean_capm<-mean(data_cap$cap_dHb[data_cap$Sex=="M"]))

tosave<-append(tosave, list(sd_cap=c(sd_cap, sd_capf, sd_capm), mean_cap=c(mean_cap, mean_capf, mean_capm)))

qqPlot(data_cap$cap_dHb, main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )
qqPlot(data_cap$cap_dHb[data_cap$Sex=="F"], main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )
qqPlot(data_cap$cap_dHb[data_cap$Sex=="M"], main="Q-Q plot of residuals capillary measurement", ylab="residuals [g/L]" )

############################
# Analyse deferrals
############################
# standard deviation of individual venous measurements
(malesd<-PostDonationResults$sdm[4]/sqrt(2))
(femalesd<-PostDonationResults$sdf[4]/sqrt(2))

data_complete$d <- qnorm(cutoffperc)*femalesd
data_complete$d[data_complete$Sex=="M"]<-qnorm(cutoffperc)*malesd
data_complete$th<-dtf
data_complete$th[data_complete$Sex=="M"]<-dtm
table(data_complete$th-data_complete$d, data_complete$Sex)


#for capillary measurements
data_complete$diff_prehb_mean <- data_complete$pre_Hb - data_complete$prevMeanHb

(malessd_pre <- sd(data_complete$diff_prehb_mean[data_complete$Sex=="M"],na.rm=T))
(malessd_pre <- sd_capm)
(femalessd_pre <- sd(data_complete$diff_prehb_mean[data_complete$Sex=="F"],na.rm=T))
(femalessd_pre <- sd_capf)


# of berekend op basis van verschil met vorig meting
sd_pre <- sd(data_cap$dHb_prescr, na.rm = T)
sd_capm
(malessd_pre <- sqrt(var(data_cap$dHb_prescr[data_cap$Sex=="M" & !is.na(data_cap$dHb_prescr)])-malesd^2))
sd_capf
(femalessd_pre <- sqrt(var(data_cap$dHb_prescr[data_cap$Sex=="F" & !is.na(data_cap$dHb_prescr)])-femalesd^2))


data_complete$pre_d<-NA
data_complete$pre_d[data_complete$Sex=="M"] <- qnorm(cutoffperc)*malessd_pre
data_complete$pre_d[data_complete$Sex=="F"] <- qnorm(cutoffperc)*femalessd_pre

tosave<-append(tosave, list(sd_cap2=c(sd_pre, femalessd_pre, malessd_pre)))

########################################
# analysis of pre-donation screenings
########################################
analysisresults_pre<-AnalysePolicyImpact_prescreening(dataframe=data_complete)
tosave<-append(tosave, list(analysisresults=analysisresults_pre$outputtableprescreening))


dataf<-data_complete[data_complete$Sex=="F",]
analysisresults_f_pre<-AnalysePolicyImpact_prescreening(dataframe=dataf)
tosave<-append(tosave, list(analysisresults_f=analysisresults_f_pre$outputtableprescreening))


datam<-data_complete[data_complete$Sex=="M",]
analysisresults_m_pre<-AnalysePolicyImpact_prescreening(dataframe=datam)
tosave<-append(tosave, list(analysisresults_m=analysisresults_m_pre$outputtableprescreening))

###################################################################
# Write tosave data to datafile
###################################################################
saveRDS(tosave, paste0(folder, "SavedDeferralData_PREDONATION_",cutoffperc,"_",Sys.Date(),".RDS"))

