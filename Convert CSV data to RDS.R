# initialize R
rm(list=ls())
while (!is.null(dev.list())) dev.off()
cat("\014")  

csvfilename<-"testdata.csv"
analysisfilename<-"analysisdata.RDS"

# read first 1500 data records to test changes to be made
data<-read.csv(csvfilename, header = T, nrows=1500)
# check out the names and classes
(classes <- sapply(data, class))

# change data as required
data[,1]<-NULL
colnames(data)<-c("KeyID","DonDate","Sex","Hb")
data$DonDate<-as.Date(data$DonDate, format = "%Y-%m-%d")
data$Sex<-as.factor(data$Sex)
head(data)

# if Ok, copy code from above and read in the full dataset
data<-read.csv(csvfilename, header = T)

data[,1]<-NULL
colnames(data)<-c("KeyID","DonDate","Sex","Hb")
data$DonDate<-as.Date(data$DonDate, format = "%Y-%m-%d")
data$Sex<-as.factor(data$Sex)
head(data)

# do some checking
table(data$Sex)
table(data$Hb,data$Sex)
min(data$DonDate)
max(data$DonDate)

# Write analysis file
saveRDS(data, file=analysisfilename)
