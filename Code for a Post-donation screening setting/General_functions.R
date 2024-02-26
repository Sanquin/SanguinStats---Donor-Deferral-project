###########################################
# General functions 
###########################################

# Set code date timestamp
generalfunctionscodedatestamp<-"20240222"

# function that plots the donation profile of an individual donor
plotdonorprofile<-function(Sel_ID, leg=F, ylim=c(0,200)) {
  # Sel_ID = KeyID of the donor to print
  # leg = include legend in the plot
  # ylim = limits of Hb levels to plot
  
  x<-data[data$KeyID%in%Sel_ID,]$DonDate
  y<-data[data$KeyID%in%Sel_ID,]$Hb
  #Add some random days to date
  x_plot <- x + sample(6,length(x),replace = TRUE) - 3
  # Add a little bit of jitter to Hb 
  y_plot <- jitter(y)
  main<-paste0("Donor ID = ",Sel_ID," (",ifelse(grepl("M", data$Sex[data$KeyID==Sel_ID]), "Male", "Female"),")")[1]
  plot(x_plot, y_plot, type="l", ylim=ylim, ylab="Haemoglobin level [g/L]", xlab="Time [Years]", main=main)
  
  # calculate updating mean values
  my <- data$meanHb[data$KeyID%in%Sel_ID]
  # calculate overall mean Hb level
  mys<- tail(my, n=1)
  # calculate upper threshold values
  lyt<-data$meanHb[data$KeyID%in%Sel_ID]+data$d[data$KeyID%in%Sel_ID]/sqrt(data$numdons[data$KeyID%in%Sel_ID])
  # calculate lower threshold values for individual measurement
  lyt2<-data$meanHb[data$KeyID%in%Sel_ID]-data$d[data$KeyID%in%Sel_ID]
  # has to be shifted to the next values
  if (length(lyt2)>1) lyt2<-c(NA,lyt2[1:(length(lyt2)-1)])
  lyt2u<-data$meanHb[data$KeyID%in%Sel_ID]+data$d[data$KeyID%in%Sel_ID]
  # has to be shifted to the next values
  if (length(lyt2u)>1) lyt2u<-c(NA,lyt2u[1:(length(lyt2u)-1)])
  # calculate extended lower threshold, accounting for mean variability as well!
  lyt3<-data$meanHb[data$KeyID%in%Sel_ID]-data$d[data$KeyID%in%Sel_ID]*(1+1/sqrt(data$numdons[data$KeyID%in%Sel_ID]))
  # has to be shifted to the next values
  if (length(lyt3)>1) lyt3<-c(NA,lyt3[1:(length(lyt3)-1)])
  
  lyt3u<-data$meanHb[data$KeyID%in%Sel_ID]+data$d[data$KeyID%in%Sel_ID]*(1+1/sqrt(data$numdons[data$KeyID%in%Sel_ID]))
  # has to be shifted to the next values
  if (length(lyt3u)>1) lyt3u<-c(NA,lyt3u[1:(length(lyt3u)-1)])
  
  # is there permanent deferral?
  def<-which(lyt<data$th[data$KeyID%in%Sel_ID]) # are there any values below the threshold?
  print(paste("Permanent deferral points:", paste(def, collapse=" ")))
  defind<-T
  if(length(def)>0) defind<-which(def>stopafter) # if so, are these later than the stopafter measurement?
  
  # plot polygons
  # for drawing a polygon data have to be added in reverse order
  # also, the reverse side of the mean needs to be calculated: 2 x mean-deferral level
  polygon(c(x,x[length(x):1]),c(lyt,2*my[length(x):1]-lyt[length(x):1]),
          col=rgb(.5, 0, 0,0.05), border=NA)
  polygon(c(x[2:length(x)],x[length(x):2]),c(lyt2[2:length(x)],2*my[length(x):2]-lyt[length(x):2]),
          col=rgb(0, .5,0, 0.05), border=NA)
  polygon(c(x[2:length(x)],x[length(x):2]),c(lyt2u[2:length(x)],lyt[length(x):2]),
          col=rgb(0, .5,0, 0.05), border=NA)
  polygon(c(x[2:length(x)],x[length(x):2]),c(lyt3[2:length(x)],lyt2[length(x):2]),
          col=rgb(0.1, .1,0.1, 0.05), border=NA)
  polygon(c(x[2:length(x)],x[length(x):2]),c(lyt3u[2:length(x)],lyt2u[length(x):2]),
          col=rgb(0.1, .1,0.1, 0.05), border=NA)
  
  #plot the lines
  abline(h=data$th[data$KeyID==Sel_ID][1], col=1, lwd=2)
  abline(h=data$th[data$KeyID==Sel_ID][1]-data$d[data$KeyID==Sel_ID][1], col=1, lwd=1, lty=4)
  lines(x, my       , lty=2, col=4)
  lines(x, lyt      , lty=3, col=2)
  lines(x, 2*my-lyt , lty=3, col=2)
  lines(x, lyt2     , lty=4, col=3)
  lines(x, lyt2u    , lty=4, col=3)
  lines(x, lyt3     , lty=5, col=8)
  lines(x, lyt3u    , lty=5, col=8)
  if(length(def)>0 & is.integer(defind[1])) abline(v=x[def[defind[1]]], col=8) # if so, print
  
  # plot points
  # mark donations
  sp<- y>=data$th[data$KeyID%in%Sel_ID]
  sp[1]<-F
  points(x_plot[sp], y_plot[sp], pch=16)
  # mark unnecessary donations
  sp<- y<data$th[data$KeyID%in%Sel_ID] & y>=lyt2
  sp[1]<-F
  points(x_plot[sp], y_plot[sp])
  # mark deferrals 
  sp<- y<lyt2
  #    points(x_plot[sp], y[sp], pch=15)
  points(x_plot[sp], y_plot[sp], pch=0)										
  
  #  points(x_plot[1], x_plot[1], col=0, pch=16)
  points(x_plot[1], y_plot[1], pch=8)
  
  if (leg) {
    cutoffperc2<-1-(1-cutoffperc)*2
    legend("bottomright",
           c("Mean Hb level",
             paste0("Donor mean ",cutoffperc2*100,"% confidence interval (CI)"),
             paste0("Donation ",cutoffperc2*100,"% confidence interval (relative to mean Hb level)"),
             paste0("Donation ",cutoffperc2*100,"% confidence interval (accounting for mean Hb variability)"),
             "Legal deferral threshold",
             "Donation deferral threshold"
           ),
           lty=c( 2, 3, 4, 5, 1, 4),
           pch=c(NA,NA,NA,NA,NA,NA),
           col=c( 4, 2, 3, 8, 1, 1),
           lwd=c( 1, 1, 1, 1, 2, 1), cex=1)
    legend("bottomleft",
           c("New donor intake",
             "Donation performed",
             "Donor deferred",
             "Donor deferred (outlier observation)"),
           lty=c(NA,NA,NA,NA),
           #pch=c( 8,16, 1,15),
           pch=c( 8,16, 1, 0),
           col=c( 1, 1, 1, 1),
           lwd=c( 1, 1, 1, 1), cex=1)
  }
  print(paste0("Donor ", Sel_ID,": ", round(sum(y[2:length(y)]<data$th[data$KeyID%in%Sel_ID][1])/(length(y)-1),2)*100,"% of attempts deferred"))
  print("")
}

# function to set gender specific limits when ploting donorprofiles
setlims<-function(l,y){
  if (missing(y)) {
    my<-c(110,160)
    fy<-c(95,145)
    my<-c(90,180)
    fy<-c(80,170)
  } else {
    my<-c(y[1],y[2])
    fy<-c(y[3],y[4])
  }
  if(Sex[KeyID==Sel_ID]=="M") return(my)
  else return(fy)
}

# function to plot multiple selected donors in a matrix
# the optional seedvalue allows changing the selection to show
plotmatrix<-function(selID, maxplots, ylim=c(80,180), seedvalue=1){
  if (length(selID)<=maxplots^2){ # only if the number of donors selected is limited
    lid<-round(sqrt(length(selID)))
    par(mfrow=c(lid,ceiling(length(selID)/lid)))
    for(l in selID) plotdonorprofile(l, ylim=ylim)
    par(mfrow=c(1,1))
  } else {
    set.seed(seedvalue)
    par(mfrow=c(maxplots,maxplots))
    for(l in sample(selID,maxplots^2,replace=F)) plotdonorprofile(l, ylim=ylim)
    par(mfrow=c(1,1))
  }
}

FitDistributions<-function(data,variable) {
  # function that fits kernel density and smoothing spline for various ranges 
  # of nr of donations per donor. The nr of splits is preset by the vector breakpoints
  # the input data should consist of variables Hb/interval, sd and Nrdon
  # the function produces three plots; the nr of observations per splitpoint, and
  # for Hb/interval and sd per splitpoint the kde, spline and normal fits. The data
  # the counts and fit parameters as well as a list of quantiles per splitpoint are 
  # returned by the function
  breakpoints <- c(0:7, (cumsum(1:20) + 7), Inf)
  labels <-c("[0,1]", "(1,2]", "(2,3]", "(3,4]", "(4,5]", "(5,6]", "(6,7]", "(7,8]", "(8,10]", "(10,13]", "(13,17]", "(17,22]", "(22,28]", "(28,35]", "(35,43]", "(43,52]", "(52,62]", "(62,73]", "(73,85]", "(85,98]", "(98,112]", "(112,127]", "(127,143]", "(143,160]", "(160,178]", "(178,197]", "(197,217]", "(217,Inf]")

  nrq<-20 # nr of quantiles to return per splitpoint
    
  data$cutted <- cut(data$Nrdon, breaks = breakpoints, include.lowest = TRUE)
  level_counts <- table(data$cutted)
  # Extract levels with non-zero counts
  axis.labels <- names(level_counts[level_counts > 0])
  data$cutted <- as.numeric(data$cutted)
  
  levelsn<-sort(unique(as.numeric(data$cutted)))
  hist(as.numeric(data$cutted), xaxt = "n", main="Number of donations per cluster", xlab="Cluster of number of donations", breaks=c(levelsn-.5, max(levelsn)+.5))
  axis(1, at = sort(unique(as.numeric(data$cutted))), labels = axis.labels)
  nrobs<-table(data$cutted, useNA="always")
  print(nrobs)
  minsubset<-min(nrobs[nrobs>0]) # set minimum subset size
  
  # set frame for plotting
  lid<-round(sqrt(length(levelsn)))
  par(mfrow=c(lid,ceiling(length(levelsn)/lid)))
  
  # distribution of variable
  par(mfrow=c(lid,ceiling(length(levelsn)/lid)))
  eval(parse(text=paste0(variable,"distr<-list(n=table(data$cutted, useNA=\"always\"))")))
  if(variable == "Hb"){
    xlab = "Hb g/L"
    xlim = c(100,200)
  } else if (variable == "interval"){
    xlab = "Days"
    xlim = c(0,1000)
  } else{
    xlab = " "
    xlim=c(0,1)
  }
  distr_table <- c()
  for (level in levelsn){
    eval(parse(text=paste0("de<-density(data$",variable,"[data$cutted==",level,"])")))
    de$s<-cumsum(de$y)/sum(de$y)
    spl <- with(de, smooth.spline(x, s, df = 25))
    eval(parse(text=paste0("normfit<-fitdist(data$",variable,"[data$cutted==",level,"], 'norm')")))
    title <- paste0("Mean ", variable, " for \n n=", labels[level])
    eval(parse(text = paste0("denscomp(normfit, xlab=xlab, addlegend=F, main=title)")))
    eval(parse(text=paste0("lines(de$x,de$",variable,")")))
    
    spl <- with(de,smooth.spline(x, s, df = 40))
    lines(predict(spl, de$x, deriv = 1), col = "blue")
    eval(parse(text=paste0(variable, "distr<-append(",variable,"distr,list(de",level,"=de))")))
    eval(parse(text=paste0(variable, "distr<-append(",variable,"distr,list(spl",level,"=spl))")))
    eval(parse(text=paste0("distr_table <- rbind(distr_table, c(",level,", labels[",level,"],length(data$",variable,"[data$cutted==",level,"]), mean(data$",variable,"[data$cutted==",level,"], na.rm=T), sd(data$",variable,"[data$cutted==",level,"], na.rm=T), median(data$",variable,"[data$cutted==",level,"], na.rm=T),    
       quantile(data$",variable,"[data$cutted==",level,"], (1:nrq)/nrq-1/2/nrq) ))")))
  }
  colnames(distr_table)<-c("group", "label", "n", "mean", "sd", "median", paste0(((1:nrq)/nrq-1/2/nrq)*100,"%"))

  # distribution of sd estimates
  lid<-round(sqrt(length(levelsn)))
  par(mfrow=c(lid,ceiling(length(levelsn)/lid)))
  eval(parse(text=paste0(variable,"sddistr<-list(n=table(data$cutted, useNA=\"always\"))")))
  distr_table_sd <- c()
  for (level in levelsn){
    eval(parse(text=paste0("dat<-data$sd[data$cutted==",level,"]")))
    dat<-dat[!is.na(dat)]
    if (length(dat)>1) {
      if(length(dat)<minsubset) minsubset<-length(dat)
      de<-density(dat)
      de$s<-cumsum(de$y)/sum(de$y)
      spl <- with(de,smooth.spline(x, s, df = 25))
      distr_table_sd <- rbind(distr_table_sd, c(level, labels[level], length(dat), mean(dat, na.rm=T), sd(dat, na.rm=T), median(dat, na.rm=T), quantile(dat, (1:nrq)/nrq-1/2/nrq)))
      normfit<-fitdist(dat, 'norm')
      title <- paste0("Sd for \n n=", labels[level])
      eval(parse(text=paste0("denscomp(normfit, xlab=xlab, addlegend=F, main=title)")))
      lines(de$x,de$y)
      spl <- with(de,smooth.spline(x, s, df = 40))
      lines(predict(spl, de$x, deriv = 1), col = "blue")
      eval(parse(text=paste0(variable, "sddistr<-append(",variable,"sddistr,list(de",level,"=de))")))
      eval(parse(text=paste0(variable, "sddistr<-append(",variable,"sddistr,list(spl",level,"=spl))")))
      eval(parse(text=paste0(variable, "sddistr<-append(",variable,"sddistr,list(n",level,"=length(dat)))")))
    }
  }
  colnames(distr_table_sd)<-c("group", "label", "n", "mean", "sd", "median", paste0(((1:nrq)/nrq-1/2/nrq)*100,"%"))
  par(mfrow=c(1,1))
  print(paste("minimum subset size:", minsubset))
  eval(parse(text=paste0("return(list(",variable,"distr=",variable,"distr, ",variable,"sddistr=",variable,"sddistr,minsubset=minsubset, distr_table=distr_table, distr_table_sd = distr_table_sd))")))
}


AnalysePolicyImpact<-function(dataframe){
  # function that analyses the policy impact for a dataset datt containing aggregated donation and Hb data per donor per donation
  dataset <- dataframe
  # define an output array with 8 columns, one row per subsequent donation
  # the items that are stored in per column are explained below

  outputsummarytable<-as.data.frame(matrix(0,maxDons,8))
  colnames(outputsummarytable)<-c("Deferred", "Non-deferred", "ShouldNotDeferred", "ShouldNotDonate", 
                                  "Should_not_have_donated","Missed_by_stopped_donor", "ShouldNotDeferred2", "RequiresReview")
  
  stopped_KeyID <<- NA # indicator for whether a donors has stopped or not
  # is set when the mean Hb level was demonstrably below 
  # the eligibility threshold at previous donation 
  stopped2_KeyID <<- NA # indicator for whether a donors has stopped or not
  # is set when the mean Hb level is below the 
  # eligibility threshold at previous donation
  stopafter <<- 3# stop donating after significant evidence only after stopafter donations have been made, is required globally
  maxDons <<- max(data$numdons)
  
  for (i in 1:maxDons ){
    calculate <<- dataset[dataset$numdons==i,]
    # 1 - Deferred donors
    # The number of donors deferred at step i are those with a Hb value that is  
    # below the deferral threshold
    outputsummarytable[i,1] <- sum(!is.na(calculate$Hb) & calculate$Hb<calculate$th)
    
    # 2 - non-Deferred donors
    # The number of donors not deferred at step i are those with a Hb value that   
    # is equal or larger than the deferral threshold
    outputsummarytable[i,2] <- sum(!is.na(calculate$Hb) & calculate$Hb >= calculate$th)
    
    # 3 - Donors that should not have been deferred as the Hb deviation relative to
    #     their mean Hb value does not provide sufficient evidence against donation
    # only count events from second donation onwards
    if(i>1) outputsummarytable[i,3] <- sum(!is.na(calculate$Hb) & calculate$Hb < calculate$th & calculate$Hb >= calculate$prevMeanHb-calculate$d, na.rm=T)
    
    # 4 - Donors that should not donate as their Hb is demonstrably below the eligibility threshold
    outputsummarytable[i,4]<-sum(!is.na(calculate$Hb) & calculate$meanHb < calculate$th - calculate$d/sqrt(i), na.rm=T)
    
    # 5 - Donors that should not have donated
    outputsummarytable[i,5] <- sum(!is.na(calculate$Hb) & calculate$prevMeanHb < calculate$th-calculate$d/sqrt(i-1) & calculate$Hb >= calculate$th, na.rm=T)
    
    # set new index for (previously) stopped donors
    # index stopped indicates that the mean Hb level was demonstrably below the eligibility threshold at previous donation 
    if(i>stopafter) stopped_KeyID <<- c(stopped_KeyID, calculate$KeyID[!is.na(calculate$Hb) & calculate$prevMeanHb < calculate$th - calculate$d/sqrt(i-1)])
    # index stopped2 indicates that the mean Hb level is below the eligibility threshold at previous donation
    if (i>stopafter) stopped2_KeyID <<- c(stopped2_KeyID, calculate$KeyID[!is.na(calculate$Hb) & calculate$prevMeanHb < calculate$th])
    
    # 6 - donations missed as a result of new deferral rule
    if(i>stopafter) outputsummarytable[i,6] <- sum(!is.na(calculate$Hb) & calculate$Hb >= calculate$th & calculate$KeyID %in% stopped_KeyID, na.rm=T)
    
    
    # 7 - Donors that should not have been deferred as the Hb deviation relative to 
    #     the absolute Hb threshold is insufficient (see also evaluation 3 above)
    # this basically presumes that anyone with a Hb level over th-d may donate
    outputsummarytable[i,7]<- sum(!is.na(calculate$Hb) & calculate$Hb< calculate$th & calculate$Hb>=calculate$th-calculate$d, na.rm=T)
    # 8 - Identified as outlier, but not deferred
    if(i>1) outputsummarytable[i,8] <- sum(!is.na(calculate$Hb) & calculate$Hb < calculate$prevMeanHb - calculate$d & calculate$Hb >= calculate$th, na.rm=T)
  }
  stopped <<- length(unique(stopped_KeyID)) 
  stopped2 <<- length(unique(stopped2_KeyID))
  
  length(unique(stopped_KeyID))/length(stopped_KeyID) #proportion of stopped donors
  length(unique(stopped2_KeyID))/length(stopped2_KeyID) #proportion of stopped2 donors
  
  outputsummarytable$defprop<-outputsummarytable$Deferred/(outputsummarytable$Deferred+outputsummarytable$`Non-deferred`)
  outputsummarytable$nondefprop<-outputsummarytable$ShouldNotDeferred/outputsummarytable$Deferred
  outputsummarytable
  return(list(outputsummarytable=outputsummarytable, stopped=stopped, stopped2=stopped2))
  
}

AnalysePolicyImpact_prescreening <- function(dataframe){
  data_complete <- dataframe
  
  maxDons <- max(data_complete$numdons)
  outputtableprescreening <- as.data.frame(matrix(0,maxDons,4))
  colnames(outputtableprescreening) <- c("Pre-donation screenings conducted", "Deferral during pre-donation screening", "Should not deferred", "Should not donate")
  
  for(i in 1:maxDons){
    calculate <- data_complete[data_complete$numdons==i,]
    #1 - number of pre-screenings conducted
    outputtableprescreening[i,1] <- sum(!is.na(calculate$pre_Hb))
    
    #2 - number of pre-screenings that were a deferral
    outputtableprescreening[i,2] <- sum(calculate$predef==1)
    
    #3 - deferrals at pre-screening that were not necessary, because the pre-screening Hb was within expected variation of the mean
    outputtableprescreening[i,3] <- sum(!is.na(calculate$pre_Hb) & calculate$pre_Hb < calculate$th & calculate$pre_Hb >= calculate$prevMeanHb-calculate$pre_d, na.rm=T)
    
    #4 - donations done after the pre-screening measurement, that should not have been done because the pre-screening Hb was not sufficient
    outputtableprescreening[i,4]<- sum(!is.na(calculate$pre_Hb) & ((calculate$pre_Hb + calculate$prevMeanHb)/i) < calculate$th - calculate$pre_d/sqrt(i), na.rm=T)
  }
  return(list(outputtableprescreening=outputtableprescreening))
}