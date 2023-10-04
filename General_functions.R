###########################################
# General functions 
###########################################

# Set code date timestamp
generalfunctionscodedatestamp<-"20231004"

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
  main<-paste0("Donor ID = ",Sel_ID," (",ifelse(Sex[KeyID==Sel_ID]=="M", "Male", "Female"),")")
  plot(x_plot, y_plot, type="l", ylim=ylim, ylab="Haemoglobin level [g/L]", xlab="Time [Years]", main=main)
  
  # calculate updating mean values
  my<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("my<-c(my,MeanHb",i,"[KeyID%in%Sel_ID])")))
  # calculate overall mean Hb level
  mys<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("mys<-c(mys,MeanHb",length(x),"[KeyID%in%Sel_ID])")))
  # calculate upper threshold values
  lyt<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("lyt<-c(lyt, MeanHb",i,"[KeyID%in%Sel_ID]+d[KeyID%in%Sel_ID]/sqrt(",i,"))")))
  # calculate lower threshold values for individual measurement
  lyt2<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("lyt2<-c(lyt2, MeanHb",i,"[KeyID%in%Sel_ID]-d[KeyID%in%Sel_ID])")))
  # has to be shifted to the next values
  if (length(lyt2)>1) lyt2<-c(NA,lyt2[1:(length(lyt2)-1)])
  lyt2u<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("lyt2u<-c(lyt2u, MeanHb",i,"[KeyID%in%Sel_ID]+d[KeyID%in%Sel_ID])")))
  # has to be shifted to the next values
  if (length(lyt2u)>1) lyt2u<-c(NA,lyt2u[1:(length(lyt2u)-1)])
  # calculate extended lower threshold, accounting for mean variability as well!
  lyt3<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("lyt3<-c(lyt3, MeanHb",i,"[KeyID%in%Sel_ID]-d[KeyID%in%Sel_ID]*(1+1/sqrt(",i,")))")))
  # has to be shifted to the next values
  if (length(lyt3)>1) lyt3<-c(NA,lyt3[1:(length(lyt3)-1)])
  lyt3u<-c()
  for(i in 1:length(x)) eval(parse(text=paste0("lyt3u<-c(lyt3u, MeanHb",i,"[KeyID%in%Sel_ID]+d[KeyID%in%Sel_ID]*(1+1/sqrt(",i,")))")))
  # has to be shifted to the next values
  if (length(lyt3u)>1) lyt3u<-c(NA,lyt3u[1:(length(lyt3u)-1)])
  
  # is there permanent deferral?
  def<-which(lyt<th[KeyID%in%Sel_ID]) # are there any values below the threshold?
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
  abline(h=th[KeyID==Sel_ID], col=1, lwd=2)
  abline(h=th[KeyID==Sel_ID]-d[KeyID==Sel_ID], col=1, lwd=1, lty=4)
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
  sp<- y>=th[KeyID%in%Sel_ID]
  sp[1]<-F
  points(x_plot[sp], y_plot[sp], pch=16)
  # mark unnecessary donations
  sp<- y<th[KeyID%in%Sel_ID] & y>=lyt2
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
  print(paste0("Donor ", Sel_ID,": ", round(sum(y[2:length(y)]<th[KeyID%in%Sel_ID])/(length(y)-1),2)*100,"% of attempts deferred"))
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

fitHbdistributions<-function(data, nrofquantiles=20) {
  # function that fits kernel density and smoothing spline for various ranges 
  # of nr of donations per donor. The nr of splits is determined by 
  # the parameter nrofquantiles
  # the input data should consist of variables Hb, sd and Nrdon
  # the function produces three plots; the nr of observations per splitpoint, and
  # for Hb and sd per splitpoint the kde, spline and normal fits with the data
  # the counts and fit parameters are returned by the function
  
  quantiles <- quantile(data$Nrdon, prob = seq(0, 1, length = nrofquantiles+1), type = 5)
  data$cutted <- cut2(data$Nrdon, cuts = unique(as.numeric(quantiles)))
  
  levelsn<-sort(unique(as.numeric(data$cutted)))
  nrsplits<-length(levels(data$cutted))
  hist(as.numeric(data$cutted), xaxt = "n", main="Number of donations per cluster", xlab="Cluster of number of donations", breaks=c(levelsn-.5, max(levelsn)+.5))
  axis(1, at = sort(unique(as.numeric(data$cutted))), labels = levels(data$cutted))
  sum(table(data$Nrdon))
  nrobs<-table(data$cutted, useNA="always")
  print(nrobs)
  correctforone<-ifelse(as.numeric(dimnames(nrobs)[[1]][1])==1,1,0)
  minsubset<-min(nrobs[nrobs>0]) # set minimum subset size
  
  # set frame for plotting
  lid<-round(sqrt(nrsplits))
  par(mfrow=c(lid,ceiling(nrsplits/lid)))
  
  # distribution of Hb levels
  par(mfrow=c(lid,ceiling(nrsplits/lid)))
  Hbdistr<-list(n=table(data$cutted, useNA="always"))
  for (i in 1:length(levels(data$cutted))){
    de<-density(data$Hb[data$cutted==levels(data$cutted)[i]])
    de$s<-cumsum(de$y)/sum(de$y)
    spl <- with(de, smooth.spline(x, s, df = 25))
    
    normfit<-fitdist(data$Hb[data$cutted==levels(data$cutted)[i]], 'norm')
    denscomp(normfit, xlab="Hb g/L", main=paste0("Mean Hb for n=", levels(data$cutted)[i]))
    lines(de$x,de$Hb)
    spl <- with(de,smooth.spline(x, s, df = 40))
    lines(predict(spl, de$x, deriv = 1), col = "blue")
    eval(parse(text=paste0("Hbdistr<-append(Hbdistr,list(de",i,"=de))")))
    eval(parse(text=paste0("Hbdistr<-append(Hbdistr,list(spl",i,"=spl))")))
  }
  # distribution of Hb sd estimates
  lid<-round(sqrt(nrsplits-correctforone))
  par(mfrow=c(lid,ceiling((nrsplits-correctforone)/lid)))
  Hbsddistr<-list()
  for (i in 1:length(levels(data$cutted))){
    dat<-data$sd[data$cutted==levels(data$cutted)[i]]
    dat<-dat[!is.na(dat)]
    if (length(dat)>1 & str_trim(dimnames(nrobs)[[1]][i])!="1") {
      if(length(dat)<minsubset) minsubset<-length(dat)
      de<-density(dat)
      de$s<-cumsum(de$y)/sum(de$y)
      spl <- with(de,smooth.spline(x, s, df = 25))
      
      normfit<-fitdist(dat, 'norm')
      denscomp(normfit, xlab="Hb g/L", main=paste0("Hb Sd for n=", levels(data$cutted)[i]))
      lines(de$x,de$y)
      spl <- with(de,smooth.spline(x, s, df = 40))
      lines(predict(spl, de$x, deriv = 1), col = "blue")
      eval(parse(text=paste0("Hbsddistr<-append(Hbsddistr,list(de",i,"=de))")))
      eval(parse(text=paste0("Hbsddistr<-append(Hbsddistr,list(spl",i,"=spl))")))
      eval(parse(text=paste0("Hbsddistr<-append(Hbsddistr,list(n",i,"=length(dat)))")))
    }
  }
  par(mfrow=c(1,1))
  print(paste("minimum subset size:", minsubset))
  return(list(Hbdistr=Hbdistr, Hbsddistr=Hbsddistr,minsubset=minsubset))
}

AnalysePolicyImpact<-function(){
  # function that analyses the policy impact for a dataset datt containing aggregated donation and Hb data per donor per donation
  
  # define an output array with 8 columns, one row per subsequent donation
  # the items that are stored in per column are explained below

  outputsummarytable<-as.data.frame(matrix(0,maxDons,8))
  colnames(outputsummarytable)<-c("Deferred", "Non-deferred", "ShouldNotDeferred", "ShouldNotDonate", 
                                  "Should_not_have_donated","Missed_by_stopped_donor", "ShouldNotDeferred2", "RequiresReview")
  
  stopped<-rep(F, length(Hb1)) # indicator for whether a donors has stopped or not
  # is set when the mean Hb level was demonstrably below 
  # the eligibility threshold at previous donation 
  stopped2<-rep(F, length(Hb1)) # indicator for whether a donors has stopped or not
  # is set when the mean Hb level is below the 
  # eligibility threshold at previous donation
  stopafter<<-3 # stop donating after significant evidence only after stopafter donations have been made, is required globally
  
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
  return(list(outputsummarytable=outputsummarytable, stopped=stopped, stopped2=stopped2))
  
}
