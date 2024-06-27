####################################################
# Exploratory re-analysis of AB MPB R values data  #
# published as Carroll et al 2017                  #
# Barry J. Cooke                                   #
# February 21, 2023                                #
# Major revision March 25                          #
# New coding June 18, 2024                         #
####################################################

library(gstat)  # inverse distance weighted, Kriging
library(fields)
library(dplyr)
library(sf)
library(raster)
library(wesanderson)

#This is the original data file for the 2017 FRI report, which had years pre-merged
#abr <- read.csv("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\rvaluerawdata20062015Aug2016.csv",header=T)

#This input file has jack pine introgression coefficients, Q, merged
infile<-paste0(getwd(),"/data/FRI/rvaluesQvalues.csv")
abr <- read.csv(infile, header = T,na.strings=".")

#quick check of lat/lon
win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr$PLOT_LON,abr$PLOT_LAT)

#quick check of some variables
win.graph(height=6,width=6)
par(mfrow=c(2,2))
hist(abr$ELEV)
hist(abr$DBH)
hist(log10(abr$HT_PITCH_T))
hist(abr$NBR_INFEST)

sum(na.omit(abr$NBR_INFEST==0))

#FRI report has mostly non-zero infestation levels for 2006, 2007, 2009
win.graph(height=4,width=4)
boxplot(log10(abr$NBR_INFEST+1)~abr$Beetle_YEAR)

############################################
# There appear to be some database errors  #
# Rather than adjust source data, here are #
# suggested fixes for DBH's and pitch tube #
# heights that are too large               #
############################################

#if DBH>1000, divide by 100
abr$DBH[abr$DBH>1000]<-abr$DBH[abr$DBH>1000]/100

#if DBH>100, divide by 10
abr$DBH[abr$DBH>100]<-abr$DBH[abr$DBH>100]/10

#if YEAR==2006, divide HT_PITCH_T by 100
abr$HT_PITCH_T[abr$Beetle_YEAR==2006]<-abr$HT_PITCH_T[abr$Beetle_YEAR==2006]/100

#if HT_PITCH_T>26, divide HT_PITCH_T by 10
switch<-(abr$HT_PITCH_T>26)&(!is.na(abr$HT_PITCH))

abr$HT_PITCH_T[switch]<-abr$HT_PITCH_T[switch]/10

#########################################
# Also, ELEV cannot equal 0, must be NA #
#########################################

abr$ELEV[abr$ELEV==0]<-NA

###############################################################################################################
# This is the new set of data files post March 20, which have years de-merged, and columns named differently. #
# A long set of steps are required to merge and QC these data.                                                #
###############################################################################################################

########
# 2006 #
########

abr.2006 <- read.csv("data/FRI/R_value_BY2006.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2006$plot_long_dd,abr.2006$plot_lat_dd)

sum(abr.2006$plot_long_dd==0)
#correct zero longitude based on looking at spreadsheet neighbouring record
abr.2006$plot_long_dd[abr.2006$plot_long_dd==0]<--119.507

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2006$plot_long_dd,abr.2006$plot_lat_dd)

#compare to older dataset
use<-abr$Beetle_YEAR==2006
win.graph(height=8,width=5)
par(mfrow=c(1,2))
plot(abr$PLOT_LON[use],abr$PLOT_LAT[use],main="2006, old dataset",ylim=c(49,60))
plot(abr.2006$plot_long_dd,abr.2006$plot_lat_dd,main="2006, new dataset",ylim=c(49,60))
length(abr$PLOT_LON[use])
length(abr.2006$plot_long_dd)

par(mfrow=c(2,3))
hist(abr.2006$R_VALUE)
hist(abr.2006$NBR_TREES)
hist(abr.2006$dbh)
hist(abr.2006$ht_pitch_tube)
hist(abr.2006$CF)
hist(abr.2006$SSI)

#several dbhs's too high
replace<-(abr.2006$dbh>80)&(!is.na(abr.2006$dbh))
sum(replace)
abr.2006$dbh[replace]
abr.2006$dbh[replace]<-abr.2006$dbh[replace]/10

#several pitch tube beights too high
replace<-(abr.2006$ht_pitch_tube>22)&(!is.na(abr.2006$ht_pitch_tube))
sum(replace)
abr.2006$ht_pitch_tube[replace]
abr.2006$ht_pitch_tube[replace]<-abr.2006$ht_pitch_tube[replace]/10
#This takes one value of 325 and coverts it to 32.5. That's probably still too high.

win.graph(width=12,height=7)
par(mfrow=c(2,3))
hist(log10(abr.2006$R_VALUE+1))
hist(abr.2006$NBR_TREES)
hist(abr.2006$dbh)
hist(abr.2006$ht_pitch_tube)
hist(abr.2006$CF)
hist(abr.2006$SSI)

abr.lm<-lm(log10(abr.2006$R_VALUE+1)~abr.2006$NBR_TREES+abr.2006$dbh+abr.2006$ht_pitch_tube+abr.2006$CF+abr.2006$SSI)
summary(abr.lm)

########
# 2007 #
########

abr.2007 <- read.csv("data/FRI/R_value_BY2007.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2007$LONG_DD,abr.2007$LAT_DD)

sum(abr.2007$LONG_DD==0)
#correct zero longitude : delete entire record
abr.2007$OBJECTID[abr.2007$LONG_DD==0] # 3 records are junk: 93, 602, 818
abr.2007<-abr.2007[-c(abr.2007$OBJECTID[abr.2007$LONG_DD==0]),]
plot(abr.2007$LONG_DD,abr.2007$LAT_DD)

#several too far West
sum(abr.2007$LONG_DD<(-121))
abr.2007$LONG_DD[abr.2007$LONG_DD<(-121)]

#can't assess exact values, so replace garbage values in order from lowest to highest
abr.2007$LONG_DD[abr.2007$LONG_DD<(-449)]<-(-119.8281)
abr.2007$LONG_DD[abr.2007$LONG_DD<(-199)]<-(-119.7627)
abr.2007$LONG_DD[abr.2007$LONG_DD<(-192)]<-(-119.4212)

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2007$LONG_DD,abr.2007$LAT_DD)

#one is still too far west; an error in the minutes >60
abr.2007$LONG_DD[abr.2007$LONG_DD<(-120.5)]
abr.2007$OBJECTID[abr.2007$LONG_DD<(-120.5)]
abr.2007$LONG_DD[abr.2007$LONG_DD<(-120.5)]<- (-119-51.743/60)

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2007$LONG_DD,abr.2007$LAT_DD)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2007$R_VALUE+1))
hist(abr.2007$SAMPLE_TRE)
hist(abr.2007$DBH)
hist(abr.2007$HT_PITCH_T)
hist(abr.2007$CF)
hist(abr.2007$SSI)

#several dbhs's too high
replace<-(abr.2007$DBH>80)&(!is.na(abr.2007$DBH))
sum(replace)
abr.2007$DBH[replace]
abr.2007$DBH[replace]<-abr.2007$DBH[replace]/100

#several pitch tube beights too high
replace<-(abr.2007$HT_PITCH_T>22)&(!is.na(abr.2007$HT_PITCH_T))
sum(replace)
abr.2007$HT_PITCH_T[replace]
abr.2007$HT_PITCH_T[replace]<-abr.2007$HT_PITCH_T[replace]/100

win.graph(width=12,height=7)
par(mfrow=c(2,3))
hist(log10(abr.2007$R_VALUE+1))
hist(abr.2007$SAMPLE_TRE)
hist(abr.2007$DBH)
hist(abr.2007$HT_PITCH_T)
hist(abr.2007$CF)
hist(abr.2007$SSI)

abr.lm<-lm(log10(abr.2007$R_VALUE+1)~abr.2007$SAMPLE_TRE+abr.2007$DBH+abr.2007$HT_PITCH_T+abr.2007$CF+abr.2007$SSI)
summary(abr.lm)

########
# 2008 #
########

abr.2008 <- read.csv("data/FRI/R_value_BY2008.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2008$plot_long_dd,abr.2008$plot_lat_dd)

#Two -999 for r value
sum(abr.2008$R_VALUE==-999)

abr.2008$R_VALUE[abr.2008$R_VALUE==-999]<-NA

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2008$R_VALUE+1))
hist(log10(abr.2008$NBR_INFEST+1))
hist(abr.2008$dbh)
hist(abr.2008$ht_pitch_tube)
hist(abr.2008$CF)
hist(abr.2008$SSI)

abr.2008$BEETLE_YR[abr.2008$BEETLE_YR==2208]<-2008

abr.lm<-lm(log10(abr.2008$R_VALUE+1)~log10(abr.2008$NBR_INFEST+1)+abr.2008$dbh+abr.2008$ht_pitch_tube+abr.2008$CF+abr.2008$SSI)
summary(abr.lm)

########
# 2009 #
########

abr.2009 <- read.csv("data/FRI/R_value_BY2009.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2009$plot_long_dd,abr.2009$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2009$r_value+1))
hist(log10(abr.2009$nbr_infest+1))
hist(abr.2009$dbh)
hist(abr.2009$ht_pitch_tube)
hist(abr.2009$CF)
hist(abr.2009$SSI)

#nbr_infest for 2009 is almost all zero. that can't be right

#There are several records where year = 0
abr.2009$beetle_yr[abr.2009$beetle_yr==0]<-2009

#There are  records where year = 2003
abr.2009$beetle_yr[abr.2009$beetle_yr==2003]<-2009

abr.lm<-lm(log10(abr.2009$r_value+1)~log10(abr.2009$nbr_infest+1)+abr.2009$dbh+abr.2009$ht_pitch_tube+abr.2009$CF+abr.2009$SSI)
summary(abr.lm)

########
# 2010 #
########

abr.2010 <- read.csv("data/FRI/R_value_BY2010.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2010$plot_long_dd,abr.2010$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2010$R_VALUE+1))
hist(log10(abr.2010$NBR_INFEST+1))
hist(abr.2010$dbh)
hist(abr.2010$ht_pitch_tube)
hist(abr.2010$CF)
hist(abr.2010$SSI)

abr.lm<-lm(log10(abr.2010$R_VALUE+1)~log10(abr.2010$NBR_INFEST+1)+abr.2010$dbh+abr.2010$ht_pitch_tube+abr.2010$CF+abr.2010$SSI)
summary(abr.lm)

########
# 2011 #
########

abr.2011 <- read.csv("data/FRI/R_value_BY2011.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2011$plot_long_dd,abr.2011$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2011$r_value+1))
hist(log10(abr.2011$nbr_infest+1))
hist(abr.2011$dbh)
hist(abr.2011$ht_pitch_tube)
hist(abr.2011$CF)
hist(abr.2011$SSI)

abr.lm<-lm(log10(abr.2011$r_value+1)~log10(abr.2011$nbr_infest+1)+abr.2011$dbh+abr.2011$ht_pitch_tube+abr.2011$CF+abr.2011$SSI)
summary(abr.lm)

########
# 2012 #
########

abr.2012 <- read.csv("data/FRI/R_value_BY2012.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2012$plot_long_dd,abr.2012$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2012$r_value+1))
hist(log10(abr.2012$nbr_infest+1))
hist(abr.2012$dbh)
hist(abr.2012$ht_pitch_tube)
hist(abr.2012$CF)
hist(abr.2012$SSI)

abr.lm<-lm(log10(abr.2012$r_value+1)~log10(abr.2012$nbr_infest+1)+abr.2012$dbh+abr.2012$ht_pitch_tube+abr.2012$CF+abr.2012$SSI)
summary(abr.lm)

########
# 2013 #
########

abr.2013 <- read.csv("data/FRI/R_value_BY2013.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2013$plot_long_dd,abr.2013$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2013$r_value+1))
hist(log10(abr.2013$nbr_infest+1))
hist(abr.2013$dbh)
hist(abr.2013$ht_pitch_tube)
hist(abr.2013$CF)
hist(abr.2013$SSI)

abr.lm<-lm(log10(abr.2013$r_value+1)~log10(abr.2013$nbr_infest+1)+abr.2013$dbh+abr.2013$ht_pitch_tube+abr.2013$CF+abr.2013$SSI)
summary(abr.lm)

########
# 2014 #
########

abr.2014 <- read.csv("data/FRI/R_value_BY2014.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2014$plot_long_dd,abr.2014$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2014$r_value+1))
hist(log10(abr.2014$nbr_infest+1))
hist(abr.2014$dbh)
hist(abr.2014$ht_pitch_tube)
hist(abr.2014$CF)
hist(abr.2014$SSI)

#too many nbr_infest == 0

#There are several records where year = 0
abr.2014$beetle_yr[abr.2014$beetle_yr==0]<-2014

abr.lm<-lm(log10(abr.2014$r_value+1)~log10(abr.2014$nbr_infest+1)+abr.2014$dbh+abr.2014$ht_pitch_tube+abr.2014$CF+abr.2014$SSI)
summary(abr.lm)

########
# 2015 #
########

abr.2015 <- read.csv("data/FRI/R_value_BY2015.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2015$plot_long_dd,abr.2015$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2015$r_value+1))
hist(log10(abr.2015$nbr_infest+1))
hist(abr.2015$dbh)
hist(abr.2015$ht_pitch_tube)
hist(abr.2015$CF)
hist(abr.2015$SSI)

abr.lm<-lm(log10(abr.2015$r_value+1)~log10(abr.2015$nbr_infest+1)+abr.2015$dbh+abr.2015$ht_pitch_tube+abr.2015$CF+abr.2015$SSI)
summary(abr.lm)

########
# 2016 #
########

abr.2016 <- read.csv("data/FRI/R_value_BY2016.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2016$plot_long_dd,abr.2016$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2016$r_value+1))
hist(log10(abr.2016$nbr_infest+1))
hist(abr.2016$dbh)
hist(abr.2016$ht_pitch_tube)
hist(abr.2016$CF)
hist(abr.2016$SSI)

abr.lm<-lm(log10(abr.2016$r_value+1)~log10(abr.2016$nbr_infest+1)+abr.2016$dbh+abr.2016$ht_pitch_tube+abr.2016$CF+abr.2016$SSI)
summary(abr.lm)

########
# 2017 #
########

abr.2017 <- read.csv("data/FRI/R_value_BY2017.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2017$plot_long_dd,abr.2017$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2017$r_value+1))
hist(log10(abr.2017$nbr_infest+1))
hist(abr.2017$dbh)
hist(abr.2017$ht_pitch_tube)
hist(abr.2017$CF)
hist(abr.2017$SSI)

abr.lm<-lm(log10(abr.2017$r_value+1)~log10(abr.2017$nbr_infest+1)+abr.2017$dbh+abr.2017$ht_pitch_tube+abr.2017$CF+abr.2017$SSI)
summary(abr.lm)

########
# 2018 #
########

abr.2018 <- read.csv("data/FRI/R_value_BY2018.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2018$plot_long_dd,abr.2018$plot_lat_dd)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2018$r_value+1))
hist(log10(abr.2018$nbr_infest+1))
hist(abr.2018$dbh)
hist(abr.2018$ht_pitch_tube)
hist(abr.2018$CF)
hist(abr.2018$SSI)

abr.lm<-lm(log10(abr.2018$r_value+1)~log10(abr.2018$nbr_infest+1)+abr.2018$dbh+abr.2018$ht_pitch_tube+abr.2018$CF+abr.2018$SSI)
summary(abr.lm)

########
# 2019 #
########

abr.2019 <- read.csv("data/FRI/R_value_BY2019.csv", header = T,na.strings=".")

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.2019$plot_long_dd,abr.2019$plot_lat_dd)

#There are some 2020 beetle years that need deleting
abr.2019 <- abr.2019[-c(abr.2019$siteID[abr.2019$beetle_yr==2020]),]

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.2019$r_value+1))
hist(log10(abr.2019$nbr_infest+1))
hist(abr.2019$dbh)
hist(abr.2019$ht_pitch_tube)
hist(abr.2019$CF)
hist(abr.2019$SSI)

abr.lm<-lm(log10(abr.2019$r_value+1)~log10(abr.2019$nbr_infest+1)+abr.2019$dbh+abr.2019$ht_pitch_tube+abr.2019$CF+abr.2019$SSI)
summary(abr.lm)

###########################
# Merge 2006 through 2019 #
###########################

l.2007<-length(abr.2007[,1])
w.2007<-ncol(abr.2007)
abr.2007[,w.2007+1]<-rep(2007,l.2007)

# yr, lat, lon, infest, r, dbh,  ht, cf, ssi : Column numbers varies across years
abr.new<-as.data.frame(rbind(
  as.matrix(abr.2006[,c(8,14,15,11,19,20,21,22,23)]),
  as.matrix(abr.2007[,c(25,7, 8, 9,21,11,12,23,24)]),
  as.matrix(abr.2008[,c(8,15,16,14,27,29,30,31,32)]),
  as.matrix(abr.2009[,c(8,21,22,14,33,35,36,37,38)]),
  as.matrix(abr.2010[,c(9,16,17,15,28,30,31,32,33)]),
  as.matrix(abr.2011[,c(8,21,22,14,33,38,39,40,41)]),
  as.matrix(abr.2012[,c(8,21,22,14,33,35,36,37,38)]),
  as.matrix(abr.2013[,c(8,21,22,14,33,35,36,37,38)]),
  as.matrix(abr.2014[,c(8,21,22,14,33,35,36,37,38)]),
  as.matrix(abr.2015[,c(8,21,22,14,35,37,38,39,40)]),
  as.matrix(abr.2016[,c(7,20,21,13,34,37,38,39,40)]),
  as.matrix(abr.2017[,c(7,20,21,13,34,37,38,39,40)]),
  as.matrix(abr.2018[,c(7,20,21,13,34,37,38,39,40)]),
  as.matrix(abr.2019[,c(7,20,21,13,34,36,37,38,39)])
))
colnames(abr.new)<-c("yr", "lat", "lon", "infest", "r", "dbh", "ht", "cf", "ssi")
rownames(abr.new)<-c(1:length(abr.new[,1]))

win.graph(height=8,width=5)
par(mfrow=c(1,1))
plot(abr.new$lon,abr.new$lat)

sum(abr.new$r==-999)
abr.new$yr[abr.new$r==-999]

win.graph(height=8,width=12)
par(mfrow=c(2,3))
hist(log10(abr.new$r+1))
hist(log10(abr.new$infest+1))
hist(abr.new$dbh)
hist(abr.new$ht)
hist(abr.new$cf)
hist(abr.new$ssi)

win.graph(height=8,width=12)
par(mfrow=c(2,3))
boxplot(log10(abr.new$r+1)~abr.new$yr)
boxplot(log10(abr.new$infest+1)~abr.new$yr)
boxplot(abr.new$dbh~abr.new$yr)
boxplot(abr.new$ht~abr.new$yr)
boxplot(abr.new$cf~abr.new$yr)
boxplot(abr.new$ssi~abr.new$yr)

#verify there are no zero years
sum(abr.new$yr==0)
abr.new[(abr.new$yr==0),]

sum(abr.new$yr==2003)
abr.new[(abr.new$yr==2003),]

abr.lm<-lm(log10(abr.new$r+1)~log10(abr.new$infest+1)+abr.new$dbh+abr.new$ht+abr.new$cf+abr.new$ssi)
summary(abr.lm)

abr.lm<-lm(log10(abr.new$r+1)~
             abr.new$yr+abr.new$lat+abr.new$lon+
             log10(abr.new$infest+1)+abr.new$dbh+abr.new$ht+abr.new$cf+abr.new$ssi)
summary(abr.lm)

abr.early<-abr.new[(abr.new$yr<2016),]
abr.late<-abr.new[(abr.new$yr>2014),]

abr.lm.early<-lm(log10(abr.early$r+1)~
                   abr.early$yr+abr.early$lat+abr.early$lon+
                   log10(abr.early$infest+1)+abr.early$dbh+abr.early$ht+abr.early$cf+abr.early$ssi)
summary(abr.lm.early)

abr.lm.late<-lm(log10(abr.late$r+1)~
                  abr.late$yr+abr.late$lat+abr.late$lon+
                   log10(abr.late$infest+1)+abr.late$dbh+abr.late$ht+abr.late$cf+abr.late$ssi)
summary(abr.lm.late)

####################
# end new analysis #
####################

#Re-plot adjusted values
win.graph(height=5,width=12)
par(mfrow=c(1,3))
hist(abr$ELEV)
hist(abr$DBH)
hist(abr$HT_PITCH_T)

win.graph(height=5,width=12)
par(mfrow=c(1,3))
hist(abr$Beetle_YEAR)
hist(log10(abr$R_VALUE+1))
hist(log10(abr$NBR_INFEST+1))

#how many zeroes if 2006 is removed?
win.graph(height=5,width=12)
par(mfrow=c(1,1))
hist(log10(abr$R_VALUE[abr$Beetle_YEAR!=2006]+1))

#PITCH TUBE HT not used in 2017 report. Is it redundant with INFESTED TREES?
win.graph(height=5,width=12)
par(mfrow=c(1,1))
plot(log10(abr$NBR_INFEST+1),abr$HT_PITCH_T)
cor(log10(abr$NBR_INFEST+1),abr$HT_PITCH_T,use="pairwise.complete.obs")

#Check forest variables
win.graph(height=5,width=12)
par(mfrow=c(2,1))
hist(abr$SSI)
hist(abr$CF)

#Check forest variables, without 2006
win.graph(height=5,width=12)
par(mfrow=c(2,1))
hist(abr$SSI[abr$Beetle_YEAR!=2006])
hist(abr$CF[abr$Beetle_YEAR!=2006])

#Does SSI correlate with other forest variables when SSI is 0? (What does SSI=0 mean?)
win.graph(height=5,width=12)
par(mfrow=c(1,1))
plot(abr$SSI,abr$PINE_PCT)
cor(abr$SSI,abr$PINE_PCT,use="pairwise.complete.obs")

hist(abr$SSI[abr$Beetle_YEAR!=2006])

plot(log10(abr$NBR_INFEST+1),abr$HT_PITCH_T)
cor(log10(abr$NBR_INFEST+1),abr$HT_PITCH_T,use="pairwise.complete.obs")

#build linear response model of R values; use stem level variables and
#complementary cluster or stand level variables
# stand or cluster:   SSI, PINE_PCT, NUMBER INFESTED, CF (climatic factor)
# stem:               PITCH_TUBE_HT, DBH
#
#R-VALUE = SSI, PINE_PCT, NUMBER INFESTED, PITCH_TUBE_HT, DBH
R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$PINE_PCT+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH)
summary(R.lm)

#Try including interactions
R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI*abr$PINE_PCT*abr$CF*log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T*abr$DBH)
summary(R.lm)

#Try including quadratic on density
R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$PINE_PCT+abr$CF+log10(abr$NBR_INFEST+1)+(log10(abr$NBR_INFEST+1))^2+abr$HT_PITCH_T+abr$DBH)
summary(R.lm)

y<-log10(abr$R_VALUE+1)
x<-log10(abr$NBR_INFEST+1)

R.lm2<-lm(y~x+x^2)
summary(R.lm2)
win.graph(height=5,width=12)
par(mfrow=c(1,1))
plot(x,y)
cor(x,y,use="pairwise.complete.obs")
boxplot(y~x)

y<-log10(abr$R_VALUE+1)
x<-abr$DBH

R.lm2<-lm(y~x+x^2)
summary(R.lm2)
win.graph(height=5,width=12)
par(mfrow=c(1,1))
plot(x,y)
cor(x,y,use="pairwise.complete.obs")
boxplot(y~x)

#Develop theme of implosive collapse. Plot key attributes over time
win.graph(height=5,width=12)
par(mfrow=c(2,2))
boxplot(log10(abr$NBR_INFEST+1)~abr$Beetle_YEAR,xlab="year")
boxplot(log10(abr$R_VALUE+1)~abr$Beetle_YEAR,xlab="year")
boxplot(abr$HT_PITCH_T~abr$Beetle_YEAR,xlab="year")
boxplot(sqrt(abr$Qvalue)~abr$Beetle_YEAR,xlab="year")

win.graph(height=5,width=12)
par(mfrow=c(1,1))
boxplot(abr$DBH~abr$Beetle_YEAR,xlab="year")


#File RTC are red tree counts from Mike Undershultz Feb 22, 2023
rtc<-read.table("data/ab/RedTreeCounts.txt",header=T)

win.graph(height=8,width=7)
par(mfrow=c(2,1))
par(mar = c(5, 4, 4, 4) + 0.3)

hist(abr$Beetle_YEAR)

plot(rtc$Year,log10(rtc$RedTrees),xlab="year",ylab="Log10 Red Trees Detected",type="l",col="red",ylim=c(3,7))
points(rtc$Year,log10(rtc$RedTrees),pch=19,col="red",cex=1.2)
axis(side=2,col="red")
par(new=TRUE)
plot(rtc$Year,log10(rtc$TreesControlled),axes=F,type="l",xlab="",ylab="",lwd=2,col="darkgreen",ylim=c(3,7))
points(rtc$Year,log10(rtc$TreesControlled),pch=15,col="darkgreen",cex=1.5)
axis(side=4,col="darkgreen",at=pretty(range(log10(na.omit(rtc$TreesControlled)))))
mtext("Log10 Green Trees Controlled", side = 4, line = 2)

#March 6, 2023 - adding Q values of contorta-pinus introgression
win.graph(height=5,width=12)
par(mfrow=c(1,1))
hist(abr$Qvalue)

R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$PINE_PCT+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH+abr$Qvalue)
summary(R.lm)

#Try including interactions
R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI*abr$PINE_PCT*abr$CF*log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T*abr$DBH*abr$Qvalue)
summary(R.lm)

min(abr$Qvalue)

use<-(abr$Qvalue<0.1)|is.na(abr$Qvalue) #restricts the data to the lodgepole-dominated gene set
sum(use)

length(abr$Qvalue)

R.lm<-lm(log10(abr$R_VALUE[use]+1)~abr$SSI[use]+abr$PINE_PCT[use]+abr$CF[use]+log10(abr$NBR_INFEST[use]+1)+abr$HT_PITCH_T[use]+abr$DBH[use]+abr$Qvalue[use])
summary(R.lm)

R.lm<-lm(log10(abr$R_VALUE[use]+1)~abr$SSI[use]+abr$PINE_PCT[use]+abr$CF[use]+log10(abr$NBR_INFEST[use]+1)+abr$HT_PITCH_T[use]+abr$DBH[use])
summary(R.lm)

R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$PINE_PCT+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH+abr$Qvalue)
summary(R.lm)

R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$PINE_PCT+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T)
summary(R.lm)

R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$PINE_PCT+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH)
summary(R.lm)

R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI*abr$PINE_PCT*abr$CF*log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T*abr$DBH*abr$Qvalue)
summary(R.lm)

R.lm<-lm(log10(abr$R_VALUE[use]+1)~abr$SSI[use]*abr$PINE_PCT[use]*abr$CF[use]*log10(abr$NBR_INFEST[use]+1)*abr$HT_PITCH_T[use]*abr$DBH[use])
summary(R.lm)

win.graph()
par(mfrow=c(1,1))
plot(sqrt(abr$Qvalue),abr$SSI)

########################################
# interpolate spatial data for mapping #
########################################

library(gstat)  # inverse distance weighted, Kriging
library(fields)
library(dplyr)
library(sf)
library(raster)
library(wesanderson)
library(rgdal)

canada<-getData("GADM",country="Canada",level=1)
ab<- canada[match(toupper("Alberta"),toupper(canada$NAME_1)),]

win.graph(height=12,width=20)
par(mfrow=c(2,5))
for (yr in 2006:2015) {

  use<-abr$Beetle_YEAR==yr
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]

  plot(lon,lat,main=paste("lat/lon",yr),xlim=c(min(abr$PLOT_LON),max(abr$PLOT_LON)),ylim=c(min(abr$PLOT_LAT),max(abr$PLOT_LAT)))
}

win.graph(height=12,width=20)
par(mfrow=c(5,2))
for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$R_VALUE))
  hist(sqrt(abr$R_VALUE+1)[use],main=paste("R value",yr))
}

#Plot annual maps of R values
#win.graph(height=12,width=20)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\R.pdf",height=20,width=9)
par(mfrow=c(5,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$R_VALUE))
z<-sqrt(abr$R_VALUE+1)[use]
zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$R_VALUE))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-sqrt(abr$R_VALUE+1)[use]

  #Create z map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on z map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab) #mask IDW to Alberta
  plot(map.ab,col=pal[(round(interp_IDW@data@min):round(interp_IDW@data@max)-interp_IDW@data@min+1)],
       xlab="Longitude", ylab="Latitude",main=paste("R value",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

#Plot annual maps of InfestCluster values
#win.graph(height=12,width=20)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\InfestCluster.pdf",height=17,width=9)
par(mfrow=c(4,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$NBR_INFEST))
z<-10*log10(abr$NBR_INFEST[use]+1)
zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2013) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$NBR_INFEST))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-10*log10(abr$NBR_INFEST[use]+1)

  #Create z map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on z map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab) #mask IDW to Alberta
  plot(map.ab,col=pal[(round(interp_IDW@data@min):round(interp_IDW@data@max)-interp_IDW@data@min+1)],
       xlab="Longitude", ylab="Latitude",main=paste("Trees per cluster",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

#Plot annual maps of HT_PITCH_T
#win.graph(height=12,width=20)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\HTPITCHTUBE.pdf",height=20,width=9)
par(mfrow=c(5,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$HT_PITCH_T))
z<-(abr$HT_PITCH_T+1)[use]
zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$HT_PITCH_T))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-(abr$HT_PITCH_T+1)[use]

  #Create z map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on forest map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab) #mask IDW to Alberta
  plot(map.ab,col=pal[(round(interp_IDW@data@min):round(interp_IDW@data@max)-interp_IDW@data@min+1)],
       xlab="Longitude", ylab="Latitude",main=paste("Pitch Tube Ht.",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

win.graph(height=12,width=20)
par(mfrow=c(2,5))
for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$DBH)&(abr$DBH!=0))
  hist(sqrt(abr$DBH)[use],main=paste("DBH",yr))
}

#Plot annual maps of DBH
#win.graph(height=12,width=20)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\DBH.pdf",height=20,width=9)
par(mfrow=c(5,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$DBH)&(abr$DBH!=0))
z<-(abr$DBH)[use]
zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$DBH)&(abr$DBH!=0))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-(abr$DBH)[use]

  #Create z map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on forest map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab) #mask IDW to Alberta
  plot(map.ab,col=pal[(round(interp_IDW@data@min):round(interp_IDW@data@max)-interp_IDW@data@min+1)],
       xlab="Longitude", ylab="Latitude",main=paste("DBH",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

win.graph(height=12,width=20)
par(mfrow=c(2,5))
for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$Qvalue))
  hist(round(100*abr$Qvalue[use])+1,main=paste("Q value",yr))
}

#Plot annual maps of Q-value
#win.graph(height=12,width=20)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\Q.pdf",height=20,width=9)
par(mfrow=c(5,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$Qvalue))
z<-round(100*(abr$Qvalue[use]))+1
zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2015) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$Qvalue))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-round(100*(abr$Qvalue[use]))+1

  #Create forest z data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 7.5' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on forest map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab) #mask IDW to Alberta
  plot(map.ab,col=pal[(round(interp_IDW@data@min):round(interp_IDW@data@max)-interp_IDW@data@min+1)],
       xlab="Longitude", ylab="Latitude",main=paste("Q value",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

win.graph(height=12,width=20)
par(mfrow=c(2,5))
for (yr in 2006:2013) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$CF))
  hist(abr$CF[use],main=paste("CF",yr))
}

#Plot annual maps of CF
#win.graph(height=12,width=20)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\CF.pdf",height=17,width=9)
par(mfrow=c(4,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$CF))
z<-round(100*(abr$CF[use]))+1
zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2013) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$CF))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-round(100*(abr$CF[use]))+1

  #Create forest map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on forest map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab) #mask IDW to Alberta
  plot(map.ab,col=pal[(round(interp_IDW@data@min):round(interp_IDW@data@max)-interp_IDW@data@min+1)],
       xlab="Longitude", ylab="Latitude",main=paste("CF",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

win.graph(height=12,width=20)
par(mfrow=c(2,5))
for (yr in 2006:2013) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$SSI))
  hist(abr$SSI[use],main=paste("SSI",yr))
}

#Plot annual maps of SSI
#win.graph(height=12,width=20)

pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\SSI.pdf",height=17,width=9)
par(mfrow=c(4,2))
par(mar=c(5,5,3,8))

use<-(!is.na(abr$SSI)&(abr$SSI!=0))
z<-abr$SSI[use]

zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

for (yr in 2006:2013) {
  use<-(abr$Beetle_YEAR==yr)&(!is.na(abr$SSI)&(abr$SSI!=0))
  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]
  z<-abr$SSI[use]

  #Create forest map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on forest map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab)
  plot(map.ab,col=pal[round(interp_IDW@data@min):round(interp_IDW@data@max)],
       xlab="Longitude", ylab="Latitude",main=paste("SSI",yr),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
}
dev.off()

R.lm<-lm(log10(abr$R_VALUE+1)~abr$SSI+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH,na.action=na.exclude)
summary(R.lm)

R.lm.full<-lm(log10(abr$R_VALUE+1)~abr$SSI*abr$CF*log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T*abr$DBH,na.action=na.exclude)
summary(R.lm.full)

R.lm.st<-lm(log10(abr$R_VALUE+1)~abr$Beetle_YEAR+abr$PLOT_LONG+abr$PLOT_LAT+abr$Beetle_YEAR*abr$PLOT_LONG+abr$Beetle_YEAR*abr$PLOT_LAT+abr$PLOT_LONG*abr$PLOT_LAT+abr$SSI+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH,na.action=na.exclude)
summary(R.lm.st)

R.lm.st<-lm(log10(abr$R_VALUE+1)~abr$Beetle_YEAR+abr$PLOT_LONG+abr$PLOT_LAT+abr$Beetle_YEAR*abr$PLOT_LONG+abr$Beetle_YEAR*abr$PLOT_LAT+abr$PLOT_LONG*abr$PLOT_LAT+
              abr$Qvalue+abr$SSI+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH,na.action=na.exclude)
summary(R.lm.st)

R.lm.st<-lm(log10(abr$R_VALUE+1)~abr$Beetle_YEAR+abr$PLOT_LONG+abr$PLOT_LAT+
              abr$Beetle_YEAR*abr$PLOT_LONG+abr$Beetle_YEAR*abr$PLOT_LAT+abr$PLOT_LONG*abr$PLOT_LAT+
            abr$Qvalue+abr$SSI+abr$CF+log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH+
              abr$Beetle_YEAR*abr$PLOT_LONG*abr$HT_PITCH_T+abr$Beetle_YEAR*abr$PLOT_LAT*abr$HT_PITCH_T+abr$PLOT_LONG*abr$PLOT_LAT*abr$HT_PITCH_T
            ,na.action=na.exclude)
summary(R.lm.st)

R.lm.st<-lm(log10(abr$R_VALUE+1)~abr$Beetle_YEAR+abr$PLOT_LONG+abr$PLOT_LAT+
              abr$Beetle_YEAR*abr$PLOT_LONG+abr$Beetle_YEAR*abr$PLOT_LAT+abr$PLOT_LONG*abr$PLOT_LAT+
              abr$Qvalue+abr$SSI+abr$CF+
              log10(abr$NBR_INFEST+1)+abr$HT_PITCH_T+abr$DBH+
              log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T+abr$DBH*abr$HT_PITCH_T+log10(abr$NBR_INFEST+1)*abr$DBH+
              log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T*abr$DBH
            ,na.action=na.exclude)
summary(R.lm.st)

R.lm.st<-lm(log10(abr$R_VALUE+1)~abr$Beetle_YEAR*abr$PLOT_LONG*abr$PLOT_LAT+
              abr$Qvalue*abr$SSI*abr$CF*log10(abr$NBR_INFEST+1)*abr$HT_PITCH_T*abr$DBH
            ,na.action=na.exclude)
summary(R.lm.st)

#plot time-series of residuals
#win.graph(height=5,width=5)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\residuals_ST_ts.pdf",height=8,width=8)
par(mfrow=c(1,1))
boxplot(residuals(R.lm.st)~abr$Beetle_YEAR)
dev.off()


#plot time-series of residuals
win.graph(height=5,width=5)
pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\residuals_ts.pdf",height=8,width=8)
par(mfrow=c(1,1))
boxplot(residuals(R.lm)~abr$Beetle_YEAR)
dev.off()

#Plot a map of residuals
#win.graph(height=12,width=20)

pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\residuals_ST_map.pdf",height=8,width=8)
par(mar=c(5,5,3,8))

use<-(!is.na(residuals(R.lm.st)))

z<-10*(residuals(R.lm.st)[use]+1)+3 #Get residuals in range on color palette

zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

  lat<-abr$PLOT_LAT[use]
  lon<-abr$PLOT_LON[use]

  #Create forest map data structure to be IDW'd
  z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                     class = "data.frame", row.names = c(NA, length(z)))
  sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

  #create Alberta bounding box: these will be the bounds of the IDW raster
  bbox <- c(
    "xmin" = -120,
    "ymin" = 49,
    "xmax" = -110,
    "ymax" = 60
  )

  #create the extrapolation grid
  grd_template <- expand.grid(
    X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
    Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
  )

  #convert grid to raster
  grd_template_raster <- grd_template %>%
    dplyr::mutate(Z = 0) %>%
    raster::rasterFromXYZ(
      crs = 26915) #North American UTM

  #fill raster with IDW based on forest map
  fit_IDW <- gstat::gstat(
    formula = z ~ 1,
    data = as(sf_z.map, "Spatial"),
    nmax = 10, nmin = 3,
    set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
  )
  interp_IDW <- interpolate(grd_template_raster, fit_IDW)

  #output forest raster with plot data overlaid
  map.ab<-mask(interp_IDW,ab)
  plot(map.ab,col=pal[round(interp_IDW@data@min):round(interp_IDW@data@max)],
       xlab="Longitude", ylab="Latitude",main=paste("residuals"),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
  plot(ab,lwd=4,add=T)
  points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
  dev.off()

#Plot a map of residuals
#win.graph(height=12,width=20)

pdf("C:\\Barry\\Data\\MPB\\fRIreportrvalueanalysis\\residuals_map.pdf",height=8,width=8)
par(mar=c(5,5,3,8))

use<-(!is.na(residuals(R.lm)))

z<-10*(residuals(R.lm)[use]+1) #Get residuals in range on color palette

zmin<-min(na.omit(z))
zmax<-max(na.omit(z))
zrange<-zmax-zmin+1

#set color palette range, based on z values
pal <- rev(wes_palette("Zissou1", zrange+1, type = "continuous")) #This scaling of zrange categories is based on min/max z values

lat<-abr$PLOT_LAT[use]
lon<-abr$PLOT_LON[use]

#Create forest map data structure to be IDW'd
z.map <- structure(list(longitude = lon, latitude = lat, data = z), .Names = c("X","Y","Z"),
                   class = "data.frame", row.names = c(NA, length(z)))
sf_z.map <- st_as_sf(z.map, coords = c("X", "Y"), crs = 26915) #North American UTM code, data to match IDW target grid

#create Alberta bounding box: these will be the bounds of the IDW raster
bbox <- c(
  "xmin" = -120,
  "ymin" = 49,
  "xmax" = -110,
  "ymax" = 60
)

#create the extrapolation grid
grd_template <- expand.grid(
  X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.125),
  Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.125) # 30' resolution
)

#convert grid to raster
grd_template_raster <- grd_template %>%
  dplyr::mutate(Z = 0) %>%
  raster::rasterFromXYZ(
    crs = 26915) #North American UTM

#fill raster with IDW based on forest map
fit_IDW <- gstat::gstat(
  formula = z ~ 1,
  data = as(sf_z.map, "Spatial"),
  nmax = 10, nmin = 3,
  set = list(idp = 1.5) # inverse distance power; higher power = more localized contagion; lower power = dispersed
)
interp_IDW <- interpolate(grd_template_raster, fit_IDW)

#output forest raster with plot data overlaid
map.ab<-mask(interp_IDW,ab)
plot(map.ab,col=pal[round(interp_IDW@data@min):round(interp_IDW@data@max)],
     xlab="Longitude", ylab="Latitude",main=paste("residuals"),xlim=c(-120.1,-109),ylim=c(48.9,60.1))
plot(ab,lwd=4,add=T)
points(lon,lat,pch=21,cex=2,col="black",bg=pal[c(z)])
dev.off()

################################################
# spatial autocorrelation in model residuals, z
################################################

library(ncf)

sample.universe<-length(z)
cutoff.high<-6 # no more than 6 degree of lat or lon
cutoff.low<-0.2 # no fewer than 0.2 degree of lat or lon (i.e. exclude ze)

#determine correlogram structure, to prepare for bootstrap
sample.size<-400
random.sample<-round(runif(n=sample.size,min=1,max=sample.universe))
resid.s.correlog<-correlog(lon[random.sample],lat[random.sample],z[random.sample],increment=0.5)
use.dist<-(resid.s.correlog$mean.of.class < cutoff.high)&(resid.s.correlog$mean.of.class > cutoff.low)

distance.classes<-sum(use.dist==TRUE)
nreps<-10
residuals.sum.correlog<-as.data.frame(matrix(0,nrow=distance.classes,ncol=nreps+1))
residuals.sum.correlog[,1]<-resid.s.correlog$mean.of.class[use.dist]

for (reps in 1:nreps) { #bootstrap random sample
 random.sample<-round(runif(n=sample.size,min=1,max=sample.universe))
 resid.s.correlog<-correlog(lon[random.sample],lat[random.sample],z[random.sample],increment=0.5)
 use.dist<-(resid.s.correlog$mean.of.class < cutoff.high)&(resid.s.correlog$mean.of.class > cutoff.low)

 #store correlogram in sum object
 residuals.sum.correlog[,reps+1]<-resid.s.correlog$correlation[use.dist]
}

residuals.sum.correlog.mean<-apply(residuals.sum.correlog[,2:reps+1],1,mean)

win.graph(height=12,width=12)
plot(residuals.sum.correlog[,1],residuals.sum.correlog.mean,ylim=c(-0.5,0.5),xlab="distance (degrees of lat/lon)",
    ylab="spatial autocorrelation in model residuals",main=paste("Correlogram"))
lines(residuals.sum.correlog[,1],residuals.sum.correlog.mean)
points(residuals.sum.correlog[,1],residuals.sum.correlog.mean,pch=21)
abline(h=0)

