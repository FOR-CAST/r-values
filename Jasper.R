####################################################
# Exploratory analysis of Jasper NP r values data  #
# published as Cooke et al 2024                    #
# Barry J. Cooke                                   #
# June 18, 2024                                    #
# last change: June 25, 2024                       #
####################################################

# r-value: springtime recruitment rate
# R-value: interannual rate of change in infestation rate (# trees or area)

##########################
#   area infested data   #
##########################

infile<-"data/Brett/JasperData/AreaInfested20132021.txt"
Jasper.AreaInfested<-read.table(infile,header=T)
win.graph(height=5,width=5)
plot(Jasper.AreaInfested$year, Jasper.AreaInfested$area,type="l",xlab="year",ylab="area infested (ha)",main="Jasper National Park")
points(Jasper.AreaInfested$year, Jasper.AreaInfested$area,pch=19)

#plot interannual change in area infested, R
Jasper.R<-Jasper.AreaInfested$area[2:length(Jasper.AreaInfested$area)]/Jasper.AreaInfested$area[1:(length(Jasper.AreaInfested$area)-1)]
win.graph(height=5,width=5)
hist(Jasper.R)

#plot on log scale
win.graph(height=6,width=6)
plot(2014:2021,log10(Jasper.R),type="l",xlab="year",ylab="log10 (Area t/Area t-1)",main="Jasper National Park")
