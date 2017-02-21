##tutorial de bias layer para SDMs
##fonte: https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/

library(dismo) # interface with MaxEnt
library(raster) # spatial data manipulation
library(MASS) # for 2D kernel density function
library(magrittr) # for piping functionality, i.e., %>%
library(maptools) # reading shapefiles

for(i in 1:nrow(birds)){
temp<-na.omit(temp[,c("lon","lat")])
write.csv(temp, paste0("H:/Species occurrences/birds/",birds$Genus[i], " ", birds$Species[i]," occurrence.csv"),row.names=F)
}

occurdat<-list.files("H:/Species occurrences/Birds",pattern=".csv$",full=T)
locations<-read.csv(occurdat[1],colClasses="numeric")
 
for(i in 2:length(occurdat)){
temp<-read.csv(occurdat[i],colClasses="numeric")
locations<-rbind(locations,temp)
}

climdat<-brick("H:/Bioclim_brick.tif")
occur.ras<-rasterize(locations,climdat,1)
plot(occur.ras)

states<-readShapePoly("H:/Shapefiles/US_States/states.shp",proj4string=occur.ras@crs)
occur.states<-mask(occur.ras,states) %>% crop(.,states)

presences<-which(values(occur.states)==1)
pres.locs<-coordinates(occur.states)[presences,]
 
dens<-kde2d(pres.locs[,1],pres.locs[,2],n=c(nrow(occur.states),ncol(occur.states)))
dens.ras<-raster(dens)
plot(dens.ras)

writeRaster(dens.ras,"H:/Species Occurrences/Bird bias file.tif")

occurrences<-read.csv(occurdat[1])
mod1<-maxent(climdat, occurrences, biasfile=dens.ras)


