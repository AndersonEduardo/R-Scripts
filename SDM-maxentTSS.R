library(raster)

options(java.parameters = "-Xmx1g" )
library(dismo)
library(SDMTools)
library(sm)



rm(list=ls())
setwd("/Users/pabloriul/Desktop/new2/")	
dir.create("data/models")
dir.create("data/models/biased")
dir.create("data/models/unbiased")

predictors <- stack(list.files(path = "data/myExpl", pattern='asc',full.names=T))
########## creates object with all .csv occurence files####################################################################
occ.sps <- list.files('data/points',pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))


#######Sampling bias surface Grid following Fitzpatrick et al. 2013#########

pts.old <- read.csv('data/SDMinput.csv',h=T)
mask <- predictors[[1]]>-1000
bias <- cellFromXY(mask, pts.old[,-1])
cells <- unique(sort(bias))
kernelXY <- xyFromCell(mask, cells)
samps <- as.numeric(table(bias))

# code to make KDE raster
KDEsur <- sm.density(kernelXY, weights=samps, display="none", nbins=0)
KDErast=SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
KDErast = SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, length(KDEsur$estimate))))
KDErast <- raster(KDErast)
KDErast <- resample(KDErast, mask)
KDErast <- KDErast*mask
names(KDErast) <- "layer"
KDEpts <- rasterToPoints(KDErast)



#######Modeling

nmodels <- 100


Sys.setenv(NOAWT=TRUE)
bg=data.frame(KDEpts[sample(seq(1:nrow(KDEpts)), size=1000, replace=F),1:2])
size <- nlayers(predictors)+2
bg[,3:size] <- extract(predictors,bg)
u <- as.data.frame(predictors,xy=T)
names(bg) <- names(u)
vars <- names(bg)[3:size]

write.csv(bg,"/Users/pabloriul/Desktop/new2/bg-biased.csv",row.names=F)


for (i in 1:length(occ.sps)){	
	

  setwd("/Users/pabloriul/Desktop/new2/")	
  sp.file <- read.csv(paste('data/points/', occ.sps[i],sep=""),h=T)
  sp.occ <- sp.file[,2:3]
  sp.occ[,3:size] <- extract(predictors,sp.occ)
  names(sp.occ) <- names(bg)
  train <- rbind(sp.occ, bg)[,vars]
  p <- c(rep(1, nrow(sp.occ)), rep(0, nrow(bg)))

  me <- maxent(x=train, p=p, args=c("-p", "-a", "nothreshold","-m", 1000, "writebackgroundpredictions", "threads=4"
  ,"randomtestpoints=25", paste("replicates=",nmodels,sep=""),"replicatetype=bootstrap", "randomseed"), 
  path=paste(getwd(),'/data/models/biased/', sp.file$Species[1], "/", sep="")) 
  
  setwd(paste(getwd(),'/data/models/biased/', sp.file$Species[1], "/", sep=""))
  pred <- predict(me, predictors) #filename=paste(sp.file$Species[1],".asc",sep="")
  mpred <- mean(pred)
  sdpred <- calc(pred,sd)
  
######### GET EACH MODEL TEST AUC AND CALCULATE TSS

	mres <- read.csv("maxentResults.csv",h=T)

	for (j in 0:paste(nmodels-1)) {

	trasam <- mres$X.Training.samples[j+1]
	tessam <- mres$X.Test.samples[j+1]
	tesAUC <- mres$Test.AUC[j+1]
	threshold <- mres$X10.percentile.training.presence.logistic.threshold[j+1] #change here the threshold

	spred <- read.csv(paste("species_",j,"_samplePredictions.csv",sep=""),h=T)

	spred <- spred[which(spred[,3]=="test"),]

	a <- nrow(spred[which(spred[,6]>threshold),])
	b <- nrow(spred[which(spred[,6]<threshold),])

	bgpred <- read.csv(paste("species_",j,"_backgroundPredictions.csv",sep=""),h=T)

	c <- nrow(bgpred[which(bgpred[,5]>threshold),])
	d <- nrow(bgpred[which(bgpred[,5]<threshold),])

	#sen <- a/(a+b)
	#spe <- d/(c+d)

	TSS <- (a/(a+b))+(d/(c+d))-1

		if (j == 0) {
		 df <- data.frame(trasam,tessam,tesAUC,TSS)
		 } else {
		 df2 <- data.frame(trasam,tessam,tesAUC,TSS)
		 df <- rbind(df,df2)
		 }



}
  
## Making predctions to a raster based in the full model  
  
 
#  me2 <- maxent(x=train, p=p, args=c("-q", "-p", "-a", "nothreshold","-m", 1000, "writebackgroundpredictions", "threads=4"
#  ,"randomtestpoints=0", "replicates=1","replicatetype=bootstrap", "randomseed"),)
  
#  crs <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
#  r <- predict(me2,predictors,crs=crs)
#  names(r) <- splist[i]
  setwd("/Users/pabloriul/Desktop/new2/data/models/biased")
#  writeRaster(r,filename=paste(sp.file$Species[1], ".asc", sep=""),overwrite=T)
  write.csv(df,paste(sp.file$Species[1], "_AUCTSS.csv", sep=""),row.names=F)
  writeRaster(mpred,filename=paste(sp.file$Species[1], "-mean.asc", sep=""),overwrite=T)
  writeRaster(sdpred,filename=paste(sp.file$Species[1], "-sd.asc", sep=""),overwrite=T)  
  rm(me);rm(pred);rm(mpred);rm(sdpred)
  cat("Modeled", i ,"species of", paste(length(occ.sps)))
  cat("\n")
}



#TGBackground




Sys.setenv(NOAWT=TRUE)
bg=data.frame(KDEpts[sample(seq(1:nrow(KDEpts)), size=1000, replace=F,prob=KDEpts[,"layer"]),1:2])
size <- nlayers(predictors)+2
bg[,3:size] <- extract(predictors,bg)
u <- as.data.frame(predictors,xy=T)
names(bg) <- names(u)
vars <- names(bg)[3:size]

write.csv(bg,"/Users/pabloriul/Desktop/new2/bg-unbiased.csv",row.names=F)


for (i in 1:length(occ.sps)){	
	

  setwd("/Users/pabloriul/Desktop/new2/")	
  sp.file <- read.csv(paste('data/points/', occ.sps[i],sep=""),h=T)
  sp.occ <- sp.file[,2:3]
  sp.occ[,3:size] <- extract(predictors,sp.occ)
  names(sp.occ) <- names(bg)
  train <- rbind(sp.occ, bg)[,vars]
  p <- c(rep(1, nrow(sp.occ)), rep(0, nrow(bg)))

  me <- maxent(x=train, p=p, args=c("-p", "-a", "nothreshold","-m", 1000, "writebackgroundpredictions", "threads=4"
  ,"randomtestpoints=25", paste("replicates=",nmodels,sep=""),"replicatetype=bootstrap", "randomseed"), 
  path=paste(getwd(),'/data/models/unbiased/', sp.file$Species[1], "/", sep="")) 
  
  setwd(paste(getwd(),'/data/models/unbiased/', sp.file$Species[1], "/", sep=""))
  pred <- predict(me, predictors) #filename=paste(sp.file$Species[1],".asc",sep="")
  mpred <- mean(pred)
  sdpred <- calc(pred,sd)
  
######### GET EACH MODEL TEST AUC AND CALCULATE TSS

	mres <- read.csv("maxentResults.csv",h=T)

	for (j in 0:paste(nmodels-1)) {

	trasam <- mres$X.Training.samples[j+1]
	tessam <- mres$X.Test.samples[j+1]
	tesAUC <- mres$Test.AUC[j+1]
	threshold <- mres$X10.percentile.training.presence.logistic.threshold[j+1] #change here the threshold

	spred <- read.csv(paste("species_",j,"_samplePredictions.csv",sep=""),h=T)

	spred <- spred[which(spred[,3]=="test"),]

	a <- nrow(spred[which(spred[,6]>threshold),])
	b <- nrow(spred[which(spred[,6]<threshold),])

	bgpred <- read.csv(paste("species_",j,"_backgroundPredictions.csv",sep=""),h=T)

	c <- nrow(bgpred[which(bgpred[,5]>threshold),])
	d <- nrow(bgpred[which(bgpred[,5]<threshold),])

	#sen <- a/(a+b)
	#spe <- d/(c+d)

	TSS <- (a/(a+b))+(d/(c+d))-1

		if (j == 0) {
		 df <- data.frame(trasam,tessam,tesAUC,TSS)
		 } else {
		 df2 <- data.frame(trasam,tessam,tesAUC,TSS)
		 df <- rbind(df,df2)
		 }



}
  
## Making predctions to a raster based in the full model  
  
 
 # me2 <- maxent(x=train, p=p, args=c("-q", "-p", "-a", "nothreshold","-m", 1000, "writebackgroundpredictions", "threads=4","randomtestpoints=0", "replicates=1","replicatetype=bootstrap", "randomseed"),)
  
 # crs <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
 # r <- predict(me2,predictors,crs=crs)
 # names(r) <- splist[i]
  setwd("/Users/pabloriul/Desktop/new2/data/models/unbiased")
 # writeRaster(r,filename=paste(sp.file$Species[1], "unb.asc", sep=""),overwrite=T)
  write.csv(df,paste(sp.file$Species[1], "_unbiased_AUCTSS.csv", sep=""),row.names=F)
  writeRaster(mpred,filename=paste(sp.file$Species[1], "-mean.asc", sep=""),overwrite=T)
  writeRaster(sdpred,filename=paste(sp.file$Species[1], "-sd.asc", sep=""),overwrite=T)  
  rm(me);rm(pred);rm(mpred);rm(sdpred)
  cat("Modeled", i ,"species of", paste(length(occ.sps)))
  cat("\n")

}

#T test to compare models

setwd("/Users/pabloriul/Desktop/new2/")	

occ.sps <- list.files('data/models/biased',pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

occ.sps2 <- list.files('data/models/unbiased',pattern="csv")
splist2 <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

setwd("/Users/pabloriul/Desktop/new2/")	

for (l in 1:length(occ.sps)) {

biased <- read.csv(paste('data/models/biased/', occ.sps[l],sep=""),h=T)
unbiased <- read.csv(paste('data/models/unbiased/', occ.sps2[l],sep=""),h=T)
tauc <- wilcox.test(biased[,3],unbiased[,3])
ttss <- wilcox.test(biased[,4],unbiased[,4])

if (l == 1) {
	
		 # table <- data.frame(paste(splist[l]),mean(biased$trasam),mean(biased$tessam),mean(biased$tesAUC),sd(biased$tesAUC),mean(unbiased$tesAUC),sd(unbiased$tesAUC),tauc		 		 $statistic,tauc$parameter,tauc$p.value,mean(biased$TSS),sd(biased$TSS),mean(unbiased$TSS),sd(unbiased$TSS),ttss$statistic,ttss$parameter,ttss$p.value)
		 # names(table) <-c("Species","TrainingSamps","TestingSamps","B_AUC","B_AUCsd","UB_AUC","UB_AUCsd","T_val","Df","P_val","B_TSS","B_TSSsd","UB_TSS","UB_TSSsd","T_val","Df","P_val")
	
	
		 table <- data.frame(paste(splist[l]),mean(biased$trasam),mean(biased$tessam),mean(biased$tesAUC),sd(biased$tesAUC),mean(unbiased$tesAUC),sd(unbiased$tesAUC),tauc		 		 		 $statistic,tauc$p.value,mean(biased$TSS),sd(biased$TSS),mean(unbiased$TSS),sd(unbiased$TSS),ttss$statistic,ttss$p.value)
		 names(table) <- c("Species","TrainingSamps","TestingSamps","B_AUC","B_AUCsd","UB_AUC","UB_AUCsd","W_val","P_val","B_TSS","B_TSSsd","UB_TSS","UB_TSSsd","W_val","P_val")

		 } else {
		
		 # table2 <- data.frame(paste(splist[l]),mean(biased$trasam),mean(biased$tessam),mean(biased$tesAUC),sd(biased$tesAUC),mean(unbiased$tesAUC),sd(unbiased$tesAUC),tauc		 		 $statistic,tauc$parameter,tauc$p.value,mean(biased$TSS),sd(biased$TSS),mean(unbiased$TSS),sd(unbiased$TSS),ttss$statistic,ttss$parameter,ttss$p.value)
		 # names(table2) <- c("Species","TrainingSamps","TestingSamps","B_AUC","B_AUCsd","UB_AUC","UB_AUCsd","T_val","Df","P_val","B_TSS","B_TSSsd","UB_TSS","UB_TSSsd","T_val","Df","P_val")
	
	
		 table2 <- data.frame(paste(splist[l]),mean(biased$trasam),mean(biased$tessam),mean(biased$tesAUC),sd(biased$tesAUC),mean(unbiased$tesAUC),sd(unbiased$tesAUC),tauc		 		 		 $statistic,tauc$p.value,mean(biased$TSS),sd(biased$TSS),mean(unbiased$TSS),sd(unbiased$TSS),ttss$statistic,ttss$p.value)
		 names(table2) <- c("Species","TrainingSamps","TestingSamps","B_AUC","B_AUCsd","UB_AUC","UB_AUCsd","W_val","P_val","B_TSS","B_TSSsd","UB_TSS","UB_TSSsd","W_val","P_val")

	
		 table <- rbind(table,table2)
		 }
		 
		 tableFinal <- cbind(table[,1],round(table[,2:15],3))
		 names(tableFinal) <- c("Species","TrainingSamps","TestingSamps","B_AUC","B_AUCsd","UB_AUC","UB_AUCsd","W_val","P_val","B_TSS","B_TSSsd","UB_TSS","UB_TSSsd","W_val","P_val")

		 write.csv(tableFinal,"modelevaluations.csv",row.names=F)
}


#Creating splist files
library(raster)


setwd("/Users/pabloriul/Desktop/new2/") 



# PREPARING ZONATION FILES

splist <- data.frame(list.files(path = "data/models/biased", pattern='-mean.asc',full.names=T))


c1 <- rep(1,nrow(splist))
c2 <- rep(1,nrow(splist))
c3 <- rep(1,nrow(splist))
c4 <- rep(1,nrow(splist))
c5 <- rep(0.25,nrow(splist))

splist <- data.frame(list.files(path = "data/models/biased", pattern='-mean.asc',full.names=T))
L <- cbind(c1,c2,c3,c4,c5,splist)
write.table(L,file="data/models/splistBiased.spp",row.names=F,col.names=F,sep="\t",quote=F)
write.table("call zig3.exe -r data/models/set-biased.dat data/models/splistBiased.spp results/biased/biased.txt 0.0 0 1.0 1",file="do_biased.bat",quote=F,row.names=F,col.names=F)

splist2 <- data.frame(list.files(path = "data/models/biased", pattern='-sd.asc',full.names=T))
L <- cbind(c1,splist2)
write.table(L,file="data/models/splistBiased2.spp",row.names=F,col.names=F,sep="\t",quote=F)




splist <- data.frame(list.files(path = "data/models/unbiased", pattern='-mean.asc',full.names=T))
L <- cbind(c1,c2,c3,c4,c5,splist)
write.table(L,file="data/models/splistUnbiased.spp",row.names=F,col.names=F,sep="\t",quote=F)
write.table("call zig3.exe -r data/models/set-unbiased.dat data/models/splistUnbiased.spp results/unbiased/unbiased.txt 0.0 0 1.0 1",file="do_unbiased.bat",quote=F,row.names=F,col.names=F)

splist2 <- data.frame(list.files(path = "data/models/unbiased", pattern='-sd.asc',full.names=T))
L <- cbind(c1,splist2)
write.table(L,file="data/models/splistUnbiased2.spp",row.names=F,col.names=F,sep="\t",quote=F)


dir.create("results")
dir.create("results/biased")
dir.create("results/unbiased")
