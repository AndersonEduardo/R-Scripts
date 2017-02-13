#selecionando a pasta dos outputs do Maxent
#setwd('/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Maxent/spHW/hinge/1')

# GET EACH MODEL TEST AUC AND CALCULATE TSS

#mres <- read.csv("maxentResults.csv",h=T)

TSSmaxent = function(maxentOutputFolder){
    #setwd(paste(maxentOutputFolder,sep=''))
    mres <- read.csv(paste(maxentOutputFolder,"/maxentResults.csv",sep=''),h=T)
    #mres = read.csv(paste(maxentFolder,spsTypes[i],"/maxentResults.csv",sep=''),h=TRUE)
    for (j in 0:(nrow(mres)-1)) { #loop em que cada iteracao eh para uma linha da planilha maxentResults.csv
        
        trasam <- mres$X.Training.samples[j+1]
        tessam <- mres$X.Test.samples[j+1]
        tesAUC <- mres$Test.AUC[j+1]
        threshold <- mres$X10.percentile.training.presence.logistic.threshold[j+1] #change here the threshold
        
        spred <- read.csv(paste(maxentOutputFolder,"/species_",j,"_samplePredictions.csv",sep=""),h=T)
        
        spred <- spred[which(spred[,3]=="test"),]
        
        a <- nrow(spred[which(spred[,6]>threshold),])
        b <- nrow(spred[which(spred[,6]<threshold),])
        
        bgpred <- read.csv(paste(maxentOutputFolder,"/species_",j,"_backgroundPredictions.csv",sep=""),h=T)
        
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
    return(df)
}	
