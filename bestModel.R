## Funcao para selecao do(s) modelo(s) com maior especificidade, 
## a partir do output do biomod2 
## [obtido por get_evaluations(myBiomodModelOut)]
## Anderson A. Eduardo
## 30/jan/2018

bestModel = function(outputDataRaw, myBiomodModelOut){
    
    ## parametros locais
    evaluationScores = outputDataRaw
    ## statResults = statResults
    myBiomodModelOut = myBiomodModelOut

    ## ss_i = ss_i
    ## sdmType = sdmType

    ## subsetting 
    outputRawTSS = base::subset( x=evaluationScores, subset= c(Eval.metric=='TSS')  )
    outputRawAUC = base::subset( x=evaluationScores, subset=c(Eval.metric=='ROC')  )

    bestModelTSSraw = outputRawTSS[which(outputRawTSS$Specificity == max(outputRawTSS$Specificity)),] 
    bestModelAUCraw = outputRawAUC[which(outputRawAUC$Specificity == max(outputRawAUC$Specificity)),]

    ##modelNames = c(as.vector(bestModelAUCraw$Model.name), as.vector(bestModelTSSraw$Model.name) )
    namePatternsRaw = as.vector(bestModelTSSraw$Model.name) #nao usando mais os melhores modelos pelo AUC    

    ## ## valores de especificidade (maximizando TSS e AUC)
    ## outputRawTSSspec = evaluationScores['TSS','Specificity',,,]
    ## outputRawTSSspec = outputRawTSSspec[complete.cases(outputRawTSSspec),]
    ## outputRawAUCspec = evaluationScores['ROC','Specificity',,,]
    ## outputRawAUCspec = outputRawAUCspec[complete.cases(outputRawAUCspec),]
    
    ## ## maior valor de especificidade para cada modelo (pela maximizacao do TSS e pela maximizacao do AUC)	
    ## TSSspec = as.numeric(apply(outputRawTSSspec, 1, max, na.rm=TRUE))
    ## AUCspec = as.numeric(apply(outputRawAUCspec, 1, max, na.rm=TRUE))
    
    ## ## maior valor de especificidade entre os modelos implementados
    ## tssMax = max(TSSspec)
    ## aucMax = max(AUCspec)
    
    ## ## formacao de vetor com nome e RUN do melhor modelo (especificidade que maximiza o TSS)
    ## bestAlgorithmTSS = statResults[which(statResults$maxTSSspecificity==tssMax),]$model
    ## bestRunTSS = statResults[which(statResults$maxTSSspecificity==tssMax),]$bestModelTSS
    ## patternsTSS = paste(bestRunTSS,bestAlgorithmTSS,sep='_')
    
    ## ## formacao de vetor com nome e RUN do melhor modelo (especificidade que maximiza o AUC)
    ## bestAlgorithmAUC = statResults[which(statResults$maxAUCspecificity==aucMax),]$model
    ## bestRunAUC = statResults[which(statResults$maxAUCspecificity==aucMax),]$bestModelAUC
    ## patternsAUC = paste(bestRunAUC,bestAlgorithmAUC,sep='_')
    
    ## ## nomes dos melhores modelos
    namePatternsRawList = vector()
    namePatternsRaw = gsub(pattern='_PA1',replacement='',x=namePatternsRaw)
    for(modName_i in 1:length(namePatternsRaw)){
        namePatternsRawSplit = strsplit(x=namePatternsRaw[[modName_i]], split='_')
        namePatternsRawList = append(x=namePatternsRawList, values=paste(namePatternsRawSplit[[1]][2],namePatternsRawSplit[[1]][1],sep='_'))
    }
        
    biomodModelsNames = grep(pattern=paste(namePatternsRawList,collapse='|'), x=myBiomodModelOut@models.computed, value=TRUE)
    otherSDMs = grep(pattern=paste(myBiomodModelOut@models.computed,collapse='|'), x=namePatternsRawList , value=TRUE, invert=TRUE)

    modelNames = c(biomodModelsNames, otherSDMs)
    
    ## output da funcao
    return(modelNames)
    
}
