## Funcao para selecao do(s) modelo(s) com maior especificidade, 
## a partir do output do biomod2 
## [obtido por get_evaluations(myBiomodModelOut)]
## Anderson A. Eduardo
## 30/jan/2018
 
bestModel = function(outputDataRaw, statResults, myBiomodModelOut){
  
	## parametros locais
	evaluationScores = outputDataRaw
	statResults = statResults
	myBiomodModelOut = myBiomodModelOut
  
	## valores de especificidade (maximizando TSS e AUC)
	outputRawTSSspec = evaluationScores['TSS','Specificity',,,]
	outputRawTSSspec = outputRawTSSspec[complete.cases(outputRawTSSspec),]
	outputRawAUCspec = evaluationScores['ROC','Specificity',,,]
  outputRawAUCspec = outputRawAUCspec[complete.cases(outputRawAUCspec),]
  
	## maior valor de especificidade para cada modelo (pela maximizacao do TSS e pela maximizacao do AUC)	
  	TSSspec = as.numeric(apply(outputRawTSSspec, 1, max, na.rm=TRUE))
  	AUCspec = as.numeric(apply(outputRawAUCspec, 1, max, na.rm=TRUE))
  
  	## maior valor de especificidade entre os modelos implementados
  	tssMax = max(TSSspec)
  	aucMax = max(AUCspec)
  
  	## formacao de vetor com nome e RUN do melhor modelo (TSS)
  	bestAlgorithmTSS = statResults[which(statResults$maxTSSspecificity==tssMax),]$model
  	bestRunTSS = statResults[which(statResults$maxTSSspecificity==tssMax),]$bestModelTSS
  	patternsTSS = paste(bestRunTSS,bestAlgorithmTSS,sep='_')
  
  	## formacao de vetor com nome e RUN do melhor modelo (AUC)
  	bestAlgorithmAUC = statResults[which(statResults$maxAUCspecificity==aucMax),]$model
  	bestRunAUC = statResults[which(statResults$maxAUCspecificity==aucMax),]$bestModelAUC
  	patternsAUC = paste(bestRunAUC,bestAlgorithmAUC,sep='_')
  
  	## nomes dos melhores modelos
  	modelNames = grep(pattern=paste(c(patternsTSS,patternsAUC),collapse='|'), x=myBiomodModelOut@models.computed, value=TRUE) 
  
	## output da funcao
  	return(modelNames)
  
}
