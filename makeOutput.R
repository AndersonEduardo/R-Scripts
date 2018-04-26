## ## Funcao para produzir uma tabela de dados
## ## a partir do output do biomod2
## ## [obtido por get_evaluations(myBiomodModelOut)]
## ## Anderson A. Eduardo
## ## 30/jan/2018

## makeOutput = function(outputDataRaw, data_frame, i_index,j_index, sdm_type, sample_size){

##     ## parametros locais  
##     evaluationScores = outputDataRaw
##     statResults = data_frame
##     i = i_index
##     j = j_index
##     SDM = sdm_type
##     sampleSize = sample_size


##     evaluationScoresTSS = subset( x=evaluationScores, subset=c(Eval.metric=='TSS' ) )
##     evaluationScoresAUC = subset( x=evaluationScores, subset=c(Eval.metric=='ROC' ) )

    

##     maxTSSspecificity = evaluationScoresTSS[ evaluationScoresTSS$Specificity == max(evaluationScoresTSS$Specificity), ] #linha com maior especificidade TSS
##     maxAUCspecificity = evaluationScoresAUC[ evaluationScoresAUC$Specificity == max(evaluationScoresAUC$Specificity), ] #linha com maior especificidade AUC
##     TSSvalue_bestModel = evaluationScoresTSS[ evaluationScoresTSS$Testing.data == max(evaluationScoresTSS$Testing.data), ] #linha com maior valor de TSS
##     AUCvalue_bestModel = evaluationScoresAUC[ evaluationScoresAUC$Testing.data == max(evaluationScoresAUC$Testing.data), ] #linha com maior valor de AUC

##     bestModelAUC = as.character(AUCvalue_bestModel$Model.name)
##     bestModelTSS = as.character(TSSvalue_bestModel$Model.name)

##     ##gsub(pattern=c('MAXENT.Phillips_|GLM_|_PA1'), replacement='', x=as.character(maxTSSspecificity$Model.name))




##     aggregate(x=evaluationScoresTSS$Testing.data, by=list(evaluationScoresTSS$Model.name), FUN=mean)

    
    
##     ##especificidade    
##     outputRawTSSspec = subset( x=evaluationScores, subset=c(Eval.metric=='TSS' ), select=Specificity ) #evaluationScores['TSS','Specificity',,,]
##     outputRawTSSspec = outputRawTSSspec[complete.cases(outputRawTSSspec),]
##     outputRawAUCspec = subset( x=evaluationScores, subset=c(Eval.metric=='ROC' ), select=Specificity ) #evaluationScores['ROC','Specificity',,,]
##     outputRawAUCspec = outputRawAUCspec[complete.cases(outputRawAUCspec),]
##     ##auc e tss
##     outputRawTSSvalue = subset( x=evaluationScores, subset=c(Eval.metric=='TSS' ), select=Testing.data ) #evaluationScores['TSS','Testing.data',,,]
##     outputRawTSSvalue = outputRawTSSvalue[complete.cases(outputRawTSSvalue),]
##     outputRawAUCvalue = subset( x=evaluationScores, subset=c(Eval.metric=='ROC' ), select=Testing.data ) #evaluationScores['ROC','Testing.data',,,]
##     outputRawAUCvalue = outputRawAUCvalue[complete.cases(outputRawAUCvalue),]
    
##     ##maior valor de especificidade de cada algoritmo implementado (tanto para TSS qto para AUC)
##     TSSspec = max(outputRawTSSspec, na.rm=TRUE) #as.numeric(apply(outputRawTSSspec, 1, max, na.rm=TRUE))
##     AUCspec = max(outputRawAUCspec, na.rm=TRUE) #as.numeric(apply(outputRawAUCspec, 1, max, na.rm=TRUE))
    
##     ##tabela auxiliar para obtencao das informacoes do melhor modelo
##     tabBestScoresTSS = data.frame(outputRawTSSspec, bestvalue=TSSspec)
##     tabBestScoresAUC = data.frame(outputRawAUCspec, bestvalue=AUCspec)
    
##     ##vetores vazios (para os nomes dos melhores modelos)
##     bestModelRunTSS = vector()
##     bestModelRunAUC = vector()
##     TSSvalues = vector()
##     AUCvalues = vector()
    
##     for (l in 1:nrow(tabBestScoresTSS)){
##     	bestRunNameTSS = names(tabBestScoresTSS)[match(tabBestScoresTSS$bestvalue[l], tabBestScoresTSS[l, 1:ncol(tabBestScoresTSS)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
##     	TSSvalueBestModel = outputRawTSSvalue[l,bestRunNameTSS] #pegando o TSS do melhor modelo
##     	bestModelRunTSS = append(bestModelRunTSS, bestRunNameTSS) #empilhando os nomes dos melhores modelos
##     	TSSvalues = append(TSSvalues, TSSvalueBestModel)
##     	##
##     	bestRunNameAUC = names(tabBestScoresAUC)[match(tabBestScoresAUC$bestvalue[l], tabBestScoresAUC[l, 1:ncol(tabBestScoresAUC)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
##     	AUCvalueBestModel = outputRawAUCvalue[l,bestRunNameAUC] #pegando o AUC do melhor modelo
##     	bestModelRunAUC = append(bestModelRunAUC, bestRunNameAUC)
##     	AUCvalues = append(AUCvalues, AUCvalueBestModel)
##     }
    
##     ##medias
##     meanTSSValue = as.numeric(apply(outputRawTSSvalue, 1, mean, na.rm=TRUE))
##     meanTSSspecificity = as.numeric(apply(outputRawTSSspec, 1, mean, na.rm=TRUE))
##     meanAUCValue = as.numeric(apply(outputRawAUCvalue, 1, mean, na.rm=TRUE))
##     meanAUCspecificity = as.numeric(apply(outputRawAUCspec, 1, mean, na.rm=TRUE))
    
##     ##gravando estatisticas basicas do melhor modelo de cada algoritmo
##     statResults = rbind(statResults,
##                         data.frame(SDM = paste(SDM),
##                                    sampleSize = sampleSize,
##                                    sp = i,
##                                    ##tss
##                                    model = rownames(tabBestScoresTSS),
##                                    meanTSSValue = meanTSSValue,
##                                    meanTSSspecificity = meanTSSspecificity,
##                                    maxTSSspecificity = TSSspec,
##                                    bestModelTSS = bestModelRunTSS,
##                                    TSSvalue_bestModel= TSSvalues,
##                                    ##auc
##                                    meanAUCValue = meanAUCValue,
##                                    meanAUCspecificity = meanAUCspecificity,
##                                    maxAUCspecificity = AUCspec,
##                                    bestModelAUC = bestModelRunAUC,
##                                    AUCvalue_bestModel= AUCvalues,
##                                    stringsAsFactors = FALSE)
##                         )
    
##     return(statResults)
    
## }
