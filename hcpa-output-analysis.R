## Graficos e analise dos resultados do HCPA_script
## Anderson A. Eduardo
## 07/fev;2019


##ambiente para salvar os graficos
setwd(projectFolder)

##planilhas de dados
statResultsSDMimproved = read.csv(paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), header=TRUE)
statResultsSDMnormal = read.csv(paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), header=TRUE)

# statResultsSDMimproved = read.csv(paste(projectFolder,'/improved','.csv',sep=''), header=TRUE)
# statResultsSDMnormal = read.csv(paste(projectFolder,'/normal','.csv',sep=''), header=TRUE)

## deixando so os cenarios que rodaram para os dois SDMs
dim(statResultsSDMnormal)
dim(statResultsSDMimproved)

## Verificando atraves do grafico de barras simples ##
##igualando os nomes na coluna modelos (estam diferentes)
statResultsSDMimproved$Model.name = gsub(pattern='_AllData|_PA1' , replacement='' , x=as.character(statResultsSDMimproved$Model.name))
statResultsSDMnormal$Model.name = gsub(pattern='_PA1' , replacement='' , x=as.character(statResultsSDMnormal$Model.name))

##contando (pela coluna que optar para ser discriminada)
xxImp = table(statResultsSDMimproved$sampleSize)
xxNor = table(statResultsSDMnormal$sampleSize)

pooledData = cbind(xxImp,xxNor) ##juntando numa matrix (necessario para o barplot)

barplot( t(pooledData), beside=T, las=2, mar=c(12,5,5,5), legend.text=c('HCPA','Normal')) #grafico
####


##obs: usar o merge() para recortar a planilha que for maior
##exemplo:
##newTab=merge(tabB,tabA,by=c('infoA','infoB'))[,-3] ##a terceira coluna e da planilha menor

statResultsSDMnormal = merge(statResultsSDMimproved,
                             statResultsSDMnormal, 
                             all.y=TRUE,
                             by=c('sampleSize','biaslevel','Model.name','Eval.metric'))[,c('SDM.y',
                                                                                           'sampleSize',
                                                                                           'biaslevel',
                                                                                           'Model.name', 
                                                                                           'Eval.metric', 
                                                                                           'Testing.data.y',
                                                                                           "Evaluating.data.y",
                                                                                           "Cutoff.y",
                                                                                           "Sensitivity.y",
                                                                                           "Specificity.y")]


##ajustando so nomes
names(statResultsSDMnormal) = names(statResultsSDMimproved)

modelNames = unique( gsub(pattern='_.*' , replacement='' , x=as.character(statResultsSDMnormal$Model.name)) )


##chacagem rapida
for(i in modelNames){
  print(paste(i, 
              ' ==>> SDMnormal = ', length(statResultsSDMnormal$Model.name[ grep(modelNames[1], statResultsSDMnormal$Model.name) ]),
              ' // ',
              'SDMimproved = ',length(statResultsSDMimproved$Model.name[ grep(modelNames[1], statResultsSDMimproved$Model.name) ]),
              sep=''))
}


##boxplot AUC
jpeg(filename='boxplotAUC.jpeg')
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data, statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='AUC')
dev.off()

idxNorm = which(is.finite(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data))
idxImp = which(is.finite(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data))
jpeg(filename='densidadeAUC.jpeg')
plot(density( statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data[idxNorm] ), ylim=c(0,40),col='blue', lwd=3, xlab='AUC values', main="")
lines(density( statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data[idxImp] ),col='red', lwd=3)
legend(x='topleft', legend=c('SDM normal', 'SDM improved'), lty=1, lwd=3, col=c('blue', 'red'))
dev.off()

##boxplot TSS
jpeg(filename='boxplotTSS.jpeg')
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data, statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS')
dev.off()

idxNorm = which(is.finite(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data))
idxImp = which(is.finite(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data))
jpeg(filename='densidadeTSS.jpeg')
plot(density( statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data[idxNorm] ), ylim=c(0,40),col='blue', lwd=3, xlab='TSS values', main="")
lines(density( statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data[idxImp] ),col='red', lwd=3)
legend(x='topleft', legend=c('SDM normal', 'SDM improved'), lty=1, lwd=3, col=c('blue', 'red'))
dev.off()

##AUC x sample size
jpeg(filename='AUC_&_sampleSize.jpeg')
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='AUC')
tendenciaSDMnormalAUC = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$sampleSize)
abline(tendenciaSDMnormalAUC)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=1.5, pch=20, col=rgb(1,0,0,0.5))
tendenciaSDMimprovedAUC = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$sampleSize)
abline(tendenciaSDMimprovedAUC, col='red')
legend('bottomleft',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(1,0,0,0.5)))
dev.off()

##TSS x sample size
jpeg(filename='TSS_&_sampleSize.jpeg')
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='TSS')
tendenciaSDMnormalTSS = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$sampleSize)
abline(tendenciaSDMnormalTSS)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=1.5, pch=20, col=rgb(1,0,0,0.5))
tendenciaSDMimprovedTSS = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$sampleSize)
abline(tendenciaSDMimprovedTSS, col='red')
legend('bottomright',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(1,0,0,0.5)))
dev.off()


##teste de significancia

wilcox.test(statResultsSDMimproved$AUCvalue_bestModel,statResultsSDMnormal$AUCvalue_bestModel) #resultado: p<0.05

wilcox.test(statResultsSDMimproved$TSSvalue_bestModel,statResultsSDMnormal$TSSvalue_bestModel)  #resultado:p<<0.05


## TSS e AUC por algoritmo

jpeg(filename='TSS_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Model.name, ylim=c(0,1), ylab='TSS', main='SDM normal')
boxplot(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Model.name, ylim=c(0,1), ylab='TSS', main='SDM improved')
dev.off()

jpeg(filename='AUC_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(12,5,5,1))
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Model.name, ylim=c(0,1), ylab='AUC', main='SDM normal')
boxplot(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Model.name, ylim=c(0,1), ylab='AUC', main='SDM improved')
dev.off()


##teste TSS
kruskal.test(TSSvalue_bestModel ~ model, data=statResultsSDMnormal) #resultado: p>0.05
kruskal.test(TSSvalue_bestModel ~ model, data=statResultsSDMimproved) #resultado: p<<0.05

##teste AUC
kruskal.test(AUCvalue_bestModel ~ model, data=statResultsSDMnormal) #resultado: p<0.05
kruskal.test(AUCvalue_bestModel ~ model, data=statResultsSDMimproved) #resultado: p<<0.05


## especificidade (escolha do melhor modelo)

jpeg(filename='especificidade_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Specificity ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Model.name, ylab='Specificity', main='Maximization of AUC')
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Specificity ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Model.name, ylab='Specificity', main='Maximization of TSS')
dev.off()

#teste especificade
kruskal.test(maxTSSspecificity ~ model, data=statResultsSDMnormal) #rsultado: p>0.05
kruskal.test(maxAUCspecificity ~ model, data=statResultsSDMnormal) #resultado:p>0.05


## niveis de vies

jpeg(filename='aucXbiasleval.jpeg', width=800) #AUC
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$biaslevel, ylim=c(0,1), xlab='Bias level', ylab='AUC', pch=19, col=rgb(0,0,0,0.5))
tendenciaNorm = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$biaslevel)
abline(tendenciaNorm)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$biaslevel, ylim=c(0,1), pch=19, col=rgb(1,0,0,0.5))
tendenciaImprov = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$biaslevel)
abline(tendenciaImprov, col='red')
dev.off()

jpeg(filename='tssXbiasleval.jpeg', width=800) #TSS
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$biaslevel, ylim=c(0,1), xlab='Bias level', ylab='TSS', pch=19, col=rgb(0,0,0,0.5))
tendenciaNorm = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$biaslevel)
abline(tendenciaNorm)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$biaslevel, ylim=c(0,1), pch=19, col=rgb(1,0,0,0.5))
tendenciaImprov = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$biaslevel)
abline(tendenciaImprov, col='red')
dev.off()
