## Script para analise do nicho de uma especie invasora, comparando nicho do ambiente nativo e do ambiente invadido.
## Anderson A. Eduardo
## novembro/2017


##abrindo pacotes necessarios
library(ecospat)

##definindo variaveis e parametros
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
#
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
predictors = predictors[[c('bioclim_01','bioclim_12')]] #selecionando as variaveis usadas
projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
#
spNat = read.csv(paste(mainSampleFolder,'???.csv',sep='')) #dados de ocorrencia ambiente nativo
spNatBG = dismo::randomPoints(??? ,1000) #pontos de fundo para ambiente nativo (background points)
spNatPredictors = extract(x=predictors, y=data.frame(rbind(spNat,spNatBG)), na.rm=TRUE)
spNatData = data.frame(spNat, pres=c(rep(1,length(spNat)), resp(0,spNatBG)), spNatPredictors)
#
spInv = read.csv(paste(mainSampleFolder,'???.csv',sep='')) #dados de ocorrencia ambiente invadido
spInvBG = dismo::randomPoints(??? ,1000) #pontos de fundo para ambiente invadido (background points)
spInvPredictors = extract(x=predictors, y=data.frame(rbind(spInv,spInvBG)), na.rm=TRUE)
spInvData = data.frame(spInv, pres=c(rep(1,length(spInv)), resp(0,spInvBG)), spInvPredictors)
                   

## Analise de nicho

##The PCA is calibrated on all the sites of the study area
pca.env <- dudi.pca(rbind(spNatData,spInvData)[,c('bioclim_01','bioclim_12')],scannf=F,nf=2)
#ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

##PCA scores for the whole study area
scores.globclim <- pca.env$li

##PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,spNatData[which(spNatData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

##PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,spInvData[which(spInvData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

##PCA scores for the whole native study area
scores.clim.nat <-suprow(pca.env,spNatData[,c('bioclim_01','bioclim_12')])$li

##PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,spInvData[,c('bioclim_01','bioclim_12')])$li

##gridding the native niche
grid.clim.nat <-ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.nat, sp=scores.sp.nat, R=100, th.sp=0)

##gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.inv, sp=scores.sp.inv, R=100, th.sp=0)

##equivalencia de nicho
##OBS: Niche equivalency test H1: Is the overlap between the native and invaded niche higher than two random niches
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep=10, alternative="greater")

Dobs = eq.test$obs$D #indice D observado
Iobs = eq.test$obs$I #indice I observado
DpValue = eq.test$p.D #p-valor indice D
IpValue = eq.test$p.I #p-valor indice I


