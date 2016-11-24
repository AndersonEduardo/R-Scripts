####################################################################################################
##ALGORITMO PARA CALCULAR TSS (TRUE SKILL STATISTICS) A PARTIR DE MAPAS BINARIOS NO FORMATO RASTER##
####################################################################################################

TSSfunction = function(binaryMap,spOccPoints){

binMap = binaryMap #raster('/home/anderson/PosDoc/teste/Maxent/Raster Layers/Caiman yacareBINARIO.asc')
binMap.data = rasterToPoints(binMap)
#CyacareBIN.data = round(CyacareBIN.data,digits=2)

spOccPoints = spOccPoints #read.csv("/home/anderson/PosDoc/dados_ocorrencia/PO_unique/Caiman yacare.csv",header=T)
#CyacarePoints = round(CyacarePoints[,2:3],digits=2)

#linesID = which(CyacareBIN.data[,1] %in% c(CyacarePoints[,1]) & CyacareBIN.data[,2] %in% c(CyacarePoints[,2])) #pegando as linhas de CyacareBin.data (que tem TODAS as coordenadas e dados de ocorrencia do mapa - raster binario output do modelo) em que estao os pontos de ocorrencia (CyaracePoints - dados de ocorrencia da especie)
#presencesPoints = rep(0,nrow(CyacareBIN.data)) #cria um vetor de zeros
#presencesPoints[linesID] = 1 #preenche com "1" as linhas com as coordenadas dos pontos de ocorrencia da especie
#dataSet = data.frame(CyacareBIN.data,CyacarePoints=presencesPoints)
#names(dataSet) = c('lon','lat','forecast','observed')


## dataSet = data.frame(CyacareBIN.data,CyacarePoints=dataVetor)
## names(dataSet) = c('lon','lat','forecast','observed')

##para dataSet ja salvo
dataSet = read.csv('dataSet_TSS.csv',header=T)

##matriz de confusao (MC)
##a = sum(dataSet$observed[1:nrow(dataSet)]==1 & dataSet$forecast[1:nrow(dataSet)]==1)
##b = sum(dataSet$observed[1:nrow(dataSet)]==0 & dataSet$forecast[1:nrow(dataSet)]==1)
##c = sum(dataSet$observed[1:nrow(dataSet)]==1 & dataSet$forecast[1:nrow(dataSet)]==0)
##d = sum(dataSet$observed[1:nrow(dataSet)]==0 & dataSet$forecast[1:nrow(dataSet)]==0)

a = sum(extract(binMap,spOccPoints[,2:3]))
b = sum(binMap.data[,3])-sum(extract(binMap,spOccPoints[,2:3]))
c = nrow(spOccPoints)-sum(extract(binMap,spOccPoints[,2:3]))
d = nrow(binMap.data)-(a+b+c)


MC = matrix(c(a,b,c,d),nrow=2,ncol=2,byrow=TRUE) #matriz de confusao

TSS = (a*d-b*c)/((a+c)*(b+d)) #true skill statistic

return(TSS)
}
