library(adehabitatHS)
library(raster)

projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
evaluation = list()
TSSvector = data.frame()
sampleSizes = c(5,15,25,35,45,55,65,75,85,95)
NumRep = 10 #numero de replicas (de cada cenario amostral)

i=2;j=55;k=1;l=1

occPointsRaw = read.csv(paste(mainSampleFolder,spsTypes[i],'/occ',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
occPointsRounded = round(occPointsRaw,digits=2)
occPoints = data.frame(occPointsRounded[,1:2],data=seq(nrow(occPointsRounded)))
occSPDS = SpatialPointsDataFrame(coords = occPoints[,1:2],data=occPoints,proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))

predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
projection(predictors) =  CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
#crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
predSPDF = as(predictors,"SpatialPixelsDataFrame")

pc = dudi.pca(predSPDF,scannf = FALSE)

pr <- slot(count.points(occSPDS, predSPDF), "data")[,1]
#pres = data.frame(pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD

(enfa1 <- enfa(pc,pr,scannf = FALSE))


####
data(lynxjura)
map <- lynxjura$map
locs <- lynxjura$locs
locs <- locs[slot(locs, "data")[,2]!="D",]
slot(map,"data")[,4] <- sqrt(slot(map,"data")[,4])
tab <- slot(map, "data")
pr <- slot(count.points(locs, map), "data")[,1]

pc <- dudi.pca(tab, scannf = FALSE)

(enfa1 <- enfa(pc, pr,scannf = FALSE))

hist(enfa1)
hist(enfa1, scores = FALSE, type = "l")
scatter(enfa1)

## randomization test
## Not run: 
(renfa <- randtest(enfa1))
plot(renfa)

