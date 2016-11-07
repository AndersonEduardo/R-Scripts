#########################################################################################
####ALGORITMO PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################

library(virtualspecies)
library(maptools)
library(dismo)
library(raster)

###PRIMEIRA PARTE: criando sps virtuais###

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
############################

for (i in 1:length(caminhosCamadasTemp)){
    predictors = stack(paste(caminhosCamadasTemp[i],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[i],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
    predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
    nameScenario = basename(caminhosCamadasTemp[i])
    
    #Especie de clima quente e umido
    parametersHW <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=295,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=2000,p2=3500,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    #Especie de clima quente e seco
    parametersHD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=260,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
 
    #Especie de clima frio e seco
    parametersCD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=50,p2=220,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    spHW <- generateSpFromFun(predictors, parametersHW) #criando a especie artifical (clima quente e umido)
    spHD <- generateSpFromFun(predictors, parametersHD) #criando a especie artifical (clima quente e seco)
    spCD <- generateSpFromFun(predictors, parametersCD) #criando a especie artifical (clima frio e seco)

    auxVector=stack(c(spHW$suitab.raster,spHD$suitab.raster,spCD$suitab.raster))
    names(auxVector) = c('spHW', 'spHD', 'spCD')
    
    for(j in 1:dim(auxVector)[3]){
        projection(auxVector[[j]]) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
        writeRaster(auxVector[[j]], filename=paste(projectFolder,'NichoReal/',names(auxVector[[j]]),'/',nameScenario,'.asc',sep=""), overwrite=T,prj=T) #salvando o raster do mapa da sp
        }
    }

###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo para fazer o pooled niche model###

###Parametros necessarios###
Npass = 1:100 #numero de pontos a serem amostrados para camadas do passado (pensando em pontos fosseis)
Npres = c(10,100,200,400,800) #numero de pontos a serem amostrados para camadas do presente
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
############################

for (h in 1:length(spsTypes)){
    nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[h],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder, full.names=T, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da sp
    for (j in 1:length(Npres)){
        for (k in 1:length(Npass)){
            amostra = data.frame()
            for (i in 1:length(nicheRealPath[1:24])){
                realNicheMap = raster(nicheRealPath[i]) #abrindo o mapa de occ da sp
                #nameScenario = names(realNicheMap)
                if (names(realNicheMap) != 'X000'){ #verificando se e o mapa da ocorrencia no presente
                    amostra_i = sampleOccurrences(realNicheMap, Npass[k],plot=FALSE) #realizando a amostragem
                }
                else{
                    amostra_i = sampleOccurrences(realNicheMap, Npres[j],plot=FALSE) #realizando a amostragem
                }
                amostra = rbind(amostra,amostra_i$sample.points)
            }
            write.csv(amostra,file=paste(mainSampleFolder,spsTypes[h],'/ptsPres',Npres[j],'ptsPass',Npass[k],'.csv',sep=''),row.names=FALSE)#salvando a planilha com os dados da amostra
        }
    }
}

###TERCEIRA PARTE: SDM usando de pontos de ocorrencia em diferentes camadas de tempo (atual a 120 kyr BP)###

#######################################################
##################### MAXENT ##########################
#######################################################

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Amostras/' #caminho para pasta onde a planilha spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
maxentFolder = '/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies

#criando as variaveis e seus niveis
l1=c(rep('=FALSE',0),rep('=TRUE',5))
l2=c(rep('=FALSE',1),rep('=TRUE',4))
l3=c(rep('=FALSE',2),rep('=TRUE',3))
l4=c(rep('=FALSE',3),rep('=TRUE',2))
l5=c(rep('=FALSE',4),rep('=TRUE',1))
feature = data.frame(rbind(l1,l2,l3,l4,l5))
names(feature) = c('hinge','threshold','product','quadratic','linear')
betamultiplier = seq(1,15,0.5)
#

for (i in 1:ncol(feature){
    for (j in 1:length(betamultiplier)){
        for (k in 1:length(spsTypes)){
            sampleFolder = paste(mainSampleFolder,spsTypes[k],sep='') #pasta com os mapas de nicho real da sp
            samplePaths = list.files(path=sampleFolder, full.names=T, pattern='.csv') #lista com os enderecos dos mapas de distribuicao da sp
            sp.occ = read.csv(paste(samplePaths[1],sep=''),header=TRUE) #abrindo a planilha de pontos de occ amostrados
            sp.occ = sp.occ[,1:2] #pegando somente as coorqdenadas na planilha
            names(sp.occ) = c('lon','lat')
            predictors = stack(paste(envVarPaths[1],'/bioclim_01.asc',sep=''),paste(envVarPaths[1],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais do presente (para calibrar o modelo com dados da atalidade)
            predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
            #ajuste do modelo
            me <- maxent(predictors,sp.occ, args=c("responsecurves=TRUE", "outputformat=logistic","randomseed=TRUE","randomtestpoints=25","replicates=10","replicatetype=subsample","outputgrids=FALSE","maximumiterations=5000",'removeduplicates=TRUE','writeclampgrid=TRUE','writemess=TRUE','threads=4','writebackgroundpredictions=TRUE',paste('betamultiplier=',betamultiplier[j],sep=''),paste('linear',feature[i,1],sep=''),paste('quadratic',feature[i,2],sep=''),paste('product',feature[i,3],sep=''),paste('threshold',feature[i,4],sep=''),paste('hinge',feature[i,5],sep='')),path=paste(maxentFolder,spsTypes[k],'/',names(feature)[i],'/',betamultiplier[j],sep='')) # run maxent model with raw output
            #projecoes do modelo
            for (l in 1:length(envVarPaths)){
                nameScenario = basename(envVarPaths[l])
                predictors = stack(paste(envVarPaths[l],'/bioclim_01.asc',sep=''),paste(envVarPaths[l],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
                predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
                crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
                proj <- predict(me@models[[1]],predictors,crs=crs) #rodando a projecao espacial do modelo
                proj = mask(proj,AmSulShape) #recortando as variaveis ambientais
                names(proj) <- nameScenario #ajustando nome
                writeRaster(proj, filename=paste(maxentFolder,spsTypes[k],'/Projections/',names(feature)[i],'/',betamultiplier[j],'/',nameScenario,'.asc',sep=''), overwrite=T) #salvando o mapa de suitability projetado
            }
        }
    }
}        
###############################################
###############################################
###############################################

###QUARTA PARTE: calculando a correlacao entre projecao do SDM e a distribuicao espacial real do nicho da sp###

projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
l1=c(rep('=FALSE',0),rep('=TRUE',5))
l2=c(rep('=FALSE',1),rep('=TRUE',4))
l3=c(rep('=FALSE',2),rep('=TRUE',3))
l4=c(rep('=FALSE',3),rep('=TRUE',2))
l5=c(rep('=FALSE',4),rep('=TRUE',1))
feature = data.frame(rbind(l1,l2,l3,l4,l5))
names(feature) = c('hinge','threshold','product','quadratic','linear')
betamultiplier = seq(1,15,0.5)

for (i in 1:length(spsTypes)){
    nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder,pattern='asc',full.names=T) #lista com os enderecos dos mapas de distribuicao da sp
    for (j in 1:length(betamultiplier)){
        for (k in 1:ncol(feature)){
            projectionsFolder = paste('/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Maxent/',spsTypes[i],'/Projections/',names(feature)[k],'/',betamultiplier[j],sep='') #pasta com as projecoes do cenario
            projectionsPath = list.files(path=projectionsFolder, pattern='asc',full.names=T) #caminhos para os .asc na paste do cenario
            outputData = data.frame(spsTypes=character(0),feature=character(0),betamultiplier=numeric(0),correlation=numeric(0),stringsAsFactors=FALSE) 
            for (l in 1:length(nicheRealPath)){
                realNiche = raster(nicheRealPath[l]) #nicho real
                sdmNiche = raster(projectionsPath[l]) #mapa de suitability gerado por SDM
                layersToCompare = stack(c(realNiche,sdmNiche)) #empilhando os arquivos .asc
                names(layersToCompare) = c('realNiche','sdmNiche') #atualizando os nomes 
                valuesToTest = getValues(layersToCompare) #pegando os dados numericos no .asc
                cor.matrix = as.data.frame(cor(valuesToTest, use="complete.obs")) #teste de correlacao espacial
                dataScenar = c(spsTypes[i],names(feature)[k],betamultiplier[j],cor.matrix[2,1]) #gerando uma linha com as info
                outputData[l,] = dataScenar #gravando numa planilha de output
            }
            write.csv(cor.matrix, file=paste(projectionsFolder,"/correlacao.csv",sep=""),row.names=TRUE) #salvando os dados do cenario
        }
    }
}
