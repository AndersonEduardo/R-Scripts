library(virtualspecies)

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = maptools::readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul

###variaveis ambientais da area de estudos###
predictors = stack(paste(caminhosCamadasTemp[1],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[1],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais

#tempRange = seq( predictors$bioclim_01@data@min, predictors$bioclim_01@data@max )
#precRange = seq( predictors$bioclim_12@data@min, predictors$bioclim_12@data@max )

###funcoes para 'localizacao' e 'largura' do nicho (uma funcao para cada variavel ambiental)###
larguraDistribTemp = abs(rnorm(100,100,100))
localDistribTemp = rnorm(100,200,20)

larguraDistribPrecip = abs(rnorm(100,10,500))
localDistribPrecip = 4000-abs(rnorm(100,0,1000))

###criando as especies artificiais###
distList = list() #objeto para armazenar uma lista com as especies

for (i in 1:5){ #loop para o numero de especies a serem criadas
    parameterSp <- formatFunctions(predictors,bioclim_01=c(fun='dnorm',mean=localDistribTemp[i],sd=larguraDistribTemp[i]),bioclim_12=c(fun='dnorm',mean=localDistribPrecip[i],sd=larguraDistribPrecip[i])) #criando as respostas da especie às variaveis ambientais
    sp_i <- generateSpFromFun(predictors, parameterSp) #criando a especie artifical (clima quente e umido)
    distList = append(distList,sp_i$suitab.raster)
}

###empilhando e tranformando em binario os mapas das distribuicoes 
disListBIN = raster::stack(distList) > 0.1

plot(sum(disListBIN))