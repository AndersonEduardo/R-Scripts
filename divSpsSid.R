library(virtualspecies)

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = maptools::readShapePoly("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul

###variaveis ambientais da area de estudos###
predictors = stack(paste(caminhosCamadasTemp[1],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[1],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais

#tempRange = seq( predictors$bioclim_01@data@min, predictors$bioclim_01@data@max )
#precRange = seq( predictors$bioclim_12@data@min, predictors$bioclim_12@data@max )

###funcoes para 'localizacao' e 'largura' do nicho (uma funcao para cada variavel ambiental)###
larguraDistribTemp = abs(rnorm(50,500,100))
localDistribTemp = rnorm(50,250,50)
larguraDistribPrecip = abs(rnorm(100,1000,500))
localDistribPrecip = abs(rnorm(100,4000,1000))

###criando as especies artificiais###
distList = list() #objeto para armazenar uma lista com as especies

for (i in 1:5){ #loop para o numero de especies a serem criadas
    
    parameterSp <- formatFunctions(
        predictors,
        bioclim_01=c(fun='dnorm',
                     mean=localDistribTemp[i],
                     sd=larguraDistribTemp[i]),
        bioclim_12=c(fun='dnorm',
                     mean=localDistribPrecip[i],
                     sd=larguraDistribPrecip[i]) )#criando as respostas da especie às variaveis ambientais
    
    sp_i <- generateSpFromFun(raster.
                              stack=predictors,
                              parameters=parameterSp) #criando a especie artifical (clima quente e umido)
    
    distList = append(distList,sp_i$suitab.raster)
    
}

###empilhando e tranformando em binario os mapas das distribuicoes 
disListBIN = raster::stack(distList) > 0.1

#dev.off()

jpeg('/home/anderson/Documentos/Projetos/divSpsSid/div100.jpg',width=1200,height=1200)
plot(sum(disListBIN));grid()
plot(AmSulShape,add=T)
dev.off()

mapaRiq = sum(disListBIN)
writeRaster(mapaRiq,'/home/anderson/Documentos/Projetos/divSpsSid/div100.asc',overwrite=TRUE)


x1 = getValues(predictors[[1]]) #extract(predictors[[1]],disListBIN,mean)
x2 = getValues(predictors[[2]]) #extract(predictors[[1]],disListBIN,mean)
y1 = getValues(mapaRiq)

jpeg('/home/anderson/Documentos/Projetos/divSpsSid/corRiqVars.jpg',width=900,height=400)
par(mfrow=c(1,2))
plot(y1 ~ x1, xlab='temperatura média anual',ylab='riqueza de especies',pch=20,col=rgb(0,0,0,0.5))
plot(y1 ~ x2, xlab='precipitação anual',ylab='riqueza de especies',pch=20,col=rgb(0,0,0,0.5))
dev.off()

mod = lm(y1 ~ x1)
summary(mod)

###########################
### Framework da matriz ###
###########################

##abrindo pacotes
library(raster)

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul


##TESTANDO COM VARIAVEL AMBIENTAL SIMULADA##

##criando gridfile da variavel ambiental simlada

r <- raster(ncols=10, nrows=10) #criando a estrutura de um raster
r[] <- runif(ncell(r)) #adicionando valores ao raster
crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84" #ajustando a projecao
plot(r) #inspecao visual

datMat = as.data.frame(r,xy=TRUE) #transformando raster em data.frame
names(datMat) = c('lon','lat','var') #ajustando os nomes das colunas do data.frame

##equacoes para as dimensoes do nicho das especies

Nsp = 10 #numero de especies a serem criadas
beta = runif(n=Nsp,min=8,max=12) #parametro para cada equacao de cada especie
alpha = runif(n=Nsp,min=0.4,max=0.6) #parametro para cada equacao de cada especie
x1 = datMat$var #variavel ambiental

##solucao numerica para a equacao de cada especie

for (i in 1:Nsp){
    fSp_i = as.integer(1/(1+exp(-lambda[i]*(x1-alpha[i]))) > 0.5) #solucao da equacao com output binario
    ##fSp_i = 1/(1+exp(-beta[i]*(x1-alpha[i]))) #solucao da equacao com output continuo ("suitability")
    datMat = data.frame(cbind(datMat,fSp_i=fSp_i)) #adicionando ao data.frame
    names(datMat)[ncol(datMat)] = paste('sp',i,sep='') #ajustando os nomes das especies no data.farme
    ##salvando graficos das equacoes de cada especie
    jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','function_sp',i,'.jpeg',sep='')) 
    plot(fSp_i~x1,xlab='var 1',ylab='Suitability')
    dev.off()
}

##criando a coluna de riqueza de especies##

datMat = cbind(datMat, Richness=rowSums(datMat[,grep('sp',names(datMat))]))

##raster da distribuicao modelada

for(i in 1:Nsp){
    SpDist = datMat[,c('lon','lat',paste('sp',i,sep=''))] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
    coordinates(SpDist) = ~lon + lat #definindo colunas das coordenadas
    gridded(SpDist) = TRUE #definindo gridded
    proj4string(SpDist) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definindo proj
    rasterSpDist = raster(SpDist) #criando objeto raster
    ##criando imagem da distribuicao de cada especie
    jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','sp',i,'.jpeg',sep=''))
    plot(rasterSpDist)
    dev.off()
}

##rater da riqueza

SpRic = datMat[,c('lon','lat','Richness')] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
coordinates(SpRic) = ~lon+lat  #definindo colunas das coordenadas
gridded(SpRic) = TRUE #definindo gridded
proj4string(SpRic) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definind proj
rasterSpRic = raster(SpRic) #criando objeto raster
##criando imagem da distribuicao da riqueza
jpeg(filename='/home/anderson/Documentos/Projetos/divSpsSid/richness.jpeg')
plot(rasterSpRic)
dev.off()


##TESTANDO COM VARIAVEIS REAIS##

###variaveis ambientais da area de estudos###
predictors = stack(paste(caminhosCamadasTemp[1],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[1],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais

datMat = as.data.frame(predictors,xy=TRUE,na.rm=TRUE) #transformando raster em data.frame
names(datMat) = c('lon','lat','bio1','bio12') #ajustando os nomes das colunas do data.frame

##equacoes para as dimensoes do nicho das especies
Nsp = 5 #numero de especies a serem criadas

## betaBio1 = runif(n=Nsp,min=0.1,max=1) #parametro para cada equacao de cada especie
## betaBio12 = runif(n=Nsp,min=0.001,max=0.01) #parametro para cada equacao de cada especie
## alphaBio1 = runif(n=Nsp,min=quantile(x=varBio1,probs=0.5,na.rm=TRUE),max=quantile(x=varBio1,probs=1,na.rm=TRUE)) #parametro para cada equacao de cada especie
## alphaBio12 = runif(n=Nsp,min=quantile(x=varBio12,probs=0.5,na.rm=TRUE),max=quantile(x=varBio12,probs=1,na.rm=TRUE)) #parametro para cada equacao de cada especie

betaBio1 = abs(rnorm(n=Nsp,mean=0.1,sd=0.1)) #parametro para cada equacao de cada especie
betaBio12 = abs(rnorm(n=Nsp,mean=0.001,sd=0.1)) #parametro para cada equacao de cada especie
alphaBio1 = abs(rnorm(n=Nsp,mean=quantile(x=varBio1,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
alphaBio12 = abs(rnorm(n=Nsp,mean=quantile(x=varBio12,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
varBio1 = datMat$bio1 #variavel ambiental bioclim01
varBio12 = datMat$bio12 #variavel ambiental bioclim12

##solucao numerica para a equacao de cada especie
for (i in 1:Nsp){
    ##equacoes do nicho da especie
    fBio1Sp_i = as.integer( 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    fBio12Sp_i = as.integer( 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    #fBio1Sp_i = 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) #solucao da equacao com output continuo ("suitability")
    #fBio12Sp_i = 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) #solucao da equacao com output continuo ("suitability")
    datMat = data.frame(cbind(datMat,fSp=fBio1Sp_i*fBio12Sp_i)) #adicionando ao data.frame
    names(datMat)[ncol(datMat)] = paste('sp',i,sep='') #ajustando os nomes das especies no data.farme
    ##salvando graficos das equacoes de cada especie
    jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','functions_sp',i,'.jpeg',sep=''))
    par(mfrow=c(1,2))
    plot(fBio1Sp_i~varBio1,xlab='Bioclim 01',ylab='Suitability',ylim=c(0,1))
    plot(fBio12Sp_i~varBio12,xlab='Bioclim 12',ylab='Suitability',ylim=c(0,1))
    dev.off()
}

##criando a coluna de riqueza de especies##
datMat = cbind(datMat, Richness=rowSums(datMat[,grep('sp',names(datMat))]))

##raster da distribuicao modelada
for(i in 1:Nsp){
    SpDist = datMat[,c('lon','lat',paste('sp',i,sep=''))] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
    coordinates(SpDist) = ~lon+lat #definindo colunas das coordenadas
    gridded(SpDist) = TRUE #definindo gridded
    proj4string(SpDist) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definindo proj
    rasterSpDist = raster(SpDist) #criando objeto raster
    ##criando imagem da distribuicao de cada especie
    jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','sp',i,'.jpeg',sep=''))
    plot(rasterSpDist)
    dev.off()
}

##rater da riqueza
SpRic = datMat[,c('lon','lat','Richness')] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
coordinates(SpRic) = ~lon+lat  #definindo colunas das coordenadas
gridded(SpRic) = TRUE #definindo gridded
proj4string(SpRic) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definind proj
rasterSpRic = raster(SpRic) #criando objeto raster

##criando imagem da distribuicao da riqueza
jpeg(filename='/home/anderson/Documentos/Projetos/divSpsSid/richness.jpeg')
plot(rasterSpRic)
dev.off()

##salvando a matriz em um arquivo .csv

write.csv(x=datMat,file='/home/anderson/Documentos/Projetos/divSpsSid/datMat.csv',row.names=FALSE)

#################################
###ajusando funcao assintotica###
#################################

dados = data.frame(x=1:10,y=c(2,9,12,15,16,16,17,17,19,21))

modelo = formula(y ~ a - b/x)
#modeloSig = formula(y ~ a/(1+exp(-k*(x-0)))) #logistico, ainda dando erro
modeloLin = formula(y ~ x)

modFit = nls(formula=modelo, data=dados, start=list(a=0,b=0))
#modFitSig = nls(formula=modeloSig, data=dados, start=list(a=0,k=0.1)) #logistico, ainda dando erro :(
modFitLin = lm(formula=modeloLin, data=dados)

##avalaindo
summary(modFit)
#summary(modFitSig)

AIC(modFit,modFitLin)

anova(modFit,modFitLin,test='Chisq')

plot(dados$y ~ dados$x, pch=19,type='b')
points(fitted(modFit) ~ dados$x,pch='+',col='blue',type='b')
points(fitted(modFitLin) ~ dados$x,pch='x',col='red',type='b')

##projetando pontos fora dos dados, para visualisar assintota

novosX = 11:20

predY = predict(modFit,newdata=data.frame(x=novosX))

jpeg('/home/anderson/Documentos/Projetos/divSpsSid/assint.jpg')
plot(fitted(modFit)~dados$x,xlim=c(0,21),ylim=c(0,40),type='b')
points(predY~novosX,col='blue',pch=20,type='b')
points(dados$y~dados$x,pch=19,cex=0.7)
dev.off()
