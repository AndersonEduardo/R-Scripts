library(raster)

projFolder = '/home/anderson/Documentos/Projetos/transferability'

environ1 = raster(ncol=100,nrow=100)
var1vector = rep(x=dnorm(x=1:ncol(environ1),mean=50,sd=20), times=nrow(environ1))
var2vector = rep(x=dnorm(x=1:nrow(environ1),mean=50,sd=20), times=ncol(environ1))
#
var1 = matrix(var1vector,nrow=nrow(environ1),ncol=ncol(environ1),byrow=FALSE) * 500
var2 = matrix(var2vector,nrow=nrow(environ1),ncol=ncol(environ1),byrow=TRUE) * 500
#
environ1[] = (var1 * var2) 
#
plot(environ1)

environ2 = raster(ncol=100,nrow=100)
var1vector = rep(x=dnorm(x=1:ncol(environ2),mean=50,sd=20), times=nrow(environ2)) * 500
var2vector = rep(x=dnorm(x=1:nrow(environ2),mean=75,sd=20), times=ncol(environ2)) * 500
#
var1 = matrix(var1vector,nrow=nrow(environ2),ncol=ncol(environ2),byrow=FALSE)
var2 = matrix(var2vector,nrow=nrow(environ2),ncol=ncol(environ2),byrow=TRUE)
#
environ2[] = (var1 * var2) * 1.5
#
plot(environ2)

#varificando climas nao analogos entre os dois ambintes
pontos = randomPoints(mask=environ1,n=100,prob=TRUE)
#pontos = data.frame(x=sample(-180:180,100),y=sample(-90:90,100))
v = extract(environ1, pontos)
messTest =  mess(environ2,v)

par(mfrow=c(1,2))
plot(environ1)
points(pontos)
plot(messTest)

##superficies para probabilidade de amostragem (para simulacao do vies amostral)

for (i in seq(0,1,0.2)){
    mu = i*50 + 50
    ##
    probSamp = raster(ncol=100,nrow=100)
    var1vector = rep(x=dnorm(x=1:ncol(probSamp),mean=mu,sd=20), times=nrow(probSamp)) * 500
    var2vector = rep(x=dnorm(x=1:nrow(probSamp),mean=mu,sd=20), times=ncol(probSamp)) * 500
    ##
    var1 = matrix(var1vector,nrow=nrow(probSamp),ncol=ncol(probSamp),byrow=FALSE)
    var2 = matrix(var2vector,nrow=nrow(probSamp),ncol=ncol(probSamp),byrow=TRUE)
    ##
    probSamp[] = (var1 * var2) * 1.5
    ##
    writeRaster(x=probSamp,filename=paste(projFolder,'/probSamp',mu,'.asc',sep=''))
    
    plot(probSamp)
}


###
hii = raster(x='/home/anderson/Downloads/hii-s-america-geo-grid/hii_s_america_grid/hii_s_amer')
dado = raster(x='/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_01.asc')

##crop hii
areaDado = extent(dado)
hiiSA = crop(hii,areaDado)

##hii ajustado
raster1 = res(dado)
raster2 = res(raster2)
fator = round(raster1[1]/raster2[1])  ##raster maior / raster menor
hiiSA2 = aggregate(hiiSA,fator)
##verificando
res(hiiSA2)
res(dado)
