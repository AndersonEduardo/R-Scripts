library(raster)
library(rgdal)
library(maptools)

var=raster("/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_01.asc",pattern='asc')
var21=raster("/home/anderson/PosDoc/dados_ambientais/dados_projeto/021/bioclim_01.asc",pattern='asc')
sp.data = read.table("/home/anderson/R/bergmman/CulpeoBergmann.txt",header=T)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

pontos=sp.data[,c(2,1)]
temp=extract(var,pontos,method='bilinear')
temp21=extract(var21,pontos,method='bilinear')
tam=sp.data[,3]
dados=as.data.frame(cbind(pontos,tam,temp,temp21))
#ausencias=randomPoints(mask=var,n=nrow(pontos))
#temp.ausencias=extract(var,ausencias)
#dados.ausencias=as.data.frame(cbind(ausencias,temp.ausencias,pres=0))
#names(dados.ausencias)=c("Long","Lat","temp","pres")
#dados.totais=as.data.frame(rbind(dados,dados.ausencias))

#ajuste  do modelo (a partir dos dados reais)
modelo=glm(tam~temp,family="gaussian",data=dados)

##grafico principal
plot(tam~temp,cex=0.5,col='black',xlab=c('Temperatura'),ylab=c("Tamanho"),grid=T)
abline(modelo)
#points(fitted(modelo)~temp,pch=20)
#points(fitted(modelo)~temp70,pch=20,col="red")
legend(150,25.5,pch=c(20,20),col=c("black","red"),c("Presente","2070"))
grid()

#usando o modelo ajustado para projetar no futuro
projecao=predict.glm(modelo,as.data.frame(temp21),type='response',se.fit=T) #calculando a projecao
points(projecao[[1]],pch=20,col='red') #plotando os tamanhos estimados a partir do modelo
modeloProj = glm(projecao[[1]]~temp21,family='gaussian') #criando uma equacao para a linha dos pontos projetados
abline(modeloProj,col='red') #plotando a linha no grafico

#comparando as temperaturas
qqplot(x=temp,y=temp21,xlab=c("Temperatura atual"),ylab=c("Temperatura a 21 kyr BP"))
points(sort(temp21),sort(temp21),type='l',col='red')
boxplot(temp21,temp,names=c('21 kyr BP','Presente'))
plot(temp21~temp,pch=20,xlab=c("Temperatura atual"),ylab=c("Temperatura a 21 kyr BP"))
points(temp21,temp21,type='l',col='red')

##comparando dados empiricos 
qqplot(x=tam,y=projecao[[1]],xlab=c("Tamanho observado"),ylab=c("Tamanho projetado (2070)"))
points(sort(tam),sort(tam),type='l',col='red')

corTest=cor(as.data.frame(cbind(tam,projecao[[1]])))
plot(projecao[[1]]~tam,pch=20,cex=0.8,xlab=c("Tamanho observado"),ylab=c("Tamanho projetado (2070)"))
abline(lm(projecao[[1]]~tam))
grid()

dif=projecao[[1]]-tam
#
tmin=min(tam)
tmax=max(tam)
rge=tmax-tmin
colfun=((tam-tmin)/(tmax-tmin))
#
difmin=min(abs(dif))
difmax=max(abs(dif))
difrge=difmax-difmin
difcolfun=((abs(dif)-difmin)/(difmax-difmin))
#

plot(AmSulShape)
points(pontos,pch=20,col=ifelse(dif>0,rgb(difcolfun,0,0,alpha=0.5),rgb(0,0,difcolfun,alpha=0.5)),cex=difcolfun*7)

plot(dif~pontos$Lat,xlab=c("Latitude"),ylab=c("Tamanho projetado - Tamanho observado"),col=ifelse(dif>0,rgb(difcolfun,0,0,alpha=0.7),rgb(0,0,difcolfun,alpha=0.7)),pch=20,cex=2)





