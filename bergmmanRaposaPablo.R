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
plot(dados$tam~dados$temp,cex=0.5,col='black',xlab=c('Temperatura'),ylab=c("Tamanho"),grid=T)
abline(modelo)
#points(fitted(modelo)~temp,pch=20)
#points(fitted(modelo)~temp70,pch=20,col="red")
legend(150,25.5,pch=c(20,20),col=c("black","red"),c("Presente","2070"))
grid()

#usando o modelo ajustado para projetar no futuro
projecao=predict.glm(modelo,newdata=data.frame(temp=dados$temp21),type='response',se.fit=T) #calculando a projecao
points(projecao[[1]],pch=20,col='red') #plotando os tamanhos estimados a partir do modelo
modeloProj = glm(projecao[[1]]~dados$temp21,family='gaussian') #criando uma equacao para a linha dos pontos projetados
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

library(ggplot2)
#ggplot(dados,aes(x=tam,fill='red'))+geom_histogram(alpha=0.5) #??

hist(dados$tam,breaks=15,xlim=c(16,28),ylim=c(0,90),col=rgb(0,0,1,0.7),main="",xlab="Tamanho")
par(new=TRUE)
hist(projecao[[1]],xlim=c(16,28),ylim=c(0,90),breaks=5,col=rgb(0,1,1,0.4),main="",xlab="",ylab="")
legend('topright',fill=c('blue','red'))
#plot(density(dados$tam),xlim=c(15,30),ylim=c(0,1))
#lines(density(projecao[[1]],xlim=c(15,30)))


corTest=cor(as.data.frame(cbind(tam,projecao[[1]])))
plot(projecao[[1]]~tam,pch=20,cex=0.8,xlab=c("Tamanho observado"),ylab=c("Tamanho projetado (2070)"))
abline(lm(projecao[[1]]~tam))
grid()

##diferencas entre dados e projecao
dif=projecao[[1]]-tam
#
tmin=min(tam)
tmax=max(tam)
colfun=((tam-tmin)/(tmax-tmin))
#
difmin=min(abs(dif))
difmax=max(abs(dif))
difcolfun=((abs(dif)-difmin)/(difmax-difmin))
#

plot(AmSulShape)
points(pontos,pch=20,col=ifelse(dif>0,rgb(difcolfun,0,0,alpha=0.5),rgb(0,0,difcolfun,alpha=0.5)),cex=difcolfun*7)
legend('left',legend=c("Tamanho menor","","Sem mudança","","Tamanho maior"),pch=20,cex=0.8,col=c('darkblue','blue','black','red','darkred',alpha=0.5),pt.cex=c(4,3,1.5,3,4),x.intersp=2,y.intersp=1.8, bg="white",title=c("Mudança nos tamanhos:"))

plot(dif~pontos$Lat,xlab=c("Latitude"),ylab=c("Tamanho projetado - Tamanho observado"),col=ifelse(dif>0,rgb(difcolfun,0,0,alpha=0.7),rgb(0,0,difcolfun,alpha=0.7)),pch=20,cex=2)





#########
x=10:30
new.x=1:21
#new.x2=c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,27,30,31,34,36,38,40)
y=10+2*x+rnorm(length(x),0,4)
dados=data.frame(cbind(varx=x,nvarx=new.x,vary=y))
plot(dados$vary~dados$varx)

mod = glm(vary~varx,family=c('gaussian'),data=dados)
summary(mod)

proj = predict.glm(mod,newdata=data.frame(varx=dados$nvarx),type='response')

fitted(mod)[1:5]
proj[1:5]


plot(dados$vary~dados$varx,xlim=c(0,35),ylim=c(0,80))
points(fitted(mod)~dados$varx,col='blue',pch=20)
points(proj~dados$nvarx,col='red',pch=1,cex=3)
points(fitted(mod)~dados$nvarx,col='green')
