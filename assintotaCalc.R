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
