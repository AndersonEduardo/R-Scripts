# algoritmo para criar aleatoriamente parametros para funcao resposta sigmuiodal de uma especies 
# aas variaveis bio1 e bio12
# Anderson A. Eduardo
# 29/jan/2019

makeRespFunc = function(x){
  datMatCurrent = x
  names(datMatCurrent) = c('lon','lat','bio1','bio12','fSp')
  
  ##condicao para nao permitir distribuicoes vazias (i.e. inexistente) ou tbm sobre a Am. Sul toda. Condicao: distribuicao > 1% ou <95% da america do sul
  #while( (sum(datMatCurrent[,paste('sp',i,sep='')]) < 0.05*(nrow(datMatCurrent))) | (sum(datMatCurrent[,paste('sp',i,sep='')]) > 0.5*(nrow(datMatCurrent))) ){
  #while( (sum(datMatCurrent$fSp, na.rm=TRUE) < 0.05*(nrow(datMatCurrent))) | (sum(datMatCurrent$fSp, na.rm=TRUE) > 0.5*(nrow(datMatCurrent))) ){
    
    ##patametros
    betaBio1 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
    betaBio12 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
    ##
    alphaBio1 = runif(n=1, min=quantile(datMatCurrent$bio1, probs=0.25, na.rm=TRUE), max=quantile(datMatCurrent$bio1, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
    alphaBio12 = runif(n=1, min=quantile(datMatCurrent$bio12, probs=0.25, na.rm=TRUE), max=quantile(datMatCurrent$bio12, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
  #}
  
  output = data.frame(betaBio1=betaBio1, betaBio12=betaBio12, alphaBio1=alphaBio1, alphaBio12=alphaBio12)
  
  class(output) = "respFuncObject"
  
  return(output)
  
}