## Anderson A. Eduardo
## 25/jan/2019

iterateFunc = function(x){
  
  datMatCurrent = x
  names(datMatCurrent) = c('lon','lat','bio1','bio12','fSp')
  
  ##condicao para nao permitir distribuicoes vazias (i.e. inexistente) ou tbm sobre a Am. Sul toda. Condicao: distribuicao > 1% ou <95% da america do sul
  #while( (sum(datMatCurrent[,paste('sp',i,sep='')]) < 0.05*(nrow(datMatCurrent))) | (sum(datMatCurrent[,paste('sp',i,sep='')]) > 0.5*(nrow(datMatCurrent))) ){
  while( (sum(datMatCurrent$fSp, na.rm=TRUE) < 0.05*(nrow(datMatCurrent))) | (sum(datMatCurrent$fSp, na.rm=TRUE) > 0.5*(nrow(datMatCurrent))) ){
    
    ##equacoes para as dimensoes do nicho das especies
    betaBio1 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
    betaBio12 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
    betaElev = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
    ##
    alphaBio1 = runif(n=1, min=quantile(datMatCurrent$bio1, probs=0.25, na.rm=TRUE), max=quantile(datMatCurrent$bio1, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
    alphaBio12 = runif(n=1, min=quantile(datMatCurrent$bio12, probs=0.25, na.rm=TRUE), max=quantile(datMatCurrent$bio12, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
    ## alphaElev = runif(n=1, min=quantile(datMatCurrent$elevation, probs=0.25, na.rm=TRUE), max=quantile(datMatCurrent$elevation, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
    
    ## betaBio1 = abs(rnorm(n=Nsp,mean=0.1,sd=0.1)) #parametro para cada equacao de cada especie
    ## betaBio12 = abs(rnorm(n=Nsp,mean=0.001,sd=0.1)) #parametro para cada equacao de cada especie
    ## alphaBio1 = abs(rnorm(n=Nsp,mean=quantile(x=varBio1,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
    ## alphaBio12 = abs(rnorm(n=Nsp,mean=quantile(x=varBio12,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
    varBio1 = datMatCurrent$bio1 #variavel ambiental bioclim01
    varBio12 = datMatCurrent$bio12 #variavel ambiental bioclim12
    ## varElev = datMatCurrent$elevation
    
    ##solucao numerica para a equacoes do nicho de cada especie
    fBio1Sp_i = 1/(1+exp(betaBio1*(varBio1-alphaBio1))) #solucao da equacao com output binario ("suitability")
    fBio12Sp_i = 1/(1+exp(-betaBio12*(varBio12-alphaBio12))) #solucao da equacao com output binario ("suitability")
    ## fElevSp_i = as.integer( 1/(1+exp(-betaElev*(varElev-alphaElev))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    
    ## fBio1Sp_i = 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) #solucao da equacao com output continuo ("suitability")
    ## fBio12Sp_i = 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) #solucao da equacao com output continuo ("suitability")
    ## fBio1Sp_i = as.integer( exp(-((varBio1-alphaBio1)^2/(2*betaBio1^2))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    ## fBio12Sp_i = as.integer( exp(-((varBio12-alphaBio12)^2/(2*betaBio12^2))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    
    ##datMatCurrent = data.frame(cbind(datMatCurrent,fSp=fBio1Sp_i*fBio12Sp_i)) #adicionando ao data.frame
    ##names(datMatCurrent)[ncol(datMatCurrent)] = paste('sp',i,sep='') #ajustando os nomes das especies no data.farme
    datMatCurrent[,'fSp'] = fBio1Sp_i*fBio12Sp_i ##*fElevSp_i
    
    ##salvando graficos das equacoes de cada especie
    ##jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','functions_sp',i,'.jpeg',sep=''))
    ##par(mfrow=c(1,2))
    ## plot(fBio1Sp_i~varBio1,xlab='Bioclim 01',ylab='Suitability',ylim=c(0,1))
    ## plot(fBio12Sp_i~varBio12,xlab='Bioclim 12',ylab='Suitability',ylim=c(0,1))
    ## plot(fElevSp_i~varElev,xlab='Elevation',ylab='Suitability',ylim=c(0,1)) 
    ##dev.off()
    
  }
  
  ##raster da distribuicao de adequabilidade climatica modelada
  SpDistribution = datMatCurrent[,c('lon','lat','fSp')] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
  coordinates(SpDistribution) = ~lon+lat #definindo colunas das coordenadas
  gridded(SpDistribution) = TRUE #definindo gridded
  proj4string(SpDistribution) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definindo proj
  rasterSpDistribution = raster(SpDistribution) #criando objeto raster
  
  return(rasterSpDistribution)
  
}