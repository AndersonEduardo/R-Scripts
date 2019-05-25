## funcao para Uncorrelated Spatio-Temporal Predictor Variables - USTPV. 
## O output eh um data.frame cujo conteudo indica quais variaveis estao pouco/não correlacionadas (numeros '1')
## e foram mantidas no dataset (usando a selecao de variaveis pela funcao 'vifcor' do  pacote usdm) quais 
## variaveis estão muito correlacionadas (numeros '0') e foram retiradas do dataset.
## Anderson A. Eduardo
## 11/fev/2019

require(raster)
require(usdm)

ustpv = function(path, ageMin, ageMax){
  
  ##parametros
  #path = caminho para a pasta do conjunto de variaveis ambientais
  #ageMin = idade minima do range temporal analisado
  #ageMax = idade maxima do range temporal analisado
  
  ##verificacoes de seguranca
  if(!is.numeric(ageMin) & !is.numeric(ageMax)){
    stop("Os parâmetros 'ageMin' e 'ageMax' precisam ser numeros.")
  }
  if(!is.character(path)){
    stop("O parâmetro 'path' precisa ser o caminho até a pasta com as variaveis ambientais.")
  }
  ##
  
  envFolder = path #'/home/anderson/gridfiles/dados_projeto'
  envVarPaths = list.dirs(envFolder, full.names=TRUE, recursive = FALSE)
  outputDF = data.frame()
  age = seq(ageMin, ageMax) #seq(0,22) #range de idades, em kyr BP
  
  for (i in age){
    tryCatch({
      
      ##definindo variaveis e parametros internos
      predictors = stack(list.files(path = paste(envFolder, '/', age[i+1], sep=''), 
                                    full.names = TRUE, 
                                    pattern = '.asc')) #predictors com todas as variaveis (presente)
      crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
      
      predictorsForVif = predictors
      vif(predictorsForVif)
      predictorsVif1 = vifcor(predictorsForVif, th=0.7)
      
      idx = grep(paste(as.character(predictorsVif1@results$Variables), collapse = '|'), names(predictors))
      
      outputDF = rbind( outputDF, c(rep(0, length(names(predictors))), age[i+1]) )
      names(outputDF) = c(names(predictors), 'age')
      
      outputDF[,idx] = 1
      
    }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  
  return(outputDF)
  
}
