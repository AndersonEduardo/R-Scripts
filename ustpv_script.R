####script/tutorial para selecao de variaveis preditoras 0kyr - 22kyr BP####

library(raster)
library(usdm)

source('/home/anderson/R-Scripts/ustpv.R')

envFolder = '/home/anderson/gridfiles/dados_projeto' #caminho da pasta com o conjunto de variaveis preditoras
ageMin = 0 #idade minima para o range temporal, em kyr BP
ageMax = 22 #idade maxima para o range temporal, em kyr BP

ustpvCalc = ustpv(path=envFolder, ageMin=ageMin, ageMax=ageMax) #rodando o ustpv

#verificacoes...
View(ustpvCalc)
apply(outputDF, 2, sum)

#conjunto de variaveis selecionadas
ustpv = c('bioclim_04', 'bioclim_10', 'bioclim_15', 'bioclim_16', 'bioclim_17')  #minimum set of Uncorrelated Spatio-Temporal Predictor Variables
#OBS.: pressuposto deste exemplo: toda America do Sul e ocupavel por todas as especies estudadas
