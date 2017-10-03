##SCRIPT PARA CORTAR UMA AREA DE INTERESSE EM UM ARQUIVO DO TIPO RASTER##

#abrindo os pacotes necessarios
library(dismo)
library(raster)

#abrir a amada a ser cortada
camada = raster("/home/anderson/PosDoc/dados_ambientais/bcmidbi_2-5m _asc_America_Sul/bcmidbi1.grd") 

#definir a area de interesse
areaCorte = c(-81.58333,-34.04167,-57.125,13) #manter essas coordenadas para os barbeiros

#cortanto a camada
camadaCortada = crop(camada,areaCorte) 

#salvando no computador
writeRaster(camadaCortada,filename="/home/anderson/bio1-AmericaSul.asc",overwrite=TRUE) 
