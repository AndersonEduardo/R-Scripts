##codigo basico para reduzir a resolucao de um raster (i.e. aumentar o tamanho das quadriculas com informacao no raster)
library(raster)
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
