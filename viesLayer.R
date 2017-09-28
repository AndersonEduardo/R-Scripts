##script para produzir um grid file com as probabilidades de encontrar um registro de ocorrencia de especie usando conjunto de informacoes previas sobre outros registros ja realizados
##setembro/2017

##GRID FILE PARA VIES AMOSTRAL FOSSIL AMERICA DO SUL

##obtendo pontos de ocorrencia (para o PASSADO)
fossilDataRaw = paleobioDB::pbdb_occurrences(limit="all", base_name="vertebrata",show=c("coords", "phylo", "ident")) #dados do PBDB

##shape file da america do sul e extent para o grid file final
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp")
SOAextent = extent(-84.69856, -31.52813, -57.29995, 15.15035)

##inspecionando pontos fosseis
plot(AmSulShape)
points(fossilDataRaw[,c('lng','lat')])

##kernel density estimation com o pacote MASS
dens = MASS::kde2d(fossilDataRaw[,c('lng')],fossilDataRaw[,c('lat')],n=100)
densRas = raster(dens)
densRas = crop(densRas,SOAextent)

##salvado raster
writeRaster(x=densRas, file='/home/anderson/Documentos/Projetos/Sps artificiais/biasLayer.grd',overwrite=TRUE)

##inspecionando
plot(densRas)
plot(AmSulShape,add=TRUE)

points(fossilDataRaw[,c('lng','lat')],cex=0.5)
