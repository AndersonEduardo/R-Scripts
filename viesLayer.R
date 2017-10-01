##script para produzir um grid file com as probabilidades de encontrar um registro de ocorrencia de especie usando conjunto de informacoes previas sobre outros registros ja realizados
##setembro/2017


##GRID FILE PARA VIES AMOSTRAL FOSSIL AMERICA DO SUL

##obtendo pontos de ocorrencia (para o PASSADO)
#fossilDataRaw = paleobioDB::pbdb_occurrences(limit="all", base_name="vertebrata",continent='SOA', show=c("coords", "phylo", "ident")) #dados do PBDB
#write.csv(fossilDataRaw, '/home/anderson/Documentos/Projetos/Sps artificiais/vertFosseis.csv',row.names=FALSE)

##abrindo arquivo salvo
fossilDataRawData = read.csv(file='/home/anderson/Documentos/Projetos/Sps artificiais/vertFosseis.csv', header=TRUE)
fossilDataRaw = fossilDataRawData[,c('lng','lat')]
fossilDataRaw = unique(fossilDataRaw)

##convertendo de data.frame para staialPointsDataFrame
coordinates(fossilDataRaw) <- ~lng+lat
proj4string(fossilDataRaw) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

##abrindo o shapefile que 'cortara' os pontos
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp",p4s=proj4string(fossilDataRaw))

##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
fossilDataRaw <- spTransform(fossilDataRaw, CRS(proj4string(AmSulShape)))
fossilData <- fossilDataRaw[AmSulShape, ]

##extent para o grid file final
#SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)

##inspecionando pontos fosseis
plot(fossilDataRaw)
points(fossilData, col='red', cex=2)
plot(AmSulShape, add=TRUE)

##convertendo spatialPoints em data.frame para o KDE2D
fossilDataDF = as.data.frame(fossilData)

##kernel density estimation com o pacote MASS
dens = MASS::kde2d(x=fossilDataDF$lng, y=fossilDataDF$lat, n=100)
densRas = raster(dens)
densRas = mask(x=densRas, mask=AmSulShape)

#densRasLarger = extend(densRas,extent(-180,0,-90,0)) ##coordenadas do hemisferio sul ocidental


##salvado raster
writeRaster(x=densRas, file='/home/anderson/Documentos/Projetos/Sps artificiais/biasLayer.grd',overwrite=TRUE)

##inspecionando
plot(densRas)
plot(AmSulShape,add=TRUE)
points(fossilDataDF[,c('lng','lat')],cex=0.5)
