##script para retirar pontos que estao fora de um poligono

##pacotes
library(maptools)
library(raster)
data(wrld_simpl)

##coordinates
test.points <- data.frame(lon = c(-5,0,10,19,19,19), lat = c(45,45,45,50,55,60))

##plot (shape-file + points)
plot(wrld_simpl, axes = TRUE, col = "grey", xlim = c(-15,+45), ylim = c(+35,+60)); box()
points(test.points$lon, test.points$lat, col = "red", pch = 20, cex = 1)

##convertendo de data.frame para staialPointsDataFrame
coordinates(test.points) <- ~lon+lat
proj4string(test.points) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

##retirando os pontos externos
insidePoints = over(test.points,wrld_simpl)

##plotando os dados limpos
plot(wrld_simpl, axes = TRUE, col = "grey",xlim = c(-15,+45),ylim = c(+35,+60)); box()
points(insidePoints$LON, insidePoints$LAT,col = "red", pch = 1, cex = 2)
