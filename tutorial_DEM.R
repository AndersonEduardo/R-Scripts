## tutorial para dowload e dados de altitude com o pacote 'elevatr'
## fonte: https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html
## acesso: 19/jan/2018

install.packages('elevatr')
library(elevatr)
cat("mapzen_key=mapzen-kGT6SUr\n", file = file.path(normalizePath("~/"), ".Renviron"), append = TRUE)

##Using get_elev_point() to Access The Mapzen Elevation Service

# Create an example data.frame
set.seed(65.7)
examp_df <- data.frame(x = runif(10, min = -73, max = -71), y = runif(10, min = 41, 
    max = 45))
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Create and example data.frame with additional columns
cats <- data.frame(category = c("H", "H", "L", "L", "L", "M", "H", "L", "M", "M"))

examp_df2 <- data.frame(examp_df, cats)

# Create an example SpatialPoints
examp_sp <- SpatialPoints(examp_df, proj4string = CRS(prj_dd))

# Create an example SpatialPointsDataFrame
examp_spdf <- SpatialPointsDataFrame(examp_sp, proj4string = CRS(prj_dd), data = cats)

# Example using data.frame with longitude and latitude
df_elev <- get_elev_point(examp_df, prj = prj_dd, src = "mapzen")

# Compare
examp_df
data.frame(df_elev)
# Example using data.frame with longitude, latitude and an additional column
df2_elev <- get_elev_point(examp_df2, prj = prj_dd, src = "mapzen")

# Compare
examp_df2
data.frame(df2_elev)


# Example using SpatialPointsDataFrame prj is taken from the
# SpatialPointsDataFrame object api_key is taken from environment variable
# mapzen_key
spdf_elev <- get_elev_point(examp_spdf)

# Compare
examp_spdf
spdf_elev


##Mapzen Terrain Tile Service

# SpatialPolygonsDataFrame example
data(lake)
elevation <- get_elev_raster(lake, z=9, src="aws")

plot(elevation)
plot(lake, add = TRUE)

# data.frame example
elevation_df <- get_elev_raster(examp_df, prj = prj_dd, z = 5)

plot(elevation_df)
plot(examp_sp, add = TRUE)

# Bounding box on edge
elev_edge <- get_elev_raster(lake, z = 10, api_key = key)

plot(elev_edge)
plot(lake, add = TRUE)

# Use expand to grab additional tiles
elev_expand <- get_elev_raster(lake, z = 10, expand = 1500, api_key = key)

plot(elev_expand)
plot(lake, add = TRUE)


