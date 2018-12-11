##function for cross-temporal background points

paleobg = function(x, colNames=names(x), envFolder, n=10000){
  
  ### parameters ###
  ##x = data.frame with (at least) occurrences (lon, lat) and age
  ##colNames = names of the columns in x for lon, lat and age
  ##envFolder = path to predictor variables' folder
  ##n = desired number of background points
  ################
  
  colNames = colNames
  envFolder = envFolder
  occTable = x[,colNames]
  names(occTable) = c('lon','lat','age')
  ages = unique(occTable$age)
  ages = sample(ages, n, replace=TRUE)
  bgData = data.frame()
  
  for (age_i in unique(ages)){
    sampleSize = length(grep(age_i, ages))
    current_occPts = occTable[ match(age_i, occTable$age), c('lon','lat')]
    
    if(age_i %in% as.numeric(list.files(envFolder))){
      envData = list.files(paste(envFolder,'/',age_i,sep=''), full.names=TRUE)
    }
    if( length(agrep("maxent", basename(envData), value = T)) > 0 ){
      envData = envData[-c(agrep("maxent", basename(envData)))]
    }
    
    envData = stack(envData)
    crs(envData) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
    
    bgData_i = dismo::randomPoints(mask=envData[[1]], n=sampleSize, p=current_occPts)
    bgData_i = data.frame(bgData_i)
    bgData_i$age = age_i
    names(bgData_i) = c('lon','lat','age')
    bgData_i = data.frame(bgData_i, extract(envData, bgData_i[,c('lon','lat')]))
    bgData = rbind(bgData, bgData_i)
  }
  return(bgData)
}