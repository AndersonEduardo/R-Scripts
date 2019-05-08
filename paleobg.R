##function for cross-temporal background points

require(dismo)
require(raster)

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
  vecNA = vector()
  
  for ( age_i in unique(ages) ){
    
    sampleSize = length(grep(age_i, ages))
    current_occPts = occTable[ match(age_i, occTable$age), c('lon','lat')]
    
    if(age_i %in% as.numeric(list.files(envFolder))){
      envData = list.files(paste(envFolder,'/',age_i,sep=''), full.names=TRUE)
      varNames = basename(list.files(paste(envFolder,'/',age_i,sep=''), full.names=TRUE))
      varNames = gsub("\\..*", "", varNames)
    }else{
      warning("Atenção: algumas idades não possuem dados ambientais no directório. NAs produzidos.")
      vecNA = append( vecNA, age_i )
      next
    }
    if( (length(agrep("maxent", basename(envData), value = T)) > 0) ){
      envData = envData[-c(agrep("maxent", basename(envData)))]
    }
    
    envData = stack(envData)
    crs(envData) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
    
    bgData_i = dismo::randomPoints(mask=envData[[1]], n=sampleSize, p=current_occPts)
    bgData_i = data.frame(bgData_i)
    bgData_i$age = age_i
    names(bgData_i) = c('lon','lat','age')
    bgData_i = round(bgData_i, 2)
    bgData_i = unique(bgData_i)
    bgData_i_vars = extract(envData, bgData_i[,c('lon','lat')])
    bgData_i_vars = data.frame(bgData_i_vars)
    names(bgData_i_vars) = varNames
    bgData_i = data.frame(bgData_i,  bgData_i_vars)
    bgData = rbind(bgData, bgData_i)
  }
  
  #dealing with unexistnt ages
  if ( length(vecNA) > 0 ){
    matrixNA = matrix( data=rep(NA, length(vecNA)*max(2,ncol(bgData))), nrow=length(vecNA), ncol=max(2,ncol(bgData)) )
    matrixNA = data.frame(matrixNA)
    names(matrixNA) = unlist( list('data', names(bgData) )[ifelse(length(names(bgData))==0, 1, 2)] )
    idxForAge = ifelse( ncol(matrixNA)==2, 1, grep('age', names(matrixNA)))
    matrixNA[ ,idxForAge] = vecNA
    bgData = rbind(bgData, matrixNA)
  }
  
  #output
  return(bgData)
}