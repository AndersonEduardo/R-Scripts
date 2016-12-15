####################################################################################################
##ALGORITMO PARA CALCULAR TSS (TRUE SKILL STATISTICS) A PARTIR DE MAPAS BINARIOS NO FORMATO RASTER##
####################################################################################################

TSSfunction = function(binaryMap,spOccPoints){
    if ((binaryMap@data@max !=1) | (binaryMap@data@min !=0) | (ncol(spOccPoints) != 2)){
        print("ERROR: Input data must be: (a) a binary raster; (b) a two column (long & lat) data.frame or matrix.")
        }else{

            ##abrindo os dados
            binMap = binaryMap 
            binMap.data = rasterToPoints(binMap)            
            spOccPoints = spOccPoints 

            ##extraindo informacoes para matriz de confusao
            a = sum(extract(binMap,spOccPoints),na.rm=TRUE) ##observed: presente - forecast: presente
            b = sum(binMap.data[,3])-sum(extract(binMap,spOccPoints),na.rm=TRUE) ##observed: ausente - forecast: presente
            c = nrow(spOccPoints)-sum(extract(binMap,spOccPoints),na.rm=TRUE) ##observed: presente - forecast: ausente
            d = nrow(binMap.data)-(a+b+c)  ##observed: ausente - forecast: ausente
            
            ##matriz de confusao
            MC = matrix(c(a,b,c,d),nrow=2,ncol=2,byrow=TRUE) #matriz de confusao

            ##calculo do TSS
            TSS = a/(a+c) + d/(b+d) - 1 #true skill statistic (Alouche et al. 2006)

            ##retornando o resultado do algoritmo
            return(TSS)
        }
}
