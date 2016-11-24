library(raster)
##TSS

CyacareBIN = raster('/home/anderson/PosDoc/teste/Maxent/Raster Layers/Caiman yacareBINARIO.asc')
CyacareBIN.data = rasterToPoints(CyacareBIN)

CyacarePoints = read.csv("/home/anderson/PosDoc/dados_ocorrencia/PO_unique/Caiman yacare.csv",header=T)
CyacarePoints.data = extract(CyacareBIN,CyacarePoints[,2:3])

dataVetor=numeric()
for (i in 1:nrow(CyacareBIN.data)){
    binaryIndicator = 0
    for (j in 1:nrow(CyacarePoints[,2:3])){
        if ( (round(CyacareBIN.data[i,1],digits=2) == round(CyacarePoints[j,2],digits=2)) && (round(CyacareBIN.data[i,2],digits=2) == round(CyacarePoints[j,3],digits=2)) ){
            binaryIndicator = 1
        }
    }
    if (binaryIndicator == 0){
        dataVetor = append(dataVetor,0)
    }else{
        dataVetor = append(dataVetor,1)
    }
}


dataSet = data.frame(CyacareBIN.data,CyacarePoints=dataVetor)
names(dataSet) = c('lon','lat','forecast','observed')

##matriz de confusao (MC)
a=b=c=d=0
for (i in 1:nrow(dataSet)){
    if ( (dataSet$forecast[i]==1) && (dataSet$observed[i]==1) ){
        a = a+1
    }
    ##
    if ( (dataSet$forecast[i]==1) && (dataSet$observed[i]==0) ){
        b = b+1
    }
    ##
    if ( (dataSet$forecast[i]==0) && (dataSet$observed[i]==1) ){
        c = c+1
    }
    ##
    if ( (dataSet$forecast[i]==0) && (dataSet$observed[i]==0) ){
        d = d+1
    }
}

MC = matrix(c(a,b,c,d),nrow=2,ncol=2,byrow=TRUE)

TSS = (a*d-b*c)/((a+c)*(b+d))

