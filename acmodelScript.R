## Script exemplo para acmodel (modelo de simulacao com automatos celulares)

library(lattice)
source('/home/anderson/R-Scripts/acmodel.R')

var1 = 1:10
var2 = 1:10

aSp1=5
aSp2=8

nichoSp1 = expand.grid(var1=var1, var2=var2)
nichoSp1$suitability = ( 1/(1+exp(nichoSp1$var1 - aSp1)) )*( 1/(1+exp(nichoSp1$var2 - aSp1)) )
nichoSp1Raster = matrix(nichoSp1$suitability, length(var1), length(var2))
nichoSp1Raster = raster(nichoSp1Raster)

wireframe(suitability~var1+var2, nichoSp1)

nichoSp2 = expand.grid(var1=var1, var2=var2)
nichoSp2$suitability = ( 1/(1+exp(nichoSp2$var1 - aSp2)) )*( 1/(1+exp(nichoSp2$var2 - aSp2)) )
nichoSp2Raster = matrix(nichoSp2$suitability, length(var1), length(var2))
nichoSp2Raster = raster(nichoSp2Raster)

wireframe(suitability~var1+var2, nichoSp2)


distSp1 = acmodel(lattice=nichoSp1Raster, seeds=10, iteractions=30) #modelo sp1

distSp2 = acmodel(lattice=nichoSp2Raster, seeds=10, iteractions=30) #modelo sp2



### teste para interacao entre especies ##

nichoSp1Raster = nichoSp1Raster>0

interacao = c(1,1) #[0,1]=competicao; [1,0]=mutualismo; [1,1]=neutralismo; [0,0]=aniquilação

interact()

interact = function(){

    for (iteracao in 1:10){

        bg = nichoSp1Raster*0

        if (iteracao == 1){
            env1 = nichoSp1Raster
            env2 = nichoSp2Raster
        }else{
            distSp1Raster = raster::rasterize(distSp1$distribution, bg)
            distSp2Raster = raster::rasterize(distSp2$distribution, bg)

            distSp1Raster[distSp1Raster>0] = interacao[1] #0=competicao, 1=inexistente
            distSp1Raster[is.na(distSp1Raster)] = interacao[2] #0=mutualismo, 1=inexistente

            distSp2Raster[distSp2Raster>0] = interacao[1] #0=competicao, 1=inexistente
            distSp2Raster[is.na(distSp2Raster)] = interacao[2] #0=mutualismo, 1=inexistente
            
            env1 = distSp2Raster*nichoSp1Raster
            env2 = distSp1Raster*nichoSp2Raster
        }

        if (iteracao == 1){
            seeds1 = seeds2 = 50
        }else{
            seeds1 = distSp1$popSize[length(distSp1$popSize)]
            seeds2 = distSp2$popSize[length(distSp2$popSize)]
        }
        
        distSp1 =  acmodel(lattice=env1, seeds=seeds1, iteractions=1, plot=FALSE)
        distSp2 = acmodel(lattice=env2, seeds=seeds2, iteractions=1, plot=FALSE)


        plot(bg, col='lightgrey')
        plot(distSp1$distribution, lwd=2, cex=2, col='blue', add=T)
        plot(distSp2$distribution, lwd=2, col='red', pch=2, add=T)
        Sys.sleep(0.2)

    }
}

