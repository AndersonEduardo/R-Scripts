library(raster)
set.seed(8354)

popSize = vector()

##lattice
nli = 10
nco = 10
r <- raster(matrix(ncol=nco,nrow=nli,abs(round(rnorm(nco*nli,5,5),0))))
e <- extent(c(0,nli,0,nco))
extent(r) <- e
plot(r)

## pontos iniciais
pts = data.frame(x=c(nli/2,2), y=c(nco/2,2))
cells_i = cellFromXY(r, pts)

for (iteraction in 1:10){

    if(iteraction == 1){
        probExt = 0 #nao extinguir na primeira iteracao (sementes)
    }else{
        probExt = 0.5 #probabildiade de extincao
    }
    
    ##extincao local
    extEvents = sample(c(TRUE,FALSE), size=length(cells_i), replace=TRUE, prob=c(probExt, 1-probExt)) #extincao=TRUE; nao extincao=FALSE
    cells_i = cells_i[!extEvents] #persistem aquelas localidades que NAO EXIBIRAM eventos de extincao

    ##migracao para localidades vizinhas
    probOcu = 0.5 #porbabilidade de migrar
    ##buffer
    buf = adjacent(r, cells_i, 4)[,'to'] 
    migEvents = sample(c(TRUE,FALSE), size=length(buf), replace=TRUE, prob=c(probOcu, 1-probOcu))
    buf = buf[migEvents]
    
    ##juntando as celulas novas com as originais
    cells_i = append(cells_i, buf)
    cells_i = unique(cells_i)

    ##contagem do N ao longo das iteracoes
    popSize = append(popSize, length(cells_i))
    
}

##visualizacao temporal
plot(popSize, type='b')

##visualizcao espacial
cells = xyFromCell(r, cells_i)
cells = data.frame(cells)
coordinates(cells) = ~x+y
gridded(cells) = TRUE
plot(r)
plot(cells,add=T, col='blue',lwd=4)
