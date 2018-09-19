acmodel = function(lattice, seeds=1, probs="suitability", iteractions=10, plot=TRUE){

    ##variaveis locais
    popSize = c(seeds)
    ##lattice
    r = lattice
    nli = dim(r)[1] #10
    nco = dim(r)[2] #10
    ## r <- raster(matrix(ncol=nco,nrow=nli,abs(round(rnorm(nco*nli,5,5),0))))
    ## e <- extent(c(0,nli,0,nco))
    ## extent(r) <- e
    ##plot(r)
    ##pontos iniciais
    if (class(seeds) == 'numeric' &  length(seeds)==1){
        pts = data.frame( x = res(r)[1] * sample(seq(nli), seeds, replace=TRUE),
                         y= res(r)[2] * sample(seq(nco), seeds, replace=TRUE) )  #data.frame(x=c(nli/2,2), y=c(nco/2,2))
        cellsStart = cellFromXY(r, pts)
    } else if ( class(seeds) %in% c('integer','numeric') & length(seeds) > 1 ) {
        pts = xyFromCell(r, seeds)
        cellsStart = cellFromXY(r, pts)
    }else{
        cat("/n Problema: seeds precisa ser um numero ou um vetor de numeros inteiros (correspondente ao numero das celulas a serem colocadas as seeds). \n")
    }

    if (probs == "suitability"){
        probOcu = getValues(lattice)
        probExt = 1 - probOcu #probabildiade de extincao
    }else{
        probOcu = probs[1] #porbabilidade de migrar
        probExt = probs[2] #probabildiade de extincao
    }

    for (iteraction in 1:iteractions){

        if(iteraction == 1){
            cellsPop = cellsStart
        }

        ##migracao para localidades vizinhas
        ##buffer
        buf = adjacent(r, cellsPop, 4)[,'to']
        if (probs == "suitability"){
            currentProbOcu = probOcu[buf]
            migEvents = sapply(seq(length(buf)), function(x)  sample(c(TRUE,FALSE), size=1, replace=TRUE, prob=c(currentProbOcu[x], 1-currentProbOcu[x])) )
        }else{
            currentProbOcu = probOcu
            migEvents = sample(c(TRUE,FALSE), size=length(buf), replace=TRUE, prob=c(probOcu, 1-probOcu))
        }

        colonizedCells = buf[migEvents]

        ##atualizacao da pos pos-migracao
        cellsPop = unique(c(cellsPop, colonizedCells))
        
        ##extincao local
        if (probs == "suitability"){
            currentProbExt = probExt[buf]
            extEvents = sapply(seq(length(cellsPop)), function(x)  sample(c(TRUE,FALSE), size=1, replace=TRUE, prob=c(currentProbExt[x], 1-currentProbExt[x])) )
        }else{
            extEvents = sample(c(TRUE,FALSE), size=length(cellsPop), replace=TRUE, prob=c(probExt, 1-probExt)) #extincao=TRUE; nao extincao=FALSE
        }

        survivingCells = cellsPop[!extEvents] #persistem aquelas localidades que NAO EXIBIRAM eventos de extincao
        
        ##atualizacao da pop pos-extincoes locais
        cellsPop = unique(survivingCells)

        ##contagem do N ao longo das iteracoes
        popSize = append(popSize, length(cellsPop))

        ##checando extincao
        if (length(cellsPop) == 0){
            stop('\n Populacao extinta! \n')
        }
        
        ##distribuicao espacial
        cells = xyFromCell(r, cellsPop)
        cells = data.frame(cells)
        coordinates(cells) = ~x+y
        ##gridded(cells) = TRUE

        ##visualizacao
        if (plot == TRUE){
            plot(r)
            plot(cells, add=TRUE, col='blue', lwd=4, pch=13, cex=2)
            Sys.sleep(0.5)
        }        
    }

    ##output    
    output = list(popSize=popSize, distribution=cells)
    
    return(output)
    
}
