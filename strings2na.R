##funcao transforma strigs em NAs em datasets de registros de ocorrencia de especies.
##06/Setembro/2018
##Anderson Aires Eduardo


strings2na = function(x, safeCols){

    if (class(x) != "data.frame"){
        stop("O dataset de entrada precisa ser um objeto da classe 'data.frame'.")
    }
    
    if (class(safeCols) != "character"){
        stop("O nome da coluna do nome das especies precisa ser um objeto da classe 'character'.")
    }

    ##variaveis locais
    ptsCoordsRaw = x
    safeCols = safeCols

    ##manipulacao
    ptsCoords = apply(ptsCoordsRaw[, which(names(ptsCoordsRaw) != safeCols) ],1, as.numeric) #transforma informacao de texto em NA (ex.: gruta azul -> NA)
    ptsCoords = data.frame(t(ptsCoords)) #transformando em data.frame
    ptsCoords = data.frame(ptsCoordsRaw[, safeCols], ptsCoords)
    names(ptsCoords) = names(ptsCoordsRaw)
    ##ptsCoords = data.frame(sp=pts$Species, lon=ptsCoords[,1], lat=ptsCoords[,2], mean=pts$Cal..Mean, min=pts$Min, max=pts$Max) #consolidando os dados

    ##saida da funcao
    return(ptsCoords)
}
