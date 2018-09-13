##funcao que define a idade de registros de ocorrencia a partir do intervalo de erro das suas idades.
##06/Setembro/2018
##Anderson Aires Eduardo


dataInstance = function(x, col_names, tempRes=1000, n=2*nrow(x)){

    if (class(x) != "data.frame"){
        stop("O dataset de entrada precisa ser um objeto da classe 'data.frame'.")
    }
    
    if (class(col_names) != "character"){
        stop("O nome da coluna do nome das especies precisa ser um objeto da classe 'character', contendo os nomes das colunas com a idade media, idades minima e maxima dos resgistros.")
    }
    
    if (class(tempRes) != "numeric"){
        stop("A resolucao temporal (tempRes) precisa ser um objeto da classe 'numeric'. Esse dado refere-se a resolucao temporal dos dados ambintais (o default é de 1000 anos).")
    }

    if (class(n) != "numeric"){
        stop("O numeri de iteracoes (n) precisa ser um objeto da classe 'numeric'. Esse dado refere-se ao numero de instancias de dados a serem criadas.")
    }


    ##variaveis locais
    ptsAgeRaw = x
    col_names = col_names
    output = list()

    for (i in seq(n)){
        
        ##manipulacao para substituicao de strings e NAs
        ptsAge = apply(ptsAgeRaw[,col_names], 1, as.numeric) #transforma informacao de texto em NA (ex.: pleistocene -> NA)
        ptsAge = data.frame(t(ptsAge)) #consolidando os dados
        ptsAgeMeanNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[1]), mean(x[2:3]), x[1]) ) #se media=NA, obtem a partir do intervalo (max e min)
        ptsAgeMinNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[2]) & is.na(x[3]), x[1], x[2]) ) #se min=NA, entao min=mean (i.e. data altamente precisa)
        ptsAgeMaxNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[2]) & is.na(x[3]), x[1], x[3]) ) #se max=NA, entao max=mean (i.e. data altamente precisa)
        ptsAge = data.frame(cbind(ptsAgeMeanNA,ptsAgeMinNA,ptsAgeMaxNA)) #consolidando os dados de idade
        ptsAge = data.frame(ptsAgeRaw[,names(ptsAgeRaw)[col_names != names(ptsAgeRaw)]], ptsAge)
        names(ptsAge) = names(ptsAgeRaw)

        ##manipulacao para obtencao uma isntancia para as idades
        inxMin = agrep(pattern='min', x=col_names, max.distance=1)
        inxMax = agrep(pattern='max', x=col_names, max.distance=1)
        
        if(inxMin == inxMax | any(is.na(c(inxMin,inxMax)))){
            stop("Falha na identificacao das colunas de idade media, minima e maxima. Nomeie essas colunas preferencialmente como 'mean', 'min' e 'max'.")
        }else{
            minmaxCols = col_names[c(inxMin, inxMax)]
        }
        
        ptsAge$age = apply( ptsAge[,minmaxCols], 1, function(x) runif(1, min=x[1], max=x[2]) ) #sorteia uma idade detro do intervalo
        ptsAge$age = round(ptsAge$age/tempRes) #deixando idades como numeros inteiros (importante para corresponder ao nome das pastas)

        ptsAge = ptsAge[,-match(col_names, names(ptsAge)) ]

        ##adicionando instancia ao output
        output[[i]] = ptsAge

    }
    
    ##saida da funcao
    return(output)

}
