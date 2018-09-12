##Anderson A. Eduardo
##12/Set/2018

cleanData = function(x, cols, p=2){

    if( !class(x) %in% c("list","data.frame") ){
        stop("O dataset de entrada deve ser um objeto da classe 'list' ou da classe 'data.fame.'")
    }

    if( class(x) != "list"){
        x = list(x)
    }
    
    if( any(sapply( seq(length(x)), function(i) ncol(x[[i]]) ) < 3) ){
        stop("O dataset de entrada deve ser um objeto da classe data.frame, com no mínimo três colunas: 'longitude', 'latitude', 'idade'.")
    }

    if ( class(cols) != "character" ){
        stop("'cols' deve ser um objeto da classe 'character', contendo os nomes das colunas referentes aos dados de 'longitude', 'latitude', 'idade' no dataset de entrada.")
    }
    
    ##variaveis locais
    dataset = x
    cols = cols
    precision = p

    #manipulando os dados
    dataset = lapply( seq(length(dataset)), function(x) data.frame(lon=round(dataset[[x]][,cols[1]], precision), lat=round(dataset[[x]][,cols[2]], precision), age=dataset[[x]][,cols[3]]) ) #arredondando para duas casas decimais
    dataset = lapply( seq(length(dataset)), function(x) dataset[[x]][complete.cases(dataset[[x]]),] ) #retirando dados incompletos
    dataset = lapply( seq(length(dataset)), function(x) unique(dataset[[x]]) ) #retirando dados duplicados

    ##ouput da funcao
    return(dataset)

}
