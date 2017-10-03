##(exercicio) algoritmo para converter uma matriz quadrada em uma matriz triangular superior por eliminacao de Gauss
##Anderson Eduardo
##30/set/2017

elimGauss = function(mat){

    mat = mat
    
    if(nrow(mat)!=ncol(mat)){
        print('ERRO: a matriz precisa ser quadrada.')
        break
    }

    for (i in 1:c(ncol(mat)-1)){
        piv = mat[i,i]
        linhaPivo = mat[i,]
        for(j in i:(nrow(mat)-1)){
            k = mat[j+1,i]/piv
            linhaAtual = mat[j+1,]
            mat[j+1,] =  linhaAtual - linhaPivo*k
        }
    }

    mat = round(mat)

    cat(" >> Algoritmo processado com sucesso!  <<  \n")
    
    return(mat)
    
}

mat = matrix(sample(10*10),10,10)
matU = elimGauss(mat)
matU
