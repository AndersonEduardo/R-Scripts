## uncertest - funcao para analise da sensibilidade de um modelo estatistico (e.g., GLM) ao erro nas medidas das variaveis (preditoras e resposta). 

sensinc = function(dados, vary_name, varx_names, errorPars, modelo){

    ##objetos locais
    dataSet = dados
    vary_name = vary_name
    varx_names =varx_names
    errorPars = errorPars
    modelo = modelo
    ss = data.frame()
    params = data.frame()
    PCCData = data.frame()
    output = list()

    ##1000 iteracoes
    for (iter in seq(1000)) {
        
        dataError_i = list() #lista vazia para o dataset

        ##construcao do dataset com variacoes nas medidas (erro)
        dataError_i = lapply(X = 1:length(errorPars),
                             FUN=function(x) if(errorPars[[x]][1]=='normal'){
                                                 dataError_i[[x]] = dataSet[,x] + rnorm(n=nrow(dataSet), mean=0, sd=as.numeric(errorPars[[x]][2]))
                                             }else{
                                                 dataError_i[[x]] = dataSet[,x] + runif( n=nrow(dataSet), min=-1*as.numeric(errorPars[[x]][2]), max=as.numeric(errorPars[[x]][2])) } )

        ##contruindo dataframe para o dataset construido
        dataError_i = data.frame( matrix(data=unlist(dataError_i), nrow=length(unlist(dataError_i))/length(errorPars), byrow=FALSE) )
        names(dataError_i) = names(dataSet) ##ajustando os nomes

        ##consolidando o dataset construido
        ##dataUncert = data.frame(vary = dataSet[,vary_name], dataError_i)
        dataUncert = dataError_i

        ##quadrado dos desvios (desvio em relacao aos dados originais)
        ss = rbind(ss,
                   colSums( (dataUncert - dataSet)^2, na.rm=TRUE))

        ##rodando o modelo com o novo dataset
        model_i = update(modelo, . ~ ., data=dataUncert)

        params = rbind(params,
                       coef(model_i))

    }

    ##ajustanado os nomes
    names(ss) = paste('SS_',names(dataSet),sep='')
    names(params) = names(coef(modelo))

    ##consolidando o output
    PCCData = rbind(PCCData,
                      cbind(params, ss))

    ##PCC
    numberOfParms = length(coef(modelo))

    for (i in seq(numberOfParms)){

        currentParm = names(coef(modelo))[i]
        PCCDataSS = PCCData[,-c(seq(numberOfParms))]

        PCCData_i = cbind( PCCData[,i], PCCDataSS )
        names(PCCData_i) = c( names(PCCData)[i], names(PCCDataSS) )

        output[[names(coef(modelo))[i]]] = pcc(X=PCCData_i[,-1], y=PCCData_i[,1], nboot=1000)
        
    }
    
    ##saida da funcao
    return(output)
}
