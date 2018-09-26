
##install.packages('xlsx') precisa desse pacote
library(xlsx)


## 1) abrindo dados 'brutos'
dadosBrutos = read.xlsx('/home/anderson/Downloads/DADOS.xlsx', sheetIndex=1, header=FALSE) 
dados = dadosBrutos #so pra garantir que pra qualquer problema tem um backup interno


## 2) manipulacao da matrix de presenca/ausencia
occMat =  data.frame(apply(dados[-c(1:3),-1], c(1,2),  as.numeric)) #matrix de presenca (tirando as primeiras linhas do dataset)
occMat = cbind(sps=as.character(dados[-c(1:3),1]), occMat) #juntando os nomes das sps
names(occMat) = c( 'sps', as.character( unlist(dados[3,-1]) ) ) #ajuste de nomes das colunas


## 3) construcao da planilha no formato 'base de dados'
base_de_dados = data.frame() #planilha que serao salvo o output

for (col_i in 2:ncol(occMat)){ #loop sobre as colunas do dataset 

    if ( sum(occMat[,col_i], na.rm=TRUE) > 0 ){ #procura por colunas com registros

        infoColeta = dados[1:3,col_i] #coletando informacao do ponto de coleta
        sps_i = occMat[ which(occMat[,col_i]>0), 'sps'] #coletando nomes das especies

        ##montando o output
        base_de_dados = rbind(base_de_dados,
                              data.frame(CAMPANHA = infoColeta[1],
                                         REPLICA = infoColeta[2],
                                         PONTO = infoColeta[3],
                                         ESPECIE = sps_i,
                                         QTD = occMat[which(occMat[,col_i]>0), col_i]) )
    }
}

## 4) salvando o output
write.csv(base_de_dados, 'base_de_dados.csv', row.names=FALSE) #verifiar qual a pasta de trabalho antes de salvar
