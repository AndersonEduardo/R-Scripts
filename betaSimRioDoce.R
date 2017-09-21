##aleatorizacao de pres-aus em comunidades, com ajuste das matrizes para calculo de beta diversidade
##set-2017


##PARTE 1: abrindo dados e funcoes necessarias

##hora do inicio
starTime = Sys.time()
print('Rodando script...')

##funcoes
source('/home/anderson/Documentos/Projetos/Rio Doce Isac/FD/__online_2017_02_01/functions/multidimFbetaD.R')
source('/home/anderson/Documentos/Projetos/Rio Doce Isac/FD/__online_2017_02_01/functions/quality_funct_space.R')

##dados
matFunc = read.table('/home/anderson/Documentos/Projetos/Rio Doce Isac/Projeto_Rio_Doce/doce_funcional_imputacao.txt') #matriz traits funcionais
matPres = read.table('/home/anderson/Documentos/Projetos/Rio Doce Isac/Projeto_Rio_Doce/doce.comunidades.txt') #matriz pres-aus

##crindo objeto para armazenar os resultados e com os cenarios
outputBetaDiv = list()
scenarios = c(0.05, 0.25, 0.50, 0.75, 1.00)


##PARTE 2: realizando as iteracoes dos cenarios e o calculo da beta div para cada um


for (s in scenarios){

    ##PARTE 2A: criando cenario de extincoes

    ##criando a comundade do cenario atual    
    doceComm = matPres[row.names(matPres)=='DOCE',] #comunidade Rio Doce
    doceCommScen = rep(0, length(doceComm)) #vetor de zeros para o Rio Doce
    doceCommScen[sample(which(doceComm==1))[1:round(s*length(doceComm[which(doceComm==1)]))]] = 1 #eliminacao aleatoria de especies para o cenario atual

        for (i in 1:100){ #for (i in 1:100){
        
        ##PARTE 2B: aleatorizacao

        ##sistema de diretorios para a iteracao atual
        mainDir = "/home/anderson/Documentos/Projetos/Rio Doce Isac" #diretorio principal
        scenarioDir = paste('outputs/','scenario_',(1-s)*100,sep='') #diretorio dos resultados da s-esimo cenario
        iterationDir = paste('iteration_',i,sep='') #diretorio dos resultados da i-esima iteracao
        dir.create(path=file.path(mainDir, scenarioDir, iterationDir), showWarnings=FALSE, recursive=TRUE) #criando diretorio (se ele nao existe)
        setwd(dir=file.path(mainDir, scenarioDir, iterationDir)) #mudando diretorio

        ##aleatorizando dados da comundiade para a iteracao atual
        presSimulated = sample(doceCommScen)  #especies aleatorias para Rio Doce
        matPresSimulated = rbind('DOCEsim'=presSimulated, matPres[row.names(matPres)!='DOCE',]) #juntando RioDoce simulado com as demais
        matPresSimulatedNonZero =  matPresSimulated[,which(colSums(matPresSimulated) != 0)] #retirando colunas == 0 (sps que nao aparecem em comunidade alguma)

        ##rm(matPresSimulated)

        ##PARTE 2C: calculo das netricas de beta div a partir dos cenarios

        ##matPresSimulated = matPres #matriz de pres-aus da simulacao atual

        matFuncNonZero = matFunc[match(names(matPresSimulatedNonZero),row.names(matFunc)),] #matriz atributos funcionais iteracao atual

        ##matriz distancia funcional
        matDisFunc = quality_funct_space(mat_funct=matFuncNonZero, nbdim=3,metric='Gower', dendro=FALSE)$details_funct_space$mat_coord 

        ##metricas da beta div
        indices = match(colnames(as.matrix(matPresSimulatedNonZero)),rownames(matDisFunc)) #indices (para ajustar as duas matrizes)
        statsBeta = multidimFbetaD(coord=matDisFunc[indices,], occ=as.matrix(matPresSimulatedNonZero), nm_asb_plot=row.names(matPresSimulatedNonZero)) #calculo das metricas de beta-div

        ##salvando as metricas calculadas
        #outputBetaDiv[[i]] = statsBeta
        save(statsBeta, file=paste('cenario_',(1-s)*100,'_iteracao_',i,sep='')) #salvando objetos em um .Rdata

    }
}

##registro do tempo gasto
print('Concluido com sucesso!')         
timeTaken = Sys.time() - starTime
round(timeTaken,2)

