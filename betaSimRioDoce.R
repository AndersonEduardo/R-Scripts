##aleatorizacao de pres-aus em comunidades, com ajuste das matrizes para calculo de beta diversidade
##set-2017


##PARTE 1: abrindo dados e funcoes necessarias


##hora do inicio
starTime = Sys.time()
print('Rodando script...')

##carregando funcoes
source('/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac/FD/__online_2017_02_01/functions/multidimFbetaD.R')
source('/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac/FD/__online_2017_02_01/functions/quality_funct_space.R')

##dados
matFunc = read.table('/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac/Projeto_Rio_Doce/doce_funcional_imputacao.txt') #matriz traits funcionais
matPres = read.table('/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac/Projeto_Rio_Doce/doce.comunidades.txt') #matriz pres-aus

##crindo objeto para armazenar os resultados e com os cenarios
#outputBetaDiv = list()
scenarios = c(0.05, 0.25, 0.50, 0.75, 1.00)


##PARTE 2: realizando as iteracoes dos cenarios e o calculo da beta div para cada um


for (s in scenarios){

    ##PARTE 2A: criando cenario de extincoes

    ##criando a comundade do cenario atual    
    doceComm = matPres[row.names(matPres)=='DOCE',] #comunidade Rio Doce
    doceCommScen = rep(0, length(doceComm)) #vetor de zeros para o Rio Doce
    doceCommScen[sample(which(doceComm==1))[1:round(s*length(doceComm[which(doceComm==1)]))]] = 1 #eliminacao aleatoria de especies para o cenario atual

        for (i in 1:100){
        
        ##PARTE 2B: aleatorizacao

        ##sistema de diretorios para a iteracao atual
        mainDir = "/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac" #diretorio principal
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


##PARTE 3: graficos


rm(list=ls())

##definicao de variaveis e parametros
nestVec = array()
betaVec = array()
turnVec = array()
pturnVec = array()
nestDF = data.frame()
betaDF = data.frame()
turnDF = data.frame()
pturnDF = data.frame()
scenarios = c(0.05, 0.25, 0.50, 0.75, 1.00)

##algoritmo recursivo para extracao e organizacao dos dados armazenados no hd
##loop sobre os cenarios
for (s in scenarios){

    ##loop sobre as iteracoes
    for (i in 1:100){
        
        ##sistema de diretorios para a iteracao atual
        mainDir = "/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac" #diretorio principal
        scenarioDir = paste('outputs/','scenario_',(1-s)*100,sep='') #diretorio dos resultados da s-esimo cenario
        iterationDir = paste('iteration_',i,sep='') #diretorio dos resultados da i-esima iteracao
        output_i = file.path(mainDir, scenarioDir, iterationDir) #diretorio com os resultados

        ##carregando os dados salvos
        rm(statsBeta) #garatindo que nao ocorra confusao...

        if(file.exists(paste(output_i,'/cenario_',(1-s)*100,'_iteracao_',i,sep=''))){
            load(paste(output_i,'/cenario_',(1-s)*100,'_iteracao_',i,sep='')) #abrindo
        }else{
            next
        }
        
        ##transformando as matrizes em vetores (i.e. pegando todas as metricas da matriz e colocando juntas num vetor)
        nestVec = append(x=nestVec, values=as.vector(statsBeta$F_nest_res))
        betaVec = append(x=betaVec, values=as.vector(statsBeta$F_beta))
        turnVec = append(x=turnVec, values=as.vector(statsBeta$F_turn))
        pturnVec = append(x=turnVec, values=as.vector(statsBeta$F_pturn))

    }
    
    nestDF = rbind(nestDF, data.frame(scenario=rep(paste('Nes',(1-s)*100,sep=''),length(nestVec)), Nestedness=nestVec))
    betaDF = rbind(betaDF, data.frame(scenario=rep(paste('Beta',(1-s)*100,sep=''),length(betaVec)), Beta=betaVec))
    turnDF = rbind(turnDF, data.frame(scenario=rep(paste('Tur',(1-s)*100,sep=''),length(turnVec)), Turnover=turnVec))
    pturnDF = rbind(pturnDF, data.frame(scenario=rep(paste('Ptur',(1-s)*100,sep=''),length(pturnVec)), P_turnover=pturnVec))
    
}

##exporando os data.frames em arquivos .csv
diretorio = '/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac' #especificar diretorio onde salvar

write.csv(nestDF, file=file.path(diretorio,'nestData.csv'), row.names=FALSE)
write.csv(betaVec, file=file.path(diretorio,'betaData.csv'), row.names=FALSE)
write.csv(turnDF, file=file.path(diretorio,'turnData.csv'), row.names=FALSE)
write.csv(pturnDF, file=file.path(diretorio,'pturnData.csv'), row.names=FALSE)

##importando arquivos csv
nestDF = read.csv(file=file.path(diretorio,'nestData.csv'), header = TRUE)
betaVec = read.csv(file=file.path(diretorio,'betaData.csv'), header = TRUE)
turnDF = read.csv(file=file.path(diretorio,'turnData.csv'), header = TRUE)
pturnDF = read.csv(file=file.path(diretorio,'pturnData.csv'), header = TRUE)


##criando os graficos
diretorio = '/home/anderson/Documentos/Projetos/Diversidade beta Rio Doce Isaac' #especificar diretorio onde salvar os graficos

##boxplots
jpeg(file.path(diretorio,'boxplotNest.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5)
nestDF$scenario = factor(nestDF$scenario,levels=c('Nes0','Nes25','Nes50','Nes75','Nes95')) #controle da ordem dos fatores no eixo x
boxplot(nestDF$Nestedness~nestDF$scenario, xlim=c(0,6), ylim=c(0,1.1), outcol=rgb(0,0,0,0.5)) + text(x=0.2, y=1.05, labels='A', cex=3)
dev.off()

jpeg(file.path(diretorio,'boxplotTurn.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5)
turnDF$scenario = factor(turnDF$scenario,levels=c('Tur0','Tur25','Tur50','Tur75','Tur95')) #controle da ordem dos fatores no eixo x
boxplot(turnDF$Turnover~turnDF$scenario, xlim=c(0,6), ylim=c(0,1.1), outcol=rgb(0,0,0,0.5)) + text(x=0.2, y=1.05, labels='B', cex=3)
dev.off()

jpeg(file.path(diretorio,'boxplotBeta.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5)
betaDF$scenario = factor(betaDF$scenario,levels=c('Beta0','Beta25','Beta50','Beta75','Beta95')) #controle da ordem dos fatores no eixo x
boxplot(betaDF$Beta~betaDF$scenario, xlim=c(0,6), ylim=c(0,1.1), outcol=rgb(0,0,0,0.5)) + text(x=0.2, y=1.05, labels='C', cex=3)
dev.off()

###densidades
nestDens95 = density(nestDF[complete.cases(nestDF$Nestedness) & nestDF$scenario=='Nes95','Nestedness'])
nestDens75 = density(nestDF[complete.cases(nestDF$Nestedness) & nestDF$scenario=='Nes75','Nestedness'])
nestDens50 = density(nestDF[complete.cases(nestDF$Nestedness) & nestDF$scenario=='Nes50','Nestedness'])
nestDens25 = density(nestDF[complete.cases(nestDF$Nestedness) & nestDF$scenario=='Nes25','Nestedness'])
nestDens0 = density(nestDF[complete.cases(nestDF$Nestedness) & nestDF$scenario=='Nes0','Nestedness'])

betaDens95 = density(betaDF[complete.cases(betaDF$Beta) & betaDF$scenario=='Beta95','Beta'])
betaDens75 = density(betaDF[complete.cases(betaDF$Beta) & betaDF$scenario=='Beta75','Beta'])
betaDens50 = density(betaDF[complete.cases(betaDF$Beta) & betaDF$scenario=='Beta50','Beta'])
betaDens25 = density(betaDF[complete.cases(betaDF$Beta) & betaDF$scenario=='Beta25','Beta'])
betaDens0 = density(betaDF[complete.cases(betaDF$Beta) & betaDF$scenario=='Beta0','Beta'])

turnDens95 = density(turnDF[complete.cases(turnDF$Turn) & turnDF$scenario=='Tur95','Turnover'])
turnDens75 = density(turnDF[complete.cases(turnDF$Turn) & turnDF$scenario=='Tur75','Turnover'])
turnDens50 = density(turnDF[complete.cases(turnDF$Turn) & turnDF$scenario=='Tur50','Turnover'])
turnDens25 = density(turnDF[complete.cases(turnDF$Turn) & turnDF$scenario=='Tur25','Turnover'])
turnDens0 = density(turnDF[complete.cases(turnDF$Turn) & turnDF$scenario=='Tur0','Turnover'])

##graficos juntos
jpeg(file.path(diretorio,'densidadeUni.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5, cex.lab=1.5)
plot(nestDens95, lty=3, lwd=2, col='yellow', ylim=c(0,10), main = '', ylab='Density', xlab='Beta Diversity')
lines(nestDens75, lty=3, lwd=2, col='purple')
lines(nestDens50, lty=3,lwd=2, col='red')
lines(nestDens25, lty=3, lwd=2, col='blue')
lines(nestDens0,lty=3, lwd=2, col='black')
#
lines(betaDens95, lty=1, lwd=2, col='yellow')
lines(betaDens75, lty=1, lwd=2, col='purple')
lines(betaDens50, lty=1, lwd=2, col='red')
lines(betaDens25, lty=1, lwd=2, col='blue')
lines(betaDens0, lty=1, lwd=2, col='black')
#
lines(turnDens95, lty=4, lwd=2, col='yellow')
lines(turnDens75, lty=4, lwd=2, col='purple')
lines(turnDens50, lty=4, lwd=2, col='red')
lines(turnDens25, lty=4, lwd=2, col='blue')
lines(turnDens0, lty=4, lwd=2, col='black')
#
legend(x='topright',legend=c('0% of extinction','25% of extinction','50% of extinction','75% of extinction','95% of extinction'), pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','red','purple','yellow'), cex=1.5)
dev.off()

##graficos separados
jpeg(file.path(diretorio,'densidadeSep_Nest.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5, cex.lab=1.5)
plot(nestDens95, lty=3, lwd=2, col='yellow', ylim=c(0,5), main = '', ylab='Density', xlab='Beta Diversity')
lines(nestDens75, lty=3, lwd=2, col='purple')
lines(nestDens50, lty=3,lwd=2, col='red')
lines(nestDens25, lty=3, lwd=2, col='blue')
lines(nestDens0,lty=3, lwd=2, col='black')
legend(x='topright',legend=c('0% of extinction','25% of extinction','50% of extinction','75% of extinction','95% of extinction'), pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','red','purple','yellow'), cex=1.5)
dev.off()

jpeg(file.path(diretorio,'densidadeSep_Beta.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5, cex.lab=1.5)
plot(betaDens95, lty=1, lwd=2, col='yellow', ylim=c(0,6), main = '', ylab='Density', xlab='Beta Diversity')
lines(betaDens75, lty=1, lwd=2, col='purple')
lines(betaDens50, lty=1, lwd=2, col='red')
lines(betaDens25, lty=1, lwd=2, col='blue')
lines(betaDens0, lty=1, lwd=2, col='black')
legend(x='topright',legend=c('0% of extinction','25% of extinction','50% of extinction','75% of extinction','95% of extinction'), pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','red','purple','yellow'), cex=1.5)
dev.off()

jpeg(file.path(diretorio,'densidadeSep_Turn.jpg'), width=600, height=600)
par(family='times', cex.axis=1.5, cex.lab=1.5)
plot(turnDens95, lty=4, lwd=2, col='yellow', ylim=c(0,10), main = '', ylab='Density', xlab='Beta Diversity')
lines(turnDens75, lty=4, lwd=2, col='purple')
lines(turnDens50, lty=4, lwd=2, col='red')
lines(turnDens25, lty=4, lwd=2, col='blue')
lines(turnDens0, lty=4, lwd=2, col='black')
legend(x='topright',legend=c('0% of extinction','25% of extinction','50% of extinction','75% of extinction','95% of extinction'), pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','red','purple','yellow'), cex=1.5)
dev.off()
