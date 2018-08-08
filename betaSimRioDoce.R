##aleatorizacao de pres-aus em comunidades, com ajuste das matrizes para calculo de beta diversidade
##set-2017


##PARTE 1: abrindo dados e funcoes necessarias


##hora do inicio
starTime = Sys.time()
print('Rodando script...')

##carregando funcoes
source('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/FD/__online_2017_02_01/functions/multidimFbetaD.R')
source('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/FD/__online_2017_02_01/functions/quality_funct_space.R')
source("http://villeger.sebastien.free.fr/R%20scripts/GFD_matcomm.R"); GFD<-GFD_matcomm

##dados
matFunc = read.table('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/Projeto_Rio_Doce/doce_funcional_imputacao.txt') #matriz traits funcionais
matPres = read.table('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/Projeto_Rio_Doce/doce.comunidades.txt') #matriz pres-aus

##crindo objeto para armazenar os resultados e com os cenarios
##outputBetaDiv = list()
scenarios = c(0.05, 0.25, 0.50, 0.75, 1.00) ##proporcao que permancera


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
        mainDir = "/home/anderson/Projetos/Isaac - diversidade beta Rio Doce" #diretorio principal
        scenarioDir = paste('outputs/','scenario_',(1-s)*100,sep='') #diretorio dos resultados da s-esimo cenario
        iterationDir = paste('iteration_',i,sep='') #diretorio dos resultados da i-esima iteracao
        dir.create(path=file.path(mainDir, scenarioDir, iterationDir), showWarnings=FALSE, recursive=TRUE) #criando diretorio (se ele nao existe)
        setwd(dir=file.path(mainDir, scenarioDir, iterationDir)) #mudando diretorio

        ##aleatorizando dados da comundiade para a iteracao atual
        presSimulated = sample(doceCommScen)  #especies aleatorias para Rio Doce
        matPresSimulated = rbind('DOCEsim'=presSimulated, matPres[row.names(matPres)!='DOCE',]) #juntando RioDoce simulado com as demais
        matPresSimulatedNonZero =  matPresSimulated[,which(colSums(matPresSimulated) != 0)] #retirando colunas == 0 (sps que nao aparecem em comunidade alguma)

        ##PARTE 2C: calculo das netricas de beta div a partir dos cenarios

        ##matPresSimulated = matPres #matriz de pres-aus da simulacao atual

        matFuncNonZero = matFunc[match(names(matPresSimulatedNonZero),row.names(matFunc)),] #matriz atributos funcionais iteracao atual

        ##matriz distancia funcional
        matDisFunc = quality_funct_space(mat_funct=matFuncNonZero, nbdim=3,metric='Gower', dendro=FALSE)$details_funct_space$mat_coord 

        ##metricas da beta div
        indices = match(colnames(as.matrix(matPresSimulatedNonZero)),rownames(matDisFunc)) #indices (para ajustar as duas matrizes)
        statsBeta = multidimFbetaD(coord=matDisFunc[indices,], occ=as.matrix(matPresSimulatedNonZero), nm_asb_plot=row.names(matPresSimulatedNonZero)) #calculo das metricas de beta-div

        ##salvando as metricas calculadas
        ##outputBetaDiv[[i]] = statsBeta
        save(statsBeta, file=paste('cenario_',(1-s)*100,'_iteracao_',i,sep='')) #salvando objetos em um .Rdata
        rm(list=c('matPresSimulated','matPresSimulatedNonZero','matFuncNonZero','matDisFunc','statsBeta'))
        gc()

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
        mainDir = "/home/anderson/Projetos/Isaac - diversidade beta Rio Doce" #diretorio principal
        scenarioDir = paste('outputs/','scenario_',(1-s)*100,sep='') #diretorio dos resultados da s-esimo cenario
        iterationDir = paste('iteration_',i,sep='') #diretorio dos resultados da i-esima iteracao
        output_i = file.path(mainDir, scenarioDir, iterationDir) #diretorio com os resultados

        ##carregando os dados salvos
        ## if(exists('statsBeta')){
        ##     rm(statsBeta) #garatindo que nao ocorra confusao...
        ## }else{
        ##     next
        ## }

        ##        if(file.exists(paste(output_i,'/cenario_',(1-s)*100,'_iteracao_',i,sep=''))){
        load(paste(output_i,'/cenario_',(1-s)*100,'_iteracao_',i,sep='')) #abrindo
        ##        }else{
        ##            next
        ##        }
        
        ##transformando as matrizes em vetores (i.e. pegando todas as metricas da matriz e colocando juntas num vetor)
        nestVec = append(x=nestVec, values=as.vector(statsBeta$F_nest_res))
        betaVec = append(x=betaVec, values=as.vector(statsBeta$F_beta))
        turnVec = append(x=turnVec, values=as.vector(statsBeta$F_turn))
        pturnVec = append(x=turnVec, values=as.vector(statsBeta$F_pturn))

        rm(statsBeta) #garatindo que nao ocorra confusao...
        
    }
    
    nestDF = rbind(nestDF, data.frame(scenario=rep(paste('Nes',(1-s)*100,sep=''),length(nestVec)), Nestedness=nestVec))
    betaDF = rbind(betaDF, data.frame(scenario=rep(paste('Beta',(1-s)*100,sep=''),length(betaVec)), Beta=betaVec))
    turnDF = rbind(turnDF, data.frame(scenario=rep(paste('Tur',(1-s)*100,sep=''),length(turnVec)), Turnover=turnVec))
    pturnDF = rbind(pturnDF, data.frame(scenario=rep(paste('Ptur',(1-s)*100,sep=''),length(pturnVec)), P_turnover=pturnVec))
    
}

##exportando os data.frames em arquivos .csv
##diretorio = '/home/anderson/Projetos/Isaac - diversidade beta Rio Doce' #especificar diretorio onde salvar

## write.csv(nestDF, file=file.path(diretorio,'nestData.csv'), row.names=FALSE)
## write.csv(betaDF, file=file.path(diretorio,'betaData.csv'), row.names=FALSE)
## write.csv(turnDF, file=file.path(diretorio,'turnData.csv'), row.names=FALSE)
## write.csv(pturnDF, file=file.path(diretorio,'pturnData.csv'), row.names=FALSE)

##importando arquivos csv
diretorio = '/home/anderson/Projetos/Isaac - diversidade beta Rio Doce' #especificar diretorio onde salvar

nestDF = read.csv(file=file.path(diretorio,'nestData.csv'), header = TRUE)
betaDF = read.csv(file=file.path(diretorio,'betaData.csv'), header = TRUE)
turnDF = read.csv(file=file.path(diretorio,'turnData.csv'), header = TRUE)
pturnDF = read.csv(file=file.path(diretorio,'pturnData.csv'), header = TRUE)

##boxplots
jpeg(file.path(diretorio,'boxplotNest.jpg'), width=600, height=600)
par(family='times', mar=c(5,5,3,2), cex.axis=2, cex.lab=2.5)
nestDF$scenario = factor(nestDF$scenario,
                         levels=c('Nes0','Nes25','Nes50','Nes75','Nes95')) #controle da ordem dos fatores no eixo x
boxplot(nestDF$Nestedness~nestDF$scenario,
        xlim=c(0,6),
        ylim=c(0,1.1),
        outcol=rgb(0,0,0,0.2),
        outcex=0.8,
        outlwd=0.8,
        lwd=2,
        names=c('0%','25%','50%','75%','95%'),
        ylab='Nestedness',
        xlab='Extinction')
text(x=0.2, y=1.05, labels='A', cex=3)
dev.off()

jpeg(file.path(diretorio,'boxplotTurn.jpg'), width=600, height=600)
par(family='times', mar=c(5,5,3,2), cex.axis=2, cex.lab=2.5)
turnDF$scenario = factor(turnDF$scenario,
                         levels=c('Tur0','Tur25','Tur50','Tur75','Tur95')) #controle da ordem dos fatores no eixo x
boxplot(turnDF$Turnover~turnDF$scenario,
        xlim=c(0,6),
        ylim=c(0,1.1),
        outcol=rgb(0,0,0,0.2),
        outcex=0.8,
        outlwd=0.8,
        lwd=2,
        names=c('0%','25%','50%','75%','95%'),
        ylab='Turnover', xlab='Extinction')
text(x=0.2, y=1.05, labels='B', cex=3)
dev.off()

jpeg(file.path(diretorio,'boxplotBeta.jpg'), width=600, height=600)
par(family='times', mar=c(5,5,3,2), cex.axis=2, cex.lab=2.5)
betaDF$scenario = factor(betaDF$scenario,
                         levels=c('Beta0','Beta25','Beta50','Beta75','Beta95')) #controle da ordem dos fatores no eixo x
boxplot(betaDF$Beta~betaDF$scenario,
        xlim=c(0,6),
        ylim=c(0,1.1),
        outcol=rgb(0,0,0,0.2),
        outcex=0.8,
        outlwd=0.8,
        lwd=2,
        names=c('0%','25%','50%','75%','95%'),
        ylab='Beta diversity',
        xlab='Extinction')
text(x=0.2, y=1.05, labels='C', cex=3)
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
par(family='times', mar=c(4.5,5.5,2,2), cex.axis=2, cex.lab=2.5)
plot(nestDens95, lty=3, lwd=2.5, col='red', ylim=c(0,10), main = '', ylab='Density', xlab='Beta Diversity')
lines(nestDens75, lty=3, lwd=2.5, col='purple')
lines(nestDens50, lty=3,lwd=2.5, col='green')
lines(nestDens25, lty=3, lwd=2.5, col='blue')
lines(nestDens0,lty=3, lwd=2.5, col='black')
                                        #
lines(betaDens95, lty=1, lwd=2.5, col='red')
lines(betaDens75, lty=1, lwd=2.5, col='purple')
lines(betaDens50, lty=1, lwd=2.5, col='green')
lines(betaDens25, lty=1, lwd=2.5, col='blue')
lines(betaDens0, lty=1, lwd=2.5, col='black')
                                        #
lines(turnDens95, lty=4, lwd=2.5, col='red')
lines(turnDens75, lty=4, lwd=2.5, col='purple')
lines(turnDens50, lty=4, lwd=2.5, col='green')
lines(turnDens25, lty=4, lwd=2.5, col='blue')
lines(turnDens0, lty=4, lwd=2.5, col='black')
                                        #
legend(x='topright',legend=c('0%','25%','50%','75%','95%','Nestedness', 'Beta diversity', 'Turnover'), title='Extinction', pch=c(22,22,22,22,22,NA,NA,NA), col=c(0,0,0,0,0,1,1,1) , pt.bg=c('black','blue','green','purple','red','black','black','black'), lty=c(NA,NA,NA,NA,NA,3,2,4), lwd=c(NA,NA,NA,NA,NA,2,2,2), cex=1.5, ncol=1)
dev.off()

##graficos separados
jpeg(file.path(diretorio,'densidadeSep_Nest.jpg'), width=600, height=600)
par(family='times', mar=c(4.5,5.5,2,2), cex.axis=2, cex.lab=2.5)
plot(nestDens95, lty=1, lwd=2.5, col='red', ylim=c(0,5), main = '', ylab='Density', xlab='Nestedness')
lines(nestDens75, lty=1, lwd=2.5, col='purple')
lines(nestDens50, lty=1,lwd=2.5, col='green')
lines(nestDens25, lty=1, lwd=2.5, col='blue')
lines(nestDens0,lty=1, lwd=2.5, col='black')
legend(x='topright',legend=c('0%','25%','50%','75%','95%'), title='Extinction', pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','green','purple','red'), cex=2, ncol=2)
dev.off()

jpeg(file.path(diretorio,'densidadeSep_Beta.jpg'), width=600, height=600)
par(family='times', mar=c(4.5,5.5,2,2), cex.axis=2, cex.lab=2.5)
plot(betaDens95, lty=1, lwd=2.5, col='red', ylim=c(0,6), main = '', ylab='Density', xlab='Beta Diversity')
lines(betaDens75, lty=1, lwd=2.5, col='purple')
lines(betaDens50, lty=1, lwd=2.5, col='green')
lines(betaDens25, lty=1, lwd=2.5, col='blue')
lines(betaDens0, lty=1, lwd=2.5, col='black')
legend(x='topright',legend=c('0%','25%','50%','75%','95%'), title='Extinction', pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','green','purple','red'), cex=2, ncol=2)
dev.off()

jpeg(file.path(diretorio,'densidadeSep_Turn.jpg'), width=600, height=600)
par(family='times', mar=c(4.5,5.5,2,2), cex.axis=2, cex.lab=2.5)
plot(turnDens95, lty=1, lwd=2.5, col='red', ylim=c(0,10), main = '', ylab='Density', xlab='Turnover')
lines(turnDens75, lty=1, lwd=2.5, col='purple')
lines(turnDens50, lty=1, lwd=2.5, col='green')
lines(turnDens25, lty=1, lwd=2.5, col='blue')
lines(turnDens0, lty=1, lwd=2.5, col='black')
legend(x='topright',legend=c('0%','25%','50%','75%','95%'), title='Extinction', pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','green','purple','red'), cex=2, ncol=2)
dev.off()


## graficos usando Standarlized Size Effect (SES), ou algo assim : ) ##


##SES = (Xobs - mean (Xnull)) / SD (Xnull)
##
##The Xobs represents the observed values for TβD and FβD, Xnull are the values found from
##the null models. Positive SES values indicates that the observed dissimilarity is higher
##than what is expected by chance. Negative SES values indicates observed dissimilarity lower
##than by chance (Swenson, 2014).
##

##importando arquivos csv
diretorio = '/home/anderson/Projetos/Isaac - diversidade beta Rio Doce' #especificar diretorio onde salvar

nestDF = read.csv(file=file.path(diretorio,'nestData.csv'), header = TRUE)
betaDF = read.csv(file=file.path(diretorio,'betaData.csv'), header = TRUE)
turnDF = read.csv(file=file.path(diretorio,'turnData.csv'), header = TRUE)
pturnDF = read.csv(file=file.path(diretorio,'pturnData.csv'), header = TRUE)


##SES nestedness##
##95% de extincao
nestNull95 = nestDF[which(nestDF$scenario == 'Nes95'),]
nestNull95Mean = mean(nestNull95$Nestedness, na.rm=TRUE)
nestNull95SD = sd(nestNull95$Nestedness, na.rm=TRUE)
SESnest95 = ( nestDF[which(nestDF$scenario == 'Nes0'),'Nestedness'] - nestNull95Mean ) / nestNull95SD

jpeg(paste(diretorio, '/SESnest95.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESnest95, main='', xlab='Nestedness', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##75% de extincao
nestNull75 = nestDF[which(nestDF$scenario == 'Nes75'),]
nestNull75Mean = mean(nestNull75$Nestedness, na.rm=TRUE)
nestNull75SD = sd(nestNull75$Nestedness, na.rm=TRUE)
SESnest75 = ( nestDF[which(nestDF$scenario == 'Nes0'),'Nestedness'] - nestNull75Mean ) / nestNull75SD

jpeg(paste(diretorio, '/SESnest75.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESnest75, main='', xlab='Nestedness', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##50% de extincao
nestNull50 = nestDF[which(nestDF$scenario == 'Nes50'),]
nestNull50Mean = mean(nestNull50$Nestedness, na.rm=TRUE)
nestNull50SD = sd(nestNull50$Nestedness, na.rm=TRUE)
SESnest50 = ( nestDF[which(nestDF$scenario == 'Nes0'),'Nestedness'] - nestNull50Mean ) / nestNull50SD

jpeg(paste(diretorio, '/SESnest50.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESnest50, main='', xlab='Nestedness', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##25% de extincao
nestNull25 = nestDF[which(nestDF$scenario == 'Nes25'),]
nestNull25Mean = mean(nestNull25$Nestedness, na.rm=TRUE)
nestNull25SD = sd(nestNull25$Nestedness, na.rm=TRUE)
SESnest25 = ( nestDF[which(nestDF$scenario == 'Nes0'),'Nestedness'] - nestNull25Mean ) / nestNull25SD

jpeg(paste(diretorio, '/SESnest25.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESnest25, main='', xlab='Nestedness', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##SES Beta diversity##
##95% de extincao
betaNull95 = betaDF[which(betaDF$scenario == 'Beta95'),]
betaNull95Mean = mean(betaNull95$Beta, na.rm=TRUE)
betaNull95SD = sd(betaNull95$Beta, na.rm=TRUE)
SESbeta95 = ( betaDF[which(betaDF$scenario == 'Beta0'),'Beta'] - betaNull95Mean ) / betaNull95SD

jpeg(paste(diretorio, '/SESbeta95.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESbeta95, main='', xlab='Beta diversity', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##75% de extincao
betaNull75 = betaDF[which(betaDF$scenario == 'Beta75'),]
betaNull75Mean = mean(betaNull75$Beta, na.rm=TRUE)
betaNull75SD = sd(betaNull75$Beta, na.rm=TRUE)
SESbeta75 = ( betaDF[which(betaDF$scenario == 'Beta0'),'Beta'] - betaNull75Mean ) / betaNull75SD

jpeg(paste(diretorio, '/SESbeta75.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESbeta75, main='', xlab='Beta diversity', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##50% de extincao
betaNull50 = betaDF[which(betaDF$scenario == 'Beta50'),]
betaNull50Mean = mean(betaNull50$Beta, na.rm=TRUE)
betaNull50SD = sd(betaNull50$Beta, na.rm=TRUE)
SESbeta50 = ( betaDF[which(betaDF$scenario == 'Beta0'),'Beta'] - betaNull50Mean ) / betaNull50SD

jpeg(paste(diretorio, '/SESbeta50.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESbeta50, main='', xlab='Beta diversity', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##25% de extincao
betaNull25 = betaDF[which(betaDF$scenario == 'Beta25'),]
betaNull25Mean = mean(betaNull25$Beta, na.rm=TRUE)
betaNull25SD = sd(betaNull25$Beta, na.rm=TRUE)
SESbeta25 = ( betaDF[which(betaDF$scenario == 'Beta0'),'Beta'] - betaNull25Mean ) / betaNull25SD

jpeg(paste(diretorio, '/SESbeta25.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESbeta25, main='', xlab='Beta diversity', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##SES Turnover##
##95% de extincao
turnNull95 = turnDF[which(turnDF$scenario == 'Tur95'),]
turnNull95Mean = mean(turnNull95$Turnover, na.rm=TRUE)
turnNull95SD = sd(turnNull95$Turnover, na.rm=TRUE)
SESturn95 = ( turnDF[which(turnDF$scenario == 'Tur0'),'Turnover'] - turnNull95Mean ) / turnNull95SD

jpeg(paste(diretorio, '/SESturn95.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESturn95, main='', xlab='Turnover', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##75% de extincao
turnNull75 = turnDF[which(turnDF$scenario == 'Tur75'),]
turnNull75Mean = mean(turnNull75$Turnover, na.rm=TRUE)
turnNull75SD = sd(turnNull75$Turnover, na.rm=TRUE)
SESturn75 = ( turnDF[which(turnDF$scenario == 'Tur0'),'Turnover'] - turnNull75Mean ) / turnNull75SD

jpeg(paste(diretorio, '/SESturn75.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESturn75, main='', xlab='Turnover', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##50% de extincao
turnNull50 = turnDF[which(turnDF$scenario == 'Tur50'),]
turnNull50Mean = mean(turnNull50$Turnover, na.rm=TRUE)
turnNull50SD = sd(turnNull50$Turnover, na.rm=TRUE)
SESturn50 = ( turnDF[which(turnDF$scenario == 'Tur0'),'Turnover'] - turnNull50Mean ) / turnNull50SD

jpeg(paste(diretorio, '/SESturn50.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESturn50, main='', xlab='Turnover', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##25% de extincao
turnNull25 = turnDF[which(turnDF$scenario == 'Tur25'),]
turnNull25Mean = mean(turnNull25$Turnover, na.rm=TRUE)
turnNull25SD = sd(turnNull25$Turnover, na.rm=TRUE)
SESturn25 = ( turnDF[which(turnDF$scenario == 'Tur0'),'Turnover'] - turnNull25Mean ) / turnNull25SD

jpeg(paste(diretorio, '/SESturn25.jpg', sep=''))
par(family='times', mar=c(4,4.5,2,2))
hist(SESturn25, main='', xlab='Turnover', cex.axis=1.5, cex.lab=2,, col='black', border='white')
abline(v=0, lty=2, lwd=2, col='red')
box()
dev.off()

##tabela

SEStableBeta = data.frame( Metric = rep('Beta diversity', 4),
                          Extinction_scenario = c('25','50','75','95'),
                          Observed_mean = c(mean(betaDF$Beta, na.rm=TRUE), mean(betaDF$Beta, na.rm=TRUE), mean(betaDF$Beta, na.rm=TRUE), mean(betaDF$Beta, na.rm=TRUE)),
                          Observed_var = c(var(betaDF$Beta, na.rm=TRUE), var(betaDF$Beta, na.rm=TRUE), var(betaDF$Beta, na.rm=TRUE), var(betaDF$Beta, na.rm=TRUE)),
                          mean_SES = c(mean(SESbeta25, na.rm=TRUE), mean(SESbeta50, na.rm=TRUE), mean(SESbeta75, na.rm=TRUE), mean(SESbeta95, na.rm=TRUE)),
                          var_SES =  c(var(SESbeta25, na.rm=TRUE), var(SESbeta50, na.rm=TRUE), var(SESbeta75, na.rm=TRUE), var(SESbeta95, na.rm=TRUE))
                          )

SEStableTurn = data.frame( Metric = rep('Turnover', 4),
                          Extinction_scenario = c('25','50','75','95'),
                          Observed_mean = c(mean(turnDF$Turnover, na.rm=TRUE), mean(turnDF$Turnover, na.rm=TRUE), mean(turnDF$Turnover, na.rm=TRUE), mean(turnDF$Turnover, na.rm=TRUE)),
                          Observed_var = c(var(turnDF$Turnover, na.rm=TRUE), var(turnDF$Turnover, na.rm=TRUE), var(turnDF$Turnover, na.rm=TRUE), var(turnDF$Turnover, na.rm=TRUE)),
                          mean_SES = c(mean(SESturn25, na.rm=TRUE), mean(SESturn50, na.rm=TRUE), mean(SESturn75, na.rm=TRUE), mean(SESturn95, na.rm=TRUE)),
                          var_SES =  c(var(SESturn25, na.rm=TRUE), var(SESturn50, na.rm=TRUE), var(SESturn75, na.rm=TRUE), var(SESturn95, na.rm=TRUE))
                          )

SEStableNest = data.frame( Metric = rep('Nestedness', 4),
                          Extinction_scenario = as.factor(c("25","50","75","95")),
                          Observed_mean = c(mean(nestDF$Nestedness, na.rm=TRUE), mean(nestDF$Nestedness, na.rm=TRUE), mean(nestDF$Nestedness, na.rm=TRUE), mean(nestDF$Nestedness, na.rm=TRUE)),
                          Observed_var = c(var(nestDF$Nestedness, na.rm=TRUE), var(nestDF$Nestedness, na.rm=TRUE), var(nestDF$Nestedness, na.rm=TRUE), var(nestDF$Nestedness, na.rm=TRUE)),
                          mean_SES = c(mean(SESnest25, na.rm=TRUE), mean(SESnest50, na.rm=TRUE), mean(SESnest75, na.rm=TRUE), mean(SESnest95, na.rm=TRUE)),
                          var_SES =  c(var(SESnest25, na.rm=TRUE), var(SESnest50, na.rm=TRUE), var(SESnest75, na.rm=TRUE), var(SESnest95, na.rm=TRUE))
                          )

##juntando numa tabela unica
SEStable = rbind(SEStableBeta, SEStableTurn, SEStableNest)
names(SEStable) = c('Metric', 'Extinction scenario (in %)', 'Observed (mean)', 'Observed (Var)', 'SES (mean)', 'SES (Var)') #ajustando nomes

write.csv(SEStable, paste(diretorio,'/SEStable.csv',sep=''), row.names=TRUE)
