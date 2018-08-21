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
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
nestDF$scenario = factor(nestDF$scenario,
                         levels=c('Nes0','Nes25','Nes50','Nes75','Nes95')) #controle da ordem dos fatores no eixo x
boxplot(nestDF$Nestedness~nestDF$scenario,
        xlim=c(0,14),
        ylim=c(0,1.1),
        at=c(1,4,7,10,13),
        boxwex = 2,
        outcol=rgb(0,0,0,0.2),
        outcex=0.8,
        outlwd=0.8,
        lwd=2,
        names=c('0%','25%','50%','75%','100%'),
        ylab='Nestedness')
text(x=0.2, y=1.05, labels='A', cex=5)
dev.off()

jpeg(file.path(diretorio,'boxplotTurn.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
turnDF$scenario = factor(turnDF$scenario,
                         levels=c('Tur0','Tur25','Tur50','Tur75','Tur95')) #controle da ordem dos fatores no eixo x
boxplot(turnDF$Turnover~turnDF$scenario,
        xlim=c(0,14),
        ylim=c(0,1.1),
        at=c(1,4,7,10,13),
        boxwex = 2,
        outcol=rgb(0,0,0,0.2),
        outcex=0.8,
        outlwd=0.8,
        lwd=2,
        names=c('0%','25%','50%','75%','100%'),
        ylab='Turnover')
text(x=0.2, y=1.05, labels='B', cex=5)
dev.off()

jpeg(file.path(diretorio,'boxplotBeta.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
betaDF$scenario = factor(betaDF$scenario,
                         levels=c('Beta0','Beta25','Beta50','Beta75','Beta95')) #controle da ordem dos fatores no eixo x
boxplot(betaDF$Beta~betaDF$scenario,
        xlim=c(0,14),
        ylim=c(0,1.1),
        at=c(1,4,7,10,13),
        boxwex = 2,
        outcol=rgb(0,0,0,0.2),
        outcex=0.8,
        outlwd=0.8,
        lwd=2,
        names=c('0%','25%','50%','75%','100%'),
        ylab='Beta diversity')
text(x=0.2, y=1.05, labels='C', cex=5)
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
##par(family='times', mar=c(4.5,5.5,2,2), cex.axis=2, cex.lab=2.5)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
plot(nestDens95, lty=3, lwd=2.5, col='red', ylim=c(0,10), main = '', ylab='Density', xlab='Metric')
lines(nestDens75, lty=3, lwd=2.5, col='purple')
lines(nestDens50, lty=3,lwd=2.5, col='green')
lines(nestDens25, lty=3, lwd=2.5, col='blue')
lines(nestDens0,lty=3, lwd=2.5, col='black')
##
lines(turnDens95, lty=4, lwd=2.5, col='red')
lines(turnDens75, lty=4, lwd=2.5, col='purple')
lines(turnDens50, lty=4, lwd=2.5, col='green')
lines(turnDens25, lty=4, lwd=2.5, col='blue')
lines(turnDens0, lty=4, lwd=2.5, col='black')
##
lines(betaDens95, lty=1, lwd=2.5, col='red')
lines(betaDens75, lty=1, lwd=2.5, col='purple')
lines(betaDens50, lty=1, lwd=2.5, col='green')
lines(betaDens25, lty=1, lwd=2.5, col='blue')
lines(betaDens0, lty=1, lwd=2.5, col='black')
##
legend(x='topright',legend=c('0%','25%','50%','75%','100%','Nestedness', 'Turnover','Beta diversity'), title='Extinction', pch=c(22,22,22,22,22,NA,NA,NA), col=c(0,0,0,0,0,1,1,1) , pt.bg=c('black','blue','green','purple','red','black','black','black'), lty=c(NA,NA,NA,NA,NA,3,4,1), lwd=c(NA,NA,NA,NA,NA,2,2,2), cex=1.5, ncol=1)
dev.off()


##graficos separados
jpeg(file.path(diretorio,'densidadeSep_Nest.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
plot(nestDens95, lty=1, lwd=4, col='red', ylim=c(0,5), main = '', ylab='Density', xlab='Nestedness')
lines(nestDens75, lty=1, lwd=4, col='purple')
lines(nestDens50, lty=1,lwd=4, col='green')
lines(nestDens25, lty=1, lwd=4, col='blue')
lines(nestDens0,lty=1, lwd=4, col='black')
legend(x='topright',legend=c('0%','25%','50%','75%','100%'), title='Extinction', pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','green','purple','red'), cex=2.5, ncol=2)
text(x=-0.15, y=4.8, labels='A', cex=5)
dev.off()

jpeg(file.path(diretorio,'densidadeSep_Turn.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
plot(turnDens95, lty=1, lwd=4, col='red', ylim=c(0,10), main = '', ylab='Density', xlab='Turnover')
lines(turnDens75, lty=1, lwd=4, col='purple')
lines(turnDens50, lty=1, lwd=4, col='green')
lines(turnDens25, lty=1, lwd=4, col='blue')
lines(turnDens0, lty=1, lwd=4, col='black')
legend(x='topright',legend=c('0%','25%','50%','75%','100%'), title='Extinction', pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','green','purple','red'), cex=2.5, ncol=2)
text(x=-0.03, y=9.6, labels='B', cex=5)
dev.off()

jpeg(file.path(diretorio,'densidadeSep_Beta.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
plot(betaDens95, lty=1, lwd=4, col='red', ylim=c(0,6), main = '', ylab='Density', xlab='Beta Diversity')
lines(betaDens75, lty=1, lwd=4, col='purple')
lines(betaDens50, lty=1, lwd=4, col='green')
lines(betaDens25, lty=1, lwd=4, col='blue')
lines(betaDens0, lty=1, lwd=4, col='black')
legend(x='topright',legend=c('0%','25%','50%','75%','100%'), title='Extinction', pch=c(22,22,22,22,22), col=c(0,0,0,0,0) , pt.bg=c('black','blue','green','purple','red'), cex=2.5, ncol=2)
text(x=0, y=5.8, labels='C', cex=5)
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

nestDFraw = read.csv(file=file.path(diretorio,'nestData.csv'), header = TRUE)
betaDFraw = read.csv(file=file.path(diretorio,'betaData.csv'), header = TRUE)
turnDFraw = read.csv(file=file.path(diretorio,'turnData.csv'), header = TRUE)
pturnDFraw = read.csv(file=file.path(diretorio,'pturnData.csv'), header = TRUE)

## ##correcao do numero de cenarios - os arquivos tinham 3 vezes o numero de simulacoes que deveriam ter (21*5*100=10500)...
## idxNest = c( which(nestDFraw$scenario=='Nes0')[1:100], which(nestDFraw$scenario=='Nes25')[1:100], which(nestDFraw$scenario=='Nes50')[1:100], which(nestDFraw$scenario=='Nes75')[1:100], which(nestDFraw$scenario=='Nes95')[1:100] )
## idxBeta = c( which(betaDFraw$scenario=='Beta0')[1:100], which(betaDFraw$scenario=='Beta25')[1:100], which(betaDFraw$scenario=='Beta50')[1:100], which(betaDFraw$scenario=='Beta75')[1:100], which(betaDFraw$scenario=='Beta95')[1:100] )
## idxTurn = c( which(turnDFraw$scenario=='Tur0')[1:100], which(turnDFraw$scenario=='Tur25')[1:100], which(turnDFraw$scenario=='Tur50')[1:100], which(turnDFraw$scenario=='Tur75')[1:100], which(turnDFraw$scenario=='Tur95')[1:100] )
## idxPturn = c( which(pturnDFraw$scenario=='Ptur0')[1:100], which(pturnDFraw$scenario=='Ptur25')[1:100], which(pturnDFraw$scenario=='Ptur50')[1:100], which(pturnDFraw$scenario=='Ptur75')[1:100], which(pturnDFraw$scenario=='Ptur95')[1:100] )

## nestDF = nestDFraw[idxNest,]
## betaDF = betaDFraw[idxBeta,]
## turnDF = turnDFraw[idxTurn,]
## pturnDF = pturnDFraw[idxPturn,]

nestDF = nestDFraw
betaDF = betaDFraw
turnDF = turnDFraw
pturnDF = pturnDFraw


##valores observados
rm(statsBeta)
load(paste(diretorio,'/outputs/Observado/cenario_obseravdo',sep='')) #abrindo
nestVecObs = as.vector(statsBeta$F_nest_res)
betaVecObs = as.vector(statsBeta$F_beta)
turnVecObs = as.vector(statsBeta$F_turn)
pturnVecObs = as.vector(statsBeta$F_pturn)


##SES nestedness
nestDFmeans = aggregate(x=nestDF$Nestedness, by=list(Mean=nestDF$scenario), FUN=mean, na.rm=TRUE)
nestDFsd = aggregate(x=nestDF$Nestedness, by=list(Mean=nestDF$scenario), FUN=sd, na.rm=TRUE)
nestDFses = (mean(nestVecObs) - nestDFmeans$x) / nestDFsd$x
SESnest = data.frame(Metric = rep('Nestedness',length(nestDFses)), Extinction=c(0,25,50,75,100), Observed=mean(nestVecObs), Null_mean=nestDFmeans$x, Null_SD=nestDFsd$x, SES=nestDFses)

##SES Beta diversity
betaDFmeans = aggregate(x=betaDF$Beta, by=list(Mean=betaDF$scenario), FUN=mean, na.rm=TRUE)
betaDFsd = aggregate(x=betaDF$Beta, by=list(Mean=betaDF$scenario), FUN=sd, na.rm=TRUE)
betaDFses = (mean(betaVecObs) - betaDFmeans$x) / betaDFsd$x
SESbeta = data.frame(Metric = rep('Beta',length(betaDFses)), Extinction=c(0,25,50,75,100), Observed=mean(betaVecObs), Null_mean=betaDFmeans$x, Null_SD=betaDFsd$x, SES=betaDFses)

##SES Turnover
turnDFmeans = aggregate(x=turnDF$Turnover, by=list(Mean=turnDF$scenario), FUN=mean, na.rm=TRUE)
turnDFsd = aggregate(x=turnDF$Turnover, by=list(Mean=turnDF$scenario), FUN=sd, na.rm=TRUE)
turnDFses = (mean(turnVecObs) - turnDFmeans$x) / turnDFsd$x
SESturn = data.frame(Metric = rep('Turnover',length(turnDFses)), Extinction=c(0,25,50,75,100), Observed=mean(turnVecObs), Null_mean=turnDFmeans$x, Null_SD=turnDFsd$x, SES=turnDFses)

##SES pturnover
pturnDFmeans = aggregate(x=pturnDF$P_turnover, by=list(Mean=pturnDF$scenario), FUN=mean, na.rm=TRUE)
pturnDFsd = aggregate(x=pturnDF$P_turnover, by=list(Mean=pturnDF$scenario), FUN=sd, na.rm=TRUE)
pturnDFses = (mean(pturnVecObs) - pturnDFmeans$x) / pturnDFsd$x
SESpturn = data.frame(Metric = rep('Pturnover',length(pturnDFses)), Extinction=c(0,25,50,75,100), Observed=mean(pturnVecObs), Null_mean=pturnDFmeans$x, Null_SD=pturnDFsd$x, SES=pturnDFses)

##tabela de output
SEStab = rbind(SESnest, SESbeta, SESturn, SESpturn)
write.csv(SEStab, paste(diretorio,'/SEStable.csv',sep=''), row.names=FALSE)

##grafico de barras
SEStabGraph = data.frame(Nestedness=nestDFses, Turnover=turnDFses, Beta=betaDFses)
rownames(SEStabGraph) = c("0%","25%","50%","75%","100%")

jpeg(paste(diretorio,'/SESbarras.jpg', sep=''), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.5, cex.lab=3.3)
barplot( t(SEStabGraph), beside=TRUE, ylab="Metric", ylim=c(-1,0.5), col=c('black','darkgray','lightgray'))
legend('topleft', legend=c("Nestedness","Beta Diversity","Turnover"), fill=c('black','darkgray','lightgray'), cex=1.5)
abline(h=0, lwd=0.5)
box()
dev.off()



########################################
######### observado ####################

setwd(paste(diretorio,'/outputs/Observado',sep=''))

##carregando funcoes
source('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/FD/__online_2017_02_01/functions/multidimFbetaD.R')
source('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/FD/__online_2017_02_01/functions/quality_funct_space.R')
source("http://villeger.sebastien.free.fr/R%20scripts/GFD_matcomm.R"); GFD<-GFD_matcomm

##dados
matFunc = read.table('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/Projeto_Rio_Doce/doce_funcional_imputacao.txt') #matriz traits funcionais
matPres = read.table('/home/anderson/Projetos/Isaac - diversidade beta Rio Doce/Projeto_Rio_Doce/doce.comunidades.txt') #matriz pres-aus

##matriz distancia funcional
matDisFunc = quality_funct_space(mat_funct=matFunc, nbdim=3,metric='Gower', dendro=FALSE)$details_funct_space$mat_coord 

##metricas da beta div
indices = match(colnames(as.matrix(matPres)),rownames(matDisFunc)) #indices (para ajustar as duas matrizes)
statsBeta = multidimFbetaD(coord=matDisFunc[indices,], occ=as.matrix(matPres), nm_asb_plot=row.names(matPres)) #calculo das metricas de beta-div

##salvando as metricas calculadas
##outputBetaDiv[[i]] = statsBeta
save(statsBeta, file='cenario_obseravdo') #salvando objetos em um .Rdata
gc()


###############################################
###############################################


################################################
######### beta diversidade taxonomica ##########


diretorio = '/home/anderson/Projetos/Isaac - diversidade beta Rio Doce' #especificar diretorio onde salvar

taxdata = read.csv("/home/anderson/Downloads/taxonomic.results.csv", header=TRUE)

taxdata[,1] = NULL
names(taxdata)


jpeg(file.path(diretorio,'boxplotTaxonomicNestedness.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
boxplot(taxdata[,c("ssne0","ssne35","ssne70","ssne104","ssne139")],
        names=c('0%','25%','50%','75%','100%'),
        ylab="Nestedness",
        xlim=c(0,21),
        at=c(1,5,10,15,20),
        boxwex=3,
        lwd=2)
text(x=0.7, y=0.166, labels='A', cex=5)
dev.off()


jpeg(file.path(diretorio,'boxplotTaxonomicTurnover.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
boxplot(taxdata[,c("ssim0","ssim35","ssim70","ssim104","ssim139")],
        names=c('0%','25%','50%','75%','100%'),
        ylab="Turnover",
        xlim=c(0,21),
        at=c(1,5,10,15,20),
        boxwex=3,
        lwd=2)
text(x=1, y=0.679, labels='B', cex=5)
dev.off()

jpeg(file.path(diretorio,'boxplotTaxonomicBeta.jpg'), width=600, height=600)
par(family='times', mar=c(5,5.3,3,2), cex.axis=2.8, cex.lab=3.3)
boxplot(taxdata[,c("ssor0","ssor35","ssor70","ssor104","ssor139")],
        names=c('0%','25%','50%','75%','100%'),
        ylab="Beta diversity",
        xlim=c(0,21),
        at=c(1,5,10,15,20),
        boxwex=3,
        lwd=2)
text(x=0.7, y=0.818, labels='C', cex=5)
dev.off()


###############################################
###############################################
