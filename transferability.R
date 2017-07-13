###ATENCAO: ESSE ARQUIVO CONTEM 'PEDACOS' DE CODIGO, TODOS VOLTADOS PARA MINHA INTENCAO DE FUTURAMENTE EXPLORAR A ###
###GENERALIDADE DOS SDMs(USANDO DIFERENTES ESTRATEGIAS DE TRANSFERABILDADE DOS MODELOS)							  ###

###PEDACO DE CODIGO 01

library(raster)

projFolder = '/home/anderson/Documentos/Projetos/transferability'

environ1 = raster(ncol=100,nrow=100)
var1vector = rep(x=dnorm(x=1:ncol(environ1),mean=50,sd=20), times=nrow(environ1))
var2vector = rep(x=dnorm(x=1:nrow(environ1),mean=50,sd=20), times=ncol(environ1))
#
var1 = matrix(var1vector,nrow=nrow(environ1),ncol=ncol(environ1),byrow=FALSE) * 500
var2 = matrix(var2vector,nrow=nrow(environ1),ncol=ncol(environ1),byrow=TRUE) * 500
#
environ1[] = (var1 * var2) 
#
plot(environ1)

environ2 = raster(ncol=100,nrow=100)
var1vector = rep(x=dnorm(x=1:ncol(environ2),mean=50,sd=20), times=nrow(environ2)) * 500
var2vector = rep(x=dnorm(x=1:nrow(environ2),mean=75,sd=20), times=ncol(environ2)) * 500
#
var1 = matrix(var1vector,nrow=nrow(environ2),ncol=ncol(environ2),byrow=FALSE)
var2 = matrix(var2vector,nrow=nrow(environ2),ncol=ncol(environ2),byrow=TRUE)
#
environ2[] = (var1 * var2) * 1.5
#
plot(environ2)

#varificando climas nao analogos entre os dois ambintes
pontos = randomPoints(mask=environ1,n=100,prob=TRUE)
#pontos = data.frame(x=sample(-180:180,100),y=sample(-90:90,100))
v = extract(environ1, pontos)
messTest =  mess(environ2,v)

par(mfrow=c(1,2))
plot(environ1)
points(pontos)
plot(messTest)

##superficies para probabilidade de amostragem (para simulacao do vies amostral)

for (i in seq(0,1,0.2)){
    mu = i*50 + 50
    ##
    probSamp = raster(ncol=100,nrow=100)
    var1vector = rep(x=dnorm(x=1:ncol(probSamp),mean=mu,sd=20), times=nrow(probSamp)) * 500
    var2vector = rep(x=dnorm(x=1:nrow(probSamp),mean=mu,sd=20), times=ncol(probSamp)) * 500
    ##
    var1 = matrix(var1vector,nrow=nrow(probSamp),ncol=ncol(probSamp),byrow=FALSE)
    var2 = matrix(var2vector,nrow=nrow(probSamp),ncol=ncol(probSamp),byrow=TRUE)
    ##
    probSamp[] = (var1 * var2) * 1.5
    ##
    writeRaster(x=probSamp,filename=paste(projFolder,'/probSamp',mu,'.asc',sep=''))
    
    plot(probSamp)
}


###
hii = raster(x='/home/anderson/Downloads/hii-s-america-geo-grid/hii_s_america_grid/hii_s_amer')
dado = raster(x='/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_01.asc')

##crop hii
areaDado = extent(dado)
hiiSA = crop(hii,areaDado)

##hii ajustado
raster1 = res(dado)
raster2 = res(raster2)
fator = round(raster1[1]/raster2[1])  ##raster maior / raster menor
hiiSA2 = aggregate(hiiSA,fator)
##verificando
res(hiiSA2)
res(dado)


###PEDACO DE CODIGO 02
##Transferibilidade pelo framework de xxx(2016?)###

library(dismo)
library(phyloclim)
library(virtualspecies)

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
source("/home/anderson/R/R-Scripts/TSSmaxent.R")


##ocorrencias
sampleData = data.frame()

for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra
            for (j in 1:NumRep){ #replicas do cenario amostral
                sampleData_i = sampleOccurrences(x=raster(nicheRealPath[sAge+1]),n=sSize,plot=FALSE)$sample.points[,1:2] #amostra d ponto
                scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
                scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
                names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
                write.csv(sampleData,paste(projectFolder,'Amostras/transferability/',spsTypes[i],'/',sAge,'k',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
                sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
            }
        }
    }
}

##TESTANDO MODELOS: calibracao 'monotemporal'

spOccFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/transferability/'
##abrindo um data.frame para armazenar os resultados de AUC
resultsEvaluationMX<-data.frame()
index=0 #auxiliara na criacao do data.frame durante o loop
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
evaluation = list()
TSSvector = data.frame()
sampleSizes = c(5,15,25,35,45,55,65,75,85,95)
NumRep = 3 #numero de replicas (de cada cenario amostral)
sampledAges = c(0,21)
predictors0k = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
projection(predictors0k)=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
predictors21k = stack(list.files(path=envVarPaths[22],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (22kyr BP)
projection(predictors21k)=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
predictors = list('0k'=predictors0k,'21k'=predictors21k)
aucOutput = data.frame()
    
for (i in 1:length(spsTypes)){
    for (j in sampledAges){
        for (k in sampleSizes){
            for (l in 1:NumRep){
    
                sp.file <- read.csv(paste(spOccFolder,spsTypes[i],'/',j,'k',k,'pts',l,"rep.csv",sep=""),header=TRUE) ### read sp occurrence
                sp.occ <- sp.file[,c('lon','lat')] ## select lat-long columns
    
                ##extraindo dados da variavel climatica nos pontos de ocorrencia
                presencesVars <- sp.file[,c('bioclim_01','bioclim_12')]
    
                ##criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
                pres = rep(1, nrow(presencesVars))
    
                ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
                presencesData = data.frame(cbind(sp.occ,pres,presencesVars))
                presencesData = presencesData[complete.cases(presencesData),]
                presencesData = presencesData[complete.cases(presencesData),]
    
                ##criando pontos de background
                background1 <- randomPoints(mask=predictors0k[[1]], n=5000, p=presencesData[,c("lon","lat")], excludep=TRUE)
                background2 <- round(background1, digits=2)
                background3 <- background2[!duplicated(background2),]
                background4 <- background3[complete.cases(background3),]
                background <- data.frame(background4)
                colnames(background) <- c("lon", "lat")
                
                ##extraindo dados da variavel climatica nos pontos de background
                
                backgroundVars <- extract(predictors[[which(sampledAges==j)]], background, method='bilinear', buffer=NULL, fun=NULL)

                ##criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
                pres = rep(0, nrow(backgroundVars))
    
                ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
                backgroundData = data.frame(cbind(background,pres,backgroundVars[,c('bioclim_01','bioclim_12')]))
                backgroundData = backgroundData[complete.cases(backgroundData),]    
    
                ##planilha de dados final
                dataSet = data.frame(rbind(presencesData,backgroundData))

                ##avaliando o modelo
                evaluation = list()
    
                for (m in 1:10){
                    tryCatch({# bootstrapping with 10 replications
            
                        ##reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
                        rand = round(0.75*runif(nrow(presencesData)))
                        presencesTrain = presencesData[rand==0,]
                        backgroundTrain = backgroundData[rand==0,]
            
                        ##juntando presencas e ausencias da calibracao
                        PresBackTrainRaw <- data.frame(rbind(presencesTrain, backgroundTrain))
                        PresBackTrainRaw = PresBackTrainRaw[!duplicated(PresBackTrainRaw[,c('longitude','latitude')]),] #selecionar colunas de longitude e latitude
                        PresBackTrainRaw <- PresBackTrainRaw[complete.cases(PresBackTrainRaw),]
                        PresBackTrain = PresBackTrainRaw

                        ##CRIANDO E RODANDO O MODELO (esquema do SWD - sample with data)##    
                        MX <- maxent(x=PresBackTrain[,c('bioclim_01','bioclim_12')],p=PresBackTrain$pres,args=c(
                                                                                           'responsecurves=true',
                                                                                           'jackknife=true',
                                                                                           'randomseed=true',
                                                                                           'randomtestpoints=0',
                                                                                           'replicates=1',
                                                                                           'writebackgroundpredictions=true',
                                                                                           'linear=true',
                                                                                           'quadratic=true',
                                                                                           'product=false',
                                                                                           'threshold=false',
                                                                                           'hinge=false',
                                                                                           'maximumiterations=1000',
                                                                                           'convergencethreshold=1.0E-5',
                                                                                           'threads=2'))
                        
                        ##avaliacao INTERNA
                        presencesTest = presencesData[rand==1,]
                        presencesTest <- presencesTest[,c('lon','lat')]
                        backgroundTest = backgroundData[rand==1,]
                        backgroundTest = backgroundTest[,c('lon','lat')]
                        
                        ##rodando a avaliacao do modelo
                        internalEval = evaluate(p=presencesTest,a=backgroundTest,m=MX,x=predictors[[which(sampledAges==j)]],type='response')

                        ##avaliacao EXTERNA
                        presencesDataExt = read.csv(paste(spOccFolder,spsTypes[i],'/',sampledAges[which(sampledAges!=j)],'k',k,'pts',l,"rep.csv",sep=""),header=TRUE) ### read sp occurrence
                        presencesDataExt = presencesDataExt[,c('lon','lat')]
                        presencesDataExt = presencesDataExt[rand==1,]

                        ##rodando a avaliacao do modelo
                        externalEval = evaluate(p=presencesDataExt,a=backgroundTest,m=MX,x=predictors[[which(sampledAges!=j)]],type='response')

                        ##salvando AUC interno e externo
                        aucOutput = rbind(aucOutput,data.frame(sp=spsTypes[i], sample_size=k, replicate=l, AUC_int=internalEval@auc, AUC_ext=externalEval@auc,int_age=j,ext_age=sampledAges[which(sampledAges!=j)]))
                        
                        
                    }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
                    )}
            }
        }
    }
}

write.csv(aucOutput,paste(spOccFolder,'/aucOutputsMonotemp.csv',sep=''),row.names = FALSE)


##TESTANDO MODELOS: calibracao 'multitemporal'


for (i in 1:length(spsTypes)){
    for (j in sampledAges){
        for (k in sampleSizes){
            for (l in 1:NumRep){
    
                sp.file.intern <- read.csv(paste(spOccFolder,spsTypes[i],'/',j,'k',k,'pts',l,"rep.csv",sep=""),header=TRUE) ### read sp occurrence
                sp.occ.intern <- sp.file.intern[,c('lon','lat')] ## select lat-long columns

                sp.file.extern <- read.csv(paste(spOccFolder,spsTypes[i],'/',sampledAges[which(sampledAges!=j)],'k',k,'pts',l,"rep.csv",sep=""),header=TRUE) ### read sp occurrence
                sp.occ.extern <- sp.file.extern[,c('lon','lat')] ## select lat-long columns
                
                ##extraindo dados da variavel climatica nos pontos de ocorrencia
                presencesVarsIntern <- sp.file.intern[,c('bioclim_01','bioclim_12')]
                presencesVarsExtern <- sp.file.extern[,c('bioclim_01','bioclim_12')]
    
                ##criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
                pres = rep(1, nrow(presencesVarsIntern))
    
                ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
                presencesDataIntern = data.frame(cbind(sp.occ.intern,pres,presencesVarsIntern))
                presencesDataIntern = presencesDataIntern[complete.cases(presencesDataIntern),]
                presencesDataIntern = presencesDataIntern[complete.cases(presencesDataIntern),]
                ##
                presencesDataExtern = data.frame(cbind(sp.occ.intern,pres,presencesVarsExtern))
                presencesDataExtern = presencesDataExtern[complete.cases(presencesDataExtern),]
                presencesDataExtern = presencesDataExtern[complete.cases(presencesDataExtern),]
    
                ##criando pontos de background
                backgroundIntern1 <- randomPoints(mask=predictors[[which(sampledAges==j)]][[1]], n=5000, p=presencesDataIntern[,c("lon","lat")], excludep=TRUE)
                backgroundIntern2 <- round(background1, digits=2)
                backgroundIntern3 <- backgroundIntern2[!duplicated(backgroundIntern2),]
                backgroundIntern4 <- backgroundIntern3[complete.cases(backgroundIntern3),]
                backgroundIntern <- data.frame(backgroundIntern4)
                colnames(backgroundIntern) <- c("lon", "lat")
                ##
                backgroundExtern1 <- randomPoints(mask=predictors[[which(sampledAges!=j)]][[1]], n=5000, p=presencesDataExtern[,c("lon","lat")], excludep=TRUE)
                backgroundExtern2 <- round(background1, digits=2)
                backgroundExtern3 <- backgroundExtern2[!duplicated(backgroundExtern2),]
                backgroundExtern4 <- backgroundExtern3[complete.cases(backgroundExtern3),]
                backgroundExtern <- data.frame(backgroundExtern4)
                colnames(backgroundExtern) <- c("lon", "lat")
                
                ##extraindo dados da variavel climatica nos pontos de background
                
                backgroundVarsIntern <- extract(predictors[[which(sampledAges==j)]], backgroundIntern1, method='bilinear', buffer=NULL, fun=NULL)
                backgroundVarsExtern <- extract(predictors[[which(sampledAges!=j)]], backgroundExtern1, method='bilinear', buffer=NULL, fun=NULL)
                

                ##criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
                pres = rep(0, nrow(backgroundVarsIntern))
    
                ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
                backgroundDataIntern = data.frame(cbind(backgroundIntern,pres,backgroundVarsIntern[,c('bioclim_01','bioclim_12')]))
                backgroundDataIntern = backgroundDataIntern[complete.cases(backgroundDataIntern),]
                ##
                backgroundDataExtern = data.frame(cbind(backgroundExtern,pres,backgroundVarsExtern[,c('bioclim_01','bioclim_12')]))
                backgroundDataExtern = backgroundDataExtern[complete.cases(backgroundDataExtern),]    
    
                ##planilha de dados final
                dataSet = list()
                for (m in 1:10){
                    rand = sample(c(rep(1,round(nrow(presencesDataIntern)*0.5)),rep(0,nrow(presencesDataIntern)-round(nrow(presencesDataIntern)*0.5))))
                    presencesTrain = rbind(presencesDataIntern[rand==0,],presencesDataExtern[rand==1,])
                    backgroundTrain = rbind(backgroundDataIntern[rand==0,],backgroundDataExtern[rand==1,])
                    dataSet_i = data.frame(c(rbind(presencesTrain,backgroundTrain)))
                    dataSet[[m]] = dataSet_i
                }

                ##avaliando o modelo
                evaluation = list()
    
                for (n in 1:10){
                    tryCatch({# bootstrapping with 10 replications
            
                        ##reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
                        ##rand = round(0.75*runif(nrow(dataSet[[n]][dataSet[[n]]$pres==1,])))
                        rand = sample(c(rep(1,round(nrow(dataSet[[n]][dataSet[[n]]$pres==1,])*0.25)),rep(0,nrow(dataSet[[n]][dataSet[[n]]$pres==1,])-round(nrow(dataSet[[n]][dataSet[[n]]$pres==1,])*0.25))))
                        presencesTrain = dataSet[[n]][dataSet[[n]]$pres==1,][rand==0,]
                        backgroundTrain = dataSet[[n]][dataSet[[n]]$pres==0,][rand==0,]
            
                        ##juntando presencas e ausencias da calibracao
                        PresBackTrainRaw <- data.frame(rbind(presencesTrain, backgroundTrain))
                        PresBackTrainRaw = PresBackTrainRaw[!duplicated(PresBackTrainRaw[,c('lon','lat')]),] #selecionar colunas de longitude e latitude
                        PresBackTrainRaw <- PresBackTrainRaw[complete.cases(PresBackTrainRaw),]
                        PresBackTrain = PresBackTrainRaw

                        ##CRIANDO E RODANDO O MODELO (esquema do SWD - sample with data)##    
                        MX <- maxent(x=PresBackTrain[,c('bioclim_01','bioclim_12')],p=PresBackTrain$pres,args=c(
                                                                                           'responsecurves=true',
                                                                                           'jackknife=true',
                                                                                           'randomseed=true',
                                                                                           'randomtestpoints=0',
                                                                                           'replicates=1',
                                                                                           'writebackgroundpredictions=true',
                                                                                           'linear=true',
                                                                                           'quadratic=true',
                                                                                           'product=false',
                                                                                           'threshold=false',
                                                                                           'hinge=false',
                                                                                           'maximumiterations=1000',
                                                                                           'convergencethreshold=1.0E-5',
                                                                                           'threads=2'))
                        
                        ##avaliacao INTERNA
                        presencesTestIntern = presencesDataIntern[rand==1,]
                        presencesTestIntern <- presencesTestIntern[,c('lon','lat')]
                        backgroundTestIntern = backgroundDataIntern[rand==1,]
                        backgroundTestIntern = backgroundTestIntern[,c('lon','lat')]
                        
                        ##rodando a avaliacao do modelo
                        internalEval = evaluate(p=presencesTest,a=backgroundTest,m=MX,x=predictors[[which(sampledAges==j)]],type='response')

                        ##avaliacao EXTERNA
                        presencesTestExtern = presencesDataExtern[rand==1,]
                        presencesTestExtern <- presencesTestExtern[,c('lon','lat')]
                        backgroundTestExtern = backgroundDataExtern[rand==1,]
                        backgroundTestExtern = backgroundTestExtern[,c('lon','lat')]

                        ##rodando a avaliacao do modelo
                        externalEval = evaluate(p=presencesTestExtern,a=backgroundTestExtern,m=MX,x=predictors[[which(sampledAges!=j)]],type='response')

                        ##salvando AUC interno e externo
                        aucOutputMultitemp = rbind(aucOutputMultitemp,data.frame(sp=spsTypes[i], sample_size=k, replicate=l, AUC_int=internalEval@auc, AUC_ext=externalEval@auc,int_age=j,ext_age=sampledAges[which(sampledAges!=j)]))
                        
                        
                    }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
                    )}
            }
        }
    }
}

write.csv(aucOutput,paste(spOccFolder,'/aucOutputsMultitemp.csv',sep=''),row.names = FALSE)



