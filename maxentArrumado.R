#######################################################
####################### MAXENT ########################
#######################################################

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
evaluation = list()
TSSvector = data.frame()
sampleSizes = c(5,15,25,35,45,55,65,75,85,95)

for (i in 1:length(spsTypes)){
    for (j in sampleSizes){
        for (k in 1:5){ #loop sobre o numero de replicas #10
        
            occPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/occ',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
            backgroundPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/bg',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de background
            names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
            dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
            
            me = maxent(
                x=dataSet[,c("bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")],
                p=dataSet$pres,
                path=paste(maxentFolder,spsTypes[i],sep=''),
                args=c('responsecurves=TRUE',
                       'jackknife=TRUE',
                       'randomseed=TRUE',
                       'randomtestpoints=25',
                       'maximumbackground=5000',
                       'replicates=1',
                       'replicatetype=subsample',
                       'writebackgroundpredictions=TRUE',
                       'linear=TRUE',
                       'quadratic=TRUE',
                       'product=FALSE',
                       'threshold=FALSE',
                       'hinge=FALSE',
                       'maximumiterations=1000',
                       'convergencethreshold=1.0E-5',
                       'threads=2'
                       ))
            
            ##rodando a avaliacao do modelo
            TSSvector = rbind(TSSvector, TSSmaxent(paste(maxentFolder,spsTypes[i],sep='')))
            evaluation = append(evaluation, evaluate(p=occPoints,a=backgroundPoints,model=me))
        
            for (l in 1:length(envVarPaths[1:24])){
                predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
                predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
                crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
                proj = predict(me,predictors,crs=crs) #realizando projetacoes (para cada replica)
                writeRaster(mean(proj),paste(maxentFolder,spsTypes[i],'/projections/projection-Time',l-1,'kyrBP','-Replica',k,'-Sample',j,'.asc',sep=''),overwrite=TRUE) #salvando a projecao media
            }
        }

        ##criando um mapa binario
        thresholdValues = NULL
        aucValues = NULL
        for (m in 1:length(evaluation)){
            thresholdValues <- append(thresholdValues, threshold(evaluation[[m]],'spec_sens'))
            aucValues = append(aucValues, evaluation[[m]]@auc)
        }
        aucMean = mean(aucValues)
        thresholdMean = mean(thresholdValues)
        TSSmean = mean(TSSvector$TSS)
        statResults = data.frame(sp=spsTypes[i],AUCmean=aucMean,TSSmean=TSSmean,ThresholdMean=thresholdMean)
        write.csv(statResults,file=paste(projectFolder,'Maxent/',spsTypes[i],'/StatisticsResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
        
        ##esvaziando o vetor para a proxima especie
        evaluation = list()
        TSSvector = data.frame()
    }
}
#######################################################
#######################################################
#######################################################
