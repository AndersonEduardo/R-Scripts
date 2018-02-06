## invasao mamiferos america do sul

## pacotes necessarios
library(biomod2)
library(raster)

## parametros e variaveis globais
projectFolder = '/home/anderson/Área de Trabalho/bioinvasao_AmSul'
sampleFolder = '/home/anderson/Área de Trabalho/bioinvasao_AmSul/occ_data'
envVarFolder = '/home/anderson/Área de Trabalho/bioinvasao_AmSul/env_data'
maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java/maxent.jar' #pasta para resultados do maxent'

## variaveis ambientais
globalVars = stack(list.files('/home/anderson/PosDoc/dados_ambientais/bio_2-5m_bil', pattern='.bil', full.names=TRUE))

plot(globalVars[[1]])
##
extEuro = drawExtent()
extEuro = c(-9.121301, 145.8906, 6.600978, 75.33003)
##
extAm = drawExtent()
extAm = c(-175.3659, -20.35405, 14.55908, 82.56467)
##
extAmSul = drawExtent()
extAmSul = c(-88.49937,-28.59139,-57.7873,13.83562)

euroVars = crop(globalVars, extEuro)
amerVars = crop(globalVars, extAm)
amerSulVars = crop(globalVars, extAmSul)

## testando correcaloes
## eurasia
testEurasia = getValues(euroVars)
cor.matrix <- as.data.frame(cor(testEurasia, use="complete.obs"))
write.csv(cor.matrix,paste(envVarFolder,'/cor_Eurasia.csv',sep=''), row.names=TRUE)
## america do norte
testAmNorte = getValues(amerVars)
cor.matrix <- as.data.frame(cor(testAmNorte, use="complete.obs"))
write.csv(cor.matrix,paste(envVarFolder,'/cor_AmNorte.csv',sep=''), row.names=TRUE)

## identificando variaveis nao correlacionadas
## eurasia
corMatrix = read.csv(paste(envVarFolder,'/cor_Eurasia.csv',sep=''),header=TRUE, row.names=1)
corMatrixTri = corMatrix
corMatrixTri[upper.tri(corMatrixTri)] = 0.0
diag(corMatrixTri) = 0.0
corMatrixClean = corMatrix[,!apply(corMatrixTri,2,function(x) any(abs(x) > 0.7))]
sort(names(corMatrixClean))
## america norte
corMatrix = read.csv(paste(envVarFolder,'/cor_AmNorte.csv',sep=''),header=TRUE, row.names=1)
corMatrixTri = corMatrix
corMatrixTri[upper.tri(corMatrixTri)] = 0.0
diag(corMatrixTri) = 0.0
corMatrixClean = corMatrix[,!apply(corMatrixTri,2,function(x) any(abs(x) > 0.7))]
sort(names(corMatrixClean))

## construcao das variaveis preditoras

euroVars = euroVars[[c("bio7", "bio8", "bio9", "bio18", "bio19" )]]
amerVars = amerVars[[c("bio7", "bio8", "bio9", "bio18", "bio19" )]]
amerSulVars = amerSulVars[[c("bio7", "bio8", "bio9", "bio18", "bio19" )]]

writeRaster(x=euroVars, filename=paste(envVarFolder,euroVars,'.grd',sep=''))

#######################################################################
############# densidade da variavel temporetatura #####################
################ na area nativa e area invadida #######################
#######################################################################


##vetor com o nome dos arquivos
spsNames = list.files(paste(sampleFolder,sep=''), full.names=FALSE )
vecNames = gsub(pattern='_invadido.csv',replacement='',x=spsNames)
vecNames = gsub(pattern='_nativo.csv',replacement='',x=vecNames)
vecNames = unique(vecNames)
vecNames = vecNames[which(vecNames != "Mus_musculus")] #retirando Mus musculus


##temperaturas criticas 
termNeut = data.frame(sps=c("Sus_scrofa",
                            "Rattus_rattus",
                            "Castor_canadensis",
                            "Mus_musculus"),
                      min=c(40, 260, 0, NA),
                      max=c(220, 340, 280, NA))

## sus acrofa = 4 e 22
## rattus rattus = 26 e 34
## castor = 0 e 28

## comparacao grafica entre temperatura media nos dados para area nativa, area invadida e TNZ
for(sp_i in vecNames){
    ##vetor de nomes de arquivo
    spsVector = list.files(paste(sampleFolder,sep=''), pattern=sp_i, full.names=TRUE)
    ##area nativa
    nativoData = read.csv(spsVector[grep(pattern='nativo',x=spsVector)], header=TRUE)
    names(nativoData) = c('sp','lon','lat')
    ##area invadida
    invadidoData = read.csv(spsVector[grep(pattern='invadido',x=spsVector)], header=TRUE)
    names(invadidoData) = c('sp','lon','lat')
    ##extraindo variaveis ambientais
    if(sp_i == "Castor_canadensis"){
        nativVars = amerVars
    }else{
        nativVars = euroVars}
    ##
    tempNativa = extract(x=nativVars[[c('bio1','bio5','bio6')]],y=nativoData[,c('lon','lat')], na.rm='TRUE') #area nativa
    tempInvadido = extract(x=amerSulVars[[c('bio1','bio5','bio6')]],y=invadidoData[,c('lon','lat')], na.rm=TRUE) #area invadida
    tempNativa = tempNativa[complete.cases(tempNativa),] #limpando
    tempInvadido = tempInvadido[complete.cases(tempInvadido),] #limpando
    ##graficos
    jpeg(paste(projectFolder,'/',sp_i,'.jpeg', sep=''))
    plot(density(tempNativa[,'bio1']),
         ylim=c(0, max(max(density(tempNativa[,'bio1'])$y),max(density(tempInvadido[,'bio1'])$y))),
         xlim=c(min(c(range(density(tempNativa[,'bio1'])$x), range(density(tempInvadido[,'bio1'])$x), as.numeric(termNeut[termNeut$sps==sp_i,][,c('min','max')]))),
                max(c(range(density(tempNativa[,'bio1'])$x), range(density(tempInvadido[,'bio1'])$x), as.numeric(termNeut[termNeut$sps==sp_i,][,c('min','max')])))),
         lwd=2,
         col='blue',
         xlab='Mean annual temperature',
         main=gsub('_',' ',sp_i))
    lines(density(tempInvadido[,'bio1']),
          lwd=2,
          col='red')
    abline(v=termNeut[termNeut$sps==sp_i,][,c('min','max')], lty=2, lwd=1.5)
    legend(x='topleft', legend=c('Native area', 'Invaded area', 'Thermal neutral zone'), lty=c(1,1,2), lwd=2, col=c('blue','red','black'), cex=0.8, bg='white')
    dev.off()
}


######### modelagem de nicho #########
############# BIOMOD2 ################

##vetor com o nome dos arquivos
spsNames = list.files(paste(sampleFolder,sep=''), full.names=FALSE )
vecNames = gsub(pattern='_invadido.csv',replacement='',x=spsNames)
vecNames = gsub(pattern='_nativo.csv',replacement='',x=vecNames)
vecNames = unique(vecNames)
vecNames = vecNames[which(vecNames != "Mus_musculus")] #retirando Mus musculus

for(sp_i in vecNames){
    
    ##vetor de nomes de arquivo
    spsFiles = list.files(paste(sampleFolder,sep=''), pattern=sp_i, full.names=TRUE)
    ##area nativa
    nativoData = read.csv(spsVector[grep(pattern='nativo',x=spsVector)], header=TRUE)
    names(nativoData) = c('sp','lon','lat')
    ##area invadida

    ##extraindo variaveis ambientais
    if(sp_i == "Castor_canadensis"){
        nativVars = amerVars
    }else{
        nativVars = euroVars}
    ##    

    ##variaveis e parametros locais especificos para o biomod2

    myRespName <- sp_i
    myResp <- nativoData[,c('lon','lat')] # variavel resposta (para biomod2)
    coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
    crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints
    myExpl = stack(euroVars)

    ##ajuste de dados de entrada para biomod2
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         PA.nb.rep = 1)

          ##parametrizando os modelos
      myBiomodOption <- BIOMOD_ModelingOptions(
        MAXENT.Phillips = list(path_to_maxent.jar= maxentFolder,
                               linear=TRUE,
                               quadratic=TRUE,
                               product=FALSE,
                               threshold=FALSE,
                               hinge=FALSE,
                               maximumiterations=1000,
                               convergencethreshold=1.0E-5,
                               threads=2))
    
    ##rodando o(s) algoritmo(s) (i.e. SDMs)
    myBiomodModelOut <- BIOMOD_Modeling(
        data = myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 100,
        DataSplit = 75,
        models.eval.meth = c('TSS','ROC'),
        do.full.models = FALSE,
        modeling.id = paste(myRespName,'_area_nativa',sep=''))

    ##My output data
    evaluationScores = get_evaluations(myBiomodModelOut)

    ##rodando algortmo de projecao (i.e. rodando a projecao)
    myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExpl_invadido,
        proj.name = paste(myRespName,'_area_invadida',sep=''),
        compress = 'TRUE',
        build.clamping.mask = 'TRUE')

}


