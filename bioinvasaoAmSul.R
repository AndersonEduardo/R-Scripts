## invasao mamiferos america do sul
## Anderson A. Eduardo
## 07/fev/2018

## pacotes necessarios
library(biomod2)
library(raster)

## parametros e variaveis globais
projectFolder = '/home/anderson/Área de Trabalho/bioinvasao_AmSul'
sampleFolder = '/home/anderson/Área de Trabalho/bioinvasao_AmSul/occ_data/todasSp'
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

##para SDMs
euroVars = crop(globalVars, extEuro)
amerVars = crop(globalVars, extAm)
amerSulVars = crop(globalVars, extAmSul)

##temperatura media anual (para zona TNZ)
tempMedAmNorte = crop(globalVars[['bio1']], extAm)
tempMedEurasia = crop(globalVars[['bio1']], extEuro)
tempMedAmSul = crop(globalVars[['bio1']], extAmSul)

## salvando variaveis temperatura
writeRaster(x=tempMedAmNorte, filename=paste(envVarFolder,'/tempMedAmNorte.grd',sep=''), overwrite=TRUE)
writeRaster(x=tempMedEurasia, filename=paste(envVarFolder,'/tempMedEurasia.grd',sep=''), overwrite=TRUE)
writeRaster(x=tempMedAmSul, filename=paste(envVarFolder,'/tempMedAmSul.grd',sep=''), overwrite=TRUE)


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
euroVars = euroVars[[c("bio1", "bio7", "bio15", "bio16", "bio17" )]]
amerVars = amerVars[[c("bio1", "bio7", "bio15", "bio16", "bio17" )]]
amerSulVars = amerSulVars[[c("bio1", "bio7", "bio15", "bio16", "bio17" )]]

## salvando variaveis preditoras
writeRaster(x=euroVars, filename=paste(envVarFolder,'/euroVars.grd',sep=''), format='raster', bylayer=TRUE, suffix=names(euroVars), overwrite=TRUE)
writeRaster(x=amerVars, filename=paste(envVarFolder,'/amerVars.grd',sep=''), format='raster', bylayer=TRUE, suffix=names(amerVars), overwrite=TRUE)
writeRaster(x=amerSulVars, filename=paste(envVarFolder,'/amerSulVars.grd',sep=''), format='raster', bylayer=TRUE, suffix=names(amerSulVars), overwrite=TRUE)

#######################################################################
############# densidade da variavel temporetatura #####################
################ na area nativa e area invadida #######################
#######################################################################

##abrindo variaveis ambientais
euroVars = stack(grep(pattern='*.euroVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))
amerVars = stack(grep(pattern='*.amerVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))
amerSulVars = stack(grep(pattern='*.amerSulVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))

##abrindo variaveis temperatura
tempMedAmNorte = raster(paste(envVarFolder,'/tempMedAmNorte.grd',sep=''))
tempMedEurasia = raster(paste(envVarFolder,'/tempMedEurasia.grd',sep=''))
tempMedAmSul = raster(paste(envVarFolder,'/tempMedAmSul.grd',sep=''))

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
        nativVars = tempMedAmNorte
    }else{
        nativVars = tempMedEurasia}
    ##
    tempNativa = extract(x=nativVars,y=nativoData[,c('lon','lat')], na.rm=TRUE) #area nativa
    tempInvadido = extract(x=amerSulVars,y=invadidoData[,c('lon','lat')], na.rm=TRUE) #area invadida
    tempNativa = tempNativa[complete.cases(tempNativa)] #limpando
    tempInvadido = tempInvadido[complete.cases(tempInvadido)] #limpando
    ##graficos
    jpeg(paste(projectFolder,'/',sp_i,'.jpeg', sep=''))
    plot(density(tempNativa),
         ylim=c(0, max(max(density(tempNativa)$y),max(density(tempInvadido)$y))),
         xlim=c(min(c(range(density(tempNativa)$x), range(density(tempInvadido)$x), as.numeric(termNeut[termNeut$sps==sp_i,][,c('min','max')]))),
                max(c(range(density(tempNativa)$x), range(density(tempInvadido)$x), as.numeric(termNeut[termNeut$sps==sp_i,][,c('min','max')])))),
         lwd=2,
         col='blue',
         xlab='Mean annual temperature',
         main=gsub('_',' ',sp_i))
    lines(density(tempInvadido),
          lwd=2,
          col='red')
    #abline(v=termNeut[termNeut$sps==sp_i,][,c('min','max')], lty=2, lwd=1.5)
    legend(x='topright', legend=c('Native area', 'Invaded area', 'Thermal neutral zone'), lty=c(1,1,2), lwd=2, col=c('blue','red','black'), cex=0.8, bg='white')
    dev.off()
}


######################################
######### modelagem de nicho #########
############# BIOMOD2 ################
######################################


##alterando ambiente de trabalho (importante porque o biomod2 vai construir pastas)
if (!file.exists(paste(projectFolder,'/modelos',sep=''))){
    dir.create(paste(projectFolder,'/modelos',sep=''), recursive=TRUE)
}

setwd(paste(projectFolder,'/modelos',sep=''))

##parametros e variaveis globais
statResults = data.frame()
AUCraw = data.frame()
TSSraw = data.frame()

##abrindo variaveis ambientais
euroVars = stack(grep(pattern='*.euroVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))
amerVars = stack(grep(pattern='*.amerVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))
amerSulVars = stack(grep(pattern='*.amerSulVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))

##abrindo variaveis temperatura
tempMedAmNorte = raster(paste(envVarFolder,'/tempMedAmNorte.grd',sep=''))
tempMedEurasia = raster(paste(envVarFolder,'/tempMedEurasia.grd',sep=''))
tempMedAmSul = raster(paste(envVarFolder,'/tempMedAmSul.grd',sep=''))

##vetor com o nome dos arquivos
spsNames = list.files(paste(sampleFolder,sep=''), full.names=FALSE )
vecNames = gsub(pattern='_invadido.csv',replacement='',x=spsNames)
vecNames = gsub(pattern='_nativo.csv',replacement='',x=vecNames)
vecNames = unique(vecNames)
vecNames = vecNames[which(vecNames != "Mus_musculus")] #retirando Mus musculus

## loop sobre a lista de especies

for(sp_i in vecNames){
    
    ##vetor de nomes de arquivo
    spsFiles = list.files(paste(sampleFolder,sep=''), pattern=sp_i, full.names=TRUE)
    ##area nativa
    nativoData = read.csv(spsFiles[grep(pattern='nativo',x=spsFiles)], header=TRUE)
    names(nativoData) = c('sp','lon','lat')
    invadidoData = read.csv(spsFiles[grep(pattern='invadido',x=spsFiles)], header=TRUE)
    names(invadidoData) = c('sp','lon','lat')
    ##area invadida

    ##extraindo variaveis ambientais
    if(sp_i == "Castor_canadensis"){
        nativVars = amerVars
    }else{
        nativVars = euroVars}

    ##variaveis e parametros locais especificos para o biomod2

    myRespName <- sp_i
    myResp <- nativoData[,c('lon','lat')] # variavel resposta (para biomod2)
    coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
    crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints
    myExpl = stack(nativVars)
    myExpl_invadido = stack(amerSulVars)

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
                               maximumiterations=700,
                               convergencethreshold=1.0E-5,
                               threads=2))
    
    ##rodando o(s) algoritmo(s) (i.e. SDMs)
    myBiomodModelOut <- BIOMOD_Modeling(
        data = myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 10,
        DataSplit = 75,
        models.eval.meth = c('TSS','ROC'),
        do.full.models = FALSE,
        modeling.id = paste(myRespName,'_area_nativa',sep=''))

    ##My output data
    evaluationScores = get_evaluations(myBiomodModelOut)

    ##tabela ccompleta do AUC
    AUCraw = rbind(AUCraw, data.frame(sp=myRespName,
                                      rbind(evaluationScores['ROC','Testing.data',,,])
                                      )
                   )
    write.csv(AUCraw, paste(projectFolder,'/AUCraw.csv',sep=''), row.names=FALSE)

    ##tabela ccompleta do TSS
    TSSraw = rbind(TSSraw, data.frame(sp=myRespName,
                                      rbind(evaluationScores['ROC','Testing.data',,,])
                                      )
                   )
    write.csv(TSSraw, paste(projectFolder,'/TSSraw.csv',sep=''), row.names=FALSE)

    ##tabela ccompleta de resultados dos modelos
    statResults = rbind(statResults,
                        data.frame(sp = myRespName,
                                   meanAUC = mean(rbind(evaluationScores['ROC','Testing.data',,,])),
                                   maxAUC = max(rbind(evaluationScores['ROC','Testing.data',,,])),
                                   minAUC = min(rbind(evaluationScores['ROC','Testing.data',,,])),
                                   meanTSS = mean(rbind(evaluationScores['TSS','Testing.data',,,])),
                                   maxTSS = max(rbind(evaluationScores['TSS','Testing.data',,,])),
                                   minTSS = min(rbind(evaluationScores['TSS','Testing.data',,,]))
                                   )
                        )
    write.csv(statResults, paste(projectFolder,'/statResults.csv',sep=''), row.names=FALSE)
    
    ##rodando modelo sem repartir dados, para projecao
    rm(myBiomodModelOut)
    myBiomodModelOut <- BIOMOD_Modeling(
        data = myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 1,
        DataSplit = 100,
        models.eval.meth = c('TSS','ROC'),
        do.full.models = TRUE,
        SaveObj = FALSE,
        modeling.id = paste(myRespName,'_fullModelForProjection',sep=''))
    
    ##rodando algortmo de projecao (i.e. rodando a projecao)
    myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExpl_invadido,
        proj.name = paste(myRespName,'_area_invadida',sep=''),
        compress = 'TRUE',
        build.clamping.mask = 'TRUE')

    ## MAPAS

    suitabMap = get_predictions(myBiomodProj)
    
    ##extraindo variaveis ambientais
    if(sp_i == "Castor_canadensis"){
        nativTemp = tempMedAmNorte
    }else{
        nativTemp = tempMedEurasia}
    
    ##limpando pts Am. Sul
    ptsFilter = data.frame(invadidoData, value=extract(tempMedAmSul, invadidoData[,c('lon','lat')]))
    invadidoData = ptsFilter[complete.cases(ptsFilter),]

    ##limpando pts area nativa
    ptsFilter = data.frame(nativoData, value=extract(nativTemp, nativoData[,c('lon','lat')]))
    nativoData = ptsFilter[complete.cases(ptsFilter),]
    
    jpeg(paste(projectFolder,'/modelos/',sp_i,'_invadido.jpeg',sep=''))
    plot(suitabMap/1000, main=gsub('_', ' ', sp_i))
    plot(tempMedAmSul>termNeut[which(termNeut$sps==sp_i),'min'] & tempMedAmSul<termNeut[which(termNeut$sps==sp_i),'max'],
         col=c(rgb(1,1,1,0),rgb(1,0,0,0.5)),
         add=TRUE,
         legend=FALSE)
    points(invadidoData[,c('lon','lat')], pch=19, col=rgb(0,0,0,0.5))
    grid()
    dev.off()
    
    jpeg(paste(projectFolder,'/modelos/',sp_i,'_nativo.jpeg',sep=''))
    plot(nativTemp/100, main='Mean temperature, occurrences and TNZ')
    plot(nativTemp > termNeut[which(termNeut$sps==sp_i),'min'] &
         nativTemp < termNeut[which(termNeut$sps==sp_i),'max'],
         col=c(rgb(1,1,1,0),rgb(1,0,0,0.5)),
         add=TRUE,
         legend=FALSE)
    points(nativoData[,c('lon','lat')], pch=19, col=rgb(0,0,0,0.5))
    grid()
    dev.off()
    
}

#######################################################
### teste de significancia da comparacao dos nichos ###
#######################################################

library(ecospat)

##variaveis globais
outputData = data.frame()

##abrindo variaveis ambientais
euroVars = stack(grep(pattern='*.euroVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))
amerVars = stack(grep(pattern='*.amerVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))
amerSulVars = stack(grep(pattern='*.amerSulVars.*grd', x=list.files(path=envVarFolder, full.names=TRUE), value=TRUE))

##abrindo variaveis temperatura
tempMedAmNorte = raster(paste(envVarFolder,'/tempMedAmNorte.grd',sep=''))
tempMedEurasia = raster(paste(envVarFolder,'/tempMedEurasia.grd',sep=''))
tempMedAmSul = raster(paste(envVarFolder,'/tempMedAmSul.grd',sep=''))

##vetor com o nome dos arquivos
spsNames = list.files(paste(sampleFolder,sep=''), full.names=FALSE )
vecNames = gsub(pattern='_invadido.csv',replacement='',x=spsNames)
vecNames = gsub(pattern='_nativo.csv',replacement='',x=vecNames)
vecNames = unique(vecNames)
vecNames = vecNames[which(vecNames != "Mus_musculus")] #retirando Mus musculus

## loop sobre a lista de especies

for(sp_i in vecNames){
    
    ##vetor de nomes de arquivo
    spsFiles = list.files(paste(sampleFolder,sep=''), pattern=sp_i, full.names=TRUE)

    ##pts area nativa
    nativoData = read.csv(spsFiles[grep(pattern='nativo',x=spsFiles)], header=TRUE)
    names(nativoData) = c('sp','lon','lat')
    nativoData = data.frame(nativoData, occ=1)

    ##pts area invadida
    invadidoData = read.csv(spsFiles[grep(pattern='invadido',x=spsFiles)], header=TRUE)
    names(invadidoData) = c('sp','lon','lat')
    invadidoData = data.frame(invadidoData, occ=1)

    ##variaveis ambientais
    if(sp_i == "Castor_canadensis"){
        nativVars = amerVars
    }else{
        nativVars = euroVars}

    ##pontos de 'ausencia'
    nativAbs = dismo::randomPoints(nativVars, 1000)
    nativAbs = data.frame(sp=sp_i, nativAbs, occ=0); names(nativAbs) = c('sp','lon','lat','occ')
    invadAbs = dismo::randomPoints(amerSulVars, 1000)
    invadAbs = data.frame(sp=sp_i, invadAbs, occ=0); names(invadAbs) = c('sp','lon','lat','occ')
    
    ##agrupando pts de presencas e ausencias
    nativoData = rbind(nativoData, nativAbs)
    invadidoData = rbind(invadidoData, invadAbs)

    ##extraindo variaveis ambientais
    nativoData = data.frame( nativoData, extract(x=nativVars,y=nativoData[,c('lon','lat')]) )
    nativoData = nativoData[complete.cases(nativoData),]

    invadidoData = data.frame( invadidoData, extract(x=amerSulVars,y=invadidoData[,c('lon','lat')]) )
    invadidoData = invadidoData[complete.cases(invadidoData),]
    
    ##The PCA is calibrated on all the sites of the study area
    pca.env <- dudi.pca(rbind(nativoData, invadidoData)[,c('bio1','bio7','bio15','bio16','bio17')],scannf=F,nf=2)
    ##ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico
    
    ##PCA scores for the whole study area
    scores.globclim <- pca.env$li
    ##PCA scores for the species native distribution
    scores.sp.realNiche <- suprow(pca.env,nativoData[which(nativoData[,'occ']==1),c('bio1','bio7','bio15','bio16','bio17')])$li

    ##PCA scores for the species invasive distribution
    scores.sp.SDMniche <- suprow(pca.env,invadidoData[which(invadidoData[,'occ']==1),c('bio1','bio7','bio15','bio16','bio17')])$li

    ##PCA scores for the whole native study area
    scores.clim.realNiche <-suprow(pca.env,nativoData[,c('bio1','bio7','bio15','bio16','bio17')])$li

    ##PCA scores for the whole invaded study area
    scores.clim.SDMniche <- suprow(pca.env,invadidoData[,c('bio1','bio7','bio15','bio16','bio17')])$li

    ##gridding the native niche
    grid.clim.realNiche <-ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.realNiche,sp=scores.sp.realNiche, R=100,th.sp=0)

    ##gridding the invasive niche
    grid.clim.SDMniche <- ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.SDMniche,sp=scores.sp.SDMniche, R=100,th.sp=0)

    ##equivalencia de nicho
    ##OBS: Compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim
    ## and z2.sim, which are built from random reallocations of occurences of z1 and z2.
    ##'alternative' argument specifies if you want to test for niche conservatism (alternative = "greater", i.e.  the
    ## niche overlap is more equivalent/similar than random) or for niche divergence (alternative = "lower",
    ## i.e. the niche overlap is less equivalent/similar than random).
    eq.test <- ecospat.niche.equivalency.test(grid.clim.realNiche, grid.clim.SDMniche,rep=100, alternative = "greater")

    ##similaridade de nicho
    ##OBS: Compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches
    ## (z2.sim) as available in the range of z2 (z2$Z). z2.sim has the same pattern as z2 but the center is
    ## randomly translatated in the availabe z2$Z space and weighted by z2$Z densities. If rand.type = 1,
    ## both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
    ## 'alternative' specifies if you want to test for niche conservatism (alternative = "greater", i.e.  the
    ## niche overlap is more equivalent/similar than random) or for niche divergence (alternative = "lower",
    ## i.e. the niche overlap is less equivalent/similar than random)
    sim.test <- ecospat.niche.similarity.test(grid.clim.realNiche, grid.clim.SDMniche, rep=100, alternative = "greater")
    
    Dobs_equiv = eq.test$obs$D #indice D observado no teste de equivalencia de nicho
    Iobs_equiv = eq.test$obs$I #indice I observado no teste de equivalencia de nicho
    DpValue_equiv = eq.test$p.D #p-valor indice D no teste de equivalencia de nicho
    IpValue_equiv = eq.test$p.I #p-valor indice I no teste de equivalencia de nicho
    ##
    Dobs_simi = sim.test$obs$D #indice D observado no teste de similaridade de nicho
    Iobs_simi = sim.test$obs$I #indice I observado no teste de similaridade de nicho
    DpValue_simi = sim.test$p.D #p-valor indice D no teste de similaridade de nicho
    IpValue_simi = sim.test$p.I #p-valor indice I no teste de similaridade de nicho

    outputData = rbind(outputData,data.frame(sp = sp_i,
                                             Schoeners_D_equiv = Dobs_equiv,
                                             p_value_equiv = DpValue_equiv,
                                             Hellinger_I_equiv = Iobs_equiv,
                                             p_value_equiv = IpValue_equiv,
                                             Schoeners_D_simi = Dobs_simi,
                                             p_value_simi = DpValue_simi,
                                             Hellinger_I_simi = Iobs_simi,
                                             p_value_simi = IpValue_simi))
    
    write.csv(outputData, file=paste(projectFolder,'/output.csv',sep=''),row.names=FALSE) #salvando os dados do cenario

}
