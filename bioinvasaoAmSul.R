## invasao mmuferos america do sul

library(biomod2)
library(raster)

globalVars = stack(list.files('/home/anderson/PosDoc/dados_ambientais/bio_2-5m_bil', pattern='.bil', full.names=TRUE))

plot(globalVars[[1]])
extEuro = drawExtent()
extAm = drawExtent()
extAmSul = drawExtent()

euroVars = crop(globalVars[[c('bio1','bio5','bio6')]], extEuro)
amerVars = crop(globalVars[[c('bio1','bio5','bio6')]], extAm)
amerSulVars = crop(globalVars[[c('bio1','bio5','bio6')]], extAmSul)



sampleFolder = '/home/anderson/Área de Trabalho/sps_invasao'

spsVector = list.files(sampleFolder, pattern='nativo', full.names=TRUE )

nativoData = read.csv(spsVector[3], header=TRUE)
names(nativoData) = c('sp','lon','lat')



##variaveis e parametros locais especificos para o biomod2

occData = nativoData[,c('lon','lat')]
##
spName = gsub('.csv','',list.files(sampleFolder, pattern='nativo')[2])
spName =  gsub(' ','_',spName)
myRespName <- spName
##
myResp <- occData # variavel resposta (para biomod2)
coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints
##
myExpl = stack(euroVars)


##ajuste de dados de entrada para biomod2
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1)
      
##rodando o(s) algoritmo(s) (i.e. SDMs)
myBiomodModelOut <- BIOMOD_Modeling(
    data = myBiomodData,
    models = c('MAXENT.Phillips'),
    models.options = myBiomodOption,
    NbRunEval = 100,
    DataSplit = 75,
    models.eval.meth = c('TSS','ROC'),
    do.full.models = FALSE,
    modeling.id = myRespName)

##My output data
evaluationScores = get_evaluations(myBiomodModelOut)

##rodando algortmo de projecao (i.e. rodando a projecao)
myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl_invadido,
    proj.name = myRespName,
    compress = 'TRUE',
    build.clamping.mask = 'TRUE')


#######################################################################
############# densidade da variavel temporetatura #####################
################ na area nativa e area invadida #######################
#######################################################################


##caminhos dos arquivos de dados
spsVector = list.files(sampleFolder, full.names=TRUE )

##vetor com o nome dos arquivos
spsNames = list.files(sampleFolder)
vecNames = gsub(pattern='_invadido.csv',replacement='',x=spsNames)
vecNames = gsub(pattern='_nativo.csv',replacement='',x=vecNames)
vecNames = unique(vecNames)


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


for(sp_i in vecNames){
    ##vetor de nomes de arquivo
    spsVector = list.files(sampleFolder, pattern=sp_i, full.names=TRUE)
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
    jpeg(paste('/home/anderson/Área de Trabalho/sps_invasao/',sp_i,'.jpeg', sep=''), main=)
    plot(density(tempNativa[,'bio1']), xlim=c(-200,500), lwd=2, col='blue')
    lines(density(tempInvadido[,'bio1']), lwd=2, col='red')
    abline(v=termNeut[termNeut$sps==sp_i,][,c('min','max')])
    dev.off()
}
