##investigando determinantes da doenca de chagas em escala macroecologica##

##pacotes
library(raster)
library(rgdal)

##arquivos raster dos dados
mapaRiqueza = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico/Mapas de riqueza/mapaRiquezaPresente.asc')
mapaRisco = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc')
mapaRiquezaHii = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaPresente.asc')
mapaRiscoHii = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoPresente.asc')
mapaHII = raster('/home/anderson/PosDoc/dados_ambientais/hii-s-america-geo-grid/res2-5/hii-2-5.asc')
estados = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp',layer='BRUFE250GC_SIR')
vetorEst = estados$NM_ESTADO

##tabela de resultados
tabDados=data.frame()

##extraindo dados para cada estado do Brasil
for (i in 1:length(vetorEst)){
    est_i = estados[estados$NM_ESTADO == vetorEst[i],]
    ##riqueza
    riqMed = mean(unlist(extract(mapaRiqueza,est_i,na.rm=TRUE)),na.rm=TRUE)
    riqSD = sd(unlist(extract(mapaRiqueza,est_i,na.rm=TRUE)),na.rm=TRUE)
    ##risco
    riscMed = mean(unlist(extract(mapaRisco,est_i,na.rm=TRUE)),na.rm=TRUE)
    riscSD = sd(unlist(extract(mapaRisco,est_i,na.rm=TRUE)),na.rm=TRUE)
    ##riqueza HII
    riqMedHii = mean(unlist(extract(mapaRiquezaHii,est_i,na.rm=TRUE)),na.rm=TRUE)
    riqSDHii = sd(unlist(extract(mapaRiquezaHii,est_i,na.rm=TRUE)),na.rm=TRUE)
    ##risco HII
    riscMedHii = mean(unlist(extract(mapaRiscoHii,est_i,na.rm=TRUE)),na.rm=TRUE)
    riscSDHii = sd(unlist(extract(mapaRiscoHii,est_i,na.rm=TRUE)),na.rm=TRUE)
    ##impacto humano
    hiiMed = mean(unlist(extract(mapaHII,est_i,na.rm=TRUE)),na.rm=TRUE)
    hiiSD = sd(unlist(extract(mapaHII,est_i,na.rm=TRUE)),na.rm=TRUE)

    tabDados = rbind(tabDados,data.frame(
                              estado = vetorEst[i],
                              riqMed=riqMed,
                              riqSD=riqSD,
                              riscMed=riscMed,
                              riscSD=riscSD,
                              riqMedHii=riqMedHii,
                              riqSDHii=riqSDHii,
                              riscMedHii=riscMedHii,
                              riscSDHii=riscSDHii,
                              hiiMed=hiiMed,
                              hiiSD=hiiSD))
}

##salvando a tabela com dados das variaveis extraidas por estado
write.csv(tabDados,file="/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDados.csv",row.names=FALSE)


###tirando valores das variaveis para municipios
municipios = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_municipios/BRMUE250GC_SIR.shp',layer='BRMUE250GC_SIR')
vetorMuni = municipios$NM_MUNICIP

tabDadosMuni = data.frame()

for (i in 1:length(vetorMuni)){
    print(paste(i,'de',length(vetorMuni),'municipios...'))
    mun_i = municipios[municipios$NM_MUNICIP == vetorMuni[i],]
    ##riqueza
    riqMed = mean(unlist(extract(mapaRiqueza,mun_i,na.rm=TRUE)),na.rm=TRUE)
    riqSD = sd(unlist(extract(mapaRiqueza,mun_i,na.rm=TRUE)),na.rm=TRUE)
    ##risco
    riscMed = mean(unlist(extract(mapaRisco,mun_i,na.rm=TRUE)),na.rm=TRUE)
    riscSD = sd(unlist(extract(mapaRisco,mun_i,na.rm=TRUE)),na.rm=TRUE)
    ##riqueza HII
    riqMedHii = mean(unlist(extract(mapaRiquezaHii,mun_i,na.rm=TRUE)),na.rm=TRUE)
    riqSDHii = sd(unlist(extract(mapaRiquezaHii,mun_i,na.rm=TRUE)),na.rm=TRUE)
    ##risco HII
    riscMedHii = mean(unlist(extract(mapaRiscoHii,mun_i,na.rm=TRUE)),na.rm=TRUE)
    riscSDHii = sd(unlist(extract(mapaRiscoHii,mun_i,na.rm=TRUE)),na.rm=TRUE)
    ##impacto humano
    hiiMed = mean(unlist(extract(mapaHII,mun_i,na.rm=TRUE)),na.rm=TRUE)
    hiiSD = sd(unlist(extract(mapaHII,mun_i,na.rm=TRUE)),na.rm=TRUE)

    tabDadosMuni = rbind(tabDadosMuni,data.frame(
                              municipio = vetorMuni[i],
                              riqMed=riqMed,
                              riqSD=riqSD,
                              riscMed=riscMed,
                              riscSD=riscSD,
                              riqMedHii=riqMedHii,
                              riqSDHii=riqSDHii,
                              riscMedHii=riscMedHii,
                              riscSDHii=riscSDHii,
                              hiiMed=hiiMed,
                              hiiSD=hiiSD))
}

##salvando a tabela com dados das variaveis extraidas por estado
write.csv(tabDadosMuni,file="/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosMuni.csv",row.names=FALSE)

##retirando o distrito federal e organizando em ordem alfabetica
#tabDados = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDados.csv",header=TRUE)
tabDados = tabDados[!tabDados$estado=='DISTRITO FEDERAL',] #retirar o DF
tabDados = tabDados[order(tabDados$estado),] #ordem alfabetica 

##dados de doenca de chagas do ministerio da saude
tabBarb = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/Incidencia media anual casos agudos 2000-2013.csv",header=TRUE)
tabBarb = tabBarb[order(tabBarb$UF),]

##adicionando dados de doenca de chagas do ministerio da saude na tabela de dados para cada estado
tabDados$NumCasos = c(tabBarb$NumCasos[1:10],NA,tabBarb$NumCasos[11:13],NA,tabBarb$NumCasos[14:24])
tabDados$MediaAno = c(tabBarb$MediaAno[1:10],NA,tabBarb$MediaAno[11:13],NA,tabBarb$MediaAno[14:24])
tabDados$Incidencia = c(tabBarb$Incidencia[1:10],NA,tabBarb$Incidencia[11:13],NA,tabBarb$Incidencia[14:24])

##adicionando dados do tamanho territorial dos estados
##link: http://www.ibge.gov.br/home/geociencias/areaterritorial/principal.shtm
tabArea = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/territorio_estados.csv",header=TRUE)
tabArea = tabArea[order(tabArea$estado),] #colocando em ordem alfabetica
tabArea = tabArea[tabArea$estado!='Distrito_Federal',] #retirando o distrito federal
tabDados$area = tabArea$area

##salvando a tabela final
write.csv(tabDados,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta.csv",row.names=FALSE)

##abrindo tabela dedados completa
tabDados = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta.csv",header=TRUE)
#tabDados = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta2.csv",header=TRUE)

##residuos: retirando o efeito dos tamanhos
residRiq = glm(riqMed~area,family='gaussian',data=tabDados)$residuals
residRiqSD = glm(riqSD~area,family='gaussian',data=tabDados)$residuals
residRisc = glm(riscMed~area,family='gaussian',data=tabDados)$residuals
residRiscSD = glm(riscMed~area,family='gaussian',data=tabDados)$residuals
#
residRiqHii = glm(riqMedHii~area,family='gaussian',data=tabDados)$residuals
residRiqSDHii = glm(riqSDHii~area,family='gaussian',data=tabDados)$residuals
residRiscHii = glm(riscMedHii~area,family='gaussian',data=tabDados)$residuals
residRiscSDHii = glm(riscMedHii~area,family='gaussian',data=tabDados)$residuals
#
residHii = glm(hiiMed~area,family='gaussian',data=tabDados)$residuals
residHiiSD = glm(hiiSD~area,family='gaussian',data=tabDados)$residuals

##modelos estatisticos

##variavel rpeditora: riqueza de especies

##modelo1: riqueza media -> n de casos
modRiqNum = glm(tabDados$NumCasos~residRiq,family='poisson')
summary(modRiqNum)
plot(tabDados$NumCasos ~ residRiq)
abline(modRiqNum,col='red')

##modelo2: desv. pad. riqueza -> n de casos
modRiqNumSD = glm(tabDados$NumCasos~residRiqSD,family='poisson')
summary(modRiqNumSD)
plot(tabDados$NumCasos ~ residRiqSD)
abline(modRiqNumSD,col='red')

##modelo3: riqueza media -> incidencia
modRiqIncid = glm(tabDados$Incidencia~residRiq,family='gaussian')
summary(modRiqIncid)
plot(tabDados$Incidencia ~ residRiq)
abline(modRiqIncid,col='red')

##modelo4: desv. pad. riqueza -> incidencia
modRiqIncidSD = glm(tabDados$Incidencia~residRiqSD,family='gaussian')
summary(modRiqIncidSD)
plot(tabDados$Incidencia ~ residRiqSD)
abline(modRiqIncidSD,col='red')

##modelo5: riqueza media -> numero medio de casos por ano
modRiqMedAno = glm(tabDados$MediaAno~residRiq,family='gaussian')
summary(modRiqMedAno)
plot(tabDados$MediaAno ~ residRiq)
abline(modRiqMedAno,col='red')

##modelo6: desv. pad. riqueza -> n de casos
modRiqMedAnoSD = glm(tabDados$MediaAno~residRiqSD,family='gaussian')
summary(modRiqIncidSD)
plot(tabDados$Incidencia ~ residRiqSD)
abline(modRiqMedAnoSD,col='red')

##variavel preditora: risco

##modelo1: risco media -> n de casos
modRiscNum = glm(tabDados$NumCasos~residRisc,family='poisson')
summary(modRiscNum)
plot(tabDados$NumCasos ~ residRisc)
abline(modRiscNum,col='red')

##modelo2: desv. pad. risco -> n de casos
modRiscNumSD = glm(tabDados$NumCasos~residRiscSD,family='poisson')
summary(modRiscNumSD)
plot(tabDados$NumCasos ~ residRiscSD)
abline(modRiscNumSD,col='red')

##modelo3: risco media -> incidencia
modRiscIncid = glm(tabDados$Incidencia~residRisc,family='gaussian')
summary(modRiscIncid)
plot(tabDados$Incidencia ~ residRisc)
abline(modRiscIncid,col='red')

##modelo4: desv. pad. risco -> incidencia
modRiscIncidSD = glm(tabDados$Incidencia~residRiscSD,family='gaussian')
summary(modRiscIncidSD)
plot(tabDados$Incidencia ~ residRiscSD)
abline(modRiscIncidSD,col='red')

##modelo5: risco media -> numero medio de casos por ano
modRiscMedAno = glm(tabDados$MediaAno~residRisc,family='gaussian')
summary(modRiscMedAno)
plot(tabDados$MediaAno ~ residRisc)
abline(modRiscMedAno,col='red')

##modelo6: desv. pad. risco -> n de casos
modRiscMedAnoSD = glm(tabDados$MediaAno~residRiscSD,family='poisson')
summary(modRiscIncidSD)
plot(tabDados$Incidencia ~ residRiscSD)
abline(modRiscMedAnoSD,col='red')

##variavel rpeditora: riqueza de especies HII

##modelo1: riqueza media -> n de casos
modRiqNumHii = glm(tabDados$NumCasos~residRiqHii,family='poisson')
summary(modRiqNumHii)
plot(tabDados$NumCasos ~ residRiqHii)
abline(modRiqNumHii,col='red')

##modelo2: desv. pad. riqueza -> n de casos
modRiqNumSDHii = glm(tabDados$NumCasos~residRiqSDHii,family='poisson')
summary(modRiqNumSDHii)
plot(tabDados$NumCasos ~ residRiqSDHii)
abline(modRiqNumSDHii,col='red')

##modelo3: riqueza media -> incidencia
modRiqIncidHii = glm(tabDados$Incidencia~residRiqHii,family='gaussian')
summary(modRiqIncidHii)
plot(tabDados$Incidencia ~ residRiqHii)
abline(modRiqIncidHii,col='red')

##modelo4: desv. pad. riqueza -> incidencia
modRiqIncidSDHii = glm(tabDados$Incidencia~residRiqSDHii,family='gaussian')
summary(modRiqIncidSDHii)
plot(tabDados$Incidencia ~ residRiqSDHii)
abline(modRiqIncidSDHii,col='red')

##modelo5: riqueza media -> numero medio de casos por ano
modRiqMedAnoHii = glm(tabDados$MediaAno~residRiqHii,family='gaussian')
summary(modRiqMedAnoHii)
plot(tabDados$MediaAno ~ residRiqHii)
abline(modRiqMedAnoHii,col='red')

##modelo6: desv. pad. riqueza -> n de casos
modRiqMedAnoSDHii = glm(tabDados$MediaAno~residRiqSDHii,family='gaussian')
summary(modRiqIncidSDHii)
plot(tabDados$Incidencia ~ residRiqSDHii)
abline(modRiqMedAnoSDHii,col='red')

##variavel preditora: risco HII

##modelo1: risco media -> n de casos
modRiscNumHii = glm(tabDados$NumCasos~residRiscHii,family='poisson')
summary(modRiscNumHii)
plot(tabDados$NumCasos ~ residRiscHii)
abline(modRiscNumHii,col='red')

##modelo2: desv. pad. risco -> n de casos
modRiscNumSDHii = glm(tabDados$NumCasos~residRiscSDHii,family='poisson')
summary(modRiscNumSDHii)
plot(tabDados$NumCasos ~ residRiscSDHii)
abline(modRiscNumSDHii,col='red')

##modelo3: risco media -> incidencia
modRiscIncidHii = glm(tabDados$Incidencia~residRiscHii,family='gaussian')
summary(modRiscIncidHii)
plot(tabDados$Incidencia ~ residRiscHii)
abline(modRiscIncidHii,col='red')

##modelo4: desv. pad. risco -> incidencia
modRiscIncidSDHii = glm(tabDados$Incidencia~residRiscSDHii,family='gaussian')
summary(modRiscIncidSDHii)
plot(tabDados$Incidencia ~ residRiscSDHii)
abline(modRiscIncidSDHii,col='red')

##modelo5: risco media -> numero medio de casos por ano
modRiscMedAnoHii = glm(tabDados$MediaAno~residRiscHii,family='gaussian')
summary(modRiscMedAnoHii)
plot(tabDados$MediaAno ~ residRiscHii)
abline(modRiscMedAnoHii,col='red')

##modelo6: desv. pad. risco -> n de casos
modRiscMedAnoSDHii = glm(tabDados$MediaAno~residRiscSDHii,family='gaussian')
summary(modRiscIncidSDHii)
plot(tabDados$Incidencia ~ residRiscSDHii)
abline(modRiscMedAnoSDHii,col='red')

##variavel preditora: indice de imapcto humano

##modelo1: impacto humano medio -> numero de casos
modHiiMeanNum = glm(NumCasos~residHii,family='poisson',data=tabDados)
summary(modHiiMeanNum)
plot(tabDados$NumCasos ~ residHii)
abline(modHiiMeanNum,col='red')

##modelo2: desv. pad. impacto humano -> numero de casos
modHiiSDNum = glm(NumCasos~residHiiSD,family='poisson',data=tabDados)
summary(modHiiSDNum)
plot(tabDados$NumCasos ~ residHiiSD)
abline(modHiiSDNum,col='red')

##modelo3: impacto humano medio -> incidencia
modHiiMeanIncid = glm(Incidencia~residHii,family='gaussian',data=tabDados)
summary(modHiiMeanIncid)
plot(tabDados$Incidencia ~ residHii)
abline(modHiiMeanIncid,col='red')

##modelo4: desv. pad. impacto humano -> incidencia
modHiiSDIncid = glm(Incidencia~residHiiSD,family='gaussian',data=tabDados)
summary(modHiiSDIncid)
plot(tabDados$Incidencia ~ residHiiSD)
abline(modHiiSDIncid,col='red')

##modelo5: impacto humano medio -> numero medio de casos por ano
modHiiMeanMedAno = glm(MediaAno~residHii,family='gaussian',data=tabDados)
summary(modHiiMeanIncid)
plot(tabDados$MediaAno ~ residHii)
abline(modHiiMeanMedAno,col='red')

##modelo6: desv. pad. impacto humano -> numero medio de casos por ano
modHiiSDMedAno = glm(MediaAno~residHiiSD,family='gaussian',data=tabDados)
summary(modHiiSDMedAno)
plot(tabDados$MediaAno ~ residHii)
abline(modHiiSDMedAno,col='red')


##tabela de resultados

tabRes = data.frame(
    preditora = c(rep(c('riqueza','risco','riquezaHII','riscoHII','HII'),3)),
    resposta = c(rep('Num_de_casos',5),rep('Incidencia',5),rep('Num_med_casos/ano',5)),
    efeito = c(modRiqNum$coefficients[2],modRiscNum$coefficients[2],modRiqNumHii$coefficients[2],modRiscNumHii$coefficients[2],modHiiMeanNum$coefficients[2],
            modRiqIncid$coefficients[2],modRiscIncid$coefficients[2],modRiqIncidHii$coefficients[2],modRiscIncidHii$coefficients[2],modHiiMeanIncid$coefficients[2],
            modRiqMedAno$coefficients[2],modRiscMedAno$coefficients[2],modRiqMedAnoHii$coefficients[2],modRiscMedAnoHii$coefficients[2],modHiiMeanMedAno$coefficients[2]),
    deviance = c(modRiqNum$deviance,modRiscNum$deviance,modRiqNumHii$deviance,modRiscNumHii$devianc,modHiiMeanNum$deviance,
                 modRiqIncid$deviance,modRiscIncid$deviance,modRiqIncidHii$deviance,modRiscIncidHii$deviance,modHiiMeanIncid$deviance,
                 modRiqMedAno$deviance,modRiscMedAno$deviance,modRiqMedAnoHii$deviance,modRiscMedAnoHii$deviance,modHiiMeanMedAno$deviance),
    null_deviance = c(modRiqNum$null.deviance,modRiscNum$null.deviance,modRiqNumHii$null.deviance,modRiscNumHii$null.deviance,modHiiMeanNum$null.deviance,
                      modRiqIncid$null.deviance,modRiscIncid$null.deviance,modRiqIncidHii$null.deviance,modRiscIncidHii$null.deviance,modHiiMeanIncid$null.deviance,
                      modRiqMedAno$null.deviance,modRiscMedAno$null.deviance,modRiqMedAnoHii$null.deviance,modRiscMedAnoHii$null.deviance,modHiiMeanMedAno$null.deviance),
    aic = c(modRiqNum$aic,modRiscNum$aic,modRiqNumHii$aic,modRiscNumHii$aic,modHiiMeanNum$aic,
            modRiqIncid$aic,modRiscIncid$aic,modRiqIncidHii$aic,modRiscIncidHii$aic,modHiiMeanIncid$aic,
            modRiqMedAno$aic,modRiscMedAno$aic,modRiqMedAnoHii$aic,modRiscMedAnoHii$aic,modHiiMeanMedAno$aic),
    p_valor = c(coef(summary(modRiqNum))[8],coef(summary(modRiscNum))[8],coef(summary(modRiqNumHii))[8],coef(summary(modRiscNumHii))[8],coef(summary(modHiiMeanNum))[8],
                coef(summary(modRiqIncid))[8],coef(summary(modRiscIncid))[8],coef(summary(modRiqIncidHii))[8],coef(summary(modRiscIncidHii))[8],coef(summary(modHiiMeanIncid))[8],
                coef(summary(modRiqMedAno))[8],coef(summary(modRiscMedAno))[8],coef(summary(modRiqMedAnoHii))[8],coef(summary(modRiscMedAnoHii))[8],coef(summary(modHiiMeanMedAno))[8])
)

write.csv(tabRes,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabOutputs.csv",row.names=FALSE)


##usando dados novos do SUS

##abrindo tabela de dados (com dados novos, do dataSUS)
tabDados = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta2.csv",header=TRUE,dec='.')

##modelos estatisticos

##OBS: X2007.2014 = casos 'agudos'; X2001.2006 = casos 'cronicos'

##variavel rpeditora: riqueza de especies

##modelo1: riqueza media -> n casos agudos
modRiqAgu = glm(tabDados$X2007.2014[-13]~residRiq[-13],family='poisson')
summary(modRiqAgu)
plot(tabDados$X2007.2014[-13] ~ residRiq[-13])
abline(modRiqAgu,col='red')

##modelo2: desv. pad. riqueza ->  n casos agudos
modRiqAguSD = glm(tabDados$X2007.2014[-13]~residRiqSD[-13],family='poisson')
summary(modRiqAguSD)
plot(tabDados$X2007.2014[-13] ~ residRiqSD[-13])
abline(modRiqAguSD,col='red')

##modelo3: riqueza media ->  n casos cronicos
modRiqCro = glm(tabDados$X2001.2006~residRiq,family='poisson')
summary(modRiqCro)
plot(tabDados$X2001.2006 ~ residRiq)
abline(modRiqCro,col='red')

##modelo4: desv. pad. riqueza -> n casos cronicos
modRiqCroSD = glm(tabDados$X2001.2006~residRiqSD,family='poisson')
summary(modRiqCroSD)
plot(tabDados$X2001.2006 ~ residRiqSD)
abline(modRiqCroSD,col='red')

##variavel preditora: risco

##modelo1: risco medio -> n casos agudos
modRiscAgu = glm(tabDados$X2007.2014[-13]~residRisc[-13],family='poisson')
summary(modRiscAgu)
plot(tabDados$X2007.2014[-13] ~ residRisc[-13])
abline(modRiscAgu,col='red')

##modelo2: desv. pad. risco ->  n casos agudos
modRiscAguSD = glm(tabDados$X2007.2014[-13]~residRiscSD[-13],family='poisson')
summary(modRiscAguSD)
plot(tabDados$X2007.2014[-13] ~ residRiscSD[-13])
abline(modRiscAguSD,col='red')

##modelo3: risco medio ->  n casos cronicos
modRiscCro = glm(tabDados$X2001.2006~residRisc,family='poisson')
summary(modRiscCro)
plot(tabDados$X2001.2006 ~ residRisc)
abline(modRiscCro,col='red')

##modelo4: desv. pad. risco -> n casos cronicos
modRiscCroSD = glm(tabDados$X2001.2006~residRiscSD,family='poisson')
summary(modRiscCroSD)
plot(tabDados$X2001.2006 ~ residRiscSD)
abline(modRiscCroSD,col='red')


##variavel rpeditora: riqueza de especies HII

##modelo1: riqueza media -> n casos agudos
modRiqAguHii = glm(tabDados$X2007.2014[-13]~residRiqHii[-13],family='poisson')
summary(modRiqAguHii)
plot(tabDados$X2007.2014[-13] ~ residRiqHii[-13])
abline(modRiqAguHii,col='red')

##modelo2: desv. pad. riqueza ->  n casos agudos
modRiqAguSDHii = glm(tabDados$X2007.2014[-13]~residRiqSDHii[-13],family='poisson')
summary(modRiqAguSDHii)
plot(tabDados$X2007.2014[-13] ~ residRiqSDHii[-13])
abline(modRiqAguSDHii,col='red')

##modelo3: riqueza media ->  n casos cronicos
modRiqCroHii = glm(tabDados$X2001.2006~residRiqHii,family='poisson')
summary(modRiqCroHii)
plot(tabDados$X2001.2006 ~ residRiqHii)
abline(modRiqCroHii,col='red')

##modelo4: desv. pad. riqueza -> n casos cronicos
modRiqCroSDHii = glm(tabDados$X2001.2006~residRiqSDHii,family='poisson')
summary(modRiqCroSDHii)
plot(tabDados$X2001.2006 ~ residRiqSDHii)
abline(modRiqCroSDHii,col='red')

##variavel preditora: risco HII

##modelo1: risco medio -> n casos agudos
modRiscAguHii = glm(tabDados$X2007.2014[-13]~residRiscHii[-13],family='poisson')
summary(modRiscAguHii)
plot(tabDados$X2007.2014[-13] ~ residRiscHii[-13])
abline(modRiscAguHii,col='red')

##modelo2: desv. pad. risco ->  n casos agudos
modRiscAguSDHii = glm(tabDados$X2007.2014[-13]~residRiscSDHii[-13],family='poisson')
summary(modRiscAguSDHii)
plot(tabDados$X2007.2014[-13] ~ residRiscSD[-13])
abline(modRiscAguSDHii,col='red')

##modelo3: risco medio ->  n casos cronicos
modRiscCroHii = glm(tabDados$X2001.2006~residRiscHii,family='poisson')
summary(modRiscCroHii)
plot(tabDados$X2001.2006 ~ residRiscHii)
abline(modRiscCroHii,col='red')

##modelo4: desv. pad. risco -> n casos cronicos
modRiscCroSDHii = glm(tabDados$X2001.2006~residRiscSDHii,family='poisson')
summary(modRiscCroSDHii)
plot(tabDados$X2001.2006 ~ residRiscSDHii)
abline(modRiscCroSDHii,col='red')


##variavel preditora: indice de imapcto humano

##modelo1: impacto humano medio -> numero de casos agudos
modHiiMeanAgu = glm(tabDados$X2007.2014[-13]~residHii[-13],family='poisson')
summary(modHiiMeanAgu)
plot(tabDados$X2007.2014[-13] ~ residHii[-13])
abline(modHiiMeanAgu,col='red')

##modelo2: desv. pad. impacto humano -> numero de casos agudos
modHiiSDAgu = glm(tabDados$X2007.2014[-13]~residHiiSD[-13],family='poisson')
summary(modHiiSDAgu)
plot(tabDados$X2007.2014[-13] ~ residHiiSD[-13])
abline(modHiiSDAgu,col='red')

##modelo3: impacto humano medio -> n casos cronicos
modHiiMeanCro = glm(tabDados$X2001.2006~residHii,family='poisson')
summary(modHiiMeanCro)
plot(tabDados$X2001.2006 ~ residHii)
abline(modHiiMeanCro,col='red')

##modelo4: desv. pad. impacto humano -> n casos cronicos
modHiiSDCro = glm(tabDados$X2001.2006~residHiiSD,family='poisson')
summary(modHiiSDCro)
plot(tabDados$X2001.2006 ~ residHiiSD)
abline(modHiiSDCro,col='red')


##tabela de outputs

tabResDATASUS = data.frame(
    preditora = c(rep(c('riqueza','risco','riquezaHII','riscoHII','HII'),2)),
    resposta = c(rep('2001-2006',5),rep('2007-1014',5)),
    efeito = c(modRiqAgu$coefficients[2],modRiscAgu$coefficients[2],modRiqAguHii$coefficients[2],modRiscAguHii$coefficients[2],modHiiMeanAgu$coefficients[2],
            modRiqCro$coefficients[2],modRiscCro$coefficients[2],modRiqCroHii$coefficients[2],modRiscCroHii$coefficients[2],modHiiMeanCro$coefficients[2]),
    deviance = c(modRiqAgu$deviance,modRiscAgu$deviance,modRiqAguHii$deviance,modRiscAguHii$devianc,modHiiMeanAgu$deviance,
                 modRiqCro$deviance,modRiscCro$deviance,modRiqCroHii$deviance,modRiscCroHii$deviance,modHiiMeanCro$deviance),
    null_deviance = c(modRiqAgu$null.deviance,modRiscAgu$null.deviance,modRiqAguHii$null.deviance,modRiscAguHii$null.deviance,modHiiMeanAgu$null.deviance,
                      modRiqCro$null.deviance,modRiscCro$null.deviance,modRiqCroHii$null.deviance,modRiscCroHii$null.deviance,modHiiMeanCro$null.deviance),
    aic = c(modRiqAgu$aic,modRiscAgu$aic,modRiqAguHii$aic,modRiscAguHii$aic,modHiiMeanAgu$aic,
            modRiqCro$aic,modRiscCro$aic,modRiqCroHii$aic,modRiscCroHii$aic,modHiiMeanCro$aic),
    p_valor = c(coef(summary(modRiqAgu))[8],coef(summary(modRiscAgu))[8],coef(summary(modRiqAguHii))[8],coef(summary(modRiscAguHii))[8],coef(summary(modHiiMeanAgu))[8],
                coef(summary(modRiqCro))[8],coef(summary(modRiscCro))[8],coef(summary(modRiqCroHii))[8],coef(summary(modRiscCroHii))[8],coef(summary(modHiiMeanCro))[8]))

write.csv(tabResDATASUS,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabOutputsDATASUS.csv",row.names=FALSE)
