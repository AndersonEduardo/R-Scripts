##investigando determinantes da doenca de chagas em escala macroecologica##

##pacotes
library(raster)
library(rgdal)

##arquivos raster dos dados
mapaAbun = ##FAZER LA NO SCRIPT DO LUCAS!
    
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


##########################################
###trabalhando na escala dos MUNICIPIOS###
##########################################


municipios = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_municipios/BRMUE250GC_SIR.shp',layer='BRMUE250GC_SIR')
tabInfecbMuni = read.csv('/home/anderson/Documentos/Projetos/macroecologia_de_chagas/infeccao-municipio-2007-2014.csv',header=TRUE)
muniShapes = municipios[match(toupper(tabInfecbMuni$municipio),municipios$NM_MUNICIP),]

riqMedMuni = sapply(extract(mapaRiqueza,muniShapes,na.rm=TRUE), mean, na.rm=TRUE)
riscMedMuni = sapply(extract(mapaRisco,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
riqMedMuniHii = sapply(extract(mapaRiquezaHii,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
riscMedMuniHii = sapply(extract(mapaRiscoHii,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
HiiMedMuni = sapply(extract(mapaHII,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
areaMuni = sapply(extract(mapaRiqueza,muniShapes,na.rm=TRUE), length)

tabDadosMuni = data.frame(municipios = muniShapes$NM_MUNICIP,
                          riqMedMuni=riqMedMuni,
                          riscMedMuni=riscMedMuni,
                          riqMedMuniHii=riqMedMuniHii,
                          riscMedMuniHii=riscMedMuniHii,
                          HiiMedMuni=HiiMedMuni,
                          areaMuni=areaMuni)

##salvando a tabela com dados das variaveis extraidas por municipio
write.csv(tabDadosMuni,file="/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosMuni.csv",row.names=FALSE)

##abrindo tabela de dados municipios
tabDadosMuni = read.cssv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosMuni.csv",header=TRUE)

##dados de doenca de chagas do ministerio da saude
tabBarbMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/infeccao-municipio-2007-2014.csv",header=TRUE)

##adicionando a tabela de dados de municipios
tabDadosMuni$casos = tabBarbMuni$casos

##salvando a tabela final
write.csv(tabDadosMuni,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosMuniCompleta.csv",row.names=FALSE)

##abrindo tabela dedados completa
tabDadosMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosMuniCompleta.csv",header=TRUE)

##modelos estatisticos

residRiq = glm(riqMedMuni ~ areaMuni, family='gaussian', data=tabDadosMuni)$residuals
residRisc = glm(riscMedMuni ~ areaMuni, family='gaussian', data=tabDadosMuni)$residuals
residRiqHii = glm(riqMedMuniHii ~ areaMuni, family='gaussian', data=tabDadosMuni)$residuals
residRiscHii = glm(riqMedMuniHii ~ areaMuni, family='gaussian', data=tabDadosMuni)$residuals
residHii = glm(HiiMedMuni ~ areaMuni, family='gaussian', data=tabDadosMuni)$residuals

##variavel rpeditora: riqueza de especies



#######################################################################################################################
#######################################################################################################################
###### CONTINUAR DAQUI: SUBSTITUIR VARIAVEIS PELOS RESIDUOS ###########################################################
#######################################################################################################################
#######################################################################################################################
##modelo1: riqueza media -> n de casos



modRiq = glm(casos~residRisc,family='poisson',data=tabDadosMuni)
summary(modRiq)
plot(tabDadosMuni$casos ~ residRisc)
abline(modRiq,col='red')

##modelo2: risco medio -> n de casos
modRisc = glm(casos~riscMedMuni,family='poisson',data=tabDadosMuni)
summary(modRisc)
plot(tabDadosMuni$casos ~ tabDadosMuni$riscMedMuni)
abline(modRisc,col='red')

##modelo3: riqueza media -> n de casos
modRiqHii = glm(casos~riqMedMuniHii,family='poisson',data=tabDadosMuni)
summary(modRiqHii)
plot(tabDadosMuni$casos ~ tabDadosMuni$riqMedMuniHii)
abline(modRiqHii,col='red')

##modelo4: risco medio -> n de casos
modRiscHii = glm(casos~riscMedMuniHii,family='poisson',data=tabDadosMuni)
summary(modRiscHii)
plot(tabDadosMuni$casos ~ tabDadosMuni$riscMedMuniHii)
abline(modRiscHii,col='red')

##modelo5: HII -> n de casos
modHii = glm(casos~HiiMedMuni,family='poisson',data=tabDadosMuni)
summary(modHii)
plot(tabDadosMuni$casos ~ tabDadosMuni$HiiMedMuni)
abline(modHii,col='red')

##tabela de resultados

tabResMuni = data.frame(
    preditora = c('riqueza','risco','riquezaHII','riscoHII','HII'),
    efeito = c(modRiq$coefficients[2],modRisc$coefficients[2],modRiqHii$coefficients[2],modRiscHii$coefficients[2],modHii$coefficients[2]),
    deviance = c(modRiq$deviance,modRisc$deviance,modRiqHii$deviance,modRiscHii$devianc,modHii$deviance),
    null_deviance = c(modRiq$null.deviance,modRisc$null.deviance,modRiqHii$null.deviance,modRiscHii$null.deviance,modHii$null.deviance),
    aic = c(modRiq$aic,modRisc$aic,modRiqHii$aic,modRiscHii$aic,modHii$aic),
    p_valor = c(coef(summary(modRiq))[8],coef(summary(modRisc))[8],coef(summary(modRiqHii))[8],coef(summary(modRiscHii))[8],coef(summary(modHii))[8]))

write.csv(tabResMuni,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabMuniOutputs.csv",row.names=FALSE)
