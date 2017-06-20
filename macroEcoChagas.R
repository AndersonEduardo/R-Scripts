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
    risco
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

##retirando o distrito federal e organizando em ordem alfabetica
#read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDados.csv",header=TRUE)
tabRes = tabRes[!tabRes$estado=='DISTRITO FEDERAL',] #retirar o DF
tabRes = tabRes[order(tabRes$estado),] #ordem alfabetica 

##dados de doenca de chagas do ministerio da saude
tabBarb = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/Incidencia media anual casos agudos 2000-2013.csv",header=TRUE)
tabBarb = tabBarb[order(tabBarb$UF),]

##adicionando dados de doenca de chagas do ministerio da saude na tabela de dados para cada estado
tabRes$NumCasos = c(tabBarb$NumCasos[1:10],NA,tabBarb$NumCasos[11:13],NA,tabBarb$NumCasos[14:24])
tabRes$MediaAno = c(tabBarb$MediaAno[1:10],NA,tabBarb$MediaAno[11:13],NA,tabBarb$MediaAno[14:24])
tabRes$Incidencia = c(tabBarb$Incidencia[1:10],NA,tabBarb$Incidencia[11:13],NA,tabBarb$Incidencia[14:24])

##adicionando dados do tamanho territorial dos estados
##link: http://www.ibge.gov.br/home/geociencias/areaterritorial/principal.shtm
tabArea = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/territorio_estados.csv",header=TRUE)
tabArea = tabArea[order(tabArea$estado),] #colocando em ordem alfabetica
tabArea = tabArea[tabArea$estado!='Distrito_Federal',] #retirando o distrito federal
tabRes$area = tabArea$area

##salvando a tabela final
write.csv(tabRes,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta.csv",row.names=FALSE)

##abrindo tabela de dados completa
tabRes = read.cev("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta.csv",header=TRUE)

##residuos: retirando o efeito dos tamanhos
residRiq = glm(riqMed~area,family='poisson',data=tabRes)$residuals
residRiqSD = glm(riqSD~area,family='poisson',data=tabRes)$residuals
residHii = glm(hiiMed~area,family='gaussian',data=tabRes)$residuals
residHiiSD = glm(hiiSD~area,family='gaussian',data=tabRes)$residuals
residRisc = glm(riscMed~area,family='gaussian',data=tabRes)$residuals
residRiscSD = glm(riscMed~area,family='gaussian',data=tabRes)$residuals

##usando dados novos do SUS

##abrindo tabela de dados
#tabDados = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta.csv",header=TRUE,dec='.') #(com dados novos, do dataSUS)
tabDados = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta2.csv",header=TRUE,dec='.') #(com dados novos, do dataSUS)


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


############################################
####trabalhando na escala dos MUNICIPIOS####
############################################


municipios = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_municipios/BRMUE250GC_SIR.shp',layer='BRMUE250GC_SIR')
tabInfecbMuni = read.csv('/home/anderson/Documentos/Projetos/macroecologia_de_chagas/infeccao-municipio-2007-2014.csv',header=TRUE)

##extraindo somente para municipios com dados SUS
muniShapes = municipios[match(toupper(tabInfecbMuni$municipio),municipios$NM_MUNICIP),] 
riqMedMuni = sapply(extract(mapaRiqueza,muniShapes,na.rm=TRUE), mean, na.rm=TRUE)
riscMedMuni = sapply(extract(mapaRisco,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
riqMedMuniHii = sapply(extract(mapaRiquezaHii,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
riscMedMuniHii = sapply(extract(mapaRiscoHii,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
HiiMedMuni = sapply(extract(mapaHII,muniShapes,na.rm=TRUE), mean,na.rm=TRUE)
areaMuni = sapply(extract(mapaRiqueza,muniShapes,na.rm=TRUE), length)

##extraindo para TODOS os municipios
riqMedMuni = sapply(extract(mapaRiqueza,municipios,na.rm=TRUE), mean, na.rm=TRUE)
riscMedMuni = sapply(extract(mapaRisco,municipios,na.rm=TRUE), mean,na.rm=TRUE)
riqMedMuniHii = sapply(extract(mapaRiquezaHii,municipios,na.rm=TRUE), mean,na.rm=TRUE)
riscMedMuniHii = sapply(extract(mapaRiscoHii,municipios,na.rm=TRUE), mean,na.rm=TRUE)
HiiMedMuni = sapply(extract(mapaHII,municipios,na.rm=TRUE), mean,na.rm=TRUE)
areaMuni = sapply(extract(mapaRiqueza,municipios,na.rm=TRUE), length)


tabDadosTodosMuni = data.frame(municipios = municipios$NM_MUNICIP,
                          riqMedMuni=riqMedMuni,
                          riscMedMuni=riscMedMuni,
                          riqMedMuniHii=riqMedMuniHii,
                          riscMedMuniHii=riscMedMuniHii,
                          HiiMedMuni=HiiMedMuni,
                          areaMuni=areaMuni)


##salvando a tabela com dados das variaveis extraidas por municipio
write.csv(tabDadosTodosMuni,file="/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuni.csv",row.names=FALSE)

##abrindo tabela de dados municipios
tabDadosTodosMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuni.csv",header=TRUE)

##dados de doenca de chagas do ministerio da saude
tabBarbMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/infeccao-municipio-2007-2014.csv",header=TRUE)

##adicionando a tabela de dados de municipios
vetorCasos = as.integer( tabDadosTodosMuni$municipios %in% toupper(tabBarbMuni$municipio) ) #vetor de 0 e 1
vetorLinhas = base::match(toupper(tabBarbMuni$municipio),tabDadosTodosMuni$municipios) #encontrando as linhas que tem dados do SUS
vetorCasos[vetorLinhas] = tabBarbMuni$casos #ajustando para '0' e 'numero de casos',  em vez de '0' e '1'
tabDadosTodosMuni$casos = vetorCasos #adicionando na tabela de dados

##salvando a tabela final
write.csv(tabDadosTodosMuni,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuniCompleta.csv",row.names=FALSE)

##abrindo tabela dedados completa
tabDadosTodosMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuniCompleta.csv",header=TRUE)


##sobredispersao

##nuemro de casos
boxplot(tabDadosTodosMuni$casos)
plot(tabDadosTodosMuni$casos)
hist(tabDadosTodosMuni$casos)
boxplot(tabDadosTodosMuni$casos[tabDadosTodosMuni$casos<100])
plot(tabDadosTodosMuni$casos[tabDadosTodosMuni$casos<100])
#
tabDadosTodosMuniControled = tabDadosTodosMuni[which(tabDadosTodosMuni$casos<100),]

##riqueza
boxplot(tabDadosTodosMuniControled$riqMedMuni)
plot(tabDadosTodosMuniControled$riqMedMuni)
hist(tabDadosTodosMuniControled$riqMedMuni)

##risco
boxplot(tabDadosTodosMuniControled$riscMedMuni)
plot(tabDadosTodosMuniControled$riscMedMuni)
hist(tabDadosTodosMuniControled$riscMedMuni)

##riqHii
boxplot(tabDadosTodosMuniControled$riqMedMuniHii)
plot(tabDadosTodosMuniControled$riqMedMuniHii)
hist(tabDadosTodosMuniControled$riqMedMuniHii)

##riscHii
boxplot(tabDadosTodosMuniControled$riscMedMuniHii)
plot(tabDadosTodosMuniControled$riscMedMuniHii)
hist(tabDadosTodosMuniControled$riscMedMuniHii)

##Hii
boxplot(tabDadosTodosMuniControled$HiiMedMuni)
plot(tabDadosTodosMuniControled$HiiMedMuni)
hist(tabDadosTodosMuniControled$HiiMedMuni)


##correlacao com area do municipio
plot(tabDadosTodosMuniControled$casos ~ tabDadosTodosMuniControled$areaMuni)
plot(tabDadosTodosMuniControled$casos[tabDadosTodosMuniControled$areaMuni<600] ~ tabDadosTodosMuniControled$areaMuni[tabDadosTodosMuniControled$areaMuni<600])

plot(tabDadosTodosMuniControled$riqMedMuni ~ tabDadosTodosMuniControled$areaMuni)
plot(tabDadosTodosMuniControled$riqMedMuni[tabDadosTodosMuniControled$areaMuni<600] ~ tabDadosTodosMuniControled$areaMuni[tabDadosTodosMuniControled$areaMuni<600])

plot(tabDadosTodosMuniControled$riscMedMuni ~ tabDadosTodosMuniControled$areaMuni)
plot(tabDadosTodosMuniControled$riscMedMuni[tabDadosTodosMuniControled$areaMuni<600] ~ tabDadosTodosMuniControled$areaMuni[tabDadosTodosMuniControled$areaMuni<600]) ##FAZER MODELOS COM OS RESIDUOS

plot(tabDadosTodosMuniControled$riqMedMuniHii ~ tabDadosTodosMuniControled$areaMuni)
plot(tabDadosTodosMuniControled$riqMedMuniHii[tabDadosTodosMuniControled$areaMuni<600] ~ tabDadosTodosMuniControled$areaMuni[tabDadosTodosMuniControled$areaMuni<600]) ##FAZER OS MODELOS COM OS RESIDUOS

plot(tabDadosTodosMuniControled$riscMedMuniHii ~ tabDadosTodosMuniControled$areaMuni)
plot(tabDadosTodosMuniControled$riscMedMuniHii[tabDadosTodosMuniControled$areaMuni<600] ~ tabDadosTodosMuniControled$areaMuni[tabDadosTodosMuniControled$areaMuni<600]) ##FAZER OS MODELOS COM OS RESIDUOS

plot(log(tabDadosTodosMuniControled$HiiMedMuni) ~ tabDadosTodosMuniControled$areaMuni)
plot(log(tabDadosTodosMuniControled$HiiMedMuni[tabDadosTodosMuniControled$areaMuni<600]) ~ tabDadosTodosMuniControled$areaMuni[tabDadosTodosMuniControled$areaMuni<600]) ##FAZER OS MODELOS COM OS RESIDUOS

tabDadosTodosMuniControled = tabDadosTodosMuniControled[which(tabDadosTodosMuniControled$areaMuni<600),]


##residuos da area do municipio
residRiq = glm(log(riqMedMuni) ~ areaMuni, family='gaussian', data=tabDadosTodosMuniControled)$residuals
residRisc = glm(log(riscMedMuni) ~ areaMuni, family='gaussian', data=tabDadosTodosMuniControled)$residuals
residRiqHii = glm(riqMedMuniHii ~ areaMuni, family='gaussian', data=tabDadosTodosMuniControled)$residuals
residRiscHii = glm(riqMedMuniHii ~ areaMuni, family='gaussian', data=tabDadosTodosMuniControled)$residuals
residHii = glm(HiiMedMuni ~ areaMuni, family='gaussian', data=tabDadosTodosMuniControled)$residuals


##modelos estatisticos

##modelo1: riqueza media -> n de casos
modRiq = glm(casos[residRiq>-1.5][-1]~residRiq[residRiq>-1.5],family='poisson',data=tabDadosTodosMuniControled)
summary(modRiq)
plot(tabDadosTodosMuniControled$casos[residRiq>-1.5][-1] ~ residRiq[residRiq>-1.5])
abline(modRiq,col='red')

##modelo2: risco medio -> n de casos
modRisc = glm(casos~riscMedMuni,family='poisson',data=tabDadosTodosMuniControled)
summary(modRisc)
plot(tabDadosTodosMuniControled$casos ~ tabDadosTodosMuniControled$riscMedMuni)
abline(modRisc,col='red')

##modelo3: riqueza media -> n de casos
modRiqHii = glm(casos[-1]~residRiqHii,family='poisson',data=tabDadosTodosMuniControled)
summary(modRiqHii)
plot(tabDadosTodosMuniControled$casos[-1] ~ residRiqHii)
abline(modRiqHii,col='red')

##modelo4: risco medio -> n de casos
modRiscHii = glm(casos[-1]~residRiscHii,family='poisson',data=tabDadosTodosMuniControled)
summary(modRiscHii)
plot(tabDadosTodosMuniControled$casos[-1] ~ residRiscHii)
abline(modRiscHii,col='red')

##modelo5: HII -> n de casos
modHii = glm(casos[-1]~residHii,family='poisson',data=tabDadosTodosMuniControled)
summary(modHii)
plot(tabDadosTodosMuniControled$casos[-1] ~ residHii)
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


##################################################################
##############  UTILIZANDO MODELO BINNOMIAL  #####################
##################################################################

##abrindo tabela de dados municipios
tabDadosTodosMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuni.csv",header=TRUE)

##dados de doenca de chagas do ministerio da saude
tabBarbMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/infeccao-municipio-2007-2014.csv",header=TRUE)

##adicionando a tabela de dados de municipios
vetorCasos = as.integer( tabDadosTodosMuni$municipios %in% toupper(tabBarbMuni$municipio) ) #vetor de 0 e 1
vetorLinhas = base::match(toupper(tabBarbMuni$municipio),tabDadosTodosMuni$municipios) #encontrando as linhas que tem dados do SUS
tabDadosTodosMuni$casos = vetorCasos #adicionando na tabela de dados

##salvando a tabela final
write.csv(tabDadosTodosMuni,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuniCompletaLogist.csv",row.names=FALSE)

##abrindo tabela dedados completa
tabDadosTodosMuni = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosTodosMuniCompletaLogist.csv",header=TRUE)


##sobredispersao

##riqueza
boxplot(tabDadosTodosMuni$riqMedMuni)
plot(tabDadosTodosMuni$riqMedMuni)
hist(tabDadosTodosMuni$riqMedMuni)

##risco
boxplot(tabDadosTodosMuni$riscMedMuni)
plot(tabDadosTodosMuni$riscMedMuni)
hist(tabDadosTodosMuni$riscMedMuni)

##riqHii
boxplot(tabDadosTodosMuni$riqMedMuniHii)
plot(tabDadosTodosMuni$riqMedMuniHii)
hist(tabDadosTodosMuni$riqMedMuniHii)

##riscHii
boxplot(tabDadosTodosMuni$riscMedMuniHii)
plot(tabDadosTodosMuni$riscMedMuniHii)
hist(tabDadosTodosMuni$riscMedMuniHii)

##Hii
boxplot(tabDadosTodosMuni$HiiMedMuni)
plot(tabDadosTodosMuni$HiiMedMuni)
hist(tabDadosTodosMuni$HiiMedMuni)


##correlacao com area do municipio
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$areaMuni)
plot(tabDadosTodosMuni$casos[tabDadosTodosMuni$areaMuni<600] ~ tabDadosTodosMuni$areaMuni[tabDadosTodosMuni$areaMuni<600])

plot(tabDadosTodosMuni$riqMedMuni ~ tabDadosTodosMuni$areaMuni)
plot(tabDadosTodosMuni$riqMedMuni[tabDadosTodosMuni$areaMuni<600] ~ tabDadosTodosMuni$areaMuni[tabDadosTodosMuni$areaMuni<600])

plot(tabDadosTodosMuni$riscMedMuni ~ tabDadosTodosMuni$areaMuni)
plot(tabDadosTodosMuni$riscMedMuni[tabDadosTodosMuni$areaMuni<600] ~ tabDadosTodosMuni$areaMuni[tabDadosTodosMuni$areaMuni<600]) ##FAZER MODELOS COM OS RESIDUOS

plot(tabDadosTodosMuni$riqMedMuniHii ~ tabDadosTodosMuni$areaMuni)
plot(tabDadosTodosMuni$riqMedMuniHii[tabDadosTodosMuni$areaMuni<600] ~ tabDadosTodosMuni$areaMuni[tabDadosTodosMuni$areaMuni<600]) ##FAZER OS MODELOS COM OS RESIDUOS

plot(tabDadosTodosMuni$riscMedMuniHii ~ tabDadosTodosMuni$areaMuni)
plot(tabDadosTodosMuni$riscMedMuniHii[tabDadosTodosMuni$areaMuni<600] ~ tabDadosTodosMuni$areaMuni[tabDadosTodosMuni$areaMuni<600]) ##FAZER OS MODELOS COM OS RESIDUOS

plot(log(tabDadosTodosMuni$HiiMedMuni) ~ tabDadosTodosMuni$areaMuni)
plot(log(tabDadosTodosMuni$HiiMedMuni[tabDadosTodosMuni$areaMuni<600]) ~ tabDadosTodosMuni$areaMuni[tabDadosTodosMuni$areaMuni<600]) ##FAZER OS MODELOS COM OS RESIDUOS

tabDadosTodosMuni = tabDadosTodosMuni[which(tabDadosTodosMuni$areaMuni<600),]


##residuos da area do municipio
residRiq = glm(riqMedMuni ~ areaMuni, family='gaussian', data=tabDadosTodosMuni)$residuals
residRisc = glm(riscMedMuni ~ areaMuni, family='gaussian', data=tabDadosTodosMuni)$residuals
residRiqHii = glm(riqMedMuniHii ~ areaMuni, family='gaussian', data=tabDadosTodosMuni)$residuals
residRiscHii = glm(riqMedMuniHii ~ areaMuni, family='gaussian', data=tabDadosTodosMuni)$residuals
residHii = glm(HiiMedMuni ~ areaMuni, family='gaussian', data=tabDadosTodosMuni)$residuals


##modelos estatisticos

##modelos lineares

##modelo1: riqueza media -> n de casos
modRiq = glm(casos~riqMedMuni,family='binomial',data=tabDadosTodosMuni)
summary(modRiq)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riqMedMuni)
pred = predict(modRiq,data.frame(riqMedMuni=tabDadosTodosMuni$riqMedMuni),type='response')
points(pred ~ tabDadosTodosMuni$riqMedMuni,col='red')

##modelo2: risco medio -> n de casos
modRisc = glm(casos~riscMedMuni,family='binomial',data=tabDadosTodosMuni)
summary(modRisc)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riscMedMuni)
pred = predict(modRisc,data.frame(riscMedMuni=tabDadosTodosMuni$riscMedMuni),type='response')
points(pred ~ tabDadosTodosMuni$riscMedMuni,col='red')

##modelo3: riqueza media HII -> n de casos
modRiqHii = glm(casos ~ riqMedMuniHii,family='binomial',data=tabDadosTodosMuni)
summary(modRiqHii)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riqMedMuniHii)
pred = predict(modRiqHii,data.frame(riqMedMuniHii=tabDadosTodosMuni$riqMedMuniHii),type='response')
points(pred ~ tabDadosTodosMuni$riqMedMuniHii,col='red')

##modelo4: risco medio HII -> n de casos
modRiscHii = glm(casos ~ riscMedMuniHii,family='binomial',data=tabDadosTodosMuni)
summary(modRiscHii)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riscMedMuniHii)
pred = predict(modRiscHii,data.frame(riscMedMuniHii=tabDadosTodosMuni$riscMedMuniHii),type='response')
points(pred ~ tabDadosTodosMuni$riscMedMuniHii,col='red')

##modelo5: HII -> n de casos
modHii = glm(casos ~ HiiMedMuni,family='binomial',data=tabDadosTodosMuni)
summary(modHii)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$HiiMedMuni)
pred = predict(modHii,data.frame(HiiMedMuni=tabDadosTodosMuni$HiiMedMuni),type='response')
points(pred~tabDadosTodosMuni$HiiMedMuni,col='red')

##modelos quadraticos

##modelo1: riqueza media -> n de casos
modRiqQuad = glm(casos ~ riqMedMuni + I(riqMedMuni^2),family='binomial',data=tabDadosTodosMuni)
summary(modRiqQuad)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riqMedMuni)
pred = predict(modRiqQuad,data.frame(riqMedMuni=tabDadosTodosMuni$riqMedMuni),type='response')
points(pred ~ tabDadosTodosMuni$riqMedMuni,col='red')

##modelo2: risco medio -> n de casos
modRiscQuad = glm(casos ~ riscMedMuni + I(riscMedMuni^2),family='binomial',data=tabDadosTodosMuni)
summary(modRiscQuad)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riscMedMuni)
pred = predict(modRiscQuad,data.frame(riqMedMuni=tabDadosTodosMuni$riqMedMuni,riscMedMuni=tabDadosTodosMuni$riscMedMuni),type='response')
points(pred ~ tabDadosTodosMuni$riscMedMuni,col='red')

##modelo3: riqueza media HII -> n de casos
modRiqHiiQuad = glm(casos ~ riqMedMuniHii + I(riqMedMuniHii^2),family='binomial',data=tabDadosTodosMuni)
summary(modRiqHiiQuad)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riqMedMuniHii)
pred = predict(modRiqHiiQuad,data.frame(riqMedMuniHii=tabDadosTodosMuni$riqMedMuniHii),type='response')
points(pred ~ tabDadosTodosMuni$riqMedMuniHii,col='red')

##modelo4: risco medio HII -> n de casos
modRiscHiiQuad = glm(casos ~ riscMedMuniHii + I(riscMedMuniHii^2),family='binomial',data=tabDadosTodosMuni)
summary(modRiscHiiQuad)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$riscMedMuniHii)
pred = predict(modRiscHiiQuad,data.frame(riscMedMuniHii=tabDadosTodosMuni$riscMedMuniHii),type='response')
points(pred ~ tabDadosTodosMuni$riscMedMuniHii,col='red')

##modelo5: HII -> n de casos
modHiiQuad = glm(casos ~ HiiMedMuni + I(HiiMedMuni^2),family='binomial',data=tabDadosTodosMuni)
summary(modHiiQuad)
plot(tabDadosTodosMuni$casos ~ tabDadosTodosMuni$HiiMedMuni)
pred = predict(modHiiQuad,data.frame(HiiMedMuni=tabDadosTodosMuni$HiiMedMuni),type='response')
points(pred~tabDadosTodosMuni$HiiMedMuni,col='red')


##modelos lineares combinados
###obs: riqueza e risco sao correlacionados! nao da pra fazer modelos para os dois juntos!

##modelo1: riqueza media + HII
modRiqHiiComb = glm(casos ~ riqMedMuni + HiiMedMuni ,family='binomial',data=tabDadosTodosMuni)
summary(modRiqHiiComb)

##modelo2: risco medio + HII
modRiscHiiComb = glm(casos ~ riscMedMuni + HiiMedMuni ,family='binomial',data=tabDadosTodosMuni)
summary(modRiscHiiComb)


##modelos quadraticos combinados
###obs: riqueza e risco sao correlacionados! nao da pra fazer modelos para os dois juntos!

##modelo1: riqueza media + HII
modRiqHiiQuadComb = glm(casos ~ riqMedMuni + I(riqMedMuni^2) + HiiMedMuni + I(HiiMedMuni^2) ,family='binomial',data=tabDadosTodosMuni)
summary(modRiqHiiQuadComb)

##modelo2: risco medio + HII
modRiscHiiQuadComb = glm(casos ~ riscMedMuni + I(riscMedMuni^2) + HiiMedMuni + I(HiiMedMuni^2) ,family='binomial',data=tabDadosTodosMuni)
summary(modRiscHiiQuadComb)


##tabela de resultados

tabResMuni = data.frame(
    modelo = c(rep('linear',5),rep('quadratico',5)),
    preditora = c(rep(c('riqueza','risco','riquezaHII','riscoHII','HII'),2)),
    efeito = c(modRiq$coefficients[2],modRisc$coefficients[2],modRiqHii$coefficients[2],modRiscHii$coefficients[2],modHii$coefficients[2],
               modRiqQuad$coefficients[3],modRiscQuad$coefficients[3],modRiqHiiQuad$coefficients[3],modRiscHiiQuad$coefficients[3],modHiiQuad$coefficients[3]),
    deviance = c(modRiq$deviance,modRisc$deviance,modRiqHii$deviance,modRiscHii$devianc,modHii$deviance,
                 modRiqQuad$deviance,modRiscQuad$deviance,modRiqHiiQuad$deviance,modRiscHiiQuad$devianc,modHiiQuad$deviance),
    null_deviance = c(modRiq$null.deviance,modRisc$null.deviance,modRiqHii$null.deviance,modRiscHii$null.deviance,modHii$null.deviance,
                      modRiqQuad$null.deviance,modRiscQuad$null.deviance,modRiqHiiQuad$null.deviance,modRiscHiiQuad$null.deviance,modHiiQuad$null.deviance),
    aic = c(modRiq$aic,modRisc$aic,modRiqHii$aic,modRiscHii$aic,modHii$aic,
            modRiqQuad$aic,modRiscQuad$aic,modRiqHiiQuad$aic,modRiscHiiQuad$aic,modHiiQuad$aic),
    p_valor = c(coef(summary(modRiq))[8],coef(summary(modRisc))[8],coef(summary(modRiqHii))[8],coef(summary(modRiscHii))[8],coef(summary(modHii))[8],
                coef(summary(modRiqQuad))[12],coef(summary(modRiscQuad))[12],coef(summary(modRiqHiiQuad))[12],coef(summary(modRiscHiiQuad))[12],coef(summary(modHiiQuad))[12]))

write.csv(tabResMuni,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabMuniOutputs.csv",row.names=FALSE)


##tabela de resultados modelos combinados

tabResMuniComb = data.frame(
    modelo = c(rep('linear',2),rep('quadratico',2)),
    preditora = c(rep(c('riq&HII','risco&HII'),2)),
    efeitoRiqOuRisc = c(modRiqHiiComb$coefficients[2],modRiscHiiComb$coefficients[2],
                        modRiqHiiQuadComb$coefficients[3],modRiscHiiQuadComb$coefficients[3]),
    efeitoHII = c(modRiqHiiComb$coefficients[2],modRiscHiiComb$coefficients[2],
               modRiqHiiQuadComb$coefficients[5],modRiscHiiQuadComb$coefficients[5]),
    deviance = c(modRiqHiiComb$deviance,modRiscHiiComb$deviance,
                 modRiqHiiQuadComb$deviance,modRiscHiiQuadComb$deviance),
    null_deviance = c(modRiqHiiComb$null.deviance,modRiscHiiComb$null.deviance,
                      modRiqHiiQuadComb$null.deviance,modRiscHiiQuadComb$null.deviance),
    aic = c(modRiqHiiComb$aic,modRiscHiiComb$aic,
            modRiqHiiQuadComb$aic,modRiscHiiQuadComb$aic),
    p_valorRiqOuRisc = c(coef(summary(modRiqHiiComb))[11],coef(summary(modRiscHiiComb))[11],
                         coef(summary(modRiqHiiQuadComb))[18],coef(summary(modRiscHiiQuadComb))[18]),
    p_valorHII = c(coef(summary(modRiqHiiComb))[12],coef(summary(modRiscHiiComb))[12],
                         coef(summary(modRiqHiiQuadComb))[20],coef(summary(modRiscHiiQuadComb))[20]))

write.csv(tabResMuniComb,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabMuniModelosCombinadosOutputs.csv",row.names=FALSE)

##visualizando em ordem de AUC
x = rbind(tabResMuni,tabResMuniComb)
x[order(x$aic),]

##visualizando em ordem do slope
x = rbind(tabResMuni,tabResMuniComb)
x[order(x$efeito),]


###PROJETANDO AREAS DE RISCO EPIDEMICO

extent(mapaHII) = extent(mapaRiqueza)

preditoras = stack(mapaRiqueza,mapaHII)
names(preditoras) = c('riqMedMuni','HiiMedMuni')

preditoras = stack(mapaRisco,mapaHII)
names(preditoras) = c('riscMedMuni','HiiMedMuni')

names(mapaRiqueza) = 'riqMedMuni'

riscoEpi = predict(preditoras,modRiqHiiQuad,type='response')
plot(riscoEpi)
plot(muniShapes,add=T)

jpeg('/home/anderson/Documentos/Projetos/macroecologia_de_chagas/riscoEpi.jpeg',width=1200,height=1200)
plot(riscoEpi)
plot(muniShapes,add=T)
dev.off()


###DADOS DO PAPER DO ARTIGO "Large-scale patterns in morphological diversity, and species assembly in Neotropical Triatominae (Heteroptera: Reduviidae)"
##Link: https://figshare.com/articles/Geographical_patterns_of_Triatominae_Heteroptera_Reduviidae_their_distribution_richness_and_morphology_in_the_Neotropics/653959

tabar = read.csv('/home/anderson/Documentos/Projetos/macroecologia_de_chagas/Pasta sem título/653959/Triatomines_-_Pres__Abs_115_sp__Ricness_2736_Cels_-_Scale_1x1.csv',header=TRUE)
barRiq = tabar[,c('Longitude','Latitude','Richness.Selected.sp')]

library(sp)
library(rgdal)
library(raster)
coordinates(barRiq)=~Longitude+Latitude
proj4string(barRiq)=CRS("+init=epsg:4326") # set it to lat-long
pts = spTransform(barRiq,CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 '))
gridded(barRiq) = TRUE
r = raster(barRiq)
projection(r) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ')

AmSulShape = maptools::readShapePoly("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp")

plot(barRiq)
plot(AmSulShape,add=T)

