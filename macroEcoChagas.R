##investigando determinantes da doenca de chagas em escala macroecologica##

##pacotes
library(raster)
library(rgdal)

##arquivos raster dos dados
mapaRiqueza = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico/Mapas de riqueza/mapaRiquezaPresente.asc')
mapaRisco = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc')
mapaHII = raster('/home/anderson/PosDoc/dados_ambientais/hii-s-america-geo-grid/res2-5/hii-2-5.asc')
estados = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp',layer='BRUFE250GC_SIR')
vetorEst = estados$NM_ESTADO

##tabela de resultados
tabRes=data.frame()

##extraindo dados para cada estado do Brasil
for (i in 1:length(vetorEst)){

    est_i = estados[estados$NM_ESTADO == vetorEst[i],]
    ##riqueza
    riqMed = mean(unlist(extract(mapaRiqueza,est_i,na.rm=TRUE)),na.rm=TRUE)
    riqSD = sd(unlist(extract(mapaRiqueza,est_i,na.rm=TRUE)),na.rm=TRUE)
    ##risco
    riscMed = mean(unlist(extract(mapaRisco,est_i,na.rm=TRUE)),na.rm=TRUE)
    riscSD = sd(unlist(extract(mapaRisco,est_i,na.rm=TRUE)),na.rm=TRUE)
    ##impacto humano
    hiiMed = mean(unlist(extract(mapaHII,est_i,na.rm=TRUE)),na.rm=TRUE)
    hiiSD = sd(unlist(extract(mapaHII,est_i,na.rm=TRUE)),na.rm=TRUE)

    tabRes = rbind(tabRes,data.frame(estado = vetorEst[i],riqMed=riqMed,riqSD=riqSD,riscMed=riscMed,riscSD=riscSD,hiiMed=hiiMed,hiiSD=hiiSD))
    
}

##salvando a tabela com dados das variaveis extraidas por estado
write.csv(tabRes,file="/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDados.csv",row.names=FALSE)

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

## CONTINUAR DAQUI: ADICIONAR DADOS DE TAMANHO NA TABELA E OBTER OS RESIDUOS ##

##adicionando dados do tamanho territorial dos estados
##link: http://www.ibge.gov.br/home/geociencias/areaterritorial/principal.shtm
tabArea = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/territorio_estados.csv",header=TRUE)
tabArea = tabArea[order(tabArea$estado),] #colocando em ordem alfabetica
tabArea = tabArea[tabArea$estado!='Distrito_Federal',] #retirando o distrito federal
tabRes$area = tabArea$area

##salvando a tabela final
write.csv(tabArea,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabDadosCompleta.csv",row.names=FALSE)

##residuos: retirando o efeito dos tamanhos
residRiq = glm(riqMed~area,family='poisson',data=tabRes)$residuals
residRiqSD = glm(riqSD~area,family='poisson',data=tabRes)$residuals
residHii = glm(hiiMed~area,family='gaussian',data=tabRes)$residuals
residHiiSD = glm(hiiSD~area,family='gaussian',data=tabRes)$residuals
residRisc = glm(riscMed~area,family='gaussian',data=tabRes)$residuals
residRiscSD = glm(riscMed~area,family='gaussian',data=tabRes)$residuals


##variavel rpeditora: riqueza de especies

##modelo1: riqueza media -> n de casos
modRiqNum = glm(tabRes$NumCasos~residRiq,family='poisson')
summary(modRiqNum)
plot(tabRes$NumCasos ~ residRiq)
abline(modRiqNum,col='red')

##modelo2: desv. pad. riqueza -> n de casos
modRiqNumSD = glm(tabRes$NumCasos~residRiqSD,family='poisson')
summary(modRiqNumSD)
plot(tabRes$NumCasos ~ residRiqSD)
abline(modRiqNumSD,col='red')

##modelo3: riqueza media -> incidencia
modRiqIncid = glm(tabRes$Incidencia~residRiq,family='poisson')
summary(modRiqIncid)
plot(tabRes$Incidencia ~ residRiq)
abline(modRiqIncid,col='red')

##modelo4: desv. pad. riqueza -> incidencia
modRiqIncidSD = glm(tabRes$Incidencia~residRiqSD,family='poisson')
summary(modRiqIncidSD)
plot(tabRes$Incidencia ~ residRiqSD)
abline(modRiqIncidSD,col='red')

##modelo5: riqueza media -> numero medio de casos por ano
modRiqMedAno = glm(tabRes$MediaAno~residRiq,family='poisson')
summary(modRiqMedAno)
plot(tabRes$MediaAno ~ residRiq)
abline(modRiqMedAno,col='red')

##modelo6: desv. pad. riqueza -> n de casos
modRiqMedAnoSD = glm(tabRes$MediaAno~residRiqSD,family='poisson')
summary(modRiqIncidSD)
plot(tabRes$Incidencia ~ residRiqSD)
abline(modRiqMedAnoSD,col='red')


##variavel preditora: risco

##modelo1: risco media -> n de casos
modRiscNum = glm(tabRes$NumCasos~residRisc,family='poisson')
summary(modRiscNum)
plot(tabRes$NumCasos ~ residRisc)
abline(modRiscNum,col='red')

##modelo2: desv. pad. risco -> n de casos
modRiscNumSD = glm(tabRes$NumCasos~residRiscSD,family='poisson')
summary(modRiscNumSD)
plot(tabRes$NumCasos ~ residRiscSD)
abline(modRiscNumSD,col='red')

##modelo3: risco media -> incidencia
modRiscIncid = glm(tabRes$Incidencia~residRisc,family='poisson')
summary(modRiscIncid)
plot(tabRes$Incidencia ~ residRisc)
abline(modRiscIncid,col='red')

##modelo4: desv. pad. risco -> incidencia
modRiscIncidSD = glm(tabRes$Incidencia~residRiscSD,family='poisson')
summary(modRiscIncidSD)
plot(tabRes$Incidencia ~ residRiscSD)
abline(modRiscIncidSD,col='red')

##modelo5: risco media -> numero medio de casos por ano
modRiscMedAno = glm(tabRes$MediaAno~residRisc,family='poisson')
summary(modRiscMedAno)
plot(tabRes$MediaAno ~ residRisc)
abline(modRiscMedAno,col='red')

##modelo6: desv. pad. risco -> n de casos
modRiscMedAnoSD = glm(tabRes$MediaAno~residRiscSD,family='poisson')
summary(modRiscIncidSD)
plot(tabRes$Incidencia ~ residRiscSD)
abline(modRiscMedAnoSD,col='red')


##variavel preditora: indice de imapcto humano

##modelo1: impacto humano medio -> numero de casos
modHiiMeanNum = glm(NumCasos~residHii,family='poisson',data=tabRes)
summary(modHiiMeanNum)
plot(tabRes$NumCasos ~ residHii)
abline(modHiiNum,col='red')

##modelo2: desv. pad. impacto humano -> numero de casos
modHiiSDNum = glm(NumCasos~residHiiSD,family='poisson',data=tabRes)
summary(modHiiSDNum)
plot(tabRes$NumCasos ~ residHiiSD)
abline(modHiiSDNum,col='red')

##modelo3: impacto humano medio -> incidencia
modHiiMeanIncid = glm(Incidencia~residHii,family='gaussian',data=tabRes)
summary(modHiiMeanIncid)
plot(tabRes$Incidencia ~ residHii)
abline(modHiiMeanIncid,col='red')

##modelo4: desv. pad. impacto humano -> incidencia
modHiiSDIncid = glm(Incidencia~residHiiSD,family='gaussian',data=tabRes)
summary(modHiiSDIncid)
plot(tabRes$Incidencia ~ residHiiSD)
abline(modHiiSDIncid,col='red')

##modelo5: impacto humano medio -> numero medio de casos por ano
modHiiMeanMedAno = glm(MediaAno~residHii,family='gaussian',data=tabRes)
summary(modHiiMeanIncid)
plot(tabRes$MediaAno ~ residHii)
abline(modHiiMeanMedAno,col='red')

##modelo6: desv. pad. impacto humano -> numero medio de casos por ano
modHiiSDMedAno = glm(MediaAno~residHiiSD,family='gaussian',data=tabRes)
summary(modHiiSDMedAno)
plot(tabRes$MediaAno ~ residHii)
abline(modHiiSDMedAno,col='red')

##tabela de resultados

tabRes = data.frame(
    preditora = c(rep(c('riqueza','risco','HII'),3)),
    resposta = c(rep('Num_de_casos',3),rep('Incidencia',3),rep('Num_med_casos/ano',3)),
    efeito = c(modRiqNum$coefficients[2],modRiscNum$coefficients[2],modHiiMeanNum$coefficients[2],
            modRiqIncid$coefficients[2],modRiscIncid$coefficients[2],modHiiMeanIncid$coefficients[2],
            modRiqMedAno$coefficients[2],modRiscMedAno$coefficients[2],modHiiMeanMedAno$coefficients[2]),
    deviance = c(modRiqNum$deviance,modRiscNum$deviance,modHiiMeanNum$deviance,
               modRiqIncid$deviance,modRiscIncid$deviance,modHiiMeanIncid$deviance,
               modRiqMedAno$coefficients[2],modRiscMedAno$coefficients[2],modHiiMeanMedAno$coefficients[2]),
    aic = c(modRiqNum$aic,modRiscNum$aic,modHiiMeanNum$aic,
            modRiqIncid$aic,modRiscIncid$aic,modHiiMeanIncid$aic,
            modRiqMedAno$aic,modRiscMedAno$aic,modHiiMeanMedAno$aic),
    p_valor = c(coef(summary(modRiqNum))[8],coef(summary(modRiscNum))[8],coef(summary(modHiiMeanNum))[8],
                coef(summary(modRiqIncid))[8],coef(summary(modRiscIncid))[8],coef(summary(modHiiMeanIncid))[8],
                coef(summary(modRiqMedAno))[8],coef(summary(modRiscMedAno))[8],coef(summary(modHiiMeanMedAno))[8]
                )
)

write.csv(tabRes,"/home/anderson/Documentos/Projetos/macroecologia_de_chagas/tabRes.csv",row.names=FALSE)
