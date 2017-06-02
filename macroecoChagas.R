##investigando determinantes da doenca de chagas em escala macroecologica##

##pacotes
library(raster)
library(rgdal)

##arquivos raster dos dados
mapaRiqueza = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/Mapas de riqueza/mapaRiquezaPresente.asc')
mapaRisco = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/Mapas de risco/mapaRiscoPresente.asc')
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

##adicionando dados do tamanho territorial dos estados
##link: http://www.ibge.gov.br/home/geociencias/areaterritorial/principal.shtm
tabArea = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/territorio_estados.csv",header=TRUE)
tabRes$area = 

##residuos: retirando o efeito dos tamanhos
residRiq = glm(riqMed~area,family='gaussian',data=tabRes)

##modelo1: riqueza media -> n de casos
modRiqNum = glm(NumCasos~riqMed,family='poisson',data=tabRes)
summary(modRiqNum)
plot(tabRes$NumCasos ~ tabRes$riqMed)
abline(modRiqNum,col='red')

##modelo2: desv. pad. riqueza -> n de casos
modRiscNum = glm(NumCasos~riscMed,family='poisson',data=tabRes)
summary(modRiscNum)
plot(tabRes$NumCasos ~ tabRes$riscMed)
abline(modRiscNum,col='red')

##modelo3: impacto humano medio -> numero de casos
modHiiMeanNum = glm(NumCasos~hiiMed,family='poisson',data=tabRes)
summary(modHiiMeanNum)
plot(tabRes$NumCasos ~ tabRes$hiiMed)
abline(modHiiNum,col='red')

##modelo4: desv. pad. impacto humano -> numero de casos
modHiiSDNum = glm(NumCasos~hiiSD,family='poisson',data=tabRes)
summary(modHiiSDNum)
plot(tabRes$NumCasos ~ tabRes$hiiSD)
abline(modHiiSDNum,col='red')

##modelo5: impacto humano medio -> incidencia
modHiiMeanIncid = glm(Incidencia~hiiMed,family='gaussian',data=tabRes)
summary(modHiiMeanIncid)
plot(tabRes$Incidencia ~ tabRes$hiiMed)
abline(modHiiMeanIncid,col='red')

##modelo6: desv. pad. impacto humano -> incidencia
modHiiSDIncid = glm(Incidencia~hiiSD,family='gaussian',data=tabRes)
summary(modHiiSDIncid)
plot(tabRes$NumCasos ~ tabRes$hiiSD)
abline(modHiiSDIncid,col='red')

##modelo7: impacto humano medio -> numero medio de casos por ano
modHiiMeanMedAno = glm(MediaAno~hiiMed,family='gaussian',data=tabRes)
summary(modHiiMeanIncid)
plot(tabRes$MediaAno ~ tabRes$hiiMed)
abline(modHiiMeanMedAno,col='red')

##modelo8: desv. pad. impacto humano -> numero medio de casos por ano
modHiiSDMedAno = glm(MediaAno~hiiSD,family='gaussian',data=tabRes)
summary(modHiiSDMedAno)
plot(tabRes$MediaAno ~ tabRes$hiiSD)
abline(modHiiSDMedAno,col='red')
