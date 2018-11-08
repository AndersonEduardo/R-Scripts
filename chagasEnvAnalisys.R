##investigando determinantes ambientais da ocorrencia de doenca de chagas nos municipios brasileiros

##pacotes
library(raster)
library(rgdal)

##arquivos raster dos dados
mapaRiqueza = raster('/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/SDM outputs/SDMclimAcumSuit.asc')
mapaRisco = raster('/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc')
mapaHII = raster('/home/anderson/gridfiles/hii-s-america-geo-grid/res2-5/hii-2-5.asc')
mapaUsoSolo = raster(...)
mapaSocioeco = raster(...)
tabInfecbMuni = read.csv('/home/anderson/Downloads/A165041200_17_141_3_DADOS_LIMPOS.csv',header=TRUE)
names(tabInfecbMuni) = c('municipio','casos')
municipios = readOGR(dsn='/home/anderson/shapefiles/br_municipios/BRMUE250GC_SIR.shp',layer='BRMUE250GC_SIR')


##extraindo para TODOS os municipios
riqMedMuni = sapply(extract(mapaRiqueza, municipios,na.rm=TRUE), mean, na.rm=TRUE)
riscMedMuni = sapply(extract(mapaRisco, municipios,na.rm=TRUE), mean,na.rm=TRUE)
HiiMedMuni = sapply(extract(mapaHII, municipios,na.rm=TRUE), mean,na.rm=TRUE)
areaMuni = sapply(extract(mapaHII, municipios,na.rm=TRUE), length)


### TESTANDO ESPACIALIZACAO DOS DADOS DE OCORRENCIA ###
### Duvida: a correlacao espacial das ocorrencias consegue explicar sozinha os casos de Doenca de Chagas?
estados = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp',layer='BRUFE250GC_SIR')
HiiMedEst = sapply(extract(mapaHII, estados, na.rm=TRUE), mean, na.rm=TRUE) ##usando HII so pra testar (depois: aqui casos/municipio)
cols = (HiiMedEst - min(HiiMedEst))/(min(HiiMedEst) + max(HiiMedEst))
plot(estados, col=rgb(cols,1-cols,1-cols,0.5))

xx = raster(estados, nrows=dim(mapaHII)[1], ncols=dim(mapaHII)[2])
yy = rasterize(estados, xx, field=HiiMedEst)

### FIM ###


##tabela de dados
tabDadosTodosMuni = data.frame(municipios = municipios$NM_MUNICIP,
                               riqMedMuni=riqMedMuni,
                               riscMedMuni=riscMedMuni,
                               riqMedMuniHii=riqMedMuniHii,
                               riscMedMuniHii=riscMedMuniHii,
                               HiiMedMuni=HiiMedMuni,
                               areaMuni=areaMuni)


##adicionando a tabela de dados de municipios
vetorCasos = as.integer( tabDadosTodosMuni$municipios %in% toupper(tabInfecbMuni$municipio) ) #vetor de 0 e 1
vetorLinhas = base::match( toupper(tabInfecbMuni$municipio),tabDadosTodosMuni$municipios ) #encontrando as linhas que tem dados do SUS
vetorCasos[vetorLinhas] = tabInfecbMuni$casos #ajustando para '0' e 'numero de casos',  em vez de '0' e '1'
tabDadosTodosMuni$casos = vetorCasos #adicionando na tabela de dados

##salvando a tabela final
write.csv(tabDadosTodosMuni,"/home/anderson/Projetos/Chagas environmental analysis/tabDadosTodosMuniCompleta.csv",row.names=FALSE)

##abrindo tabela dedados completa
tabDadosTodosMuni = read.csv("/home/anderson/Projetos/Chagas environmental analysis/tabDadosTodosMuniCompleta.csv",header=TRUE)
