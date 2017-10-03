library(raster)
library(rgdal)
library(picante)

projectFolder = "/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico/"
municipios = readOGR(dsn='/home/anderson/PosDoc/shapefiles/br_municipios/BRMUE250GC_SIR.shp',layer='BRMUE250GC_SIR')

##municipios com casos
tabInfecbMuni = read.csv('/home/anderson/Documentos/Projetos/macroecologia_de_chagas/infeccao-municipio-2007-2014.csv',header=TRUE)
muniChagas = municipios[match(toupper(tabInfecbMuni$municipio),municipios$NM_MUNICIP),]

#distribuicoes
listaDist = grep(list.files(paste(projectFolder,"Projecoes",sep=""),full.names=TRUE),pattern='Otimista|Pessimista',inv=T,value=T)
listaDistBIN = grep(listaPresente,pattern='BIN.asc',value=TRUE)
camadasDist = stack(listaPresenteBIN)

extractedData = sapply(extract(camadasDist,municipios[municipios$NM_MUNICIP,],na.rm=TRUE), colSums) > 0
extractedData = data.frame(extractedData)
colnames(extractedData) = paste(municipios$NM_MUNICIP[1:10])
extractedData = t(extractedData)


#distancia Sørensen
distancia = vegdist(extractedData,binary=TRUE)

#classificacao (hierarchical clustering)
cidadesCluster <- hclust(distancia, method = "average")

plot(cidadesCluster)

rect.hclust(cidadesCluster, 2)

grp <- cutree(cidadesCluster, 2)

ord <- metaMDS(extractedData, dist = "bray")  #cca(x3)
plot(ord, display = "sites")
ordihull(ord, grp, lty = 2, col = "red")

##ordenacao

# The metaMDS function automatically transforms data and checks solution
# robustness
comm.bc.mds <- metaMDS(extractedData, dist = "bray")

# Assess goodness of ordination fit (stress plot)
stressplot(comm.bc.mds)

##We can plot the ordination results in a variety of different ways:

# plot site scores as text
ordiplot(comm.bc.mds, display = "sites", type = "text")

# plot site scores as text
ordiplot(comm.bc.mds, display = "sites", type = "text")

# automated plotting of results - tries to eliminate overlapping labels
ordipointlabel(comm.bc.mds)

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig <- ordiplot(comm.bc.mds, type = "none")
# plot just the samples, colour by habitat, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = metadata$habitat == 
    "Fescue")
points(mds.fig, "sites", pch = 19, col = "blue", select = metadata$habitat == 
    "Mixedgrass")
# add confidence ellipses around habitat types
ordiellipse(comm.bc.mds, metadata$habitat, conf = 0.95, label = TRUE)
# overlay the cluster results we calculated earlier
ordicluster(comm.bc.mds, comm.bc.clust, col = "gray")
