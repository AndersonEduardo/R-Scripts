## Script para realiacao de PNO (predicted niche occupancy)
## 06/junho/2016

##install.packages('phyloclim')
library(phyloclim)

##parametros
path_bioclim = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas/bio7.asc"
path_model = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/resultados nicho climatico/Projecoes" #cuidado, tem mtos arquivos aqui

##usando a funcao
pnoTest = pno(path_bioclim = path_bioclim, path_model = path_model)

##inspecao grafica
plot(pnoTest[,1],pnoTest[,2],type='l')
points(pnoTest[,1],pnoTest[,3],type='l', col='blue')
