#funcoes
source('/home/anderson/Documentos/Projetos/Rio Doce Isac/FD/__online_2017_02_01/functions/multidimFbetaD.R')
source('/home/anderson/Documentos/Projetos/Rio Doce Isac/FD/__online_2017_02_01/functions/quality_funct_space.R')

##dados
matFunc = read.table('/home/anderson/Documentos/Projetos/Rio Doce Isac/Projeto_Rio_Doce/doce_funcional_imputacao.txt') #matriz traits funcionais
matPres = read.table('/home/anderson/Documentos/Projetos/Rio Doce Isac/Projeto_Rio_Doce/doce.comunidades.txt') #matriz pres-aus

##matriz distancia funcional
matDisFunc = quality_funct_space(mat_funct=matFunc,nbdim=3,metric='Gower',dendro=FALSE)$details_funct_space$mat_coord 

##metricas da beta div
indices = match(colnames(as.matrix(matPres)),rownames(matDisFunc))
statsBeta = multidimFbetaD(coord=matDisFunc[indices,], occ=as.matrix(matPres))

#############################################################################################################################
##CONTINUAR: RETIRAR COLUNAS DE MATPRES E LINHAS DE MATFUNC QUANDO ALGUMA COLUNA DE MATPRES FOR IGUAL A ZERO NUMA SIMULACAO.
#############################################################################################################################
