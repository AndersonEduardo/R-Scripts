## Analises Testudinidae ##


##PARTE 1 (ambiente R): definindo pasta de trabalho e abrindo dados


setwd("J:/Flauta")
source("J:/Flauta/PyRateMaster/pyrate_utilities.R")
extant_sps = read.csv("J:/Flauta/testudinidae/testudinidae_atuais.csv",header=FALSE,stringsAsFactors = FALSE)
extant_sps = as.vector(extant_sps[,1])
extract.ages.pbdb(file="J:/Flauta/testudinidae/pbdb_data.csv",extant_species=extant_sps,replicates=10)

##esse script deve gerar os seguintes arquivos no diretorio de trabalho:
#pbdb_data_SpeciesList.txt  
#pbdb_data_PyRate.py  


##PARTE 2 (prompt de comandos): implementacao da analise com o PyRate


#no terminal de comando (Windows) verificar se está tudo Ok com os dados de entrada com o comando:
python PyRate.py  J:/Flauta/pbdb_data_PyRate.py   -data_info

#para rodar apenas uma analise simples, desconsiderando as replicas, executar no terminal:
PyRate.py J:/Flauta/testudinidae/pbdb_data_PyRate.py -mHPP -mG -j 1  -n 20000000 -s 5000 -b 0.2 -singleton 1 -thread 2 2

##para rodar para todas as replicas, ha duas formas:
## a) modo sofrido: executar 10 vezes o comando acima variando o termo -j de 1 a 10 (ou seja la qual for o numero de replicas)
## b) modo espertinho: criar um arquivo (no Bloco de Notas) com o seguinte script:
@echo off
cd J:/Flauta/PyRateMaster
for /l %%i in (1,1,10) do python PyRate.py J:/Flauta/testudinidae/pbdb_data_PyRate.py -mHPP -mG -j %%i -n 20000000 -s 5000 -b 0.2 -singleton 1 -thread 2 2
##...salvar no diretorio de trabalho (testudinidae, neste2 caso) e executa-lo, digitando seu nome completo no prompt de comando.

##Para gerar o script e grafico para resultados gerados pelo PyRate:
python PyRate.py -plot J:/Flauta/testudinidae/pyrate_mcmc_logs/ -b 200 -tag pbdb_data


##PARTE 3 (prompt de comandos): analise de correlacao


##primeiramente, transformar dos dados para o formato adequado para rodar essa analise. Executar no prompt de comando:
python PyRateContinuous.py -ginput J:/Flauta/testudinidae/pyrate_mcmc_logs/ -b 200 -tag pbdb_data

##isso criara um arquivo no diretorio de trabalho chamado (neste caso) "pbdb_data_SE_est.txt"

##para rodar a analise de correlacao (default eh exponencial, -m 0 ), executar o seguinte comando no prompt de comando:
python PyRateContinuous.py -d J:/Flauta/testudinidae/pyrate_mcmc_logs/pbdb_data_SE_est.txt -c J:/Flauta/testudinidae/co2Cz.txt -m 0

##obs: o arquivo a ser indicado no arqgumento -c <caminho do arquivo> e um .txt separado por tabulacoes. Ele deve ter 2 colunas, uma para o tempo e a outra para os dados (concentracao de CO2, no presente caso).


##PARTE 4 (ambiente R): interpretar e criar grafico dos coeficientes de correlacao gerados pela analise do pyRate


##assumindo que o R ja esta trabalhando no seu diretorio de trabalho (veja inicio deste script), abrir os resultados da analise de correlacao gerados pelo pyrate:
corPy = read.table('pbdb_data_SE_est_co2_0_exp.log',header=T)

##veja se esta tudo ok com os dados varregados:
head(corPy)

##coeficiente de correlacao da variavel ambiental com a taxa de especiacao (Gl):
CImarginGl = 1.96*sd(corPy$Gl)

jpeg('Gl.jpeg')
my_hist=hist(corPy$Gl,breaks=50,plot=F)
my_color= ifelse(my_hist$breaks > mean(corPy$Gl)+CImarginGl,rgb(0.8,0,0,0.5), ifelse(my_hist$breaks < mean(corPy$Gl)-CImarginGl, rgb(0.8,0,0,0.5), rgb(0.2,0.2,0.2,0.2) ))
plot(my_hist, col=my_color, border=FALSE, main='Coeficiente de correlação Gl', xlab='Gl', ylab='Frequência', xlim=c(-0.004,0.001))
abline(v=0)
dev.off()

##coeficiente de correlacao da variavel ambiental com a taxa de extincao (Gm):
CImarginGm = 1.96*sd(corPy$Gm)

jpeg('Gm.jpeg')
my_hist=hist(corPy$Gm,breaks=50,plot=F)
my_color= ifelse(my_hist$breaks > mean(corPy$Gm)+CImarginGm,rgb(0.8,0,0,0.5), ifelse(my_hist$breaks < mean(corPy$Gm)-CImarginGm, rgb(0.8,0,0,0.5), rgb(0.2,0.2,0.2,0.2) ))
plot(my_hist, col=my_color, border=FALSE, main='Coeficiente de correlação Gm', xlab='Gm', ylab='Frequência', xlim=c(-0.004,0.001))
abline(v=0)
dev.off()

##corelacao linear

corLinPy = read.table('pbdb_data_SE_est_co2_0_linear.log',header=T)
head(corLinPy)

##coeficiente de correlacao da variavel ambiental com a taxa de especiacao (Gl):
CImarginGl = 1.96*sd(corLinPy$Gl)

jpeg('GlCorLin.jpeg')
my_hist=hist(corLinPy$Gl,breaks=50,plot=F)
my_color= ifelse(my_hist$breaks > mean(corLinPy$Gl)+CImarginGl,rgb(0.8,0,0,0.5), ifelse(my_hist$breaks < mean(corLinPy$Gl)-CImarginGl, rgb(0.8,0,0,0.5), rgb(0.2,0.2,0.2,0.2) ))
plot(my_hist, col=my_color, border=FALSE, main='CO2 e taxa de especiação', xlab='Coeficiente de correlação', ylab='Frequência', xlim=c(-0.002,0.001))
abline(v=0)
dev.off()

##coeficiente de correlacao da variavel ambiental com a taxa de extincao (Gm):
CImarginGm = 1.96*sd(corPy$Gm)

jpeg('GmLin.jpeg')
my_hist=hist(corLinPy$Gm,breaks=50,plot=F)
my_color= ifelse(my_hist$breaks > mean(corLinPy$Gm)+CImarginGm,rgb(0.8,0,0,0.5), ifelse(my_hist$breaks < mean(corLinPy$Gm)-CImarginGm, rgb(0.8,0,0,0.5), rgb(0.2,0.2,0.2,0.2) ))
plot(my_hist, col=my_color, border=FALSE, main='CO2 e taxa de extinção', xlab='Coeficiente de extinção', ylab='Frequência', xlim=c(-0.002,0.001))
abline(v=0)
dev.off()


##PARTE 5 (ambiente R): criar grafico de dispersao entre especiacao X covariavel e extincao X covariavel


setwd('/home/anderson/PyRate/Flauta') #diretorio em que estao os arquivos

##abrindo a covariavel
covar = read.table('co2Cz.txt', header=TRUE)

head(covar) #verificando

##dados da taxa de especiacao e extincao ao longo do tempo
##OBS: os valores abaixo (L_mean e M_mean) sao gerados pelo pyrate e estao no arquivo.R gerado por ele para os graficos

L_mean=c(NA, 0.134625878596,0.134369867639,0.1341601135,0.13404832311,0.13395536363,0.133855512751,0.133726362079,0.133667572663,0.133629389942,0.133604413263,0.133594613519,0.133560628152,0.133547054015,0.133540023266,0.133564866614,0.133615898902,0.133644209379,0.133670187214,0.133728496267,0.13380029276,0.133831313181,0.133891875521,0.133914599016,0.133941750746,0.133963193099,0.133941012134,0.133956637921,0.133971723711,0.133997226642,0.134080815723,0.134100418883,0.134185888922,0.134277997608,0.134401709597,0.134533641206,0.134638040599,0.134812275292,0.134837290653,0.134853875876,0.134853977763,0.13488529581,0.134938900369,0.134968653404,0.135019259984,0.135137960374,0.135258101322,0.135433477859,0.135592119253,0.135842408775,0.13624875337,0.136891922478,0.137801740574,0.139495551291,0.141918738752,0.146114689461,0.150563518548,0.151792991109,0.152204187723,0.152668467994,0.152971348173,0.153399083817,0.153655967382)

M_mean=c(NA, 0.880398646278,0.843037987828,0.582481497748,0.151049693127,0.0802738206705,0.074940348817,0.0744441128617,0.0743418241698,0.074279512379,0.0742015840324,0.0741144589784,0.0740054844281,0.0738754124884,0.0737615049263,0.0736967866103,0.0736821306755,0.0737105151112,0.0737600346443,0.0738564386283,0.074079378024,0.0744698052286,0.0750143618948,0.0757560764573,0.0765445874217,0.0774698912068,0.0786422428222,0.0799438344736,0.0813473497762,0.0826857845648,0.0837437659928,0.0847188484718,0.0858493530629,0.0869081134691,0.0874955862646,0.0880480167799,0.0887723005725,0.0899059351354,0.0912921971094,0.092655866673,0.093975473729,0.0951741312513,0.096136489364,0.0969960791467,0.0976298876344,0.0983472099588,0.0990127119268,0.0995575309364,0.100151183004,0.100730743926,0.101559881344,0.102220187375,0.102665516683,0.103065053588,0.103728451351,0.104651981235,0.105843490096,0.107453553464,0.108926297189,0.109474733904,0.109491679093,0.109507230347,0.109513429663)

speciationRate = data.frame(age=rev(age)*-1, lambda=rev(L_mean))
extinctionRate = data.frame(age=rev(age)*-1, mu=rev(M_mean))

##combinando os dados (pareando pela idade taxas (de extincao e especiacao) com o CO2)

pairedDataCorSpeciation = data.frame()
pairedDataCorExtinction = data.frame()

##pareando CO2 com idades de especiacao
index = match(round(speciationRate$age,0), round(covar$time,0))
pairedDataCorSpeciation = rbind(pairedDataCorSpeciation, data.frame(lambda=speciationRate$lambda, co2=covar$co2[index]))

##pareando CO2 com idades de extincao
index = match(round(extinctionRate$age,0), round(covar$time,0))
pairedDataCorExtinction = rbind(pairedDataCorExtinction, data.frame(mu=extinctionRate$mu, co2=covar$co2[index]))

##agora, com os data.frames gerados cria-se os graficos de dispersao

jpeg(filename='especiacaoXco2.jpeg')

plot(pairedDataCorSpeciation$lambda ~ pairedDataCorSpeciation$co2, pch=19, col=rgb(0,0,0,0.5), xlab='CO2', ylab='Taxa de especiação')
reta = glm(pairedDataCorSpeciation$lambda ~ pairedDataCorSpeciation$co2)
abline(reta, col='red'); grid()

dev.off()

jpeg(filename='extincaoXco2.jpeg')

plot(pairedDataCorExtinction$mu ~ pairedDataCorExtinction$co2, pch=19, col=rgb(0,0,0,0.5), xlab='CO2', ylab='Taxa de extinção', xlim=c(0,2000))
reta = lm(pairedDataCorSpeciation$mu ~ pairedDataCorExtinction$co2)
abline(reta, col='red'); grid()

dev.off()



