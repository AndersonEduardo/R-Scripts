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
##...salvar no diretorio de trabalho (testudinidae, neste caso) e executa-lo, digitando seu nome completo no prompt de comando.

##Para gerar o script e grafico para resultados gerados pelo PyRate:
python PyRate.py -plot J:/Flauta/testudinidae/pyrate_mcmc_logs/ -b 200 -tag pbdb_data


##PARTE 3 (prompt de comandos): analise de correlacao


##primeiramente, transformar dos dados para o formato adequado para rodar essa analise. Executar no prompt de comando:
python PyRateContinuous.py -ginput J:/Flauta/testudinidae/pyrate_mcmc_logs/ -b 200 -tag pbdb_data

##isso criara um arquivo no diretorio de trabalho chamado (neste caso) "pbdb_data_SE_est.txt"

##para rodar a analise de correlacao, executar o seguinte comando no prompt de comando:
python PyRateContinuous.py -d J:/Flauta/testudinidae/pyrate_mcmc_logs/pbdb_data_SE_est.txt -c J:/Flauta/testudinidae/co2Cz.txt

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



