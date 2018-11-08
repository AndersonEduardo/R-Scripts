## Extraindo informacao de texto (em arquivos salvos na memoria fisica do PC)
## Anderson A. Eduardo
## 08/novembro/2018

library(stringr)

## extrai linha com a informacao de interesse do arquivo
x = grep('L_mean', readLines('/home/anderson/Projetos/Tubarões - Elisa Cravo/neoselachii/analiseOrdens/ordensDadosPBDB/pyrate_mcmc_logs/carcharhiniformes_extinct_1_G_marginal_rates_RTT.r'), value=T)

## usando o gsub para 'limpar' a linha
x = gsub('L_mean=c\\(','',x[1])
x = gsub("\\)","",x)

## quebrando a linha de texto em um vetor
x2 = strsplit(x, ',')
x2 = unlist(x2)

## convertendo em valores numericos
x2 = as.numeric(x2)
