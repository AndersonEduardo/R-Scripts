## Analise de sazonalidade usando Fourier ##
## FONTE: https://anomaly.io/detect-seasonality-using-fourier-transform-r/

# Instalar pacote necessario
install.packages("TSA")
library(TSA)

# abrindo conjunto de dados (seus dados reais aqui!)
dataSet = read.csv("pasta1/pasta2/blabla.csv")

# ## inventando dados simples para testes ##
# xs <- seq(-2*pi,2*pi,pi/100)
# dataSet <- cos(1*xs)
# plot(xs,dataSet,type="l"); title("Dados de teste"); abline(h=0,lty=3)
# #########################################

# computar a Fourier Transform (vai aparecer um grafico em que os picos sao os periodos dos potenciais ciclos)
p = periodogram(dataSet)

# ajustando dados obtidos
dd = data.frame(freq=p$freq, spec=p$spec)
order = dd[order(-dd$spec),]
top3 = head(order, 3) #pegando os 3 periodos de ciclo 'mais provaveis'

# so pra visualizar as tres maiores "power frequencies"
top3

# convertendo frequencia em periodos de tempo
time = 1/top3$f
time

# para meus dados inventados, a cada 202 'passos' (aproximadamente) o ciclo se repete. Mas sera possivel? Veja os primeiros 200 passos do dataSet:
plot(dataSet[0:202]) #um ciclo da sazonalidade em meus dados! Funcina!! : )
abline(h=0, col='red', lty=2)
