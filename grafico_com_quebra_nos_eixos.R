## Script para fazer grafico customizado, com quebra no eixo X e/ou Y
## 11/jun/2018
## Anderson A. Eduardo


##install.packages('plotrix')  ##instalar esse pacote
library(plotrix)

##inventando os dados para teste
x = 1:100 
y1 =  x[1:50] + rnorm(50, mean = 10, sd = 2) 
y2 =  x[51:100] + rnorm(50, mean = 100, sd = 2)
dataSet = data.frame(x=x, y=c(y1,y2))


##plotrix
gap.plot(dataSet$x,dataSet$y, gap=c(65,143), gap.axis="y", pch=16, col="blue", bty='n', ylab='Nome do eixo Y', xlab='Nome do eixo X') #plotando o grafico
abline(h=seq(64,68,0.2), pch= 20, col="white", lwd=2)  #escondendo as linhas entre as caixas (aqui tem que ir testando os numeros até ficar certo)
box() #corrigindo a linha de entorno do grafico
axis.break(axis=2, breakpos=65, brw=0.04) #barra da quebra do eixo
