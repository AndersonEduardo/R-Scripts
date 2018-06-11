##Script com tutorial para implementar modelos estatisticos usando GLM com proporcoes na variavel reposta
##06/juunho/2018


## criando dados para exemplo ##

##mortes
mortes_totais = runif(n=100, min=0, max=100)
mortes_totais = round(mortes_totais)

##populacao
pop_total = runif(n=100, min=100:200, max=200)
pop_total = round(pop_total)

##nomes das localidades
myFun <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

localidades = myFun(100)

##variaveis ambietais
var1 = rnorm(n=100, mean=pop_total, sd=5)
var2 = 1/(rnorm(n=100, mean=pop_total, sd=5))*100


##consolidando dataset
dataSet = data.frame(local=localidades, mortes=mortes_totais, populacao=pop_total, var1=var1, var2=var2, stringsAsFactors=FALSE)


## GLM das proporcoes de mortes ##

##veja como devem estar as colunas do dataset
View(dataSet)

##modelo linear
modelLinear = glm( cbind(mortes, populacao) ~ var1 + var2, family=binomial, data=dataSet ) 
summary(modelLinear) #outputs

##modelo quadratico
modelQuad = glm( cbind(mortes, populacao) ~ var1 + I(var1^2) + var2 + I(var2^2), family=binomial, data=dataSet ) 
summary(modelQuad) #outputs

##caso queira comparar modelos fazer:
AIC(modelLinear, modelQuad) #metodo 1
anova(modelLinear, modelQuad, test="Chisq") #metodo 2 (pra itnerpretar, olhe o p-valor para a diferenca entre os modelos)


## Para fazer graficos ##

##valores criticos para gerar as projecoes dos modelos
var1min = round(min(dataSet$var1))
var1max = round(max(dataSet$var1))
var2min = round(min(dataSet$var2))
var2max = round(max(dataSet$var2))


##plotando dados observados
plot(x = dataSet$var1,
     y = c(dataSet$mortes/dataSet$populacao),
     pch = 20, xlab = 'Variavel 1', ylab = 'Mortes') ##repetir o grafico para visualizar com outra variavel
##plotando projecoes do modelo 1
points(x = seq(var1min, var1max, length=100),
       y = predict(modelLinear, newdata=data.frame(var1=seq(var1min, var1max, length=100), var2=seq(var2min, var2max, length=100)), type='response'),
       type='l', col='blue', ylim=c(0,1))
##plotando projecoes do modelo 2
points(x = seq(var1min, var1max, length=100),
       y = predict(modelQuad, newdata=data.frame(var1=seq(var1min, var1max, length=100), var2=seq(var2min, var2max, length=100)), type='response'),
       type='l', col='red', ylim=c(0,1))
##plotando legenda
legend('topright', legend=c('Linear Model','Quadratic Model'), lty=1, col=c('blue','red'))


