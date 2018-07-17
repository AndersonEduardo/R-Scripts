
##abrindo a funcao##

source('/home/anderson/R-Scripts/sensinc.R')


##criando dados para teste##

set.seed(402)

varx1 = 1:100+rnorm(100,0,10)
varx2 = sample(varx1)
vary = 1.2*varx1 + 0.1*varx2

par(mfrow=c(1,2))
plot(vary ~ varx1)
plot(vary ~ varx2)

dataSet = data.frame(varx1=varx1, varx2=varx2, vary=vary)


##modelo de regressao##


modelo = glm(vary ~ varx1 + varx2, data=dataSet)

summary(modelo)


##analise de incerteza##


#IMPORTANTE: especificacao do erro para as medidas (tem que ser uma lista como a que segue; distribuicoes disponiveis: normal e uniforme)
errorPars =  list(vary = c(distribution='uniform', range=0.001),
                  varx1 = c(distribution='normal', sd=10),
                  varx2 = c(distribution='uniform', range=5)) 

#rodando a funcao da analise (os nomes das variaveis precisam estar na mesma ordem do dataset original)
teste = sensinc(dados=dataSet, varx_names=c("varx1","varx2"), vary_name="vary", errorPars=errorPars, modelo=modelo) 

print(teste) #inspecionando os resultados


##inspecao grafica dos resultados (as barras representam o intervalo de confianca 95%)##


#so pra verificar os nomes
names(teste) 

##efeito do 'erro' nas medidas sobre o intercepto do modelo 
plot(teste$"(Intercept)")
abline(h=0, lty=2)
title(ylab="Intercept sensitivity")

##efeito do 'erro' nas medidas sobre o efeito da variavel x1 no modelo
plot(teste$"varx1")
abline(h=0, lty=2)
title(ylab="Sensitivity of varx1 effect")

##efeito do 'erro' nas medidas sobre o efeito da variavel x2 no modelo
plot(teste$"varx2")
abline(h=0, lty=2)
title(ylab="Sensitivity of varx2 effect")
