## questao 1

l1 = c(0,2,6,10,9)
l2 = c(2,0,5,9,8)
l3 = c(6,5,0,4,5)
l4 = c(10,9,4,0,3)
l5 = c(9,8,5,3,0)

m = matrix(data = c(l1,l2,l3,l4,l5), nrow = 5, ncol = 5, byrow = TRUE )
rownames(m) = 1:5
colnames(m) = rownames(m)

mdist = as.dist(m)

hc1 <- hclust(mdist, method = "single" )

plot(hc1)

## questao 2

df = data.frame(rbind(c(1,1),c(1,2),c(2,1),c(2,2),c(5,1),c(6,1),c(5,2)))
C1 = c(3,0)
C2 = c(5,0)

clusts = kmeans(x=df, centers=rbind(C1,C2))
clusts          

plot(df, col=clusts$cluster, pch=19)

## questao 3

library(fossil)

g1 <- sample(1:2, size=10, replace=TRUE)
g2 <- sample(1:3, size=10, replace=TRUE)
rand.index(g1, g2)

C1 = c(1, 3, 4, 5)
C2 = c(2, 6, 7)

df = data.frame(id=1:7,
                P = c(1,2,1,2,2,1,2),
                C = c(1,2,1,1,1,2,2)
                )
  )
  
  
rand.index(df$P, df$C)


## questao 4

##questao conceitual. MArquei a segunda opcao.

## questao 5

library(caret)

df = read.csv('/home/anderson/Downloads/cinco.csv', h=T)

x = df[,1:4]
y = df[,5]
xtask = data.frame(Aparencia = as.factor("Ensolarado"), 
                   Temperatura = as.factor("Quente"), 
                   Umidade=as.factor("Normal"), 
                   Vento=as.factor("Verdadeiro"))

model = train(x, y, 'nb', trControl=trainControl(method='cv',number=10))

predict(model$finalModel, xtask)


## questao 6-A

df = read.csv("/home/anderson/Downloads/seis.csv", h=T, sep=";")

x = df[,1:4]
y = df[,5]
xtask = data.frame(Sexta.Sabado = as.factor("Sim"),
                   Faminto = as.factor("Sim"),
                   Clientes = as.factor("Cheio"),
                   Tipo = as.factor("Frances"))
                   
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
dtree_fit <- train(Espera ~. , data=df, method = "rpart",
                   parms = list(split = "information"),
                   trControl=trctrl,
                   tuneLength = 10)

predict(dtree_fit, xtask)


##questao 6-B

xtask = data.frame(Sexta.Sabado = as.factor("Sim"),
                   Faminto = as.factor("Sim"),
                   Clientes = as.factor("Alguns"),
                   Tipo = as.factor("Italiano"))

trctrl <- trainControl(method = "repeatedcv", repeats = 3)
grid = expand.grid(k = c(1,3,6,10))

knn_fit <- train(Espera ~., 
                 data = df, 
                 method = "knn",
                 trControl = trctrl,
                 preProcess = c("center", "scale"))

predict(knn_fit, xtask)

## questao 7

df = read.csv("/home/anderson/Downloads/regressao.csv",h=T)

trainctrl <- trainControl(method = "LOOCV")
model <- train(Y ~ ., data = df, method = "glm", trControl = trainctrl) #decision tree



### questao 8

x = df[,1:5]
y = df[,6]
xtask = data.frame(X1 = 245,
                   X2 = 4,
                   X3 = 9700,
                   X4 = 4600,
                   X5 = 1835)



fit <- rpart(Y ~ ., data = df, control = ("maxdepth = 3"))
predict(fit, xtask)


## questao 9

xtask = data.frame(X1 = 245,
                   X2 = 4,
                   X3 = 9700,
                   X4 = 4600,
                   X5 = 1835)

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

knn_fit <- train(Y ~., 
                 data = df, 
                 method = "knn",
                 preProcess = c("center", "scale"),
                 tuneLength = 10)

predict(knn_fit, xtask)


## questao 10

library(arules)

df = read.csv("/home/anderson/Downloads/dez.csv", h=T)
rules <- apriori(df)
inspect(rules)

rules.sub <- rules.sub <- subset(rules, subset = rhs %ain% 
                                   c("Manteiga=S", "Pao=S"))


inspect(rules.sub)


