library(data.tree)

##versao 2

cladeDF = data.frame(sp=as.numeric(0), start=as.numeric(0), end=as.numeric(NA), origin=as.numeric(-1), status=as.character("waiting"), pathString=-1, stringsAsFactors=FALSE)
#cladeMaxSize = 5
tmax = 100
t = 0
itMax = 10
#theta = 0.1


v = 0.8
T = tmax/2
u = 0.001

theta =  0.5 - 0.5/( (1+exp(-v*(1:tmax-T)))*u + 1 ) 
plot(theta,ylim=c(0,1))

for(it in 1:itMax){
    waitingVec = which(cladeDF$status=='waiting') 
    for (i in waitingVec){
        while (t<tmax){
            if (runif(1) < theta[t+1]){
                cladeDF[i,'status'] = 'internal'
                cladeDF[i,'end'] = t
                cladeDF = rbind(cladeDF, data.frame(
                                             sp=c(nrow(cladeDF), nrow(cladeDF)+1),
                                             start=cladeDF[i,'end'],
                                             end=NA,
                                             origin=cladeDF[i,'sp'],
                                             status='waiting',
                                             pathString=c(paste(cladeDF[i,'pathString'],cladeDF[i,'sp'],sep='/'),
                                                          paste(cladeDF[i,'pathString'],cladeDF[i,'sp'],sep='/'))
                                         )
                                )
                break
            }else{
                t = t+1
            }
            cladeDF[i,'end'] = t
            cladeDF[i,'status'] = 'leaf'
#            cladeDF[i,'pathString']=paste(cladeDF[i,'pathString'],cladeDF[i,'sp'],sep='/')
        }
        cladeDF[i,'pathString']=paste(cladeDF[i,'pathString'],cladeDF[i,'sp'],sep='/')
    }
    if (it == itMax){
        index = which(is.na(cladeDF$end))
        cladeDF$end[index] = t
        cladeDF$status[index] = 'leaf'
        cladeDF$pathString[index] = paste(cladeDF[index,'pathString'],cladeDF[index,'sp'],sep='/')
    }
}

cladeDFbruto = cladeDF
##
cladeDF = cladeDF[cladeDF$status=='leaf',]
##

## cladeDF$pathString = paste(-1,
##                            cladeDF$origin,
##                            cladeDF$sp,sep='/')

cladeTree = as.Node(cladeDF)

#print(cladeTree)

plot(cladeTree)


jpeg('diminuindo.jpeg')
plot(as.dendrogram(cladeTree),center=TRUE)
dev.off()

CladeNewick = ToNewick(cladeTree)



################################################
#### ANALISE NUMERICA DO MODELO EXPONENCIAL ####
################################################



S0 = 1
S = vector(); S[1] = S0
tMax = 100
uVec = c(0.0001, 0.00025)
mVec = c(0.01,0.003, 0.002, 0.0009)
outputTab = data.frame(time = 1:tMax)
output = list()

for (u in uVec){
    for (m in mVec){
        for (i in 1:(tMax-1)){
            S[i+1] = S[i]*(1 + u/(u+m))
        }
        outputTab = data.frame(cbind(outputTab,S))
    }
    output[[paste('u',u, sep='')]] = outputTab
    outputTab = data.frame(time = 1:tMax)
}

#u=0.0001
postscript(file='/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_mStatic_u0001.eps', onefile=FALSE, horizontal=FALSE)
par(mar=c(5,6,1,1))
plot(output[[1]][,2] ~ output[[1]]$time, type='b', pch=20, cex=3, ylim=c(0,50), ylab='S', xlab='Tempo', cex.lab=1.9, cex.axis=1.5)
points(output[[1]][,3]~output[[1]]$time, type='b', pch=20, cex=3, col='yellow')
points(output[[1]][,4]~output[[1]]$time, type='b', pch=20, cex=3, col='orange')
points(output[[1]][,5]~output[[1]]$time, type='b', pch=20, cex=3, col='red')
grid()
legend(x="topleft", legend=c('m=0.01','m=0.003','m=0.002','m=0.0009'), pch=20, cex=2, lty=1, col=c('black', 'yellow', 'orange', 'red'))
dev.off()

## jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxm__mStatic_u0001.jpg')
## plot(as.numeric(output[[1]][100,2:5]) ~ mVec, type='b', pch=20, ylab=expression(S(t[max])), xlab='m')
## grid()
## dev.off()

#u=0.00025
postscript(file='/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_mStatic_u00025.eps', onefile=FALSE, horizontal=FALSE)
par(mar=c(5,6,1,1))
plot(output[[2]][,2] ~ output[[2]]$time, type='b', pch=20, cex=3, ylim=c(0,50), ylab='S', xlab='Tempo', cex.lab=1.9, cex.axis=1.5)
points(output[[2]][,3]~output[[2]]$time, type='b', pch=20, cex=3, col='yellow')
points(output[[2]][,4]~output[[2]]$time, type='b', pch=20, cex=3, col='orange')
points(output[[2]][,5]~output[[2]]$time, type='b', pch=20, cex=3, col='red')
grid()
legend(x="topright", legend=c('m=0.01','m=0.003','m=0.002','m=0.0009'), pch=20, cex=2, lty=1, col=c('black', 'yellow', 'orange', 'red'))
dev.off()

postscript('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxm__mStatic_u0001xu00025.eps', onefile=FALSE, horizontal=FALSE)
par(mar=c(5,6,1,1))
plot(log(as.numeric(output[[1]][100,2:5])) ~ mVec, ylim=c(0,20), type='b', pch=20, cex=3, col='black', ylab=expression(log(S(t[max]))), xlab='m', cex.lab=1.9, cex.axis=1.5)
points(log(as.numeric(output[[2]][100,2:5])) ~ mVec, type='b', pch=20, cex=3, col='red')
grid()
legend(x='topright', legend=c('u=0.0001','u=0.00025'), pch=20, cex=2, lty=1, col=c('black','red'))
grid()
dev.off()

##tempode de duplicacao

tempo = vector()
mVec = c(1:10)/1000
uVec = c(0.0001, 0.00025)
outputTab = data.frame(m = mVec)
output = list()

for(u in uVec){
    for(m in mVec){
        tempo = append(tempo,
                       log(2)/log(1 + (u/(u+m)))
                       )
    }
    outputTab = data.frame(cbind(outputTab,tempo))
    output[[paste('u',u,sep='')]] = outputTab
    outputTab = data.frame(m = mVec)
    tempo = vector()
}

postscript('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/tempoDuplicacao_mStatic_u00025.eps', onefile=FALSE, horizontal=FALSE)
plot(output[[1]]$tempo~output[[1]]$m, ylim=c(0,70), xlab='m', ylab='Tempo de duplicação', type='b', pch=20, cex=3, col='black', cex.lab=1.9, cex.axis=1.5)
points(output[[2]]$tempo~output[[2]]$m, type='b', pch=20, cex=3, col='red')
grid()
legend(x='topleft', legend=c('u=0001','u=00025'), pch=20, cex=3, lty=1, col=c('black','red'))
dev.off()


######################
##VARIACAO TEMPORAL###
######################


## m' ##

tMax = 100
vVec = c(0.08, 0.12, 0.13)
T = 50
#m = 1/(1 + exp(-v*(t-T)))
u = 0.001
S0 = 1
S = vector(); S[1] = S0
outputTab = data.frame(t=c(1:tMax))
output = list()

for(v in vVec){
    for(i in 1:(tMax-1)){
#        S[i+1] = S[i] + S[i]*( (1+exp(-v*(i-T)))*u )/(  (1+exp(-v*(i-T)))*u + 1 ) 
        S[i+1] = S[i] + S[i]*( 1 - 1/( (1+exp(-v*(i-T)))*u + 1 ) ) 
    }
    outputTab = data.frame(cbind(outputTab,S))
    output[[paste('v',v,sep='')]] = outputTab
    outputTab = data.frame(t=c(1:tMax))
}

postscript('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_m1.eps', onefile=FALSE, horizontal=FALSE)
plot((output[[1]]$S) ~ output[[1]]$t, ylim=c(0,max(output[[3]]$S)), ylab='S', xlab='Tempo', type='b', pch=20, cex=3, col='black', cex.lab=1.9, cex.axis=1.5)
points((output[[2]]$S) ~ output[[2]]$t, type='b', pch=20, cex=3, col='yellow')
points((output[[3]]$S) ~ output[[3]]$t, type='b', pch=20, cex=3, col='red')
grid()
legend('topleft', legend=c('v = 0.08', 'v = 0.12', 'v = 0.13'), pch=20, cex=2, lty=1, col=c('black','yellow','orange'))
dev.off()

##tempo de diplicação

tMax = 100
vVec = c(0.1, 0.2, 1)
T = 50
td = vector()
outputTab = data.frame(t=c(1:tMax))
output = list()

for(v in vVec){
    for(i in 1:tMax){
    td = append(td,
                log(2) / log( ( (1 + exp(-v*(i-T)))*2*u + 1 ) / ((1 + exp(-v*(i-T)))*u + 1) )
                )
    }
    outputTab = data.frame(cbind(outputTab,td))
    output[[paste('v',v,sep=' = ')]] = outputTab
    td = vector()
    outputTab = data.frame(t=c(1:tMax))
}
    
postscript('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/tempoDuplicacao_m1.eps', onefile=FALSE, horizontal=FALSE)
plot(output[[1]]$td~output[[1]]$t, type='b', pch=20, cex=3, col='black', ylab='Tempo de duplicação', xlab='Tempo', cex.lab=1.9, cex.axis=1.5)
points(output[[2]]$td~output[[2]]$t, type='b', pch=20 cex=3,, col='yellow')
points(output[[3]]$td~output[[3]]$t, type='b', pch=20, cex=3, col='red')
grid()
legend(x='topleft', legend=c(names(output)), pch=20, cex=2, lty=1, col=c('black','yellow','red'))
dev.off()


## m'' ##


tMax = 100
vVec = c(0.1, 0.5, 1)
T = 50
#m = 1 - 1/(1 + exp(-v*(t-T)))
u = 0.001
S0 = 1
S = vector(); S[1] = S0
outputTab = data.frame(t=c(1:tMax))
output = list()

for(v in vVec){
    for(i in 1:(tMax-1)){
        S[i+1] = S[i] + S[i]*( (1+exp(-v*(i-T)))*u ) / ( (1+exp(-v*(i-T)))*u + exp(-v*(i-T)) )
    }
    outputTab = data.frame(cbind(outputTab,S))
    output[[paste('v = ',v,sep='')]] = outputTab
    outputTab = data.frame(t=c(1:tMax))
}

postscript(file='/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_m2.eps', onefile=FALSE, horizontal=FALSE)
plot(log(output[[1]]$S) ~ output[[1]]$t, ylim=c(0,max(log(output[[3]]$S))), ylab='log(S)', xlab='Tempo', type='b', pch=20, cex=3, col='black', cex.lab=1.9, cex.axis=1.5)
points(log(output[[2]]$S) ~ output[[2]]$t, type='b', pch=20, cex=3, col='yellow')
points(log(output[[3]]$S) ~ output[[3]]$t, type='b', pch=20, cex=3, col='red')
grid()
legend('topleft', legend=c(x=names(output)), pch=20, cex=2, lty=1, col=c('black','yellow','orange'))
dev.off()

##tempo de diplicação

tMax = 100
vVec = c(0.1, 0.2, 1)
T = 50
td = vector()
outputTab = data.frame(t=c(1:tMax))
output = list()
#
for(v in vVec){
    for(i in 1:tMax){
    td = append(td,
                log(2) / ( log( ( (1+exp(-v*(i-T)))*2*u + exp(-v*(i-T)) ) / ( (1 + exp(-v*(i-T)))*u + exp(-v*(i-T)) ) ) )
                )
    }
    outputTab = data.frame(cbind(outputTab,td))
    output[[paste('v',v,sep=' = ')]] = outputTab
    td = vector()
    outputTab = data.frame(t=c(1:tMax))
}
 #   
postscript(file='/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/tempoDuplicacao_m2.eps', onefile=FALSE, horizontal=FALSE, width = 480, height = 480)
plot((output[[3]]$td) ~ output[[3]]$t, type='b', pch=20, cex=3, col='red', ylab='Tempo de duplicação', xlab='Tempo', cex.lab=1.9, cex.axis=1.5)
points((output[[2]]$td) ~ output[[2]]$t, type='b', pch=20, cex=3, col='yellow')
points((output[[1]]$td) ~ output[[1]]$t, type='b', pch=20, cex=3, col='black')
grid()
legend(x='topright', legend=c(rev(names(output))), pch=20, cex=2, lty=1, col=c('red','yellow','black'))
dev.off()




##############################




createPop = function(dados, m){
    leftPop = names(dados[1,])[1] ##nao importa a linha
    rightPop = names(dados[1,])[ncol(dados)] ##nao importa a linha

    ##para esquerda
    if (runif(1) < m & sum(dados[,1]) > 0){
        popNew = paste('pop', as.numeric(gsub("pop",'',names(dados[1])))-1,sep='')
        dados[,popNew] = 0
        dados = dados[,c(popNew,names(dados)[-ncol(dados)])]
    }

    ##para direita
    if (runif(1) < m & sum(dados[,ncol(dados)]) > 0){
        popNew = paste('pop+', as.numeric(gsub("pop",'',names(dados)[ncol(dados)]))+1,sep='')
        dados[,popNew] = 0
    }
    return(dados)
}

positionPop = function(h,i,dados){
    if (ncol(dados)==1){
        pRight = 0
        pActual = as.numeric(dados[h,i])
        pLeft = 0
        return(c(pLeft,pActual,pRight))
    }
    if (i==1 & ncol(dados)>1){
        pRight = as.numeric(dados[h,i+1])
        pActual = as.numeric(dados[h,i])
        pLeft = 0
        return(c(pLeft,pActual,pRight))
    }
    if (i == ncol(dados) & ncol(dados)>1){
        pRight = 0
        pActual = as.numeric(dados[h,i])
        pLeft = as.numeric(dados[h,i-1])
        return(c(pLeft,pActual,pRight))
    }else{
        pLeft = as.numeric(dados[h,i-1])
        pActual = as.numeric(dados[h,i])
        pRight = as.numeric(dados[h,i+1])
        return(c(pLeft,pActual,pRight))
    }
}

mutationPop = function(dados,u){
    if(runif(1) < u){
        dados[nrow(dados)+1, sample(1:ncol(dados),1)] = 0.01
        dados[nrow(dados),is.na(dados[nrow(dados),])] = 0
    }
    return(dados)
}


dados = data.frame('pop0'=10000)
tMax = 10
m = 0.1
output = vector()
u = 0.001

for (time in 1:tMax){
    dados = createPop(dados,m)
    for (h in 1:nrow(dados)){
        for (i in 1:ncol(dados)){
            ##posicao atual e vizinhanca
            positions = positionPop(h,i,dados)
            pLeft = positions[1]
            pActual = positions[2]
            pRight = positions[3]
            ## equacao ##
            p = pActual - m*pActual +(m/2)*pLeft + (m/2)*pRight
            ##atualiza vetor
            dados[h,i] = p
            ##outputs
        }
    }
    dados = mutationPop(dados,u)
}


f = t(apply(dados, 1, function(x) x/colSums(dados)))

##verificando se soma 1
colSums(f)
plot(f[1,],t='b')
points(f[2,],t='b',col='red')
points(f[3,],t='b',col='blue')



################################################

nVec = vector()
nTemp = vector()

for(i in 1:100){

    m = i/100 #0.5
    u = 0.1
    tempo = 1:1000
    n = 2 + tempo*m
    
    rho = ((m)/(n-1)) / (u + m/(n-1)) #qunato maior, mais parecidas as pops
    ##
    plot(rho ~ tempo, ylim=c(0,1))
    
    nVec = append(nVec,  n[which(rho<0.1)][1]) #momento que atinge o limite escolhido
    nTemp = append(nTemp, tempo[which(rho<0.1)][1]) #momento que atinge o limite escolhido
}

plot(nVec ~ c(1:100/100), t='b')
plot(nTemp ~ c(1:100/100), t='b')

plot(nVec ~ nTemp, t='b')


rho=0.01 #1:100/1000
m=1:500/1000
n=(m-m*rho+rho*u)/(rho*u)
##
plot(n ~ m)

m=1:500/10000
rho=0.1
u=0.001
t=-(-m+m*rho+rho*u)/(m*rho*u)
plot(t~m,t='b')

n=2+m*t
##
plot(n~t,t='b')
