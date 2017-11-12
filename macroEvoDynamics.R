
cladeDF = data.frame(sp=as.numeric(0), start=as.numeric(0), end=as.numeric(0), origin=as.numeric(-1), status=as.character("leaf"), stringsAsFactors=FALSE)
time = 1:15
theta = 0.1
cladeSizeNow = nrow(cladeDF)

for (t in time){
    for (i in which(cladeDF$status=='leaf')){
        if (runif(1) < theta){
            
            cladeDF = rbind(cladeDF, data.frame(
                                         sp=c(nrow(cladeDF), nrow(cladeDF)+1),
                                         start=t,
                                         end=NA,
                                         origin=cladeDF$sp[i],
                                         status='leaf'))

            cladeDF[i,'status'] = 'internal'
            cladeDF[i,'end'] = t
        }
    }
    cladeDF[is.na(cladeDF$end),'end'] = t
}

cladeDF

##versao 2

cladeDF = data.frame(sp=as.numeric(0), start=as.numeric(0), end=as.numeric(NA), origin=as.numeric(-1), status=as.character("waiting"), stringsAsFactors=FALSE)
theta = 0.2
#cladeMaxSize = 5
tmax = 10
t = 0
itMax=10

for(it in 1:itMax){
    waitingVec = which(cladeDF$status=='waiting') 
    for (i in waitingVec){
#        print(cladeDF)
        while (t<tmax){
            if (runif(1) < theta){
                ##print ("Mutação!!!")
                cladeDF[i,'status'] = 'internal'
                cladeDF[i,'end'] = cladeDF$start[i] + t
                cladeDF = rbind(cladeDF, data.frame(
                                             sp=c(nrow(cladeDF), nrow(cladeDF)+1),
                                             start=cladeDF[i,'end'],
                                             end=NA,
                                             origin=cladeDF[i,'sp'],
                                             status='waiting'))
                break
            }else{
                t = t+1
                                        #                print(paste('TEMPO:', t))
            }
            cladeDF[i,'end'] = cladeDF$start[i] + t
            cladeDF[i,'status'] = 'leaf'
        }
    }
    if (it == itMax){
        index = which(is.na(cladeDF$end))
        cladeDF$end[index] = cladeDF$start[index] + t
        cladeDF$status[index] = 'leaf'
    }
}

cladeDF

library(data.tree)

cladeDF$pathString = paste(-1,
                           cladeDF$origin,
                           cladeDF$sp,sep='/')

cladeTree = as.Node(cladeDF)

print(cladeTree)

plot(cladeTree)

plot(as.dendrogram(cladeTree),center=TRUE)

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
jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_mStatic_u0001.jpg')
plot(output[[1]][,2] ~ output[[1]]$time, type='b', pch=20, ylim=c(0,50), ylab='Tempo', xlab='S')
points(output[[1]][,3]~output[[1]]$time, type='b', pch=20, col='yellow')
points(output[[1]][,4]~output[[1]]$time, type='b', pch=20, col='orange')
points(output[[1]][,5]~output[[1]]$time, type='b', pch=20, col='red')
grid()
legend(x="topleft", legend=c('m=0.01','m=0.003','m=0.002','m=0.0009'), pch=20, lty=1, col=c('black', 'yellow', 'orange', 'red'))
dev.off()

## jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxm__mStatic_u0001.jpg')
## plot(as.numeric(output[[1]][100,2:5]) ~ mVec, type='b', pch=20, ylab=expression(S(t[max])), xlab='m')
## grid()
## dev.off()

#u=0.00025
jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_mStatic_u00025.jpg')
plot(output[[2]][,2] ~ output[[2]]$time, type='b', pch=20, ylim=c(0,50), ylab='Tempo', xlab='S')
points(output[[2]][,3]~output[[2]]$time, type='b', pch=20, col='yellow')
points(output[[2]][,4]~output[[2]]$time, type='b', pch=20, col='orange')
points(output[[2]][,5]~output[[2]]$time, type='b', pch=20, col='red')
grid()
legend(x="topright", legend=c('m=0.01','m=0.003','m=0.002','m=0.0009'), pch=20, lty=1, col=c('black', 'yellow', 'orange', 'red'))
dev.off()

jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxm__mStatic_u0001xu00025.jpg')
plot(log(as.numeric(output[[1]][100,2:5])) ~ mVec, ylim=c(0,20), type='b', pch=20, cex=1, col='black', ylab=expression(log(S(t[max]))), xlab='m')
points(log(as.numeric(output[[2]][100,2:5])) ~ mVec, type='b', pch=20, cex=1, col='red')
grid()
legend(x='topright', legend=c('u=0.0001','u=0.00025'), pch=20, lty=1, col=c('black','red'))
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

jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/tempoDuplicacao_mStatic_u00025.jpg')
plot(output[[1]]$tempo~output[[1]]$m, ylim=c(0,70), xlab='m', ylab='Tempo de duplicação', type='b', pch=20, col='black')
points(output[[2]]$tempo~output[[2]]$m, type='b', pch=20, col='red')
grid()
legend(x='topleft', legend=c('u=0001','u=00025'), pch=20, lty=1, col=c('black','red'))
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
##
jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_m1.jpg')
plot((output[[1]]$S) ~ output[[1]]$t, ylim=c(0,max(output[[3]]$S)), ylab='S', xlab='Tempo', type='b', pch=20, col='black')
points((output[[2]]$S) ~ output[[2]]$t, type='b', pch=20, col='yellow')
points((output[[3]]$S) ~ output[[3]]$t, type='b', pch=20, col='red')
grid()
legend('topleft', legend=c('v = 0.08', 'v = 0.12', 'v = 0.13'), pch=20, lty=1, col=c('black','yellow','orange'))
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
    
jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/tempoDuplicacao_m1.jpg')
plot(output[[1]]$td~output[[1]]$t, type='b', pch=20, col='black', ylab='Tempo de duplicação', xlab='Tempo')
points(output[[2]]$td~output[[2]]$t, type='b', pch=20, col='yellow')
points(output[[3]]$td~output[[3]]$t, type='b', pch=20, col='red')
grid()
legend(x='topleft', legend=c(names(output)), pch=20, lty=1, col=c('black','yellow','red'))
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
        S[i+1] = S[i] + S[i]*( 1 - (exp(-v*(i-T)) - 1) / ( (1 + exp(-v*(i-T)))*(u + 1) - 1 ) )
    }
    outputTab = data.frame(cbind(outputTab,S))
    output[[paste('v = ',v,sep='')]] = outputTab
    outputTab = data.frame(t=c(1:tMax))
}

jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/Sxt_m2.jpg')
plot(log(output[[1]]$S) ~ output[[1]]$t, ylim=c(0,max(log(output[[3]]$S))), ylab='log(S)', xlab='Tempo', type='b', pch=20, col='black')
points(log(output[[2]]$S) ~ output[[2]]$t, type='b', pch=20, col='yellow')
points(log(output[[3]]$S) ~ output[[3]]$t, type='b', pch=20, col='red')
grid()
legend('topleft', legend=c(x=names(output)), pch=20, lty=1, col=c('black','yellow','orange'))
dev.off()

##tempo de diplicação

tMax = 100
vVec = c(0.1, 0.5, 1)
T = 50
td = vector()
outputTab = data.frame(t=c(1:tMax))
output = list()

for(v in vVec){
    for(i in 1:tMax){
    td = append(td,
                log(2) / log( 2 - (exp(-v*(i-T)) - 1) / ( (1 + exp(-v*(i-T)))*(u + 1) - 1 )  )
                )
    }
    outputTab = data.frame(cbind(outputTab,td))
    output[[paste('v',v,sep=' = ')]] = outputTab
    td = vector()
    outputTab = data.frame(t=c(1:tMax))
}
    
jpeg('/home/anderson/Documentos/Projetos/Environmental spatial structure and macroevolutionary dynamics/tempoDuplicacao_m2.jpg')
plot((output[[3]]$td) ~ output[[3]]$t, type='b', pch=20, col='red', ylab='Tempo de duplicação', xlab='Tempo')
points((output[[2]]$td) ~ output[[2]]$t, type='b', pch=20, col='yellow')
points((output[[1]]$td) ~ output[[1]]$t, type='b', pch=20, col='black')
grid()
legend(x='topright', legend=c(rev(names(output))), pch=20, lty=1, col=c('black','yellow','red'))
dev.off()
