
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
