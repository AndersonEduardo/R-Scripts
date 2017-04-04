##funcção para calcular uma distribuição estocastica de valores de AUC a partir de background points##

AUCrand =  function(x,p,path,args){

    sampleSize = nrow(x) #tamanho da amostra de BG
    Nocc = sum(p) #numero de registros de ocorrencias
    AUCvector = numeric()

    for (its in 1:100){
        
        randomSample = sample(sampleSize,Nocc) #amostra aleatoria  de presencas a partir da amostra de BG
        occData = x[randomSample,]
        bgData = x[-randomSample,]
        dataSet = data.frame(pres=c(rep(1,nrow(occData)),rep(0,nrow(bgData))),rbind(occData,bgData))
        
        MX = maxent(x=dataSet[,-1],p=dataSet[,1],args=args)
        
        AUC_i= numeric()
        for(i in 1:length(MX@models)){AUC_i=append(AUC_i,evaluate(p=dataSet[dataSet$pres==1,],a=dataSet[dataSet$pres==0,],model=MX@models[[i]])@auc)}
        AUCvector = append(AUCvector, mean(AUC_i))
        
    }
    return (sort(AUCvector))
}
