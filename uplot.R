##codigo para funcao uplot, que plota o resultado da analise de incerteza a partir do output da funcao 'uniche'.
##Anderson A. Eduardo
##12/Set/2018

uplot = function(x, shape=NULL){

    if(class(x) != 'uniche'){
        stop("Os dados de entrada (x) deve ser um objeto da classe 'uniche'.")
        }
        
    if( is.null(shape) ){
        points(x$dataset[,c('lon','lat')], col=rgb(abs(x$uniche$PCC$original),1-abs(x$uniche$PCC$original),0,0.9), cex=1+3*abs(x$uniche$PCC$original), pch=20)
    }else{
        plot(shape)
        points(x$dataset[,c('lon','lat')], col=rgb(abs(x$uniche$PCC$original),1-abs(x$uniche$PCC$original),0,0.9), cex=1+3*abs(x$uniche$PCC$original), pch=20)
    }
}
