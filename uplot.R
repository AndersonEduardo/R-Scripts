##codigo para funcao uplot, que plota o resultado da analise de incerteza a partir do output da funcao 'uniche'.
##Anderson A. Eduardo
##12/Set/2018

uplot = function(x, shape=NULL, niche.metric=NULL, legend=TRUE){

    if (class(x) != 'uniche'){
        stop("Os dados de entrada (x) deve ser um objeto da classe 'uniche'.")
    }

    if (is.null(shape) & is.null(niche.metric)){
        stop("Quando um shapefile não é fornecido o argumento 'niche.metric' deve ser especificado (as opções são: 'marginality' e 'volume').")
    }

    if (length(niche.metric) > 1){
        stop("Quando um shapefile não é fornecido o argumento 'niche.metric' deve ser especificado com uma das opções: 'marginility' e 'volume'.")
    }
    
    if (any(!niche.metric %in% c('marginality','volume'))){
        stop("O argumento de entrada 'niche.metric' deve ser 'marginality' ou 'volume'.")
    }

    
    if( is.null(shape) ){
        if (niche.metric == 'marginality'){

            sensitivity = rgb(abs(x$uniche.marginality$PCC$original),1-abs(x$uniche.marginality$PCC$original),0,0.9)
            tempUncert = 2.5 + 2 * apply( xx$uniche.marginality$X, 2, sd) / max(apply( xx$uniche.marginality$X, 2, sd)) 
            
            points(x$dataset[,c('lon','lat')], col=sensitivity, cex=tempUncert, pch=20)
            text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)
            title('Marginality')

            if (legend==TRUE){
                legend('topright', legend=c('High','Intermediate','Low'), fill=rgb(c(1,0.5,0),c(0,0.5,1),0,0.9), border=NA, title='Sensitivity', bty='n')}

        }else{

            sensitivity = rgb(abs(x$uniche.volume$PCC$original),1-abs(x$uniche.volume$PCC$original),0,0.9)
            tempUncert = 2.5 + 2 * apply( xx$uniche.volume$X, 2, sd) / max(apply( xx$uniche.volume$X, 2, sd)) 
            
            points(x$dataset[,c('lon','lat')], col=rgb(abs(x$uniche.volume$PCC$original),1-abs(x$uniche.volume$PCC$original),0,0.9), cex=1+3*abs(x$uniche.volume$PCC$original), pch=20)
            text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)
            
            if (legend==TRUE){
                legend('topright', legend=c('High','Intermediate','Low'), fill=rgb(c(1,0.5,0),c(0,0.5,1),0,0.9), border=NA, title='Sensitivity', bty='n')
            }
            title('Volume')
        }
        
    }else{
        
        par(mfrow=c(2,2), las=2)

        ##
        
        tempUncert = 2.5 + 2 * apply( xx$uniche.marginality$X, 2, sd) / max(apply( xx$uniche.marginality$X, 2, sd)) 
        sensitivity = rgb(abs(x$uniche.marginality$PCC$original),1-abs(x$uniche.marginality$PCC$original),0,0.9)

        ## plot(x$uniche.marginality)
        ## points(xx$uniche.marginality$PCC$original, pch=20, cex=tempUncert, col=sensitivity)
        ## abline(h=0, col='red', lty=2)
        
        plot(seq(ncol(x$uniche.marginality$X)), apply( x$uniche.marginality$PCC[,c('min. c.i.','max. c.i.')], 1, sum)/2, xaxt="n", xlab='Points', ylab='PCC', main='Marginality', ylim=c(-1,1), pch=20, cex=tempUncert-1.5, col=sensitivity)
        arrows( x0=seq(ncol(x$uniche.marginality$X)), y0=x$uniche.marginality$PCC$min, x1=seq(ncol(x$uniche.marginality$X)), y1=x$uniche.marginality$PCC$max, code=0)
        axis(1, at=seq(ncol(x$uniche.marginality$X)), las=2)
        abline(h=0, col='red', lty=2)

        plot(shape, main='Marginality')
        points(x$dataset[,c('lon','lat')], col=sensitivity, cex=tempUncert, pch=20)
        text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)

        if (legend==TRUE){
            legend('topright', legend=c('High sensitivity','Intermediate sensitivity','Low sensitivity', 'Temporal uncertainity'), fill=rgb(c(1,0.5,0,1),c(0,0.5,1,1),c(0,0,0,1),c(0.9,0.9,0.9,0)), pch=c(NA,NA,NA,1), pt.cex=2, border=NA, bty='n')
        }
        
        ##

        tempUncert = 2.5 + 2 * apply( xx$uniche.volume$X, 2, sd) / max(apply( xx$uniche.volume$X, 2, sd)) 
        sensitivity = rgb(abs(x$uniche.volume$PCC$original),1-abs(x$uniche.volume$PCC$original),0,0.9)
        
        plot(seq(ncol(x$uniche.volume$X)), apply( x$uniche.volume$PCC[,c('min. c.i.','max. c.i.')], 1, sum)/2, xaxt="n", xlab='Points', ylab='PCC', main='Volume', ylim=c(-1,1), pch=20, cex=tempUncert-1.5, col=sensitivity)
        arrows( x0=seq(ncol(x$uniche.volume$X)), y0=x$uniche.volume$PCC$min, x1=seq(ncol(x$uniche.volume$X)), y1=x$uniche.volume$PCC$max, code=0)
        axis(1, at=seq(ncol(x$uniche.volume$X)), las=2)
        abline(h=0, col='red', lty=2)
        
        plot(shape, main='Volume')
        points(x$dataset[,c('lon','lat')], col=sensitivity, cex=tempUncert, pch=20)
        text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)

        if (legend==TRUE){
            legend('topright', legend=c('High sensitivity','Intermediate sensitivity','Low sensitivity', 'Temporal uncertainity'), fill=rgb(c(1,0.5,0,1),c(0,0.5,1,1),c(0,0,0,1),c(0.9,0.9,0.9,0)), pch=c(NA,NA,NA,1), pt.cex=2, border=NA, bty='n')
        }
        
        ##
    }
}
