##codigo para funcao uplot, que plota o resultado da analise de incerteza a partir do output da funcao 'uniche'.
##Anderson A. Eduardo
##12/Set/2018

uplot2 = function(x, shape=NULL, niche.metric=NULL, legend=TRUE){
  
  if (class(x) != 'uniche'){
    stop("Os dados de entrada (x) deve ser um objeto da classe 'uniche'.")
  }
  
  if (is.null(shape) & is.null(niche.metric)){
    stop("Quando um shapefile nao e fornecido o argumento 'niche.metric' deve ser especificado (as opcoes sao: 'TSS' e 'ROC').")
  }
  
  if (length(niche.metric) > 1){
    stop("Quando um shapefile nao e fornecido o argumento 'niche.metric' deve ser especificado com uma das opcoes: 'marginility' e 'ROC'.")
  }
  
  if (any(!niche.metric %in% c('TSS','ROC'))){
    stop("O argumento de entrada 'niche.metric' deve ser 'TSS' ou 'ROC'.")
  }
  
  
  if( is.null(shape) ){
    if (niche.metric == 'TSS'){
      
      sensitivity = rgb(abs(x$uniche.TSS$PRCC$original),1-abs(x$uniche.TSS$PRCC$original),0,0.9)
      tempUncert = 2.5 + 2 * apply( x$uniche.TSS$X, 2, sd) / max(apply( x$uniche.TSS$X, 2, sd)) 
      
      points(x$dataset[,c('lon','lat')], col=sensitivity, cex=tempUncert, pch=20)
      text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)
      title('TSS')
      
      if (legend==TRUE){
        legend('topright', legend=c('High','Intermediate','Low'), fill=rgb(c(1,0.5,0),c(0,0.5,1),0,0.9), border=NA, title='Sensitivity', bty='n')}
      
    }else{
      
      sensitivity = rgb(abs(x$uniche.ROC$PRCC$original),1-abs(x$uniche.ROC$PRCC$original),0,0.9)
      tempUncert = 2.5 + 2 * apply( x$uniche.ROC$X, 2, sd) / max(apply( x$uniche.ROC$X, 2, sd)) 
      
      points(x$dataset[,c('lon','lat')], col=rgb(abs(x$uniche.ROC$PRCC$original),1-abs(x$uniche.ROC$PRCC$original),0,0.9), cex=1+3*abs(x$uniche.ROC$PRCC$original), pch=20)
      text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)
      
      if (legend==TRUE){
        legend('topright', legend=c('High','Intermediate','Low'), fill=rgb(c(1,0.5,0),c(0,0.5,1),0,0.9), border=NA, title='Sensitivity', bty='n')
      }
      title('ROC')
    }
    
  }else{
    
    par(mfrow=c(2,2), las=2)
    
    ##
    
    tempUncert = 2.5 + 2 * apply( x$uniche.TSS$X, 2, sd) / max(apply( x$uniche.TSS$X, 2, sd)) 
    sensitivity = rgb(abs(x$uniche.TSS$PRCC$original),1-abs(x$uniche.TSS$PRCC$original),0,0.9)
    xData = seq(ncol(x$uniche.TSS$X))
    yData = apply( x$uniche.TSS$PRCC[,c('min. c.i.','max. c.i.')], 1, sum)/2
    if (abs(max(yData)) > 1 | abs(min(yData)) > 1){
      yAxis = 2*range(yData)
    } else{
      yAxis = ylim=c(-1,1) 
    }
    
    plot(x = xData, y = yData, xaxt="n", xlab='', ylab='PRCC', main='TSS', ylim = yAxis, pch=20, cex=tempUncert-1.5, col=sensitivity)
    arrows( x0=seq(ncol(x$uniche.TSS$X)), y0=x$uniche.TSS$PRCC$min, x1=seq(ncol(x$uniche.TSS$X)), y1=x$uniche.TSS$PRCC$max, code=0)
    axis(1, at=seq(ncol(x$uniche.TSS$X)), labels=paste('Point ', x$dataset$id, sep=''), las=2)
    abline(h=0, col='red', lty=2)
    
    plot(shape, main='TSS')
    points(x$dataset[,c('lon','lat')], col=sensitivity, cex=tempUncert, pch=20)
    text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)
    
    if (legend==TRUE){
      legend('topright', legend=c('High sensitivity','Intermediate sensitivity','Low sensitivity', 'Temporal uncertainity'), fill=rgb(c(1,0.5,0,1),c(0,0.5,1,1),c(0,0,0,1),c(0.9,0.9,0.9,0)), pch=c(NA,NA,NA,1), pt.cex=2, border=NA, bty='n')
    }
    
    ##
    
    tempUncert = 2.5 + 2 * apply( x$uniche.ROC$X, 2, sd) / max(apply( x$uniche.ROC$X, 2, sd)) 
    sensitivity = rgb(abs(x$uniche.ROC$PRCC$original),1-abs(x$uniche.ROC$PRCC$original),0,0.9)
    xData = seq(ncol(x$uniche.ROC$X))
    yData = apply( x$uniche.ROC$PRCC[,c('min. c.i.','max. c.i.')], 1, sum)/2
    if (abs(max(yData)) > 1 | abs(min(yData)) > 1){
      yAxis = 2*range(yData)
    } else{
      yAxis = ylim=c(-1,1) 
    }
    
    plot(x = xData, y = yData, xaxt="n", xlab='', ylab='PRCC', main='ROC', ylim=yAxis, pch=20, cex=tempUncert-1.5, col=sensitivity)
    arrows( x0=seq(ncol(x$uniche.ROC$X)), y0=x$uniche.ROC$PRCC$min, x1=seq(ncol(x$uniche.ROC$X)), y1=x$uniche.ROC$PRCC$max, code=0)
    axis(1, at=seq(ncol(x$uniche.ROC$X)), labels=paste('Point ', x$dataset$id, sep=''), las=2)
    abline(h=0, col='red', lty=2)
    
    plot(shape, main='ROC')
    points(x$dataset[,c('lon','lat')], col=sensitivity, cex=tempUncert, pch=20)
    text(x$dataset[,c('lon','lat')], labels=x$dataset$id, cex= 0.7)
    
    if (legend==TRUE){
      legend('topright', legend=c('High sensitivity','Intermediate sensitivity','Low sensitivity', 'Temporal uncertainity'), fill=rgb(c(1,0.5,0,1),c(0,0.5,1,1),c(0,0,0,1),c(0.9,0.9,0.9,0)), pch=c(NA,NA,NA,1), pt.cex=2, border=NA, bty='n')
    }
    
    ##
  }
}
