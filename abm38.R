# --------------------------------------------------------------------- #

# R code for Adaptive Beta Method
# version 38, June 9, 2015
# from the paper
# Adaptive credible intervals on stratigraphic ranges 
#   when recovery potential is unknown
# by Steve Wang, Phil Everson, Heather Zhou, David Chudzicki, Dasol Park
# Paleobiology 2015
# Before using, please contact the authors for updates, bug fixes, etc.:
# Steve Wang, scwang@swarthmore.edu

# --------------------------------------------------------------------- #






# -------------------- Adaptive Beta method -------------------- #



abm38 <- function(x, distance=T, ext=T, base=NULL, conf=.9, PLOT=0)    {

  # x = vector of locations or dates of fossil occurrences
  # distance = T if measurements represent distance above a base,
  #            F if measurements represent time
  # ext = T if extinction, F if origination 
  # base = value to consider 0 if known; otherwise the min is used (if upperbound) and
  #         sample size is decreased by 1
  # conf = confidence level
  # PLOT = 1 to show plots



  # ---------- FUNCTION DEFINITIONS ---------- #

  # following are 'wrapper' functions required for R's integrate() function to work
  # note: these do not check for xmax < th; if xmax > th, an error will not be generated

  prmean <- 0
  prSD <- 2

  integrand.neglambdas <- function(L,th,x)  {
  # L = vector of lambdas,  th = theta value (scalar),  x = vector of strat. positions
    k <- length(L)    
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1-L[i])/th * 1/(1-x/th)^L[i] ) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)
  }

  integrand.poslambdas <- function(L,th,x)  {
  # L = vector of lambdas,  th = theta value (scalar),  x = vector of strat. positions
    k <- length(L)    
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1+L[i])/th * (x/th)^L[i] ) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)
  }

  integrand.thetasnegL <- function(th,L,x)  {
  # th = vector of thetas,  L = lambda value (scalar), x = vector of strat. positions
    k <- length(th) 
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1-L)/th[i] * 1/(1-x/th[i])^L ) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)  
  }

  integrand.thetasposL <- function(th,L,x)  {
  # th = vector of thetas,  L = lambda value (scalar), x = vector of strat. positions
    k <- length(th) 
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1+L)/th[i] * (x/th[i])^L) ) )
    output <- output + dnorm(L, prmean,prSD, log=T) + log(1/th)
    output <- exp(output)
    return(output)  
  }

  
  # the function below is used in likelihood calculations
  drefbeta <- function(x,L)  {
    if(L<=0)  { 
      return(dbeta(x, 1, 1-L))
    }  else  return(dbeta(x, 1+L, 1))
  }



  # ---------- PRE-PROCESS DATA ---------- #

  # If a base is specified, check that it is valid
  if(!is.null(base)) 
    if( (distance & ext)   & (base > min(x))  | 
        (distance & !ext)  & (base < max(x))  |
        (!distance & ext)  & (base < max(x))  |
        (!distance & !ext) & (base > min(x)) )
      stop("Invalid value for base")

  # Convert units relative to base of section or other zero point
  xraw <- x                            # save a copy for later use
  if((distance & ext) | (!distance & !ext))  {
    if(!is.null(base))                 # if a base is specified
      x <- x - base 
    if(is.null(base))  {               # if a base is not specified
      base <- min(x)                   # condition on smallest value and
      x <- x - base 
      x <- sort(x, decreasing=F)[-1]   #   reduce sample size by 1
    } 
  } 
  if((distance & !ext) | (!distance & ext))  {
    if(!is.null(base))                 # if a base is specified
      x <- base - x
    if(is.null(base))  {               # if a base is not specified
      base <- max(x)                   # condition on smallest value and
      x <- base - x   
      x <- sort(x, decreasing=F)[-1]   #   reduce sample size by 1
    } 
  } 

  # Scale data so theta is approx. 100 (solely for numerical stability)
  xmax <- max(x);     n <- length(x)
  simplethhat <- (n+1)/n * xmax
  scalefactor <- 100/simplethhat
  x <- x * scalefactor
  xmax <- max(x)

  # Set iteration parameters
  upperlimth <- 500;     numstepsth <- 1000
  lowerlimL <- -10;      upperlimL <-  10;      numstepsL <- 40
  Lvals <- seq(lowerlimL, upperlimL, length.out=numstepsL)
  Ldens <- rep(NA, numstepsL)
  thetavals <- seq(xmax, upperlimth, length.out=numstepsth)
  thdens <- rep(NA, numstepsth)



  # ---------- ESTIMATE LAMBDA ---------- #

  # increment lambda values, integrating over theta values for each
  for(i in 1:numstepsL)  {
    Ldens[i] <- ifelse( Lvals[i]<=0, 
                  integrate(integrand.thetasnegL, xmax,upperlimth, L=Lvals[i], x=x)$value,
                  integrate(integrand.thetasposL, xmax,upperlimth, L=Lvals[i], x=x)$value )
  }
  # normalize lambda pdf to unit area  
  Ldens <- Ldens/sum(Ldens) 

  # calculate posterior quantities
  # cutoff <- which.max( cumsum(Ldens) >= .5 )       # posterior median
  # Lmed <- Lvals[cutoff]
  Lmean <- sum(Lvals*Ldens)
  Lhat <- Lmean



  # ---------- ESTIMATE THETA ---------- #

  # increment theta values, integrating over lambda values for each
  for(i in 1:numstepsth)  {
    thdens[i] <- ( integrate(integrand.neglambdas, -Inf,0, th=thetavals[i], x=x)$value
                 + integrate(integrand.poslambdas,  0,Inf, th=thetavals[i], x=x)$value )
  }
  # normalize theta pdf to unit area  
  thdens <- thdens/sum(thdens) 

  # calculate posterior quantities
  cutoff <- which.max( cumsum(thdens) >= conf )     # upper quantile of posterior
  CIupper <- thetavals[cutoff]
  cutoff <- which.max( cumsum(thdens) >= .5 )       # posterior median
  thmed <- thetavals[cutoff] 
  # thmean <- sum(thetavals*thdens)
  thhat <- thmed



  # ---------- COLLECT RESULTS ---------- #

  # Un-scale data back to original scale
  thhat <- thhat/scalefactor
  xmax <- xmax/scalefactor
  CIupper <- CIupper/scalefactor
  thetavals <- thetavals/scalefactor

  # Un-convert units back to original units (from units relative to base)
  if((distance & ext) | (!distance & !ext))  {
    thhat <- thhat + base
    xmax <- xmax + base
    CIupper <- CIupper + base
    thetavals <- thetavals + base
  }
  if((distance & !ext) | (!distance & ext))  {
    thhat <- base - thhat 
    xmax <- base - xmax 
    CIupper <- base - CIupper 
    thetavals <- base - thetavals 
  }
  temp <- t(as.matrix(c(thhat, xmax, CIupper, Lhat)))
  colnames(temp) <- c("th-hat","xmax","CIupper","L-hat")



  # ---------- PLOT RESULTS ---------- #

  if(PLOT)  {
    par(mfrow=c(2,2), mar=c(.8,1,.6,.6)*6)

    # determine histogram limits
    xmin <- min(xraw);      xmax <- max(xraw)
    # note: xmin and xmax are min and max of observed data; xlo and xhi are plot limits
    if(distance)  {
      if(ext)  {
        xlo <- base;     xhi <- thhat
      }  else  { 
        xlo <- thhat;    xhi <- base
      }
      xbase <- xlo;      xlabel <- "distance"
    }
    if(!distance)  {
      if(ext)  {
        xlo <- base;     xhi <- thhat
      }  else  {
        xlo <- thhat;    xhi <- base
      }
      xbase <- xhi;      xlabel <- "age"
    }

    # histogram of counts and recovery function
    hist(xraw, nclass=6, prob=T, xlim=c(xlo,xhi), col="gray", xlab=xlabel, border="white", 
         yaxt="s", ylab="recovery potential", main="Estimated recovery potential")
    points(xraw, rep(0,length(xraw)), cex=.5, pch="|")
    points(thhat,0, pch=8, col="red", cex=1.1)
    if((distance & ext) | (!distance & !ext))
      curve(drefbeta((x-xbase)/abs(xhi-xlo), Lhat)/abs(xhi-xlo), from=xlo, to=xhi, add=T)
    if((!distance & ext) | (distance & !ext))
      curve(drefbeta((x-xbase)/abs(xhi-xlo), -Lhat)/abs(xhi-xlo), from=xlo, to=xhi, add=T)

    # blank space
    plot(0,0, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")  

    # determine posterior plot limits
    if(distance)  {
      if(ext)  {
        xlims <- c(xmax, max(thetavals));      subset <- thetavals <= CIupper
      }  else  {
        xlims <- c(min(thetavals), xmin);      subset <- thetavals >= CIupper
      }
      xlabel <- "distance"
    }
    if(!distance)  {
      if(ext)  {
        xlims <- c(xmin, min(thetavals));      subset <- thetavals >= CIupper
      }  else  {
        xlims <- c(max(thetavals), xmax);      subset <- thetavals <= CIupper
      }
      xlabel <- "age"
    }

    # marginal posterior for theta
    plot(thetavals, thdens, xlab="theta", ylab="probability density", type="l", xlim=xlims)
    title(main="Posterior for theta")
    polygon( c(thetavals[subset],rev(thetavals[subset])), 
             c(thdens[subset],rep(0,numstepsth)[subset]), col=gray(.8))
    abline(h=0, col=gray(.8))
    points(thhat,0, pch=8, col="red", cex=.8)

    # marginal posterior for lambda
    plot(Lvals,Ldens, xlab="lambda", ylab="probability density", type="l", xlim=c(-6,4))
    title(main="Posterior for lambda")
    abline(h=0, col=gray(.8))
    points(Lhat,0, pch=18, col="red", cex=1.2)
  }

  return(temp)
}


