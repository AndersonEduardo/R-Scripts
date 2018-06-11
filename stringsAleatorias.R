##funcao para gerar strigs aleatorias
##06/junho/2018

myFun <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

##exemplo
localidades = myFun(100)
