

#f = sqr( eta^2 / ( 1 - eta^2 ) ).

computef <- function(petasq) {
  sqrt(petasq/(1-petasq))
} 

computef(0.139)  

computed <- function(petasq) {
  
}


options(scipen=999)


actxt <- function(t, n, lambda) {
  (exp(-lambda * t) * (lambda * t)^n)/factorial(n)
}
