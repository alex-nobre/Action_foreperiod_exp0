

# From https://stats.stackexchange.com/questions/3930/are-there-default-functions-for-discrete-uniform-distributions-in-r
ddu<-function(x,k) ifelse(x>=1 & x<=k & round(x)==x,1/k,0)

pdu<-function(x,k) ifelse(x<1,0,ifelse(x<=k,floor(x)/k,1))

# Valores de x
x_dexp <- c(0.6, 1.2, 1.8)

# Taxa
meanx <- 1.2
xrate <- 1/meanx

# Criar distribuicoes
y_dexp <- dexp(x_dexp, xrate)
y_antiexp <- rev(y_dexp) #antiexponencial: valores da exponencial em ordem reversa

plot(x_dexp, y_dexp)
plot(x_dexp, y_antiexp)

# Compute cumulative probability density
y_pexp <- pexp(x_dexp, xrate)

plot(x_dexp, y_pexp)

# Compute hazard function
exphazard <- function(timepoint) {
  dexp(timepoint, xrate)/(1 - pexp(timepoint, xrate))
}

exphazard(x_dexp)

xuni <- c(1,2,3)
unihazard <- function(timepoint) {
  ddu(timepoint, 3)/(1 - pdu(timepoint, 3))
}


