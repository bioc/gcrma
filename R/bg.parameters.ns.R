left.sigma <- function(x,mu) sqrt(mean((x[x<mu]-mu)^2))


bg.parameters.ns <- function(x,affinities,order.aff=NULL,Q=.25,nbreaks=40,
                             monotonize.mu=TRUE,
                             monotonize.sigma=FALSE){
  ##ns stands for nonspecific
  n <- round(length(x)/nbreaks)
  if(is.null(order.aff)) order.aff=order(affinities)
  ##break up the scatter plot and get "mean" and "sigma"
  qs <- sapply(0:(nbreaks-1),function(k) {
    o1 <- k*n+1;
    o2 <- min((k+1)*n,length(x))
    y <- x[order.aff[o1:o2]]
    y <- log(y)
    mu <- quantile(y,Q)
    sigma <- sqrt(mean((y[y<mu]-mu)^2))
    c(affinities[order.aff[(o1+o2)/2]],mu,sigma)
  })
  if(monotonize.mu) qs[2,] <- -PAV(-qs[2,])$y
  ##fill in the blanks with approx
  bg.mu <- approx(qs[1,],qs[2,],xout=affinities,rule=2)$y 
  
  if(monotonize.sigma){
    qs[3,] <- -PAV(-qs[3,])$y
    bg.sigma <- approx(qs[1,],qs[3,],xout=affinities,rule=2)$y 
  }
  else
    ##bg.sigma <- approx(qs[1,],qs[3,],xout=affinities,rule=2)$y 
    bg.sigma <- rep(median(qs[3,]),length(affinities))
  
  return(list(bg.mu = bg.mu,bg.sigma=bg.sigma))
}
