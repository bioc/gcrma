gcrma.bg.transformation.fast <- function(x,bhat,var.y,k){
  
  x <- x - bhat ##this is an unbiased estimate
  Index <- x > 0
  x[!Index] <- 0 ##trick to not get -Inf in log
  
  ##alpha is the weith used to create weighted estimate
  ##its basically  for smoothing edge
  alpha <- rep(0,length(x))
  xplusk <- x[Index]+k  #x plus k
  logxplusk <- log(xplusk)
  alpha[Index] <- (logxplusk - log(k))*(logxplusk - log(k) + k/xplusk) /
    ( ( logxplusk - log(k) + k/xplusk)^2 + var.y[Index]/xplusk^2)
  
  y <- exp(alpha*log(x+k) + (1-alpha)*log(k))
}
