average.for.PAV <-function(y, wt = rep(1, length(y)))
{
  ## compute a weighted average of a vector, y
  if(any(is.na(wt))) stop("NA's not allowed for wt")
  if(any(wt < 0))
    stop("wt must be a vector of NON-NEGATIVE weights")
  if(length(wt) != length(y)) stop("y and wt must be vectors of the same length")
  ## if any observations have Infinite weight, return the simple
  ## (unweighted) average of only those observations (giving no
  ## weight to observations with finite weight)
  if(any(wt == Inf)) {
    wt[wt < Inf] <- 0
    wt[wt == Inf] <- 1
  }
  ## if all weights are zero, return the simple (unweighted)
  ## average of y
  if(sum(wt) == 0)
    wt <- rep(1, length(wt))
  return(sum((y * wt)/sum(wt)))
}

PAV <- function(y, wt = rep(1,length(y)))
{
  ## This is a modification of Derick's PAV program
  ##
  ## (Weighted) Pool-Adjacent-Violators (PAV) algorithm
  ## for non-parametric monotonic (decreasing) regression of y on x
  n <- length(y)
  if(n != length(wt))
    stop("y, and wt must be vectors of equal length")
  yhat <- y       # initialize while loop
  j <- count <- 1
  k <- 2
  support <- vector("numeric", n)
  support[count] <- j
  while(k <= n) {
    while(yhat[j] < yhat[k]) {
      yhat[j:k] <- average.for.PAV(y[j:k], wt[j:k])
      if(yhat[support[count]] < yhat[k]) {
        j <- support[count]
        if(count > 1)
          count <- count - 1
      }
      else {
        k <- ifelse(k == n, k, k + 1)
      }
    }
    
    count <- count + 1
    support[count] <- j
    j <- k
    k <- k + 1
  }
  return(list(y = yhat, wt))
}

