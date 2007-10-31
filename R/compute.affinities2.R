compute.affinities2 <- function (cdfname, verbose = TRUE) 
{
    if (verbose) 
        cat("Computing affinities")
    data(affinity.spline.coefs)
    affinity.basis.matrix <- ns(1:25, df = length(affinity.spline.coefs)/3)
    cleancdf <- cleancdfname(cdfname, addcdf = FALSE)
    cdfpackagename <- paste(cleancdf, "cdf", sep = "")
    probepackagename <- paste(cleancdf, "probe", sep = "")
    getCDF(cdfpackagename)
    getProbePackage(probepackagename)
    p <- get(probepackagename)
    p <- check.probes(p, cdfname)
    affinity.info <- new("AffyBatch", cdfName = cdfname)
    my.xy2i <- get("xy2i", paste("package:", cdfpackagename, sep = ""))
    all <- my.xy2i(p$x,p$y)
    all25 <- all[which(nchar(p$sequence)==25)]
    pseq <- p$sequence[which(nchar(p$sequence)==25)]
    prlen <- unique(nchar(pseq))
    stopifnot(length(prlen) == 1 & prlen[1]==25)
    A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5])
    T13 <- 0
    C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10])
    G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15])
    if (verbose) 
        cat(".")
    apm <- vector("numeric", length(pseq))
    for (i in seq(along = apm)) {
        charMtrx <- .Call("gcrma_getSeq", pseq[i], PACKAGE = "gcrma")
        A <- cbind(charMtrx[1, ] %*% affinity.basis.matrix, charMtrx[2, 
            ] %*% affinity.basis.matrix, charMtrx[3, ] %*% affinity.basis.matrix)
        apm[i] <- A %*% affinity.spline.coefs
      }
    tmp.exprs = matrix(NA, nrow = max(all),ncol = 1)
    tmp.exprs[all25] = apm
    
    exprs(affinity.info) = tmp.exprs
    if (verbose) 
        cat("Done.\n")
    return(affinity.info)
}
