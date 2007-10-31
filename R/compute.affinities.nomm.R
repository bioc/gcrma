compute.affinities.nomm <- function (cdfname, verbose = TRUE) 
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
   pmIndex <- unlist(indexProbes(affinity.info, "pm"))
    #pm probes with seq
    my.xy2i <- get("xy2i", paste("package:", cdfpackagename, sep = ""))
    index1 <- match(pmIndex,my.xy2i(p$x,p$y))
    subIndex <- index1[!is.na(index1)]
    pmseq <- p$sequence[subIndex]
    prlen <- unique(nchar(pmseq))
    stopifnot(length(prlen) == 1 & prlen[1]==25)
    A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5])
    T13 <- 0
    C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10])
    G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15])
    if (verbose) 
        cat(".")
    apm <- vector("numeric", length(pmseq))
    for (i in seq(along = apm)) {
        charMtrx <- .Call("gcrma_getSeq", pmseq[i], PACKAGE = "gcrma")
        A <- cbind(charMtrx[1, ] %*% affinity.basis.matrix, charMtrx[2, 
            ] %*% affinity.basis.matrix, charMtrx[3, ] %*% affinity.basis.matrix)
        apm[i] <- A %*% affinity.spline.coefs
      }
    tmp <- get("xy2i", paste("package:", cdfpackagename, sep = ""))
    tmp.exprs = matrix(NA, nrow = max(pmIndex),ncol = 1)
    tmp.exprs[pmIndex[match(subIndex,pmIndex)]] = apm
    exprs(affinity.info) = tmp.exprs
    if (verbose) 
        cat("Done.\n")
    return(affinity.info)
}
