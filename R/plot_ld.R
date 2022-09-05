#' Linkage Disequilibrium Decay Plotter
#'
#' @param obj Output of ldfast function.
#'
#' @export
plot_ld <- function(chromLoopDF, distinctChromPlot = FALSE) {
  p <- list()
  
  if (!distinctChromPlot) {
    # we have to combine
    
  p[[1]] = ggplot2::qplot(x = chromLoopDF$distance, y = chromLoopDF$meanr2, data = chromLoopDF$rtg) +
    ggplot2::theme_bw() + ggplot2::labs(title = paste("LDdecay Plot"))
  }
  else {
    for (i in unique(chromLoopDF$chrom)) {
      buffDF = chromLoopDF[chromLoopDF$obs == i,]
      scamp <- scam::scam(formula = meanr2~s(distance, bs = "mdcx"), data = buffDF)
      buffDF$spline <- scam::predict.scam(scamp)
      
      p[[i]] = ggplot2::qplot(x = buffDF$distance, y = buffDF$meanr2) +
        ggplot2::geom_line(ggplot2::aes(x = buffDF$distance, y = buffDF$spline)) +
        ggplot2::theme_bw() + ggplot2::labs(title = paste("LDdecay Plot for Chrom", i))
    }
  }
  return(p)
}


#' Read's in matrix and loci information and outputs two vectors ready to graph with
#'
#' @param mat ldmatrix from ldfast
#' 
#' @param loci vector of positions on the genome of different markers
#'
#' @export
readinss <- function(mat, loci, chrom) { ## this is solid to read in wanted information
  nr <- nrow(mat)
  nc <- ncol(mat)
  ov <- rep(NA, (nr - 1) * (nc - 1)) # r2 for ov
  ov2 <- ov # dist for ov2
  for (i in seq(1, nr - 1)) {
    for (j in seq(i + 1, nc)) {
      ov[nc * (i - 1) + (j - i)] = mat[i, j] ## r2 values
      ov2[nc * (i - 1) + (j - i)] = abs(loci[i] - loci[j]) ## dist b/w loci
    }
  }
  nao <- is.na(ov)
  
  ## we need it to calculate one point per distance
  ovn <- ov[!nao]
  ovn2 <- ov2[!nao]
  
  uniqueDist <- base::unique(ovn2)
  obs <- rep(NA, length(uniqueDist))
  meanr2 <- rep(NA, length(uniqueDist))
  # medianr2 <- rep(NA, length(uniqueDist))
  
  for (i in seq(1, length(uniqueDist))) {
    meanr2[i] = mean(ovn[ovn2 == uniqueDist[i]])
    obs[i] = length(ovn[ovn2 == uniqueDist[i]])
  }
  data.frame(chrom, distance = uniqueDist, meanr2, obs)
}


#' Chromosome Looper
#'
#' @param obj Output of ldfast function.
#'
#' @export
chromLoop <- function(obj) {
  stopifnot(!is.null(obj$ldmat), !is.null(obj$loc))
  allchrom = rep(1)
  nullchrom = !is.null(obj$chrom)
  if (nullchrom) {
    allchrom = base::unique(obj$chrom)
  }
  
  opframe = data.frame(chrom = c(), distance = c(), meanr2 = c(), obs = c())
  
  ## if chromosome argument is supplied
  if (nullchrom) {
    for (i in allchrom) {
      bChrom = (obj$chrom == i)
      opframe = rbind(opframe, readinss(mat = obj$ldmat[bChrom, bChrom], loci = obj$loc[bChrom], chrom = i))
    }
  }
  ## if no chrom arg is supplied
  else {
    opframe = rbind(opframe, readinss(mat = obj$ldmat, loci = obj$loc, chrom = NA))
  }
  opframe
}


