#' Linkage Disequilibrium Decay Plotter
#'
#' @param obj Output of ldfast function.
#'
#' @export
plot_ld <- function(ldout, distinctChromPlot = FALSE, ldThreshold = 0.2) {
  p <- list()
  
  if (!distinctChromPlot) {
    chromLoopD <- chromLoop(ldout)
    chromLoopDF <- ldDecPrep(chromLoopD, method = "mean")
    
    scamp <- scam::scam(formula = r2v~s(dist, bs = "mdcx"), data = chromLoopDF)
    chromLoopDF$fitted <- scam::predict.scam(scamp)
    
    p[[1]] = ggplot2::qplot(x = dist, y = r2v, data = chromLoopDF) +
      ggplot2::geom_hline(yintercept = ldThreshold, colour = "blue", alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes(x = dist, y = fitted))
  }
  
  else {
    chromLoopDF <- chromLoopDC(ldout)
    for (i in unique(chromLoopDF$chrom)) {
      buffDF = chromLoopDF[chromLoopDF[['chrom']] == i,] ## CHECK THIS WORKS
      scamp <- scam::scam(formula = r2v~s(dist, bs = "mdcx"), data = buffDF)
      buffDF$spline <- scam::predict.scam(scamp)
      
      p[[i]] = ggplot2::qplot(x = buffDF$dist, y = buffDF$r2v) +
        ggplot2::geom_line(ggplot2::aes(x = buffDF$dist, y = buffDF$spline)) +
        ggplot2::theme_bw() + ggplot2::labs(title = paste("LDdecay Plot for Chrom", i)) +
        ggplot2::geom_hline(yintercept = ldThreshold, colour = "blue", alpha = 0.5)
    }
  }
  return(p)
}


#' Read's in matrix and loci information and outputs two vectors ready to graph with
#'
#' @param ldmat ldmatrix from ldfast
#' 
#' @param loci vector of positions on the genome of different markers
#' 
#' @param chrom the chromosome the data is coming from
#'
#' @export
readinss <- function(ldmat, loci, chrom) { ## this is solid to read in wanted information
  nr <- nrow(ldmat)
  nc <- ncol(ldmat)
  ov <- rep(NA, (nr - 1) * (nc - 1)) # r2 for ov
  ov2 <- ov # dist for ov2
  for (i in seq(1, nr - 1)) {
    for (j in seq(i + 1, nc)) {
      ov[nc * (i - 1) + (j - i)] = ldmat[i, j] ## r2 values
      ov2[nc * (i - 1) + (j - i)] = abs(loci[i] - loci[j]) ## dist b/w loci
    }
  }
  
  nao <- is.na(ov)
  
  ## we need it to calculate one point per distance
  ovn <- ov[!nao]
  ovn2 <- ov2[!nao]
  
  uniqueDist <- base::unique(ovn2)
  obs <- rep(NA, length(uniqueDist))
  r2v <- rep(NA, length(uniqueDist))
  # medianr2 <- rep(NA, length(uniqueDist))
  
  for (i in seq(1, length(uniqueDist))) {
    r2v[i] = mean(ovn[ovn2 == uniqueDist[i]])
    obs[i] = length(ovn[ovn2 == uniqueDist[i]])
  }
  data.frame(chrom, dist = uniqueDist, r2v, obs)
}


#' findDist
#'
#' This function uses indices of an ldmat along with indices of a loci vector
#' to calculate the distance between two markers and returns the distance along
#' with the corresponding LD value.
#' 
#' @param ldmat LD matrix output from ldfast function.
#' 
#' @param loci Vector of marker positions corresponding with ldmat. For example 
#' the ith element of loci (loci[i]) corresponds to the element at ldmat[i,] and
#' ldmat[,i].
#' 
#' @param chrom The chromosome the current data is from.
#' 
#' @export
findDist <- function(ldmat, loci, chrom) { ## this is solid to read in wanted information
  nr <- nrow(ldmat)
  nc <- ncol(ldmat)
  ov <- rep(NA, (nr - 1) * (nc - 1)) # r2 for ov
  ov2 <- ov # dist for ov2
  for (i in seq(1, nr - 1)) {
    for (j in seq(i + 1, nc)) {
      ov[nc * (i - 1) + (j - i)] = ldmat[i, j] ## r2 values
      ov2[nc * (i - 1) + (j - i)] = abs(loci[i] - loci[j]) ## dist b/w loci in bp
    }
  }
  nao <- is.na(ov)
  ovn <- ov[!nao]
  ovn2 <- ov2[!nao]
  data.frame(chrom, dist = ovn2, r2 = ovn)
}



#' chromLoop
#' 
#' This function iterates over chromosomes and runs findDist on their loci.
#' 
#'
#' @param obj Output of ldfast function.
#'
#' @export
chromLoop <- function(obj) {
  stopifnot(!is.null(obj$ldmat), !is.null(obj$loc))
  
  allchrom = rep(1)
  chromSupplied = !is.null(obj$chrom)
  
  if (chromSupplied) {
    allchrom = base::unique(obj$chrom) ## find all chroms
  }
  opframe = data.frame(chrom = c(), distance = c(), meanr2 = c(), obs = c())
  
  
  
  if (chromSupplied) {
    for (i in allchrom) {
      bChrom = (obj$chrom == i)
      opframe = rbind(opframe, findDist(ldmat = obj$ldmat[bChrom, bChrom], loci = obj$loc[bChrom], chrom = i))
    }
  }
  else {
    opframe = rbind(opframe, findDist(ldmat = obj$ldmat, loci = obj$loc, chrom = NA))
  }
  opframe
}


#' chromLoopDC
#' 
#' This function iterates over chromosomes and runs findDist on their loci.
#' 
#'
#' @param obj Output of ldfast function.
#'
#' @export
chromLoopDC <- function(obj) {
  stopifnot(!is.null(obj$ldmat), !is.null(obj$loc))
  allchrom = rep(1)
  chromSupplied = !is.null(obj$chrom)
  if (chromSupplied) {
    allchrom = base::unique(obj$chrom) ## find all chroms
  }
  
  
  
  opframe = data.frame(chrom = c(), dist = c(), r2v = c(), obs = c())
  
  if (chromSupplied) {
    for (i in allchrom) {
      bChrom = (obj$chrom == i)
      opframe = rbind(opframe, readinss(ldmat = obj$ldmat[bChrom, bChrom], loci = obj$loc[bChrom], chrom = i))
    }
  }
  else {
    opframe = rbind(opframe, readinss(ldmat = obj$ldmat, loci = obj$loc, chrom = NA))
  }
  opframe
}



#' ldDecPrep
#' 
#' @param cldf Dataframe containing good data.
#' 
#' @export
ldDecPrep <- function(cldf, combineChrom = TRUE, method = c("mean", "median")) {
  dfnrow = nrow(cldf)
  if (combineChrom) {
    df <- data.frame(dist = unique(cldf[['dist']]), r2v = rep(NA, length(unique(cldf[['dist']]))))
    for (i in seq(1, nrow(df))) {
      if (method == "mean") {
      df[['r2v']][i] = mean(cldf[cldf[['dist']] == df[['dist']][i],][['r2']])
      }
      else {
        df[['r2v']][i] = median(cldf[cldf[['dist']] == df[['dist']][i],][['r2']])
      }
    }
  }

  else {
    print("not ready yet")
    df = NULL
  }
  df
}

