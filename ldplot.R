findDistT <- function(ldmat, loci, chrom, poolChrom = TRUE) {
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
  
  if (poolChrom) {
    opdf <- data.frame(chrom, dist = ovn2, r2 = ovn)
  }
  
  else {
    uniqueDist <- base::unique(ovn2)
    obs <- rep(NA, length(uniqueDist))
    r2v <- rep(NA, length(uniqueDist))
    for (i in seq(1, length(uniqueDist))) {
      r2v[i] = mean(ovn[ovn2 == uniqueDist[i]])
      obs[i] = length(ovn[ovn2 == uniqueDist[i]])
    }
    opdf <- data.frame(chrom, dist = uniqueDist, r2v, obs)
  }
  
  opdf
}


chromLoopN <- function(ldout, poolChrom = TRUE) {
  stopifnot(!is.null(ldout$ldmat), !is.null(ldout$loc))
  allchrom = rep(1)
  chromSupplied = !is.null(ldout$chrom)
  if (chromSupplied) {
    allchrom = base::unique(ldout$chrom) ## find all chroms
  }
  opframe = data.frame(chrom = c(), dist = c(), r2v = c(), obs = c())
  
  if (poolChrom) {
    if (chromSupplied) {
      for (i in allchrom) {
        bChrom = (ldout$chrom == i)
        opframe = rbind(opframe, findDistT(ldmat = ldout$ldmat[bChrom, bChrom], loci = ldout$loc[bChrom], chrom = i, poolChrom = TRUE))
      }
    }
    else {
      opframe = rbind(opframe, findDistT(ldmat = ldout$ldmat, loci = ldout$loc, chrom = NA, poolChrom = TRUE))
    }
  }
  
  else { ## not pool chromed
    if (chromSupplied) {
      for (i in allchrom) {
        bChrom = (ldout$chrom == i)
        opframe = rbind(opframe, findDistT(ldmat = ldout$ldmat[bChrom, bChrom], loci = ldout$loc[bChrom], chrom = i, poolChrom = FALSE))
      }
    }
    else {
      opframe = rbind(opframe, findDistT(ldmat = ldout$ldmat, loci = ldout$loc, chrom = NA, poolChrom = FALSE))
    }
  }
  
  opframe
}



ldDecPrepN <- function(cldf, combineChrom = TRUE, method = c("mean", "median")) {
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


ldPlotN <- function(ldout, combineChrom = TRUE, ldThresh = 0.2) {
  p <- list()
  
  if (!combineChrom) {
    df2 <- chromLoopN(ldout, poolChrom = FALSE)
    
    for (i in unique(df2$chrom)) {
      buffDF = df2[df2[['chrom']] == i,]
      scamp <- scam::scam(formula = r2v~s(dist, bs = "mdcx"), data = buffDF)
      buffDF$spline <- scam::predict.scam(scamp)
      
      p[[i]] = ggplot2::qplot(x = dist, y = r2v, data = buffDF) + 
        ggplot2::geom_line(ggplot2::aes(x = dist, y = spline), colour = "black") +
        ggplot2::geom_hline(yintercept = ldThresh, colour = "blue", alpha = 0.5) +
        ggplot2::theme_bw() + ggplot2::labs(title = paste("LDdecay Plot for Chrom", i))
    }
    
  }
  else {
    df1 <- chromLoopN(ldout, poolChrom = TRUE)
    df2 <- ldDecPrepN(df1, combineChrom = TRUE, method = "mean")
    scamp <- scam::scam(formula = r2v~s(dist, bs = "mdcx"), data = df2)
    df2$spline <- scam::predict.scam(scamp)
    
    p[[1]] = ggplot2::qplot(x = dist, y = r2v, data = df2) +
      ggplot2::geom_line(ggplot2::aes(x = dist, y = spline), colour = "black") +
      ggplot2::geom_hline(yintercept = ldThresh, colour = "blue", alpha = 0.5) +
      ggplot2::theme_bw() + ggplot2::labs(title = "LDdecay Plot") 
    
  }

  return(p)
}
