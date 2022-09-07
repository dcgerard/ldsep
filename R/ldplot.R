
## THIS FILE IS NO LONGER REQUIRED



#' findDistT
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

#' chromLoopN
#' 
#' This function iterates over chromosomes and runs helper functions to find 
#' distances between markers.
#' 
#'
#' @param ldout Output of ldfast function.
#' 
#' @param poolChrom Boolean to choose either combining chromosomes or not.
#'
#' @export
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


#' ldDecPrepN
#' 
#' @param cldf Dataframe containing good data.
#' 
#' @export
ldDecPrepN <- function(cldf, poolChrom = TRUE, method = c("mean", "median")) {
  dfnrow = nrow(cldf)
  if (poolChrom) {
    df <- data.frame(dist = unique(cldf[['dist']]), r2 = rep(NA, length(unique(cldf[['dist']]))))
    for (i in seq(1, nrow(df))) {
      if (method == "mean") {
        df[['r2']][i] = mean(cldf[cldf[['dist']] == df[['dist']][i],][['r2']])
      }
      else {
        df[['r2']][i] = median(cldf[cldf[['dist']] == df[['dist']][i],][['r2']])
      }
    }
  }
  
  else {
    print("not ready yet")
    df = NULL
  }
  df
}


#' ldPlotN
#' 
#' This function returns ld plots for the specified chromosome.
#' 
#' @param ldout ldfast function output.
#' 
#' @param poolChrom Boolean for selecting if you want to combine chromosomes 
#' in the final LD decay plot. If FALSE will return distinct plots for each 
#' chromosome.
#' 
#' @param ldThresh Linkage disequilibrium threashold to be displayed on the plot.
#' 
#' @param allData Fits spline on all points if TRUE. If FALSE will fit for means
#' /medians or percentiles.
#' 
#' @export
ldPlotN <- function(ldout, poolChrom = TRUE, ldThresh = 0.2, allData = TRUE) {
  p <- list()
  
  if (!poolChrom) {
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
    
    if (allData) {
      scamp <- scam::scam(formula = r2~s(dist, bs = "mdcx"), data = df1)
      df1$spline <- scam::predict.scam(scamp)
    }
    if (!allData) {
      df1 <- ldDecPrepN(df1, poolChrom = TRUE, method = "mean")
      scamp <- scam::scam(formula = r2~s(dist, bs = "mdcx"), data = df1)
      df1$spline <- scam::predict.scam(scamp)
    }
    
#    p[[1]] = ggplot2::qplot(x = dist, y = r2, data = df1) +
#      ggplot2::geom_line(ggplot2::aes(x = dist, y = spline), colour = "black") +
#      ggplot2::geom_hline(yintercept = ldThresh, colour = "blue", alpha = 0.5) +
#      ggplot2::theme_bw() + ggplot2::labs(title = "LDdecay Plot") 
    print(max(df1[['spline']]) / 2)
    p[[1]] = ggplot2::ggplot(data = df1) + ggplot2::geom_line(ggplot2::aes(x = dist, y =spline))
 ##     ggplot2::geom_hline(yintercept = ldThreshN)
    
  }
  
  return(p)
}

#' dataMaster
#' 
#' This function returns ld plots for the specified chromosome.
#' 
#' @param ldout ldfast function output.
#' 






## this function takes an ldout from {ldfast} and pairs distance to ldest values to graph ldDecay (ss stands for summary stat)
dataMaster <- function(ldout, combineDist = TRUE, combineSS = c('mean', 'median', 'quantile'), quantile = 0.9) {
  ldim <- dim(ldout$ldmat)[1] # will be square
  pm <- matrix(rep(NA, (ldim^2)), nrow = ldim)
  ## calculate distances and put them in adjacent matricies
  for (r in seq(1, ldim - 1)) {
    for (c in seq(1 + r, ldim)) {
      pm[r, c] = abs(ldout[['loc']][r] - ldout[['loc']][c])
    }
  }
  opp <- list(ld = matCorn(ldout$ldmat),d =  matCorn(pm))
  
  if (combineDist) {
    if (length(combineSS) > 1) combineSS = 'mean' ## this is sort of sketchy but works for now
    opp <- combineLoc(dMout = opp, ss = combineSS, quantile = quantile)
    if (combineSS == 'quantile') combineSS = paste(combineSS, quantile)
    opp[['ss']] = combineSS
  }
  opp
}

## this function extracts the corner elements of a matrix
matCorn <- function(sqMat) {
  ldim <- dim(sqMat)[1]
  opv <- rep(NA, sum(seq(1, ldim - 1)))
  h = 1
  for (r in seq(1, ldim - 1)) {
    for (c in seq(1 + r, ldim)) {
      opv[h] = sqMat[r, c]
      h = h + 1
    }
  }
  opv
}

## this function takes vec of d and ldest and combines distances with summary stats
combineLoc <- function(dMout, ss = c('mean', 'median', 'quantile'), quantile) {
  o <- unique(dMout[['d']])
  v <- o
  for (i in seq(1, length(o))) {
    ## this is just creating a vector of ldests of the same dist and averaging
    if (ss == 'mean') v[i] = mean(dMout[['ld']][dMout[['d']] == dMout[['d']][i]])
    else if (ss == 'median') v[i] = median(dMout[['ld']][dMout[['d']] == dMout[['d']][i]])
    else if (ss == 'quantile') v[i] = quantile(dMout[['ld']][dMout[['d']] == dMout[['d']][i]], quantile)
  }
  list(ld = v, d = o)
}


plotLD <- function(dMout, dof = 8) {
  df <- data.frame(d = dMout[['d']], ld = dMout[['ld']])

  ## fit a spline here
  scamD <- scam::scam(formula = ld~s(d, bs = "mdcx", k = dof), data = df)
  df$s <- scam::predict.scam(scamD)
  
  ggplot2::ggplot(data = df) + ggplot2::geom_point(ggplot2::aes(x = d, y = ld)) +
    ggplot2::geom_line(ggplot2::aes(x = d, y = s))
}

