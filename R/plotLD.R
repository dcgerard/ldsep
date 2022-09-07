#' plotLD
#'
#' This function takes a list containing two vectors one with ld estimates and 
#' another consisting of the length (bp) between the two markers used to 
#' calculate the corresponding ld estimate.
#'
#' @param dMout List output from dataMaster function
#' 
#' @param dof Number of degrees freedom for spline
#' 
#' @param points Boolean if you wish to see scatter plot as well as spline.
#' 
#' @export
plotLD <- function(dMout, dof = 8, points = FALSE) {
  df <- data.frame(d = dMout[['d']], ld = dMout[['ld']])
  
  ## fit a spline here
  scamD <- scam::scam(formula = ld~s(d, bs = "mdcx", k = dof), data = df)
  df$s <- scam::predict.scam(scamD)
  
  opp <- ggplot2::ggplot(data = df) +
    ggplot2::geom_line(ggplot2::aes(x = d, y = s), colour = "blue")
  
  if (points) {
    opp <- opp + ggplot2::geom_point(ggplot2::aes(x = d, y = ld))
  }
  opp
}


#' dataMaster
#' 
#' This function takes an ld output from code{ldfast} and pairs distance between
#' markers to ld estimates producing a list with an ld estimate component, a 
#' distance component and an optional summary statistic used component (if 
#' combineDist == TRUE).
#' 
#' @param ldout Output from code{ldest::ldfast()} function.
#' 
#' @param combineDist Boolean argument to combine ld estimates by distance.
#' 
#' @param combineSS Summary statistic to be used if combineDist == TRUE (default
#' is mean).
#' 
#' @param quantile Quantile to be used if combineSS == "quantile".
#' 
#' @export
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


#' combineLoc
#' 
#' This function takes a vector of distances between loci and a ld estimate 
#' vector and will calculate a summary statistic grouping by distance.
#' 
#' @param dMout list output from dataMaster function
#' 
#' @param ss Summary statistic to be calculated for each distance.
#' 
#' @param quantile Quantile to be used if ss == 'quantile'.
#' 
#' @export
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


#' matCorn
#' 
#' This function extracts the corner elements of a matrix
#' 
#' @param sqMat A square Matrix.
#' 
#' @export
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

