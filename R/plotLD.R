#' plotLD
#' 
#' This function takes an ldout from `ldfast` and a vector of genetic distances 
#' and plots LD at various distances.
#' 
#' @param ldout Output from code{ldest::ldfast()} function.
#' 
#' @param groupPoints the number of discrete distance groupings.
#' 
#' 
#' @export
plotLD <- function(ldout, groupPoints = 20) {
  dMout <- dataMaster(ldout)
  dfone <- data.frame(r2 = dMout$ld, dist = dMout$d)
  bpdata <- sumDf(dfone, ngroup = groupPoints)
  pldata <- sumStat(bpdata)

  opPlot <- ggplot2::ggplot(pldata) +
    ggplot2::geom_line(ggplot2::aes(x = md, y = tens), color = "red", linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(x = md, y = tf), color = "black", linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(x = md, y = med), color = "black") +
    ggplot2::geom_line(ggplot2::aes(x = md, y = sf), color = "black", linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(x = md, y = nine), color = "red", linetype = 2) +
    ggplot2::labs(x = "dist (bp)", y = "STAT")

  return(opPlot)
}


#' dataMaster
#' 
#' This function takes an ld output from code{ldfast} and converts the data into
#' a dataframe.
#' 
#' @param ldout Output from code{ldest::ldfast()} function.
#' 
#' 
#' @export
dataMaster <- function(ldout) {
  ldim <- dim(ldout$ldmat)[1] # will be square
  pm <- matrix(rep(NA, (ldim^2)), nrow = ldim)
  ## calculate distances and put them in adjacent matricies
  for (r in seq.int(1, ldim - 1)) {
    for (c in seq.int(1 + r, ldim)) {
      pm[r, c] = abs(ldout[['loc']][r] - ldout[['loc']][c])
    }
  }
  opp <- list(ld = ldout[['ldmat']][upper.tri(ldout[['ldmat']])],d = pm[upper.tri(pm)])
  opp
}

#' sumStat
#' 
#' This function calculates quantiles based on data from sumDF.
#' 
#' @param DFncomb a sumDF output
#' 
sumStat <- function(DFncomb) {
  group <- unique(DFncomb[['val']])
  lgroup <- length(group)
  opdf <- data.frame(group, 
                     md = rep(NA, lgroup),
                     tens = rep(NA, lgroup),
                     tf = rep(NA, lgroup),
                     med = rep(NA, lgroup),
                     sf = rep(NA, lgroup),
                     nine = rep(NA, lgroup))
  
  for (i in group) {
    opdf[['md']][i] = mean(DFncomb[DFncomb[['val']] == i,][['dist']])
    opdf[['tens']][i] = stats::quantile(DFncomb[DFncomb[['val']] == i,][['r2']], .10)
    opdf[['tf']][i] = stats::quantile(DFncomb[DFncomb[['val']] == i,][['r2']], .25)
    opdf[['med']][i] = stats::quantile(DFncomb[DFncomb[['val']] == i,][['r2']], .50)
    opdf[['sf']][i] = stats::quantile(DFncomb[DFncomb[['val']] == i,][['r2']], .75)
    opdf[['nine']][i] = stats::quantile(DFncomb[DFncomb[['val']] == i,][['r2']], .90)
  }
  return(opdf)
}


#' sumDF
#' 
#' This function groups the data according to genetic distance into n groups.
#' 
#' @param DFncomb Datamaster output.
#' 
#' @param ngroup Number of groups.
#' 
sumDf <- function(DFncomb, ngroup = 4) {
  ## scaffold for this function
  DFncomb[['val']] = NA
  rb <- seq.int(from = 0, to = 1, by = 1/ngroup)
  lenrb <- length(rb)
  for (i in seq.int(1, lenrb - 1)) {
    DFncomb[['val']][DFncomb[['dist']] >= stats::quantile(DFncomb[['dist']], rb[i])] = i
  }
  return(DFncomb)
}
