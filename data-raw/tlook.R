## code to prepare `tlook` dataset

## Implement the method of Cohen (1957) to derive lookup tables to estimate
## the parameters of the doubly truncated normal.
##
## Cohen, A. C. (1957). On the solution of estimating equations for
## truncated and censored samples from normal populations.
## Biometrika, 44(1/2), 225-236.

# x0: lower bound of truncated normal
# w: width of truncated normal
# mu: Mean of underlying normal
# sig: Standard deviation of underlying normal
# vrat: sample mean over w
# srat: sample variance over w

# Idea: We can easily get (xi1, xi2) to (vrat, srat)
# We actually want inverse function (vrat, srat) to (xi1, xi2)
# xi is a known function of mu and sig, so we can get mu and sig from xi.

z1bar <- function(xi1, xi2) {
  stats::dnorm(xi1) / (stats::pnorm(xi1, lower.tail = FALSE) - stats::pnorm(xi2, lower.tail = FALSE))
}

z2bar <- function(xi1, xi2) {
  stats::dnorm(xi2) / (stats::pnorm(xi1, lower.tail = FALSE) - stats::pnorm(xi2, lower.tail = FALSE))
}

v1 <- function(x, x0) {
  mean(x - x0)
}

v2 <- function(x, x0) {
  mean((x - x0)^2)
}

H1 <- function(xi1, xi2) {
  z1 <- z1bar(xi1 = xi1, xi2 = xi2)
  z2 <- z2bar(xi1 = xi1, xi2 = xi2)
  (z1 - z2 - xi1) / (xi2 - xi1)
}

H2 <- function(xi1, xi2) {
  z1 <- z1bar(xi1 = xi1, xi2 = xi2)
  z2 <- z2bar(xi1 = xi1, xi2 = xi2)
  (1 + xi1 * z1 - xi2 * z2 - (z1 - z2)^2) / (xi2 - xi1)^2
}

sigfromxi <- function(xi1, xi2, w) {
  w / (xi2 - xi1)
}

mufromxi <- function(x0, sigma, xi1) {
  x0 - sigma * xi1
}

griddf <- expand.grid(xi1 = seq(-10, 0.5, length.out = 101),
                      xi2 = seq(-0.45, 10.05, length.out = 101))
griddf$vrat <- NA_real_
griddf$srat <- NA_real_

for (i in seq_len(nrow(griddf))) {
  griddf$vrat[[i]] <- H1(xi1 = griddf$xi1[[i]], xi2 = griddf$xi2[[i]])
  griddf$srat[[i]] <- H2(xi1 = griddf$xi1[[i]], xi2 = griddf$xi2[[i]])
}

griddf <- griddf[griddf$vrat >= 0 &
                   griddf$vrat <= 1 &
                   griddf$srat >= 0 &
                   griddf$srat <= 1 &
                   griddf$xi2 > griddf$xi1 &
                   !is.na(griddf$xi1) &
                   !is.na(griddf$xi2) &
                   !is.na(griddf$vrat) &
                   !is.na(griddf$srat), ]

## sanity check ----
num <- 2001
w <- 1
x0 <- 0
sig <- sigfromxi(xi1 = griddf[num, "xi1"], xi2 = griddf[num, "xi2"], w = w)
mu <- mufromxi(xi1 = griddf[num, "xi1"], x0 = x0, sigma = sig)
mu
sig

samp_mean <- griddf$vrat[[num]] * w
samp_var <- griddf$srat[[num]] * w^2
ldsep:::mom_tnorm(mu = mu, sigma = sig, a = x0, b = x0 + w)
truncnorm::etruncnorm(a = x0, b = x0 + w, mean = mu, sd = sig)
sqrt(truncnorm::vtruncnorm(a = x0, b = x0 + w, mean = mu, sd = sig))
samp_mean
sqrt(samp_var)

plot(griddf$vrat, griddf$srat)

## Save data ----
tlook <- griddf
usethis::use_data(tlook, overwrite = TRUE, internal = TRUE)
