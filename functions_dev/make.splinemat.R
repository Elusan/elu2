make.splinemat <- function (nseasons, order, dtfine = 1/100) 
{
  if (dtfine == 1) {
    d <- matrix(1, 1, 1)
  }
  else {
    dtspl <- 1/nseasons
    knots <- seq(0, 1, by = dtspl)
    x <- seq(0, 1 - dtfine, by = dtfine)
    if (order > 1) {
      d <- mgcv::cSplineDes(x, knots, ord = order)
    }
    else {
      if (order < 1) {
        warning("Specified spline order (", order, ") not valid!")
        order <- 1
      }
      nx <- length(x)
      nknots <- length(knots)
      d <- matrix(0, nx, nknots - 1)
      for (i in 1:(nknots - 1)) {
        inds <- which(x >= knots[i] & x < knots[i + 1])
        d[inds, i] <- 1
      }
    }
  }
  return(d)
}