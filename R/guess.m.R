guess.m <- function (inp, all.return = FALSE) 
{
  meancatch <- mean(unlist(inp$obsC))
  flag <- FALSE
  if ("obsI" %in% names(inp)) {
    if (length(inp$obsI) > 0) {
      out <- get.catchindexoverlap(inp)
      if (!is.null(out)) {
        ty <- out$ty
        y <- out$y
        tz <- out$tz
        z <- out$z
        if (length(y) == length(z)) {
          flag <- TRUE
        }
      }
    }
  }
  if (flag) {
    x <- y/z
    mod0 <- lm(z ~ x)
    a <- mod0$coefficients[1]
    b <- mod0$coefficients[2]
    MSY <- -0.25 * a^2/b
    if (MSY <= 0 | summary(mod0)$r.squared < 0.3) {
      MSY <- meancatch
    }
    if (all.return) {
      Emsy <- -0.5 * a/b
      return(list(MSY = MSY, Emsy = Emsy, a = a, b = b, 
                  x = x, y = y, z = z, mod0 = mod0))
    }
    else {
      return(MSY)
    }
  }
  else {
    return(meancatch)
  }
}
