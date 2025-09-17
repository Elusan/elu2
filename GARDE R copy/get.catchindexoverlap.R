get.catchindexoverlap <- function (inp) 
{
  y <- inp$obsC
  ty <- inp$timeC
  if (length(inp$obsI) > 0 & inp$nseasons == 1) {
    if (class(inp$obsI) == "list") {
      z <- inp$obsI[[1]]
      tz <- inp$timeI[[1]]
    }
    else {
      z <- inp$obsI
      tz <- inp$timeI
    }
    zinds <- na.omit(match(round(ty), round(tz)))
    z <- z[zinds]
    tz <- tz[zinds]
    yinds <- na.omit(match(round(tz), round(ty)))
    y <- y[yinds]
    ty <- ty[yinds]
    return(list(ty = ty, y = y, tz = tz, z = z))
  }
  else {
    return(NULL)
  }
}
