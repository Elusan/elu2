make.ffacvec <- function (inp, ffac) 
{
  if (!"ns" %in% names(inp)) 
    stop("inp needs to be a checked input list - use 'check.inp()'!")
  inp$ffacvec <- rep(1, inp$ns)
  inp$ffac <- ffac
  ind <- which(inp$time == inp$manstart)
  inp$ffacvec[ind] <- ffac
  inp$ffacvec <- inp$ffacvec + 1e-08
  return(inp)
}