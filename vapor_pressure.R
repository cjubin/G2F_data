get.ea <- function(rhmin, rhmax, tmin, tmax){
  esmn <- get.esmn(tmin) # other fun same as James suggested
  esmx <- get.esmx(tmax) 
  ea <- ((esmn * rhmax/100) + (esmx * rhmin/100)) / 2
  return(ea)
}

get.esmn <- function(tmin){
  esmn <- .6108 * exp((17.27 * tmin) / (tmin + 237.3))
  return(esmn)
}


get.esmx <- function(tmax){
  esmx <- .6108 * exp((17.27 * tmax) / (tmax + 237.3))
  return(esmx)
}

get.es <- function(tmin, tmax){
  esmn <- get.esmn(tmin)
  esmx <- get.esmx(tmax)
  es <- (esmn + esmx)/2
  return(es)
}
