logistic_discrete_fishing_mortality = function( catch, abundance) {
  F =  catch / abundance
  F = 1.0 - F
  F = -log( F )
  F[ !is.finite(F)] = 0
  return(F)
}
