logistic_discrete_reference_points = function( r, K ) {
  if ( !identical( dim(r), dim(K) )) {
    K = array(K, dim=c(dim(K), dim(r)[3]))
  }
  MSY    = r * K / 4.0 ; # maximum height of of the latent productivity (yield)
  aMSY   = K /2.0 ; # abundance at MSY
  FMSY   = 2.0 * MSY / K ; # fishing mortality at MSY
  return(list(MSY=MSY, aMSY=aMSY, FMSY=FMSY))
}