
birth_death_fishing_extract = function( posteriors, catch ) {

  Xdim = dim(posteriors$X)
  nsim = Xdim[3]

  abundance = posteriors$X[]
  for ( i in 1:nsim ) abundance[,,i] = abundance[,,i] * posteriors$K

  fraction.fished = abundance[] * 0
  for ( i in 1:nsim ) {
    fraction.fished[,,i] = catch / abundance[,,i]
  }
  fishing.mortality = -log( 1.0 - fraction.fished )
  fishing.mortality[ !is.finite(fishing.mortality)] = 0

  K = array(posteriors$K, dim=c(dim(K), nsim))

  MSY    = r * K / 4.0 ; # maximum height of of the latent productivity (yield)
  aMSY   = K /2.0 ; # abundance at MSY
  FMSY   = 2.0 * MSY / K ; # fishing mortality at MSY
  return(list(MSY=MSY, aMSY=aMSY, FMSY=FMSY))
}
