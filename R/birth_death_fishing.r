
birth_death_fishing = function( selection="stan_code", res=NULL, vn=NULL, sppoly=NULL, poly_match=NULL, time_match=NULL, breaksat=NULL, coastLayout=NULL, catches=NULL, wgts=NULL, ... ) {

  message(" ")
  message("todo :: use carstm estimated numbers as IOA")
  message("todo :: add zero-inflation")
  message(" ")

  if (grepl("parameters", selection) ) {


    return(params)
  }


  if (grepl("stan_code", selection) ) {

    # simple discrete ode form, spatially explicit (au), birth/death separated
    # n[au,t] = n [au,t−1] + b[au] n[au,t−1] (1 − b[au,t−1]  ) − c [au,t−1]

    out = "
      data {
        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int CAT[U,N];
        int IOA[U,N];
      }

      transformed data {
        real eps = 1e-12;
        real IOAmin[U];
        real IOAmax[U];
        real IOAmed[U];
        real IOAerror[U];
        for (u in 1:U) {
          IOAmax[u] = 1e-12;
          IOAmin[u] = 1e12;
          for (n in 1:N) {
            real test = IOA[u,n] + CAT[u,n] ;
            if (test > IOAmax[u]) IOAmax[u] = test;
            if (test <= IOAmin[u]) IOAmin[u] = test;
          }
          IOAmed[u] = (IOAmax[u] + IOAmin[u])/2.0 ;
          IOAerror[u] = IOAmax[u] * 0.05 ; // 5% error..  estim of sd
        }
      }

      parameters {
        real <lower=eps, upper=1.0> X[U,N];
        real <lower=0.0> g[U,N];
        real <lower=0.0, upper=0.5> m[U];
        real <lower=1.0> K[U];
        real <lower=0.0, upper=0.5>  psd[U];  // process error
        real <lower=0.0, upper=1.0> res[U,N];
        real <lower=0.0, upper=0.5>  ressd[U];  /////
        real <lower=0.0, upper=1.0>  gsd[U]; /////
        real <lower=0.5, upper=1.5>   q[U];
        real <lower=0.0, upper=0.5>   msd; ////
        real <lower=0.0, upper=0.2>   qsd;
      }

      model {
        psd ~ normal(0.0, 0.2) ;
        qsd ~ normal(0.0, 0.05);  // ranges from 0.5 to 1.5
        ressd ~ normal(0.0, 0.05) ;
        msd ~ normal(0.0, 0.05) ;
        gsd ~ normal(0.0, 0.05) ; // 0.1
        m ~ normal(0.0, msd);
        q ~ normal(1.0, qsd );

        for (u in 1:U) {
          g[u,] ~ normal( 0.0, gsd[u] );
          K[u] ~ normal( IOAmax[u], IOAerror[u] ) ;
          X[u,1] ~ normal( 0.5, psd[u] ) ;
          for (n in 1:N) {
            res[u,n] ~ normal( q[u] * (X[u,n] - CAT[u,n]/K[u] ), ressd[u] );  // no catch observation error  .. keep separate to force positive value
            IOA[u,n] ~ poisson( K[u] * res[u,n] ); // survey index observation model
            if (n < N) {
              X[u,n+1] ~ normal( X[u,n] * (1.0 + g[u,n] - m[u]*X[u,n] ) - CAT[u,n]/K[u], psd[u] ) ; //process model
            }
          }
        }
      }
    "
    return(out)
  }


  if (grepl("stan_code_testing_hurdle", selection) ) {

    # simple discrete ode form, spatially explicit (au), birth/death separated
    # n[au,t] = n [au,t−1] + b[au] n[au,t−1] (1 − b[au,t−1]  ) − c [au,t−1]

    out = "
      functions {
        int num_zero(int[] y) { // count no of zero-values
          int nz = 0;
          for (n in 1:size(y))
            if (y[n] == 0)
              nz += 1;
          return nz;
        }
      }

      data {
        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int CAT[U,N];
        int IOA[U,N];
      }

      transformed data {
        real eps = 1e-12;
        real IOAmin[U];
        real IOAmax[U];
        real IOAmed[U];
        real IOAerror[U];
        int<lower=0 > N0[U]; // no of zero-values in each au
        int<lower=0 > Ngt0[U];  // no of non-zero values in each au
        int<lower=1> IOA_nz[U], N - num_zero(y);  // values of IOA > 0
        for (u in 1:U) {
          N0[u] = num_zero( IOA[u,] );
          Ngt0[u] = N - N0[u];
          IOAmax[u] = 1e-12;
          IOAmin[u] = 1e12;
          for (n in 1:N) {
            real test = IOA[u,n] + CAT[u,n] ;
            if (test > IOAmax[u]) IOAmax[u] = test;
            if (test <= IOAmin[u]) IOAmin[u] = test;
          }
          IOAmed[u] = (IOAmax[u] + IOAmin[u])/2.0 ;
          IOAerror[u] = IOAmax[u] * 0.05 ; // 5% error..  estim of sd
        }
      }

      parameters {
        real <lower=eps, upper=1.0> X[U,N];
        real <lower=0.0> g[U,N];
        real <lower=0.0, upper=0.5> m[U];
        real <lower=1.0> K[U];
        real <lower=0.0, upper=0.5>  psd[U];  // process error
        real <lower=0.0, upper=1.0> res[U,N];
        real <lower=0.0, upper=0.5>  ressd[U];  /////
        real <lower=0.0, upper=1.0>  gsd[U]; /////
        real <lower=0.5, upper=1.5>   q[U];
        real <lower=0.0, upper=0.5>   msd; ////
        real <lower=0.0, upper=0.2>   qsd;
      }

      model {
        psd ~ normal(0.0, 0.2) ;
        qsd ~ normal(0.0, 0.05);  // ranges from 0.5 to 1.5
        ressd ~ normal(0.0, 0.05) ;
        msd ~ normal(0.0, 0.05) ;
        gsd ~ normal(0.0, 0.05) ; // 0.1
        m ~ normal(0.0, msd);
        q ~ normal(1.0, qsd );

        for (u in 1:U) {
          g[u,] ~ normal( 0.0, gsd[u] );
          K[u] ~ normal( IOAmax[u], IOAerror[u] ) ;
          X[u,1] ~ normal( 0.5, psd[u] ) ;
          for (n in 1:N) {
            res[u,n] ~ normal( q[u] * (X[u,n] - CAT[u,n]/K[u] ), ressd[u] );  // no catch observation error  .. keep separate to force positive value
            IOA[u,n] ~ poisson( K[u] * res[u,n] ); // survey index observation model
            if (n < N) {
              X[u,n+1] ~ normal( X[u,n] * (1.0 + g[u,n] - m[u]*X[u,n] ) - CAT[u,n]/K[u], psd[u] ) ; //process model
            }
          }
        }
      }
    "
    return(out)
  }


  if (grepl("stochastic_simulation", selection) ) {

    ellp = list(...)
    attach(ellp)

    Xdim = dim(X)
    nau = Xdim[1]
    ntu = Xdim[2]
    nsim = Xdim[3]

    if (!exists("ST")) ST = c(
      "@ -> b*X  -> X" ,
      "X -> d*(X/K)*X -> @",
      "X -> f*X -> catch"
    )

    if (!exists("SC"))  SC = c( "X",  "catch" )

    if (!exists("nthreads")) nthreads = 4
    if (!exists("nprojections")) nprojections = 30
    if (!exists("istart"))  istart = ntu

    nssims = min( 10, nsim )
    iss = sample.int( nsim, nssims )

    sim = array( NA, dim=c( nssims, length(SC), nprojections, nau ) )
    tspan = 1:nprojections

    for (iau in 1:nau) {

      # posterior subset
      aup = data.frame(
        X = as.integer( trunc(abundance[iss, iau, istart]) ),
        catch = as.integer( rep(0, length(X))),
        f = f[iss, iau, istart],
        b = g[iss, iau, istart],
        m = m[iss, iau, istart],
        K = as.integer( K[iss, iau] )
      )

      # simulate over each posterior subset
      for (i in 1:nssims) {
        sim[i,,iau]  = slot( run( model=mparse(
          transitions = ST,
          compartments = SC,
          gdata = c( K=aup$K[i], g=aup$g[i], m=aup$m[i], f=aup$f[i] ),
          u0 = aup[i, SC],
          tspan = tspan
          ), threads=nthreads
        ), "U")[]
      }

    }

    detach(ellp)
    out = list( sim=sim, iss =iss)

    return(out)
  }



  if (grepl("posteriors", selection) ) {

    ellp = list(...)
    attach(ellp)

    posteriors = rstan::extract( res )  # posteriors = mcmc posteriors from STAN

    Xdim = dim(posteriors$X)
    nsim = Xdim[1]
    nau = Xdim[2]
    ntu = Xdim[3]

    posteriors$numbers = posteriors$X[] * 0
    for ( i in 1:ntu ) posteriors$numbers[,,i] = posteriors$X[,,i] * posteriors$K[]  # X is fraction of K

    posteriors$fraction.fished = posteriors$X[] * 0
    for ( i in 1:nsim ) posteriors$fraction.fished[i,,] = catches[] / posteriors$numbers[i,,]

    posteriors$fishing.mortality = -log( pmax(1e-6, 1.0 - posteriors$fraction.fished ) )
    posteriors$fishing.mortality[ !is.finite(posteriors$fishing.mortality)] = 0

    # K = array( K, dim=c(dim(K), nsim) )
    # MSY    = r * K / 4.0 ; # maximum height of of the latent productivity (yield)
    # aMSY   = K /2.0 ; # numbers at MSY
    # FMSY   = 2.0 * MSY / K ; # fishing mortality at MSY

    # add other params computed post-mcmc from conditional distributions
    # convert num to biomass
    posteriors$biomass = posteriors$numbers[] * NA  # in kg
    for (i in 1:nsim) posteriors$biomass[i,,] = posteriors$numbers[i,,] * wgts[]

    detach(ellp)
    return( posteriors )
  }



  if (grepl("summary", selection)) {

    sims = lapply( apply( posteriors$biomass, 1,
      function(biom) {
        biom = as.matrix(biom)
        biom[!is.finite(biom)] = NA
        o = list()
        o$cfaall    = colSums( biom * sppoly$au_sa_km2/ sppoly$au_sa_km2/10^6, na.rm=TRUE )
        o$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea/ sppoly$au_sa_km2/10^6, na.rm=TRUE )
        o$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea/ sppoly$au_sa_km2/10^6, na.rm=TRUE )
        o$cfa23     = colSums( biom * sppoly$cfa23_surfacearea/ sppoly$au_sa_km2/10^6, na.rm=TRUE )
        o$cfa24     = colSums( biom * sppoly$cfa24_surfacearea/ sppoly$au_sa_km2/10^6, na.rm=TRUE )
        o$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea/ sppoly$au_sa_km2/10^6, na.rm=TRUE )
        return(o)
      }
    ), data.frame )

    RES = data.frame( index = 1:dim(posteriors$biomass)[3] )

    RES$cfaall = apply( sapply( sims, function(x) x[,"cfaall"], simplify =TRUE), 1, mean)
    RES$cfaall_sd = apply( sapply( sims, function(x) x[,"cfaall"], simplify =TRUE), 1, sd )
    RES$cfaall_median = apply( sapply( sims, function(x) x[,"cfaall"], simplify =TRUE), 1, median )
    RES$cfaall_lb = apply( sapply( sims, function(x) x[,"cfaall"], simplify =TRUE), 1, quantile, probs=0.025 )
    RES$cfaall_ub = apply( sapply( sims, function(x) x[,"cfaall"], simplify =TRUE), 1, quantile, probs=0.975 )

    RES$cfanorth = apply( sapply( sims, function(x) x[,"cfanorth"], simplify =TRUE), 1, mean )
    RES$cfanorth_sd = apply( sapply( sims, function(x) x[,"cfanorth"], simplify =TRUE), 1, sd )
    RES$cfanorth_median = apply( sapply( sims, function(x) x[,"cfanorth"], simplify =TRUE), 1, median )
    RES$cfanorth_lb = apply( sapply( sims, function(x) x[,"cfanorth"], simplify =TRUE), 1, quantile, probs=0.025 )
    RES$cfanorth_ub = apply( sapply( sims, function(x) x[,"cfanorth"], simplify =TRUE), 1, quantile, probs=0.975 )

    RES$cfasouth = apply( sapply( sims, function(x) x[,"cfasouth"], simplify =TRUE), 1, mean )
    RES$cfasouth_sd = apply( sapply( sims, function(x) x[,"cfasouth"], simplify =TRUE), 1, sd )
    RES$cfasouth_median = apply( sapply( sims, function(x) x[,"cfasouth"], simplify =TRUE), 1, median )
    RES$cfasouth_lb = apply( sapply( sims, function(x) x[,"cfasouth"], simplify =TRUE), 1, quantile, probs=0.025 )
    RES$cfasouth_ub = apply( sapply( sims, function(x) x[,"cfasouth"], simplify =TRUE), 1, quantile, probs=0.975 )

    RES$cfa23 = apply( sapply( sims, function(x) x[,"cfa23"], simplify =TRUE), 1, mean )
    RES$cfa23_sd = apply( sapply( sims, function(x) x[,"cfa23"], simplify =TRUE), 1, sd )
    RES$cfa23_median = apply( sapply( sims, function(x) x[,"cfa23"], simplify =TRUE), 1, median )
    RES$cfa23_lb = apply( sapply( sims, function(x) x[,"cfa23"], simplify =TRUE), 1, quantile, probs=0.025 )
    RES$cfa23_ub = apply( sapply( sims, function(x) x[,"cfa23"], simplify =TRUE), 1, quantile, probs=0.975 )

    RES$cfa24 = apply( sapply( sims, function(x) x[,"cfa24"], simplify =TRUE), 1, mean )
    RES$cfa24_sd = apply( sapply( sims, function(x) x[,"cfa24"], simplify =TRUE), 1, sd )
    RES$cfa24_median = apply( sapply( sims, function(x) x[,"cfa24"], simplify =TRUE), 1, median )
    RES$cfa24_lb = apply( sapply( sims, function(x) x[,"cfa24"], simplify =TRUE), 1, quantile, probs=0.025 )
    RES$cfa24_ub = apply( sapply( sims, function(x) x[,"cfa24"], simplify =TRUE), 1, quantile, probs=0.975 )

    RES$cfa4x = apply( sapply( sims, function(x) x[,"cfa4x"], simplify =TRUE), 1, mean )
    RES$cfa4x_sd = apply( sapply( sims, function(x) x[,"cfa4x"], simplify =TRUE), 1, sd )
    RES$cfa4x_median = apply( sapply( sims, function(x) x[,"cfa4x"], simplify =TRUE), 1, median )
    RES$cfa4x_lb = apply( sapply( sims, function(x) x[,"cfa4x"], simplify =TRUE), 1, quantile, probs=0.025 )
    RES$cfa4x_ub = apply( sapply( sims, function(x) x[,"cfa4x"], simplify =TRUE), 1, quantile, probs=0.975 )

    return(RES)

  }


  if (grepl("spplot", selection) ) {



    #  wrapper around spplot .. based on carstm_plot
    require(sp)

    ellp =list(...)   # if plotting, all ellipsis contents are expected to be sppolt args

    if (is.null(vn)) stop("must have a  vn object")
    if (is.null(spmatrix)) stop("must have a  spmatrix object")
    if (is.null(sppoly)) stop("must have an sppoly object")

    slot(sppoly, "data")[,vn] = NA

    # first index is spatial strata

    res_dim = dim( spmatrix )

    if (is.null(time_match)) {
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.8" )
    }

    if (is.null(poly_match)) poly_match = sppoly[["AUID"]]

    i_poly = match( poly_match, sppoly[["AUID"]] )  # should match exactly but in case a subset is sent as sppoly

    data_dimensionality = length( dim(spmatrix[[vn]]) )
    if (data_dimensionality==1) {
      slot(sppoly, "data")[, vn] = spmatrix[[vn]] [ i_poly ]  # year only
    }

    if (!is.null(time_match)) {
      n_indexes = length( time_match )
      if (data_dimensionality==2) {
        if (n_indexes==1) slot(sppoly, "data")[, vn] = spmatrix[[vn]] [ i_poly, time_match[[1]] ]  # year only
      }
      if (data_dimensionality==3) {
        if (n_indexes==1) slot(sppoly, "data")[, vn] = spmatrix[[vn]] [ i_poly, time_match[[1]] , ]  # year only
        if (n_indexes==2) slot(sppoly, "data")[, vn] = spmatrix[[vn]] [ i_poly, time_match[[1]], time_match[[2]] ] # year/subyear
      }
    }

    if (length(poly_match) > 1 ) {

      dev.new();

      if ( !exists("mypalette", ellp)) mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")

      if ( is.null(breaksat)) breaksat=interval_break(X=sppoly[[vn]], n=length(mypalette), style="quantile")

      if ( !exists("main", ellp ) )  ellp[["main"]]=vn
      ellp$obj = sppoly
      ellp$zcol=vn
      ellp$col.regions=mypalette
      ellp$at=breaksat
      ellp$col="transparent"

      do.call(spplot, ellp )

    }


  }


}
