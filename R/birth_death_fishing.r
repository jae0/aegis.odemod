
birth_death_fishing = function( selection="stan_code", res=NULL, vn=NULL, sppoly=NULL, poly_match=NULL, time_match=NULL, breaksat=NULL, coastLayout=NULL, ... ) {


  if (grepl("stan_code", selection) {

    # simple discrete ode form, spatially explicit (au), theta logistic, birth/death separated
    # n[au,t] = n [au,t−1] + b[au] n[au,t−1] (1 − b[au,t−1]^theta[au] ) − c [au,t−1]

    out = "
      data {
        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int CAT[U,N];
        int IOA[U,N];
        int lag;
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
        real <lower=0.0, upper=0.2>   qsd[U];
        real <lower=0.1, upper=10>   theta[U];
        real <lower=-1.0, upper =1.0>  g_ar[U, lag];
        real g_ark[U];
      }

      model {
        psd ~ normal(0.0, 0.2) ;
        qsd ~ normal(0.0, 0.1);  // ranges from 0.5 to 1.5
        ressd ~ normal(0.0, 0.1) ;
        msd ~ normal(0.0, 0.05) ;
        gsd ~ normal(0.0, 0.1) ;

        theta ~ lognormal(0.0, 0.1) ;
        m ~ normal(0.0, msd);
        q ~ normal(1.0, qsd );

        // AR(k=lag) model for g
        g_ark ~ normal( 0.0, 0.2 ); //, shrinks towards 0
        g[,1] ~ normal( 1.0, gsd );
        for (u in 1:U) {
          g_ar[u,] ~ normal( 0.0, 0.2 ); // autoregression (AR(k=lag)) ..  shrink to 0
          for (n in 2:N) {
            real gmu = 0.0;
            if ( lag > 0 && n > (lag+1) ) {
              gmu = g_ark[u];
              for ( k in 1:lag ) {
                gmu += g_ar[u, k] * g[u,n-k];
              }
            }
            g[u,n] ~ normal( gmu, gsd[u] );
          }
        }

        for (u in 1:U) {
          K[u] ~ normal( IOAmax[u], IOAerror[u] ) ;
          X[u,1] ~ normal( 0.5, psd[u] ) ;
          for (n in 1:N) {
            res[u,n] ~ normal( q[u] * (X[u,n] - CAT[u,n]/K[u] ), ressd[u] );  // no catch observation error  .. keep separate to force positive value
            IOA[u,n] ~ poisson( K[u] * res[u,n] ); // survey index observation model
            if (n < N) {
              X[u,n+1] ~ normal( X[u,n] * (1.0 + g[u,n] - m[u]*X[u,n]^theta[u] ) - CAT[u,n]/K[u], psd[u] ) ; //process model
            }
          }
        }
      }
    "
    return(out)
  }

  if (grepl("stochastic_simulation", selection) {

    ellp = list(...)
    attach(ellp)  

    Xdim = dim(X)
    nau = Xdim[1]
    ntu = Xdim[2]
    nsim = Xdim[3]

    if (!exists("ST")) ST = c(
      "@ -> b*X  -> X" ,
      "X -> d*(X/K)^(1+theta)*X -> @",
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
        theta = theta[iss, iau, istart], 
        K = as.integer( K[iss, iau] )
      )

      # simulate over each posterior subset
      for (i in 1:nssims) {
        sim[i,,iau]  = run( model=mparse(
          transitions = ST,
          compartments = SC,
          gdata = c( K=aup$K[i], g=aup$g[i], m=aup$m[i], theta=aup$theta[i], f=aup$f[i] ),
          u0 = aup[i, SC],
          tspan = tspan
          ), threads=nthreads
        )@U[]
      }

    }
    
    detach(ellp)
    out = list( sim=sim, iss =iss)

    return(out)
  }

  if (grepl("extract", selection) {
    
    ellp = list(...)
    attach(ellp)  

    Xdim = dim(X)
    nau = Xdim[1]
    ntu = Xdim[2]
    nsim = Xdim[3]

    abundance = X[]
    for ( i in 1:nsim ) abundance[,,i] = abundance[,,i] * K

    fraction.fished = abundance[] * 0
    for ( i in 1:nsim ) {
      fraction.fished[,,i] = catch / abundance[,,i]
    }
    fishing.mortality = -log( 1.0 - fraction.fished )
    fishing.mortality[ !is.finite(fishing.mortality)] = 0

    # K = array( K, dim=c(dim(K), nsim) )
    # MSY    = r * K / 4.0 ; # maximum height of of the latent productivity (yield)
    # aMSY   = K /2.0 ; # abundance at MSY
    # FMSY   = 2.0 * MSY / K ; # fishing mortality at MSY
    
    out = list( 
      abundance=abundance, 
      fishing.mortality=fishing.mortality, 
      fraction.fished=fraction.fished 
    )

    detach(ellp)
    return( out )
  }


  if (grepl("plot", selection) {
    #  wrapper around spplot .. based on carstm_plot
    require(sp)

    ellp =list(...)   # if plotting, all ellipsis contents are expected to be sppolt args

    if (is.null(vn)) stop("must have a  vn object")
    if (is.null(res)) stop("must have a  res object")
    if (is.null(sppoly)) stop("must have an sppoly object")

    sppoly@data[,vn] = NA

    # first index is spatial strata
    data_dimensionality = length( dim(res[[vn]]) )

    if (is.null(poly_match)) poly_match = match( res$AUID, sppoly[["AUID"]] )  # should match exactly but in case a subset is sent as sppoly

    if (data_dimensionality==1) {
      sppoly@data[, vn] = res[[vn]] [ poly_match ]  # year only
    }

    if (!is.null(time_match)) {
      n_indexes = length( time_match )
      if (data_dimensionality==2) {
        if (n_indexes==1) sppoly@data[, vn] = res[[vn]] [ poly_match, time_match[[1]] ]  # year only
      }
      if (data_dimensionality==3) {
        if (n_indexes==1) sppoly@data[, vn] = res[[vn]] [ poly_match, time_match[[1]] , ]  # year only
        if (n_indexes==2) sppoly@data[, vn] = res[[vn]] [ poly_match, time_match[[1]], time_match[[2]] ] # year/subyear
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
