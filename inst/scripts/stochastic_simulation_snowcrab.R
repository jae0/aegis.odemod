

# stochastic simulation (birth-death) model of population dynamics

# use numerical abundance for small area analysis .. n<1000 preferable ..
# 03.snowcrab... or comparable must has been completed .. ie carstm model


  require(SimInf)
  require(ggplot2)
  require(rstan)

  require(aegis)
  require(bio.snowcrab)
  # require(aegis.odemod)
  # loadfunctions("aegis.odemod")


  year.assessment = 2019


  p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=2004:year.assessment )  # landings recorded by position since 2004

  crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))
  sppoly = areal_units( p=p )  # will redo if not found
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  # plot(sppoly)
  # spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )


  # loadfunctions("aegis.odemod")

  standata = numerical_abundance_catch(p, lag=1 )
  if (0) {
    subset = c(1:10)
    standata$IOA = standata$IOA[subset,]
    standata$CAT = trunc( standata$CAT[subset,] )
    standata$U = length(subset)  #
  }

  stancode_compiled = rstan::stan_model( model_code=birth_death_fishing( "stan_code" ) )

  res = rstan::sampling(
    stancode_compiled,
    data=standata,
    control = list(adapt_delta=0.925, max_treedepth=14),
    iter=6000,
    warmup=5000,
    refresh = 100,
#   seed=123,
    chains=3
  )

  posteriors = rstan::extract( res )  # posteriors = mcmc posteriors from STAN
 
  # add other params computed post-mcmc from conditional distributions
  posteriors = c( 
    posteriors, 
    birth_death_fishing( "extract", 
      X=posteriors$X, 
      K=posteriors$K, 
      catch=standata$CAT 
    )
  )

  sims = birth_death_fishing( "stochastic_simulation", 
    X=posteriors$abundance, 
    K=posteriors$K, 
    catch=standata$CAT, 
    g=posteriors$g, 
    m=posteriors$m,
    f=posteriors$fishing.mortality 
  )
  

  if (0){

    outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
    
    fn_root = paste( "birth_death_mcmc", "defaultgrid_smallarea", sep="." )

    fnres  = file.path(outdir, paste( "res", fn_root, "rdata", sep=".") )
    save( res , file=fnres, compress=TRUE )

    fnsims  = file.path(outdir, paste( "sims", fn_root, "rdata", sep=".") )
    save( sims , file=fnres, compress=TRUE )

    load( fnres )
    load( fnsims )

    # ---------------------

    au = 1
    au = 2

    hist(posteriors$K[,au], "fd")

    plot( posteriors$g ~ posteriors$m )

    hist( posteriors$g - posteriors$m - posteriors$fishing.mortality)

    ny = dim(posteriors$fishing.mortality)[3]
    plot(0,0, xlim=c(0, ny ), ylim=c(0,1.25), na.rm=TRUE, type="n")
    for (sim in 1:dim(posteriors$fishing.mortality)[1]){
      lines( posteriors$fishing.mortality[sim,au, ]  ~  c(1:ny), col=alpha(au, 0.05)  )
    }

    ny = dim(posteriors$g)[3]
    plot( 0,0, xlim=c(0, ny ), ylim=c(0,1), na.rm=TRUE, type="n")
    for (sim in 1:dim(posteriors$g)[1]){
      lines( posteriors$g[sim,au, ] ~  c(1:ny), col=alpha(au, 0.05)  )
    }

    ny = dim(posteriors$m)[3]
    plot(0,0, xlim=c(0, ny ), ylim=c(0,2), na.rm=TRUE, type="n")
    for (sim in 1:dim(posteriors$m)[1]){
      lines( posteriors$m[sim,au, ] ~  c(1:ny), col=alpha(au, 0.05)  )
    }

    ny = dim(posteriors$r)[3]
    plot(0,0, xlim=c(0, ny ), ylim=c(0,2), na.rm=TRUE, type="n")
    for (sim in 1:dim(posteriors$r)[1]){
      lines( posteriors$r[sim,au, ] ~  c(1:ny), col=alpha(au, 0.05)  )
    }

    plot(0,0, xlim=range(posteriors$K[,au]), ylim=range(posteriors$r[,au,]), na.rm=TRUE, type="n")
    for (tu in 1:dim(posteriors$r)[3]){
      points( posteriors$r[,au,tu ] ~ posteriors$K[,au], col=alpha(tu, 0.1), pch=19, cex=0.6 )
    }


    plot( standata$CAT[au,] ~ c(1:ny), ylim=range(c(posteriors$cat[,au,] * posteriors$K[,au])))
    for (i in 1:nposts) {
      lines( posteriors$cat[i,au,] * posteriors$K[i,au] ~ c(1:ny), col=alpha("orange", 0.01) )
    }


    hist(posteriors$q[,au], "fd")

    hist(posteriors$qsd[], "fd")

    hist(posteriors$rsd[,au], "fd")

    hist(posteriors$psd[,au], "fd")

    hist(posteriors$theta[,au], "fd")

    hist(posteriors$thetasd[,au], "fd")


    nposts =dim(posteriors$X)[1]
    ny = dim(posteriors$X)[3]
    plot( 0,0, xlim=c(0, ny ), ylim=c(0, 1.1), type="n")
    for (i in 1:nposts) {
      lines( posteriors$X[ i, au, ] ~ c(1:ny), col=alpha("green", 0.01) )
    }
    lines( apply( posteriors$X[ , au, ], 2, median) ~ c(1:ny), col=alpha("black", 0.99) )



    
    ny = dim(posteriors$abundance)[3]
    nposts =dim(posteriors$abundance)[1]
    plot( 0,0, xlim=c(0, ny ), ylim=range( c(posteriors$abundance[ ,au, ], posteriors$K[ ,au], standata$CAT[au,]), na.rm=TRUE), type="n")
    for (i in 1:nposts) {
      lines( posteriors$abundance[ i, au, ] ~ c(1:ny), col=alpha("green", 0.01) )
      abline( h=posteriors$K[i , au]), col=alpha("orange", 0.1) )
    }
    lines( standata$IOA[ au, ] ~ c(1:ny), col=alpha("red", 0.99), lwd=4 )
    lines( standata$CAT[ au, ] ~ c(1:ny), col=alpha("magenta", 0.99), lwd=4 )
    lines( apply( posteriors$abundance[ , au, ], 2, median) ~ c(1:ny), col=alpha("darkgrey", 0.99), lwd=4 )
    abline( h=median( posteriors$K[ , au]), col=alpha("orange", 0.9), lwd=4 )


  }


  yrs0 = p$assessment.years
  yrs.last = max(yrs0) + 0.5
  ndata = length(yrs0)
  hdat = 1:ndata

  au = 1

  istart = ndata


  vn = 1
  plot( 0,0, xlim=c(0, nprojections+1), ylim=range(sim[,vn,], na.rm=TRUE), type="n")
  for (i in 1:nsims) {
    lines( sim[i,vn,] ~ c(1:nprojections), col=alpha("green", 0.25) )
  }


  vn = 2
  plot( 0,0, xlim=c(0, nprojections+1), ylim=range(sim[,vn,], na.rm=TRUE)/10, type="n")
  for (i in 1:nsims) {
    lines( diff(sim[i,vn,]) ~ c(2:nprojections), col=alpha("green", 0.25) )
  }

  vn = 3
  plot( 0,0, xlim=c(0, nprojections+1), ylim=range(sim[,vn,], na.rm=TRUE)/10, type="n")
  for (i in 1:nsims) {
    lines( diff(sim[i,vn,]) ~ c(2:nprojections), col=alpha("green", 0.25) )
  }


  au = 1
  ny = dim(abundance)[3]
  nposts = dim(abundance)[1]

  nyp = ny + nprojections
  plot( 0,0, xlim=c(0, nyp ), ylim=range( abundance[ ,au, ], na.rm=TRUE), type="n")
  for (i in 1:nposts) {
    lines( abundance[ i, au, ] ~ c(1:ny), col=alpha("green", 0.01) )
  }
  lines( standata$IOA[ au, ] ~ c(1:ny), col=alpha("red", 0.99) )
  lines( apply( abundance[ , au, ], 2, median) ~ c(1:ny), col=alpha("black", 0.99) )
  for (i in 1:nsims) {
    lines( sim[ i, vn, ] ~ c((ny+1):nyp), col=alpha("magenta", 0.1) )
  }
  lines( apply( sim[ , vn, ], 2, median) ~ c((ny+1):nyp), col=alpha("red", 0.9), lwd=5 )
  for (i in 1:nsims) {
    abline( h=posteriors$K[i, au] , col=alpha("blue", 0.1) )
  }
  abline( h=median( posteriors$K[, au]) , col=alpha("black", 0.9) )
