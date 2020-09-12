

# stochastic simulation (birth-death) model of population dynamics

# use numerical abundance for small area analysis .. n<1000 preferable ..
# 03.snowcrab... or comparable must has been completed .. ie carstm model

--- NOTE --- add R1 and this becomes delay difference ...
--- NOTE --- simulations not complete ,.. waiting for model form to be finalized ....

  require(SimInf)
  require(ggplot2)
  require(rstan)

  require(aegis)
  require(bio.snowcrab)
  # require(aegis.odemod)
  loadfunctions("aegis.odemod")
  loadfunctions("ecomod")


  year.assessment = 2019   # NOTE: for 4X, the season 2019-2020 -> 2019

  p = bio.snowcrab::load.environment(
    year.assessment=year.assessment,
    assessment_years = 2004:year.assessment,
    vars.tomodel="R0.mass",
    modeldir = project.datadirectory("bio.snowcrab", "modelled", "testing" ),  ## <--- important: alter save location for this  .. default is "*/modelled"

    areal_units_timeperiod="default",
    areal_units_overlay="snowcrab_managementareas",

    areal_units_resolution_km = 1,
    #areal_units_resolution_km = 10,
    #areal_units_resolution_km = 25,
    # areal_units_source = "lattice",
    areal_units_source = "snowcrab_polygons_tesselation",
    areal_units_constraint_nmin= 10,
    # using alt biomass index estmates
    # carstm_model_label="production",
    # carstm_model_label = paste( "testing", areal_units_source, areal_units_resolution_km, areal_units_constraint_nmin, sep="_" ),
    # carstm_model_label="testing_lattice_10_10",
    # carstm_model_label="testing_lattice_25_3",
    # carstm_model_label="testing_snowcrab_polygons_tesselation_1_10",
    carstm_model_label="testing_snowcrab_polygons_tesselation_1_10_zeroinflated",
    # carstm_model_label="testing_snowcrab_polygons_tesselation_1_4",
    carstm_modelengine = "inla",
    libs="carstm"
  )


  crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))

  # p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot

  sppoly = areal_units( p=p )  # will redo if not found
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  # plot(sppoly)
  # spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )


  standata = numerical_abundance_catch( p  )

  if (0) {
    subset = c(1:10)
    standata$IOA = standata$IOA[subset,]
    standata$CAT = trunc( standata$CAT[subset,] )
    standata$U = length(subset)  #
  }

  stancode_compiled = rstan::stan_model( model_code=birth_death_fishing( "stan_code" ) )

  res = rstan::sampling(
    stancode_compiled,
    data=standata[c("N", "U", "IOA", "CAT")],
    control = list(adapt_delta=0.95, max_treedepth=16),
    iter=5000,
    warmup=4000,
    refresh = 100,
#   seed=123,
    chains=3
  )


  if (0){
    outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
    fn_root = paste( "birth_death_mcmc", "defaultgrid_smallarea", sep="." )
    fnres  = file.path(outdir, paste( "res", fn_root, "rdata", sep=".") )
    save( res , file=fnres, compress=TRUE )
    load( fnres )
  }

  posteriors = birth_death_fishing( selection="posteriors", res=res, wgts=standata$wgts, catches=standata$CAT )

  res_summ = birth_death_fishing( selection="summary", posterior=posteriors, sppoly=sppoly )
  res_summ$yrs = p$yrs

  vn = "cfaall"
  vn = "cfanorth"
  vn = "cfasouth"
  vn = "cfa4x"

  plot(res_summ[,vn] ~ res_summ$yrs, type="b", col="darkorange", pch=20 )
  lines(res_summ[, paste(vn, "lb", sep="_")] ~ res_summ$yrs, col="slateblue", lty="dotted")
  lines(res_summ[, paste(vn, "ub", sep="_")] ~ res_summ$yrs, col="slateblue", lty="dotted")


  iyr = which(p$yrs== 2019)
  vn = "fishing.mortality";  spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly@data[,vn] = (spmatrix[,iyr]); spplot(sppoly, vn)
  vn="biomass"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly@data[,vn] = (spmatrix[,iyr]); spplot(sppoly, vn)
  vn="X"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly@data[,vn] = (spmatrix[,iyr]); spplot(sppoly, vn)  # normalized numerical abundance
  vn="fraction.fished"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly@data[,vn] = (spmatrix[,iyr]); spplot(sppoly, vn)
  vn="g"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly@data[,vn] = (spmatrix[,iyr]); spplot(sppoly, vn)

  vn="g_ar"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly@data[,vn] = (spmatrix[,1]); spplot(sppoly, vn)


  vn="K"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly@data[,vn] = (spmatrix); spplot(sppoly, vn)
  vn="m"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly@data[,vn] = (spmatrix); spplot(sppoly, vn)
  vn="theta"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly@data[,vn] = (spmatrix); spplot(sppoly, vn)
  vn="q"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly@data[,vn] = (spmatrix); spplot(sppoly, vn)



  # time_match = list(year="2000")
  # time_match = list(year="2000", dyear="0.8" )

  spmatrix = apply( posteriors[[vn]], c(2,3), median)

  birth_death_fishing( "spplot", p=p, spmatrix=spmatrix, vn=vn, time_match=time_match, sppoly=sppoly, sp.layout=p$coastLayout )


# ------------simulations

  sims = birth_death_fishing( "stochastic_simulation",
    posteriors=posteriors$numbers,
    K=posteriors$K,
    catch=standata$CAT,
    g=posteriors$g,
    m=posteriors$m,
    f=posteriors$fishing.mortality
  )


  if (0){
    outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
    fnsims  = file.path(outdir, paste( "sims", fn_root, "rdata", sep=".") )
    save( sims , file=fnres, compress=TRUE )
    load( fnsims )
  }


  isu = 1 # X state
  itu = 1 # time units
  spmatrix = reformat_to_array(
    input=apply( sims[[,1,1]], 1, median),
    matchfrom = list( AUID=AUID, yr_factor=yr_factor),
    matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
  )

  # time_match = list(year="2000")
  # time_match = list(year="2000", dyear="0.8" )

  birth_death_fishing( "plot", p=p, spmatrix=sim_array, vn=vn, time_match=time_match, sppoly=sppoly, sp.layout=p$coastLayout )






  if (0){

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


ny = dim(posteriors$numbers)[3]
nposts =dim(posteriors$numbers)[1]
nau=dim(posteriors$numbers)[2]


for (au in 1:nau) {
  print(au)

    nposts =dim(posteriors$X)[1]
    ny = dim(posteriors$X)[3]
    plot( 0,0, xlim=c(0, ny ), ylim=c(0, 1.1), type="n")
    for (i in 1:nposts) {
      lines( posteriors$X[ i, au, ] ~ c(1:ny), col=alpha("green", 0.01) )
    }
    lines( apply( posteriors$X[ , au, ], 2, median) ~ c(1:ny), col=alpha("black", 0.99) )

  u = readline()
}


# nposts=500

for (au in 1:nau) {
  print(au)

  plot( 0,0, xlim=c(0, ny ), ylim=range( c(posteriors$numbers[ ,au, ], posteriors$K[ ,au], standata$CAT[au,]), na.rm=TRUE), type="n")
  for (i in 1:nposts) {
    lines( posteriors$numbers[ i, au, ] ~ c(1:ny), col=alpha("green", 0.25) )
    abline( h=posteriors$K[i , au], col=alpha("orange", 0.1) )
  }
  lines( standata$IOA[ au, ] ~ c(1:ny), col=alpha("red", 0.99), lwd=4 )
  lines( standata$CAT[ au, ] ~ c(1:ny), col=alpha("magenta", 0.99), lwd=4 )
  lines( apply( posteriors$numbers[ , au, ], 2, median) ~ c(1:ny), col=alpha("darkgrey", 0.99), lwd=4 )
  abline( h=median( posteriors$K[ , au]), col=alpha("blue", 0.9), lwd=4 )

   u = readline()

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
  ny = dim(numbers)[3]
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







