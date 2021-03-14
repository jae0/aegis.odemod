

# stochastic simulation (birth-death) model of population dynamics
# this continues from 03.abundance estimation.carstm.R .. but uses alt methods ..
# it would be a replbacement of 05.*

# WARNING:: if the polygons are not optimal this can take many days to complete


# use numerical abundance for small area analysis .. n<1000 preferable ..
# 03.snowcrab... or comparable must has been completed .. ie carstm model

#--- NOTE --- add R1 and this becomes delay difference ...
#--- NOTE --- simulations not complete ,.. waiting for model form to be finalized ....

  require(SimInf)
  require(ggplot2)
  require(cmdstanr)

  require(aegis)
  require(bio.snowcrab)
  # require(aegis.odemod)
  loadfunctions("aegis.odemod")
  loadfunctions("ecomod")


  if (0) {

    year.assessment = 2020

    p = bio.snowcrab::snowcrab_parameters( 
      project_class="carstm", 
      assessment.years=2000:year.assessment, 
      areal_units_type="tesselation",
      carstm_model_label = "tesselation",   # default is the name of areal_units_type  
      selection = list(type = "number")
    )
    
  )


  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))

  # p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot

  sppoly = areal_units( p=p )  # will redo if not found
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  # plot(sppoly)
  # spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )

	color_scheme_set("brightblue") 


  standata = numerical_abundance_catch( p  )

  if (0) {
    subset = c(1:10)
    standata$IOA = standata$IOA[subset,]
    standata$CAT = floor( standata$CAT[subset,] )
    standata$U = length(subset)  #
  }

  mod = stan_initialize( stan_code=birth_death_fishing( "stan_code" ) )
  mod$compile()

  mod_data = standata[c("N", "U", "IOA", "CAT")]

  # see https://mc-stan.org/cmdstanr/reference/model-method-sample.html for more options
  fit = mod$sample(
    data = mod_data,
    iter_warmup = 2000,
    iter_sampling = 500,
    seed = 123,
    chains = 3,
    parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
    max_treedepth = 14,
    adapt_delta = 0.9,
    refresh = 500
  )



      if (0) {

        fit = fishery_model( p=p, DS="fit", tag=p$areal_units_type )  # to get samples
      
        print( fit, max_rows=30 )
        # fit$summary("K", "r", "q")
        
        fit$cmdstan_diagnose()
        fit$cmdstan_summary()
  
          # (penalized) maximum likelihood estimate (MLE) 
        fit_mle = mod$optimize(data =mod_data, seed = 123)
        fit_mle$summary( c("K", "r", "q") )

        mcmc_hist(fit$draws("K")) +
          vline_at(fit_mle$mle(), size = 1.5)

        # Variational Bayes  
        fit_vb = mod$variational(data =mod_data, seed = 123, output_samples = 4000)
        fit_vb$summary(c("K", "r", "q"))

        bayesplot_grid(
          mcmc_hist(fit$draws("K"), binwidth = 0.025),
          mcmc_hist(fit_vb$draws("K"), binwidth = 0.025),
          titles = c("Posterior distribution from MCMC", "Approximate posterior from VB")
        )

      }

 

      color_scheme_set("gray")
      mcmc_dens(fit$mcmc, regex_pars="K",  facet_args = list(nrow = 3, labeller = ggplot2::label_parsed ) ) + facet_text(size = 14 )   
      # mcmc_hist( fit$draws("K"))


   

  if (0){
    outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
    fn_root = paste( "birth_death_mcmc", "defaultgrid_smallarea", sep="." )
    fnres  = file.path(outdir, paste( "fit", fn_root, "rdata", sep=".") )
    save( fit , file=fnres, compress=TRUE )
    load( fnres )
  }

  posteriors = birth_death_fishing( selection="posteriors", fit=fit, wgts=standata$wgts, catches=standata$CAT )

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
  vn = "fishing.mortality";  spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly[,vn] = (spmatrix[,iyr]); plot(sppoly[,vn] )
  vn="biomass"; spmatrix = apply( posteriors[[vn]], c(2,3), median);  sppoly[,vn] = log(spmatrix[,iyr]); plot(sppoly[,vn] )

  vn="X"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly[,vn] = (spmatrix[,iyr]); plot(sppoly[,vn] )  # normalized numerical abundance
  vn="fraction.fished"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly[,vn] = (spmatrix[,iyr]); plot(sppoly[,vn] )
  vn="g"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly[,vn] = (spmatrix[,iyr]); plot(sppoly[,vn] )

  vn="g_ar"; spmatrix = apply( posteriors[[vn]], c(2,3), median); sppoly[,vn] = (spmatrix[,1]); plot(sppoly[,vn] )


  vn="K"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly[,vn] = (spmatrix); plot(sppoly[,vn] )
  vn="m"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly[,vn] = (spmatrix); plot(sppoly[,vn] )
  vn="theta"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly[,vn] = (spmatrix); plot(sppoly[,vn] )
  vn="q"; spmatrix = apply( posteriors[[vn]], c(2), median); sppoly[,vn] = (spmatrix); plot(sppoly[,vn] )



  # time_match = list(year="2000")
  # time_match = list(year="2000", dyear="0.8" )

  vn="biomass"; spmatrix = apply( posteriors[[vn]], c(2,3), median)
  
  birth_death_fishing( "spplot", p=p, spmatrix=spmatrix, vn=vn, time_match=time_match, 
    sp.layout=p$coastLayout, 
    # at = seq(... ) ,
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=vn
  )


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







