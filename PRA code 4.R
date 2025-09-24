#Set correct folder
setwd("C:\\Sync\\Personal Sub-Projects\\Conferences\\PRA-BDT Course, September 2025\\Mark Section")

# Load packages; truncnorm needed for PRA, yarrr for transparency in plot shading
library(nimble)
library(truncnorm)
library(yarrr)

# Save graphical parameters
oldp <- par()

# Simple NIMBLE example - litters (originally from the BUGS manual)
littersCode <- nimbleCode({
  for (i in 1:G) {
     for (j in 1:N) {
        # likelihood (data model)
        r[i,j] ~ dbin(p[i,j], n[i,j])
        # latent process (random effects)
        p[i,j] ~ dbeta(a[i], b[i])
     }
     # prior for hyperparameters
     a[i] ~ dgamma(1, .001)
     b[i] ~ dgamma(1, .001)
   }
})
## data and constants as R objects
G <- 2
N <- 16
n <- matrix(c(13, 12, 12, 11, 9, 10,
              9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13,
              7, 5, 10, 7, 6, 10, 10, 10, 7), nrow = 2)
r <- matrix(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9,
     12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0),
     nrow = 2)

littersConsts <- list(G = G, N = N, n = n)
littersData <- list(r = r)
littersInits <- list( a = c(2, 2), b=c(2, 2) )

## create the NIMBLE model object
littersModel <- nimbleModel(littersCode,
          data = littersData, constants = littersConsts, inits = littersInits)

# compile the NIMBLE model object for computation
cLittersModel <- compileNimble(littersModel)
# Can view objects (variables, not constants) within the compiled model object,
# for example:
cLittersModel$p
cLittersModel$r
# Can also calculate derived quantities
cLittersModel$calculate('a')
cLittersModel$calculate('b')

# Now build the MCMC
# First, configure the MCMC
littersConf <- configureMCMC(littersModel, print = TRUE)
# We can add to the list of nodes/parameters
littersConf$addMonitors(c('a', 'b', 'p'))
# Finally, build the MCMC (i.e. implement the algorithm(s)) and compile into C++
littersMCMC <- buildMCMC(littersConf)
cLittersMCMC <- compileNimble(littersMCMC, project = littersModel)

# We can run the MCMC in R, just to show how slow it is...
niter <- 1000
nburn <- 100
set.seed(1)
inits <- function() {
      a <- runif(G, 1, 20)
      b <- runif(G, 1, 20)
      p <- rbind(rbeta(N, a[1], b[1]), rbeta(N, a[2], b[2]))
      return(list(a = a, b = b, p = p))
}
# This run will be slow (using the algorithm in R code)
print(system.time(samples.slow <- runMCMC(littersMCMC, niter = niter, nburnin = nburn,
                          inits = inits, nchains = 3, samplesAsCodaMCMC = TRUE)))
# This run should be much quicker (using the compiled C++ code)
print(system.time(samples <- runMCMC(cLittersMCMC, niter = niter, nburnin = nburn,
                          inits = inits, nchains = 3, samplesAsCodaMCMC = TRUE)))

# Look at sample traces - part of the MCMC diagnostics
plot(samples[ , "a[1]"], type="l", ylab="a[1]" )
plot(samples[ , "a[2]"], type="l", ylab="a[2]" )
plot(samples[ , "b[1]"], type="l", ylab="b[1]" )
plot(samples[ , "b[2]"], type="l", ylab="b[2]" )
# NB - plots suggest lack of convergence

# Can try many more simulations
niter <- 1000000
nburn <- 100000
print(system.time(samples <- runMCMC(cLittersMCMC, niter = niter, nburnin = nburn,
                          inits = inits, nchains = 3, samplesAsCodaMCMC = TRUE)))
plot(samples[ , "a[1]"], type="l", ylab="a[1]" )
# Better solution requires implementing block sampling, beyond scope of this session

# Can simulate from Data, but need to flag - and not during MCMC
littersModel$isData('r')
littersModel$isData('p')
littersModel$r
littersModel$p
# By default, this does not work
littersModel$simulate('r')
littersModel$simulate('p')
# No change to r
littersModel$r
# p now has values
littersModel$p
# Need to force simulation for r
littersModel$simulate('r', includeData = TRUE)



# NOW - the PRA example

# Teruel - rainfall data from Aragon
climate.data <- read.csv(file="teruel_monthly.csv")
climate.data <- climate.data[(climate.data$station=="TERUEL")|(climate.data$station=="TERUEL "),]
# Tree ring data
rings.data <- read.csv(file="gazol_dead_rw.csv")#
# Collapse tree ring data by annual means
rings.by.year <- aggregate(ring_width~year,rings.data,mean)
rings.by.year
# Collapse climate data by annual means
rain.by.year <- aggregate(rainfall~year,climate.data,mean)
rain.by.year
# Merge the data - only for years in common
combined.data <- merge(rings.by.year,rain.by.year,by="year")
combined.data

# Extract variables
Rainfall <- combined.data[,3]
Ring_Width <- combined.data[,2]

# Scatterplots, adapt graphical parameters
par(mar=c(5,4,0,2))
# Plot tree ring width versus mean monthly rainfall
plot(Rainfall,Ring_Width,xlab="Mean Monthly Rainfall (mm)",ylab="Mean Tree Ring Width (mm)",pch=16)
# Same again, but include origin on x-axis, showing risk of drought
plot(Rainfall,Ring_Width,xlab="Mean Monthly Rainfall (mm)",ylab="Mean Tree Ring Width (mm)",pch=16,xlim=c(0,45))
# Reset par
par(oldp)

# Number of Years
n <- length(Rainfall)

# NIMBLE model code
Model1.Code <- nimbleCode({
    lm.alpha  ~ dnorm ( 0, sd=100 )
    lm.beta           ~ dnorm ( 0, sd=100 )
    lm.tau            ~ dgamma( 0.01, 0.01 )
    lm.sigma         <- 1 / sqrt(lm.tau)
    xmean            <- mean( x[1:ndata] )
    for(i in 1:ndata){
      lm.mu[i] <- lm.alpha + lm.beta*x[i]
      z[i]      ~ dnorm( lm.mu[i], sd=lm.sigma )
    }
} )

# Set model constants - if covariates fixed, can include them here
Model1.Constants <- list( ndata=n, x=Rainfall )

# Set model response - "data" to be treated as random (not constant)
Model1.Data      <- list( z=Ring_Width )

# Set model as combination of three components above
Model1.Nimble    <- nimbleModel( Model1.Code, constants=Model1.Constants, data=Model1.Data )

# At this point, can investigate the node names
Model1.Nimble$getNodeNames()

# Compile the model (no link to C++ yet)
Model1.Comp      <- compileNimble( Model1.Nimble )

# Configure MCMC - sets samplers etc, we use defaults here
Model1.Conf  <- configureMCMC( Model1.Comp )
# Add some monitors, for the derived variable sigma
Model1.Conf$addMonitors( c("lm.sigma") )

# Build the MCMC -  we use defaults here
Model1.MCMC <- buildMCMC( Model1.Conf )

# Compile the MCMC into C++ - much faster than in R
Model1.MCMC.Comp <- compileNimble( Model1.MCMC )

### Run the MCMC - need to set burn-in length and subsequent chain length
nsims <- 2000 ; nburnin <- 500 ; niter <- nsims+nburnin
set.seed(1)
Model1.samples <- runMCMC( Model1.MCMC.Comp, nburnin=nburnin, niter=niter )

# A shortcut version compiles and runs the original compiled model in one call,
# although this is not sensible if you will be running the analysis multiple
# times
Model1.samples.shortcut <- nimbleMCMC( Model1.Comp, nburnin=nburnin, niter=niter,
    monitors=c("lm.alpha.centred", "lm.beta", "lm.tau", "lm.alpha", "lm.sigma") )

# You can run the MCMC without compilation - so in R, rather than C++. This is
# generally not recommended as it will be *much* slower...for illustration, use
# only a tenth of the chain length here...
Model1.samples.slow <- runMCMC( Model1.MCMC, nburnin=nburnin/10, niter=niter/10 )

# Plot the MCMC traces
par( mfrow=c(3,1), mar=c(0.1,5,0.1,0.1) )
plot(Model1.samples[ , "lm.alpha"], type="l", ylab="lm.alpha" )
plot(Model1.samples[ , "lm.beta" ], type="l", ylab="lm.beta"  )
plot(Model1.samples[ , "lm.sigma"], type="l", ylab="lm.sigma" )

# Basic summaries of the model parameters
summary(Model1.samples)

# Can further investigate samples from the posteriors:
summary(Model1.samples[ , "lm.alpha"])
par(oldp)
hist(Model1.samples[ , "lm.alpha"])
summary(Model1.samples[ , "lm.beta"])
hist(Model1.samples[ , "lm.beta"])
summary(Model1.samples[ , "lm.sigma"])
hist(Model1.samples[ , "lm.sigma"])

# Can look at cross-correlations:
cor(Model1.samples[ , "lm.alpha"],Model1.samples[ , "lm.beta"])
# Very high (negative) correlation! This explains the poor traces...

# For comparison, the standard linear regression from lm()
summary( lm(Ring_Width~Rainfall) )

# NIMBLE model code - note, centres covariate to improve MCMC sampling
Model2.Code <- nimbleCode({
    lm.alpha.centred  ~ dnorm ( 0, sd=100 )
    lm.alpha         <- lm.alpha.centred - lm.beta*xmean
    lm.beta           ~ dnorm ( 0, sd=100 )
    lm.tau            ~ dgamma( 0.01, 0.01 )
    lm.sigma         <- 1 / sqrt(lm.tau)
    xmean            <- mean( x[1:ndata] )
    for(i in 1:ndata){
      lm.mu[i] <- lm.alpha.centred + lm.beta*(x[i]-xmean)
      z[i]      ~ dnorm( lm.mu[i], sd=lm.sigma )
    }
} )
Model2.Constants <- list( ndata=n, x=Rainfall )
Model2.Data      <- list( z=Ring_Width )
Model2.Nimble    <- nimbleModel( Model2.Code, constants=Model2.Constants, data=Model2.Data )
Model2.Nimble$getNodeNames()
Model2.Comp      <- compileNimble( Model2.Nimble )
Model2.Conf  <- configureMCMC( Model2.Comp )
Model2.Conf$addMonitors( c("lm.alpha", "lm.sigma") )
Model2.MCMC <- buildMCMC( Model2.Conf )
Model2.MCMC.Comp <- compileNimble( Model2.MCMC )
nsims <- 2000 ; nburnin <- 500 ; niter <- nsims+nburnin
set.seed(2)
Model2.samples <- runMCMC( Model2.MCMC.Comp, nburnin=nburnin, niter=niter )
# Can look at cross-correlations:
cor(Model2.samples[ , "lm.alpha.centred"],Model1.samples[ , "lm.beta"])
# Much better!
par( mfrow=c(3,1), mar=c(0.1,5,0.1,0.1) )
plot(Model2.samples[ , "lm.alpha"], type="l", ylab="lm.alpha" )
plot(Model2.samples[ , "lm.beta" ], type="l", ylab="lm.beta"  )
plot(Model2.samples[ , "lm.sigma"], type="l", ylab="lm.sigma" )
summary(Model2.samples)
summary( lm(Ring_Width~Rainfall) )

# Reset graphical parameters
par(oldp)

# Now obtain V and R given the above model fit
# Set up the simulations
# How many thresholds to consider between the minimum and maximum
nthresholds  <- 100
# Number of simulations
nuncertainty <- 1000
thresholds   <- seq( 20, 45, len=nthresholds )
xmean        <- mean(Rainfall)
xsigma       <- sd(Rainfall)
zmean        <- mean(Ring_Width)
# Create empty structures for simulations
Rdiff <- Vdiff <- musimH <- musimNOTH <- xsimH <- zsimH <- xsimNOTH <- zsimNOTH <- array(NA,dim=c(nsims,nthresholds))
SDR <- LCIR <- UCIR <- SDV <- LCIV <- UCIV <- numeric(nthresholds)
calcSDzNOTH  <- calcSDzH <- numeric(nthresholds)
SDzNOTH <- LCIzNOTH <- UCIzNOTH <- SDzH <- LCIzH <- UCIzH <- EzH <- EzNOTH <- R <- V <- xweight <- numeric(nthresholds)
EzHj <- EzNOTHj <- Rj <- Vj <- numeric(nuncertainty)
lm.alpha <- Model2.samples[ , "lm.alpha" ]
lm.beta  <- Model2.samples[ , "lm.beta"  ]
lm.sigma <- Model2.samples[ , "lm.sigma" ]
# Run the uncertainty analysis for each threshold - note nsims comes from the MCMC
for(i in 1:nthresholds){
   xsimH[,i]      <- rtruncnorm( nsims, a=-Inf, b=thresholds[i], mean=xmean, sd=xsigma )
   xsimNOTH[,i]   <- rtruncnorm( nsims, a=thresholds[i], b=Inf , mean=xmean, sd=xsigma )
   xweight[i]     <- pnorm     ( thresholds[i]                 , mean=xmean, sd=xsigma )
   musimH[,i]     <- lm.alpha + lm.beta*xsimH[,i]
   musimNOTH[,i]  <- lm.alpha + lm.beta*xsimNOTH[,i]
   zsimH[,i]      <- rnorm( nsims, musimH[,i]   , sd=lm.sigma )
   zsimNOTH[,i]   <- rnorm( nsims, musimNOTH[,i], sd=lm.sigma )
   EzNOTH[i]      <- mean ( zsimNOTH[,i] )
   calcSDzNOTH[i] <- sd   ( zsimNOTH[,i] ) / sqrt( n*(1-xweight[i]) )
   EzH[i]         <- mean ( zsimH[,i]    )
   calcSDzH[i]    <- sd   ( zsimH[,i]    ) / sqrt( n*   xweight[i]  )
   V[i]           <- EzNOTH[i] - EzH[i]
   R[i]           <- EzNOTH[i] - weighted.mean( x=c(zsimH[,i],zsimNOTH[,i]), w=rep(c(xweight[i],1-xweight[i]),each=nsims) )
   # Now for the uncertainty of V and R
   for(j in 1:nuncertainty){
      xsim      <- rnorm( n, mean=xmean, sd=xsigma )
      xsimHj    <- xsim[ xsim <= thresholds[i] ]
      xsimNOTHj <- xsim[ xsim >  thresholds[i] ]
      # Deal with a computational problem if the full set n simulated values exceed a threshold, in either direction
      if( length(xsimHj)==0 ){
         xsimHj <- thresholds[i]
      }
      if( length(xsimNOTHj)==0 ){
         xsimNOTHj <- thresholds[i]
      }
      # If the sample size is larger than the number of simulations, need to sample with replacement
      if(nsims >= n){
         H.samples    <- sample( 1:nsims, length(xsimHj)    )
         NOTH.samples <- sample( 1:nsims, length(xsimNOTHj) )
      }else{
         H.samples    <- sample( 1:nsims, length(xsimHj)   , replace=TRUE )
         NOTH.samples <- sample( 1:nsims, length(xsimNOTHj), replace=TRUE )
      }
      musimHj    <- lm.alpha[ H.samples]    + lm.beta[H.samples]   *xsimHj
      musimNOTHj <- lm.alpha[ NOTH.samples] + lm.beta[NOTH.samples]*xsimNOTHj
      zsimHj     <- rnorm( length(musimHj), musimHj   , sd=lm.sigma[H.samples]    )
      zsimNOTHj  <- rnorm( length(musimNOTHj), musimNOTHj, sd=lm.sigma[NOTH.samples] )
      EzNOTHj[j] <- mean( zsimNOTHj )
      EzHj[j]    <- mean( zsimHj    )
      Vj[j]      <- EzNOTHj[j] - EzHj[j]
      Rj[j]      <- EzNOTHj[j] - mean(c(zsimHj,zsimNOTHj))
   }
   SDzNOTH[i]  <- sd( EzNOTHj )
   LCIzNOTH[i] <- quantile( EzNOTHj, 0.025 )
   UCIzNOTH[i] <- quantile( EzNOTHj, 0.975 )
   SDzH[i]     <- sd( EzHj )
   LCIzH[i]    <- quantile( EzHj   , 0.025 )
   UCIzH[i]    <- quantile( EzHj   , 0.975 )
   SDV[i]      <- sd( Vj )
   LCIV[i]     <- quantile( Vj     , 0.025 )
   UCIV[i]     <- quantile( Vj     , 0.975 )
   SDR[i]      <- sd( Rj )
   LCIR[i]     <- quantile( Rj     , 0.025 )
   UCIR[i]     <- quantile( Rj     , 0.975 )
}

# Plot the PRA results, with uncertainty
par( mfrow=c(1,2) )
# Plot conditional expectations for Ring Width as a function of Rainfall thresholds
plot( thresholds, EzH, type="l", xlab="Threshold", ylab="Cond Expectation of Z",
         ylim=range(c(LCIzNOTH,LCIzH,UCIzNOTH,UCIzH)))
# Uncertainty intervals - shading (pink for situation with hazard)
polygon( c(thresholds,thresholds[nthresholds:1]),
         c(LCIzH,UCIzH[nthresholds:1]),
         col=transparent("pink",0.5), border=NA )
# Uncertainty intervals - shading (grey for situation without hazard)
polygon( c(thresholds,thresholds[nthresholds:1]), c(LCIzNOTH,UCIzNOTH[nthresholds:1]),
         col=transparent("lightgrey",0.5), border=NA )
abline( h=zmean, col="blue", lty=3 )
lines( thresholds, EzH   , col="red"   )
lines( thresholds, EzNOTH, col="black" )
legend( "bottomright", legend=c("E[z|¬H]","E[z]","E[z|H]"),
         col=c("black","blue","red"), lty=c(1,3,1), cex=0.75 )
# Plot V and R for the PRA as functions of Rainfall thresholds
plot( thresholds, V, col="black", type="l",
         ylim=range(c(LCIV,UCIV,LCIR,UCIR)),
         xlab="Threshold", ylab="" )
polygon(c(thresholds,thresholds[nthresholds:1]),
         c(LCIV,UCIV[nthresholds:1]),
         col=transparent("lightgrey",0.5), border=NA )
polygon(c(thresholds,thresholds[nthresholds:1]),
         c(LCIR,UCIR[nthresholds:1]),
         col=transparent("pink",0.5), border=NA )
lines( thresholds, V, col="black" )
lines( thresholds, R, col="red"   )
legend( "bottomright", legend=c("V","R"), col=c("black","red"), lty=1, cex=0.75 )
