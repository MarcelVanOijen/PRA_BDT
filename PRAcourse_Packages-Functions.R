# Probabilistic Risk Analysis and Bayesian Decision Theory
# - R-packages and functions
# Marcel van Oijen [Independent researcher, Edinburgh, UK - VanOijenMarcel@gmail.com]
# Mark Brewer [BioSS Office, The James Hutton Institute, Craigiebuckler, Aberdeen AB15 8QH]

  library(geodata)
  library(MCMCpack)
  library(mvtnorm)
  library(nimble)
  library(scales)
  library(terra)

# Data

### Sparse Dataset: $\texttt{\{x\_4,z\_4\}}$
  set.seed(1)
  n_4  <- 4 ; Sigma_4 <- matrix( c(1,0.5,0.5,1), nrow=2 )
  xz_4 <- rmvnorm( n_4, c(0,0), Sigma_4 ) %*% chol( Sigma_4 )
  xz_4 <- sweep( xz_4, 2, colMeans(xz_4) )
  xz_4 <- xz_4 %*% solve( chol(cov(xz_4)) ) %*% chol(Sigma_4)
  x_4  <- xz_4[,1]       ; z_4 <- xz_4[,2]
  m_4  <- colMeans(xz_4) ; S_4 <- cov(xz_4)
  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  plot( x_4, z_4, xlim=range(xz_4), ylim=range(xz_4), xlab="x", ylab="z" )

### A Collection of Linear Datasets: $\texttt{l\_xz.L}$
  set.seed(1)
  mu  <- c(0,0) ; Sigma  <- diag(1,2)          ; Sigma[1,2] <- Sigma[2,1] <- 0.5
  n_d <- 1e2    ; l_xz.L <- vector("list",n_d) ; n <- 1e3
  for(d in 1:n_d) { l_xz.L[[d]] <- rmvnorm( n, mu, Sigma ) }
  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  plot( l_xz.L[[1]][,1], l_xz.L[[1]][,2], main="dataset 1", xlab="x", ylab="z" )
  plot( l_xz.L[[2]][,1], l_xz.L[[2]][,2], main="dataset 2", xlab="x", ylab="z" )

### A Sequence of Linear Datasets (n linearly increasing)
  set.seed(1)
  mu  <- c(0,0) ; Sigma  <- diag(1,2) ; Sigma[1,2] <- Sigma[2,1] <- 0.5
  n_d <- 1e2    ; l_xz.L2 <- vector("list",n_d)
  for(d in 1:n_d) { l_xz.L2[[d]] <- rmvnorm( d, mu, Sigma ) }
  par( mfrow=c(1,3), mar=c(4,4,1,4) )
  plot( l_xz.L2[[1  ]][,1], l_xz.L2[[1  ]][,2], main="dataset 1",
        xlab="x", ylab="z" )
  plot( l_xz.L2[[30 ]][,1], l_xz.L2[[30 ]][,2], main="dataset 30",
        xlab="x", ylab="z" )
  plot( l_xz.L2[[n_d]][,1], l_xz.L2[[n_d]][,2], main=paste("dataset",n_d),
        xlab="x", ylab="z" )

### A Sequence of Linear Datasets (n exponentially increasing)
  set.seed(1)
  mu  <- c(0,0) ; Sigma  <- diag(1,2) ; Sigma[1,2] <- Sigma[2,1] <- 0.5
  n_d <- 14     ; l_xz.L3 <- vector("list",n_d)
  for(d in 1:n_d) { l_xz.L3[[d]] <- rmvnorm( 2^d, mu, Sigma ) }
  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  plot( l_xz.L3[[1  ]][,1], l_xz.L3[[1  ]][,2], main="dataset 1",
        xlab="x", ylab="z" )
  plot( l_xz.L3[[n_d]][,1], l_xz.L3[[n_d]][,2], main=paste("dataset",n_d),
        xlab="x", ylab="z" )

### A Collection of Nonlinear Datasets: $\texttt{l\_xz.NL}$
  set.seed(1)
  n_d <- 1e2 ; l_xz.NL <- vector("list",n_d) ; n <- 1e3 ; sz <- 0.1
  for(d in 1:n_d) {
    x <- runif( n, 0, 3 ) ; ez <- rnorm( n, 0, sz) ; z <- 1-exp(-x) + ez
    l_xz.NL[[d]] <- cbind(x,z) }
  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  plot( l_xz.NL[[1]][,"x"], l_xz.NL[[1]][,"z"], main="dataset 1", xlab="x", ylab="z" )
  plot( l_xz.NL[[2]][,"x"], l_xz.NL[[2]][,"z"], main="dataset 2", xlab="x", ylab="z" )

### A Sequence of Nonlinear Datasets (n exponentially increasing)
  set.seed(1)
  n_d <- 14 ; l_xz.NL3 <- vector("list",n_d) ; sz <- 0.1
  for(d in 1:n_d) {
    x <- runif( 2^d, 0, 3 ) ; ez <- rnorm( 2^d, 0, sz) ; z <- 1-exp(-x) + ez
    l_xz.NL3[[d]] <- cbind(x,z) }
  par( mfrow=c(1,3), mar=c(5,4,2,4) )
  plot( l_xz.NL3[[1  ]][,1], l_xz.NL3[[1  ]][,2], main="dataset 1",
        xlab="x", ylab="z", ylim=c(-0.5,1.5) )
  plot( l_xz.NL3[[6  ]][,1], l_xz.NL3[[6  ]][,2], main="dataset 6",
        xlab="x", ylab="z", ylim=c(-0.5,1.5) )
  plot( l_xz.NL3[[n_d]][,1], l_xz.NL3[[n_d]][,2], main=paste("dataset",n_d),
        xlab="x", ylab="z", ylim=c(-0.5,1.5) )

### German Forestry Data: $\texttt{x\_r3, z\_Fs, z\_Q, z\_Pa, z\_Ps}$
  d.data <- "data_Germany/"
  load( paste0( d.data, "xz_Germany.RData" ) )
  par ( mfrow=c(2,2), mar=c(4,4,1,2) )
  plot( x_r3, z_Fs, main="Fs", xlab=""         , ylab="mortality (%)" )
  plot( x_r3, z_Q , main="Q" , xlab=""         , ylab=""              )
  plot( x_r3, z_Pa, main="Pa", xlab="rain (mm)", ylab="mortality (%)" )
  plot( x_r3, z_Ps, main="Ps", xlab="rain (mm)", ylab=""              )
  
### Trivariate Gaussian Dataset
  set.seed(1)
  n_G3  <- 1e3 ; Sigma_G3 <- matrix(0.5,nrow=3,ncol=3) ; diag(Sigma_G3) <- 1
  xz_G3 <- rmvnorm( n_G3, c(0,0,0), Sigma_G3 )
  pairs(xz_G3)

  PRA <- function( x, z, thr=0 ) {
    n       <- length(x)      ; H         <- which(x < thr) ; n_H <- length(H)
    Ez_H    <- mean( z[ H] )  ; s_Ez_H    <- sqrt( var(z[ H]) /    n_H  )
    Ez_NotH <- mean( z[-H] )  ; s_Ez_NotH <- sqrt( var(z[-H]) / (n-n_H) )
    pH      <- n_H / n        ; V         <- Ez_NotH - Ez_H ; R  <- pH * V
    s_pH    <- sqrt( pH*(1-pH) / n )
    s_V     <- sqrt( s_Ez_H^2 + s_Ez_NotH^2 )
    s_R     <- sqrt( s_pH^2*s_V^2 + s_pH^2*V^2 + pH^2*s_V^2 )
    return( c(pH=pH,V=V,R=R,s_pH=s_pH,s_V=s_V,s_R=s_R) )
  }

  x <- c( 200, 600 )           ; z <- c(  70,  90 )
  PRA( x, z, 500 )
  x <- c( 200, 400, 600, 800 ) ; z <- c(  70,  80,  90, 100 )
  PRA( x, z, 500 )
  
  plotPRA <- function( pra ) {
    # This function requires 'pra' to be a vector c(pH,V,R,s_pH,s_V,s_R)
    mpH <- pra[1] ; mVR <- pra[2:3] ; spH <- pra[4] ; sVR <- pra[5:6]
    par( mfrow=c(1,2) )
    bpH  <- barplot( mpH, xlim=c(0,1), ylim=range(0,mpH+2*spH),
                     col="red", width=0.5 )
    segments( bpH, mpH-spH, bpH, mpH+spH ) ; ew <- bpH / 4
    segments( bpH-ew, mpH-spH, bpH+ew, mpH-spH )
    segments( bpH-ew, mpH+spH, bpH+ew, mpH+spH )
    bVR  <- barplot( mVR, ylim=range(0,mVR-2*sVR,mVR+2*sVR), beside=T,
                     col=c("cyan","yellow") )
    segments( bVR, mVR-sVR, bVR, mVR+sVR ) ; ew <- (bVR[2,1]-bVR[1,1]) / 4
    segments( bVR-ew, mVR-sVR, bVR+ew, mVR-sVR )
    segments( bVR-ew, mVR+sVR, bVR+ew, mVR+sVR )
  }
  plotPRA( PRA( x, z, 500 ) )

  PRAm <- function( x, z, thr=-1:1 ) {
    n   <- length(x) ; n_thr <- length(thr)
    H   <- vector("list",n_thr)
    n_H <- pH <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_thr)
    H[[1]]  <- which( x < thr[1] ) ; n_H[1] <- length(H[[1]])
    for(i in 2:n_thr) { H[[i]] <- which( thr[i-1] <= x & x < thr[i])
                       n_H[i] <- length(H[[i]]) } ; n_NotH <- n - sum(n_H)
    H.all   <- which( x < thr[n_thr] )
    pH      <- n_H / n             ; s_pH      <- sqrt( pH*(1-pH) / n )
    Ez_NotH <- mean( z[-H.all] )   ; s_Ez_NotH <- sqrt( var(z[-H.all] ) / n_NotH )   
    for(i in 1:n_thr) {
      Ez_Hi <- mean( z[ H[[i]] ] ) ; s_Ez_Hi <- sqrt( var(z[ H[[i]]]) /  n_H[i] )
      V[i]  <- Ez_NotH - Ez_Hi     ; s_V[i]  <- sqrt( s_Ez_NotH^2 + s_Ez_Hi^2 ) }
    R       <- pH * V
    s_R     <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2 )
    R.sum   <- sum(R) ; pH.sum <- sum(pH) ; V.wsum <- R.sum / pH.sum
    return( list( sum = c( pH.sum=pH.sum, V.wsum=V.wsum, R.sum=R.sum ),
                  seq = cbind( thr, pH, V, R, s_pH, s_V, s_R ) ) )
  }

  x_NL    <- l_xz.NL[[1]][,1]    ; z_NL <- l_xz.NL[[1]][,2]
  thr.seq <- seq(0.2,2.8,by=0.2)
  plot(x_NL,z_NL) ; abline(v=thr.seq)
  PRAm_NL <- PRAm( x_NL, z_NL, thr.seq )$seq
  par( mfrow=c(2,3), mar=c(5,3,2,0) )
  xtxt <- "Upper bound of interval" ; nms <- thr.seq
  barplot( PRAm_NL[,"pH"], main="p[H]", xlab=xtxt, ylab="",
           ylim=c(0,0.1), names=nms )
  barplot( PRAm_NL[,"V" ], main="V"   , xlab=xtxt, ylab="",
           ylim=c(0,1  ), names=nms )
  barplot( PRAm_NL[,"R" ], main="R"   , xlab=xtxt, ylab="",
           ylim=c(0,0.1), names=nms )
  barplot( PRAm_NL[,"s_pH"], main="s_p[H]", xlab=xtxt, ylab="",
           ylim=c(0,0.1), names=nms, space=5 )
  barplot( PRAm_NL[,"s_V" ], main="s_V"   , xlab=xtxt, ylab="",
           ylim=c(0,1), names=nms, space=5 )
  barplot( PRAm_NL[,"s_R" ], main="s_R"   , xlab=xtxt, ylab="",
           ylim=c(0,0.1), names=nms, space=5 )

  PRA ( x_4, z_4,        0.5  )
  PRAm( x_4, z_4, c(-0.5,0.5) )
  PRA ( l_xz.L[[1]][,1] , l_xz.L [[1]][,2], 1 )
  PRAm( l_xz.L[[1]][,1] , l_xz.L [[1]][,2]    )
  PRA ( l_xz.NL[[1]][,1], l_xz.NL[[1]][,2], 1 )
  PRAm( l_xz.NL[[1]][,1], l_xz.NL[[1]][,2], c(0.5,1) )

  l_xz      <- l_xz.NL ; n_d <- length(l_xz)
  thr.seq   <- c( 0.5, 1.5 )
  PRA.tbl   <- t( sapply( 1:n_d, function(d){
    PRA ( l_xz[[d]][,"x"], l_xz[[d]][,"z"], thr=thr.seq[2]) } ) )
  PRAm.tbl1 <- t( sapply( 1:n_d, function(d){
    PRAm( l_xz[[d]][,"x"], l_xz[[d]][,"z"], thr=thr.seq)$seq[1,] } ) )
  PRAm.tbl2 <- t( sapply( 1:n_d, function(d){
    PRAm( l_xz[[d]][,"x"], l_xz[[d]][,"z"], thr=thr.seq)$seq[2,] } ) )
  # Standard deviation of the PRA-results over the n_d samples
  s_pH <- sd(PRA.tbl[,"pH"]) ; s_V <- sd(PRA.tbl[,"V"]) ; s_R <- sd(PRA.tbl[,"R"])
  c( s_pH , s_V , s_R )
  # Standard deviation of the PRAm-results over the n_d samples
  # the ".1" and ".2" refer to (x <= thr[1]) and (thr[1] < x <= thr[2]) resp.
  s_pH1 <- sd(PRAm.tbl1[,"pH"]) ; s_pH2 <- sd(PRAm.tbl2[,"pH"])
  s_V1  <- sd(PRAm.tbl1[,"V" ]) ; s_V2  <- sd(PRAm.tbl2[,"V" ])
  s_R1  <- sd(PRAm.tbl1[,"R" ]) ; s_R2  <- sd(PRAm.tbl2[,"R" ])
  rbind( c( thr.seq[1], s_pH1 , s_V1 , s_R1 ),
         c( thr.seq[2], s_pH2 , s_V2 , s_R2 ) )

  par( mfrow=c(3,3), mar=c(2,2,2,2) )
  range.s_pH <- range( 0, PRA.tbl[,"s_pH"], s_pH, s_pH1, s_pH2 )
  range.s_V  <- range( 0, PRA.tbl[,"s_V" ], s_V , s_V1 , s_V2  )
  range.s_R  <- range( 0, PRA.tbl[,"s_R" ], s_R , s_R1 , s_R2  )
  hist( PRA.tbl  [,"s_pH"], xlim=range.s_pH, main="sd( p[H] )", xlab="", ylab="" )
    abline( v=s_pH, col="red", lwd=1 )
  hist( PRA.tbl  [,"s_V" ], xlim=range.s_V , main="sd( V )"   , xlab="", ylab="" )
    abline( v=s_V , col="red", lwd=1 )
  hist( PRA.tbl  [,"s_R" ], xlim=range.s_R , main="sd( R )"   , xlab="", ylab="" )
    abline( v=s_R , col="red", lwd=1 )
  hist( PRAm.tbl1[,"s_pH"], xlim=range.s_pH, main="sd( p[H1] )", xlab="", ylab="" )
    abline( v=s_pH1, col="red", lwd=1 )
  hist( PRAm.tbl1[,"s_V" ], xlim=range.s_V , main="sd( V1 )"   , xlab="", ylab="" )
    abline( v=s_V1 , col="red", lwd=1 )
  hist( PRAm.tbl1[,"s_R"] , xlim=range.s_R , main="sd( R1 )"   , xlab="", ylab="" )
    abline( v=s_R1 , col="red", lwd=1 )
  hist( PRAm.tbl2[,"s_pH"], xlim=range.s_pH, main="sd( p[H2] )", xlab="", ylab="" )
    abline( v=s_pH2, col="red", lwd=1 )
  hist( PRAm.tbl2[,"s_V" ], xlim=range.s_V , main="sd( V2 )"   , xlab="", ylab="" )
    abline( v=s_V2 , col="red", lwd=1 )
  hist( PRAm.tbl2[,"s_R" ], xlim=range.s_R , main="sd( R2 )"   , xlab="", ylab="" )
    abline( v=s_R2 , col="red", lwd=1 )

  EzVz_Gauss <- function( m.=m, S.=S, thr.=thr ) {
    mx     <- m.[1]   ; mz <- m.[2]
    Vx     <- S.[1,1] ; Vz <- S.[2,2]   ; Vxz  <- S.[1,2] ; r <- Vxz/sqrt(Vx*Vz)
    pthr   <- dnorm(thr., mx, sqrt(Vx)) ; Fthr <- pnorm(thr., mx, sqrt(Vx))
    Ez_xlo <- mz - Vxz * pthr / Fthr
    Ez_xhi <- mz + Vxz * pthr / (1-Fthr)
    Vz_xlo <- Vz + r * (thr.-mx) * (Ez_xlo-mz) - (Ez_xlo-mz)^2
    Vz_xhi <- Vz + r * (thr.-mx) * (Ez_xhi-mz) - (Ez_xhi-mz)^2
    result <- c( Ez_xlo, Ez_xhi, Vz_xlo, Vz_xhi )
    names ( result ) <- c( "Ez_xlo", "Ez_xhi", "Vz_xlo", "Vz_xhi" )
    return( result ) }

  set.seed(1)
  mu  <- c(0,0) ; Sigma <- diag(1,2) ; Sigma[1,2] <- Sigma[2,1] <- 0.5
  xz  <- rmvnorm( 1e6, mu, Sigma ) ; x <- xz[,1] ; z <- xz[,2]
  thr <- -1
  EzVz_Gauss( mu, Sigma, thr )
  c( mean( z[x< thr] ), mean( z[x>=thr] ), var( z[x< thr] ), var( z[x>=thr] ) )

  PRA0_Gauss <- function( m.=m, S.=S, thr.=thr ) {
    mx <- m.[1] ; sx <- sqrt(S.[1,1]) ; Vxz <- S.[1,2]
    pH <- pnorm( thr., mx, sx )
    V  <- Vxz * dnorm(thr., mx, sx) / (pH * (1-pH))
    R  <- pH * V
    return( c( pH=pH, V=V, R=R ) ) }

  PRA_Gauss <- function( m.=m, S.=S, n.=n, thr.=thr ) {
    pH <- V <- R <- rep( NA, 1e3 )
    for(j in 1:1e3){
      S     <- riwish( n.-1, S. * (n.-1) ) ; m <- rmvnorm( 1, m., S/n. )
      PRA   <- PRA0_Gauss( m, S, thr. )
      pH[j] <- PRA["pH"] ; V[j] <- PRA["V"] ; R[j] <- PRA["R"]
    }
    return( c( pH=mean(pH),   V=mean(V),   R=mean(R),
             s_pH=sd  (pH), s_V=sd  (V), s_R=sd  (R) ) )
  }

  l_xz <- l_xz.L    ; n_d <- length(l_xz) ; n <- dim(l_xz[[1]])[1]
  xz   <- l_xz[[1]] ; m   <- colMeans(xz) ; S <- cov(xz)
  x    <- xz[,1]    ; z   <- xz[,2]
  PRA_Gauss( m, S, n, -1 )

### PRAm
  PRAm_Gauss <- function( m.=m, S.=S, n.=n, thr.=-1:0 ) {
    n_thr <- length(thr.)
    pH    <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_thr)
    for(i in 1:n_thr){
      pHj <- Vj <- Rj <- rep( NA, 1e3 )
      for(j in 1:1e3){
        S       <- riwish( n.-1, S. * (n.-1) ) ; m <- rmvnorm( 1, m., S/n. )
        mx      <- m[1] ; sx <- sqrt(S[1,1]) ; Vxz <- S[1,2]
        Ez_NotH <- EzVz_Gauss(m,S,thr.[n_thr])["Ez_xhi"]
        if(i==1){
          pHj[j] <- pnorm( thr.[i], mx, sx )
          Ez_Hi  <- EzVz_Gauss(m,S,thr.[i])["Ez_xlo"]
        }else{
          pHj[j] <- pnorm( thr.[i], mx, sx ) - pnorm( thr.[i-1], mx, sx )
          Ez_Hi  <- (
            pnorm(thr.[i  ],mx,sx) * EzVz_Gauss(m,S,thr.[i  ])["Ez_xlo"] -
            pnorm(thr.[i-1],mx,sx) * EzVz_Gauss(m,S,thr.[i-1])["Ez_xlo"] ) / pHj[j]
        }
        Vj[j] <- Ez_NotH - Ez_Hi
        Rj[j] <- pHj[j] * Vj[j]
      }
      pH[i]   <- mean(pHj) ;   V[i] <- mean(Vj) ;   R[i] <- mean(Rj)
      s_pH[i] <- sd  (pHj) ; s_V[i] <- sd  (Vj) ; s_R[i] <- sd  (Rj)
    }
    R.sum <- sum(R) ; pH.sum <- sum(pH) ; V.wsum <- R.sum / pH.sum
    return( list( sum = c( pH.sum=pH.sum, V.wsum=V.wsum, R.sum=R.sum ),
                  seq = cbind( thr., pH, V, R, s_pH, s_V, s_R ) ) )
  }
  
  PRAm_Gauss( m, S, n, -2:-1 )

### Method 5
  PRA_Gauss5 <- function( m.=m, S.=S, n.=n, thr.=thr ) {
    PRA <- PRA0_Gauss( m., S., thr. )
    pH  <- PRA["pH"] ; V <- PRA["V"] ; R <- PRA["R"]
      # VarEz_H    <- EzVz_Gauss( m., S., thr. )[3] / (n.*pH)
      # VarEz_NotH <- EzVz_Gauss( m., S., thr. )[4] / (n.*(1-pH))
      r          <-  S[1,2] / sqrt( prod(diag(S)) )
      VarEz_H    <- (S[2,2] + r*(thr.-m[1])*(R-V) - (R-V)^2) / (n.*pH)
      VarEz_NotH <- (S[2,2] + r*(thr.-m[1])* R    -  R^2   ) / (n.*(1-pH))
    s_pH <- sqrt( pH*(1-pH) / n. )
    s_V  <- sqrt( VarEz_NotH + VarEz_H )
    s_R  <- sqrt( s_pH^2*s_V^2 + s_pH^2*V^2 + pH^2*s_V^2 )
    result <- c( pH, V, R, s_pH, s_V, s_R )
    names ( result ) <- c( "pH", "V", "R", "s_pH", "s_V", "s_R" )
    return( result )
  }
  
  PRA_Gauss5( m, S, n, -1 )
  PRA( x, z, thr=-1 )
  rbind( PRA_Gauss ( m, S, n, -1 ),
         PRA_Gauss5( m, S, n, -1 ) )
  
  PRAm_Gauss( m, S, n, -2:-1 )$seq

  PRA.tbl  <- t( sapply( 1:n_d, function(d){
    PRA_Gauss ( colMeans(l_xz.L[[d]]), cov(l_xz.L[[d]]), n, -1 ) } ) )
  PRA.tbl5 <- t( sapply( 1:n_d, function(d){
    PRA_Gauss5( colMeans(l_xz.L[[d]]), cov(l_xz.L[[d]]), n, -1 ) } ) )
  # Mean PRA-results over the n_d samples
  m_pH  <- mean(PRA.tbl [,"pH"]) ; m_V  <- mean(PRA.tbl [,"V"]) ; m_R  <- mean(PRA.tbl [,"R"])
  m_pH5 <- mean(PRA.tbl5[,"pH"]) ; m_V5 <- mean(PRA.tbl5[,"V"]) ; m_R5 <- mean(PRA.tbl5[,"R"])
  # "True" mean of the PRA-results over the n_d samples
  rbind( c( m_pH , m_V , m_R  ), c( m_pH5, m_V5, m_R5 ) )
  # "True" standard deviation of the PRA-results over the n_d samples
  s_pH  <- sd(PRA.tbl [,"pH"]) ; s_V  <- sd(PRA.tbl [,"V"]) ; s_R  <- sd(PRA.tbl [,"R"])
  s_pH5 <- sd(PRA.tbl5[,"pH"]) ; s_V5 <- sd(PRA.tbl5[,"V"]) ; s_R5 <- sd(PRA.tbl5[,"R"])
  rbind( c( s_pH , s_V , s_R  ), c( s_pH5, s_V5, s_R5 ) )

  par( mfrow=c(2,3), mar=c(2,2,2,2) )
  range.pH <- range( 0, PRA.tbl [,"pH"], PRA.tbl5[,"pH"] )
  range.V  <- range( 0, PRA.tbl [,"V" ], PRA.tbl5[,"V" ] )
  range.R  <- range( 0, PRA.tbl [,"R" ], PRA.tbl5[,"R" ] )
  hist( PRA.tbl [,"pH"], xlim=range.pH, main="pH", xlab="", ylab="" )
    abline( v=m_pH, col="red", lwd=1 )
  hist( PRA.tbl [,"V" ], xlim=range.V , main="V" , xlab="", ylab="" )
    abline( v=m_V , col="red", lwd=1 )
  hist( PRA.tbl [,"R" ], xlim=range.R , main="R" , xlab="", ylab="" )
    abline( v=m_R , col="red", lwd=1 )
  hist( PRA.tbl5[,"pH"], xlim=range.pH, main="pH5", xlab="", ylab="" )
    abline( v=m_pH5, col="red", lwd=1 )
  hist( PRA.tbl5[,"V" ], xlim=range.V , main="V5" , xlab="", ylab="" )
    abline( v=m_V5 , col="red", lwd=1 )
  hist( PRA.tbl5[,"R" ], xlim=range.R , main="R5" , xlab="", ylab="" )
    abline( v=m_R5 , col="red", lwd=1 )

  par( mfrow=c(2,3), mar=c(2,2,2,2) )
  range.s_pH <- range( 0, PRA.tbl [,"s_pH"], PRA.tbl5[,"s_pH"] )
  range.s_V  <- range( 0, PRA.tbl [,"s_V" ], PRA.tbl5[,"s_V" ] )
  range.s_R  <- range( 0, PRA.tbl [,"s_R" ], PRA.tbl5[,"s_R" ] )
  hist( PRA.tbl [,"s_pH"], xlim=range.s_pH, main="sd( pH )", xlab="", ylab="" )
    abline( v=s_pH, col="red", lwd=1 )
  hist( PRA.tbl [,"s_V" ], xlim=range.s_V , main="sd( V )" , xlab="", ylab="" )
    abline( v=s_V , col="red", lwd=1 )
  hist( PRA.tbl [,"s_R" ], xlim=range.s_R , main="sd( R )" , xlab="", ylab="" )
    abline( v=s_R , col="red", lwd=1 )
  hist( PRA.tbl5[,"s_pH"], xlim=range.s_pH, main="sd( pH5 )", xlab="", ylab="" )
    abline( v=s_pH5, col="red", lwd=1 )
  hist( PRA.tbl5[,"s_V" ], xlim=range.s_V , main="sd( V5 )" , xlab="", ylab="" )
    abline( v=s_V5 , col="red", lwd=1 )
  hist( PRA.tbl5[,"s_R" ], xlim=range.s_R , main="sd( R5 )" , xlab="", ylab="" )
    abline( v=s_R5 , col="red", lwd=1 )

  n         <- 5e3
  x_U       <- runif(n)     ; z_U <- 0.5 + x_U/2             + rnorm(n,0,0.05)
  thr.seq_U <- quantile( x_U, (1:19)/20 )
  PRA.seq_U <- t( sapply( thr.seq_U, function(t){PRA(x_U,z_U,t)} ) )
  x_E       <- rexp(n)      ; z_E <- 1   - exp(-x_E)/2       + rnorm(n,0,0.05)
  thr.seq_E <- quantile( x_E, (1:19)/20 )
  PRA.seq_E <- t( sapply( thr.seq_E, function(t){PRA(x_E,z_E,t)} ) )
  x_B       <- rbeta(n,5,1) ; z_B <- 0.5 + pbeta(x_B,5,1)/2  + rnorm(n,0,0.05)
  thr.seq_B <- quantile( x_B, (1:19)/20 )
  PRA.seq_B <- t( sapply( thr.seq_B, function(t){PRA(x_B,z_B,t)} ) )
  x_t       <- rt(n,1,30)   ; z_t <- 0.5 + pt(x_t,1,30)/2    + rnorm(n,0,0.05)
  thr.seq_t <- quantile( x_t, (1:19)/20 )
  PRA.seq_t <- t( sapply( thr.seq_t, function(t){PRA(x_t,z_t,t)} ) )
  x_G       <- rnorm(n)     ; z_G <- 0.5 + pnorm(x_G)/2      + rnorm(n,0,0.05)
  thr.seq_G <- quantile( x_G, (1:19)/20 )
  PRA.seq_G <- t( sapply( thr.seq_G, function(t){PRA(x_G,z_G,t)} ) )
  x_L       <- rnorm(n)     ; z_L <- 0.5 + 0.5/(1+exp(-2*x_L)) + rnorm(n,0,0.05)
  thr.seq_L <- quantile( x_L, (1:19)/20 )
  PRA.seq_L <- t( sapply( thr.seq_L, function(t){PRA(x_L,z_L,t)} ) )

  par(mfrow=c(2,3), mar=c(4,2,2,1))
  plot( x_U, z_U, xlab="", ylim=range(z_U,PRA.seq_U),
        pch=".", main="Uniform" )
  points( thr.seq_U, PRA.seq_U[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_U, PRA.seq_U[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_U, PRA.seq_U[,"R"] , type="b", lwd=3, col="red" )
  plot( x_E, z_E, xlab="", xlim=range(thr.seq_E), ylim=range(z_E,PRA.seq_E),
        pch=".", main="Exponential" )
  points( thr.seq_E, PRA.seq_E[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_E, PRA.seq_E[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_E, PRA.seq_E[,"R"] , type="b", lwd=3, col="red" )
  plot( x_B, z_B,  xlab="", ylim=range(z_B,PRA.seq_B),
          pch=".", main="Beta" )
  points( thr.seq_B, PRA.seq_B[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_B, PRA.seq_B[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_B, PRA.seq_B[,"R"] , type="b", lwd=3, col="red" )
  plot( x_t, z_t,  xlab="x or thr", xlim=range(thr.seq_t), ylim=range(z_t,PRA.seq_t),
        pch=".", main="t" )
  points( thr.seq_t, PRA.seq_t[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_t, PRA.seq_t[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_t, PRA.seq_t[,"R"] , type="b", lwd=3, col="red" )
  plot( x_G, z_G,  xlab="x or thr", ylim=range(z_G,PRA.seq_G),
        pch=".", main="Gaussian" )
  points( thr.seq_G, PRA.seq_G[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_G, PRA.seq_G[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_G, PRA.seq_G[,"R"] , type="b", lwd=3, col="red" )
  plot( x_L, z_L,  xlab="x or thr", ylim=range(z_L,PRA.seq_L),
        pch=".", main="Gaussian-Logistic" )
  points( thr.seq_L, PRA.seq_L[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_L, PRA.seq_L[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_L, PRA.seq_L[,"R"] , type="b", lwd=3, col="red" )
  legend ( "topleft", legend=c("z","pH","V","R"),
           col=c("black","blue","green","red"),
           pch=c(1,NA,NA,NA), lty=c(NA,1,1,1), cex=0.75 )

  PRAi <- function( xz, thr=c(0,0) ) {
    x1  <- xz[,1]  ; x2 <- xz[,2] ; z <- xz[,3]
    n_c <- 2^2 - 1 ; n  <- length(x1)
    H   <- vector("list",n_c)
    n_H <- pH <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_c)
    H[[1]]  <- which(x1 <  thr[1] & x2 <  thr[2]) ; n_H[1] <- length(H[[1]])
    H[[2]]  <- which(x1 <  thr[1] & x2 >= thr[2]) ; n_H[2] <- length(H[[2]])
    H[[3]]  <- which(x1 >= thr[1] & x2 <  thr[2]) ; n_H[3] <- length(H[[3]])
    NotH    <- which(x1 >= thr[1] & x2 >= thr[2]) ; n_NotH <- length(NotH)
    pH      <- n_H / n         ; s_pH      <- sqrt( pH*(1-pH) / n )
    Ez_NotH <- mean( z[NotH] ) ; s_Ez_NotH <- sqrt( var(z[NotH] ) / n_NotH )   
    for(i in 1:n_c) {
      Ez_Hi   <-       mean( z[ H[[i]] ] )
      s_Ez_Hi <- sqrt( var ( z[ H[[i]] ] ) / n_H[i] )
      V[i]    <- Ez_NotH - Ez_Hi
      s_V[i]  <- sqrt( s_Ez_NotH^2 + s_Ez_Hi^2 ) }
    R     <- pH * V
    s_R   <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2  )
    R.sum <- sum(R) ; pH.sum <- sum(pH) ; V.wsum <- R.sum / pH.sum
    return( list( sum = c( pH.sum=pH.sum, V.wsum=V.wsum, R.sum=R.sum ),
                  cat = cbind( 1:3, pH, V, R, s_pH, s_V, s_R ) ) ) }

  round( PRAi(xz_G3,c(0,0))$cat, 3 )
  round( PRAi(xz_G3,c(0,0))$sum, 3 )

  set.seed(1)
  n  <- 1e2
  x1 <- rbeta( n, 3, 3 ) ; x2 <- rbeta( n, 3, 3 )
  z  <- as.integer( x1 >= 0.5 | x2 >= 0.5 )
  xz <- cbind( x1, x2, z )
  pairs( xz )

  thr <- c( 0.5, 0.5 )
  round( PRAi( xz, thr )$cat, 3 )

  set.seed(1)
  n_t <- 8 ; x <- z <- array( dim=c(4,4,n_t) )
  for(r in 1:4){ for(c in 1:4){ x[r,c,] <- rbeta( n_t, r+c, 10-r-c ) } }
  z   <- x^2
  plot( rast(x), nr=2, reset=T, pax=list(lab=F), main=paste0("x[",1:n_t,"]") )
  plot( rast(z), nr=2, reset=T, pax=list(lab=F), main=paste0("z[",1:n_t,"]") )

  whereQ <- function( x=array(dim=c(nlon,nlat,n_t)),
                      z=array(dim=c(nlon,nlat,n_t)), thr.=thr ) {
    ns    <- prod( dim(x)[1:2] ) ; n_t <- dim(x)[3]
    freqH <- function(x,thr.=thr){ sum(x<thr.) }
    n_tH  <- apply(x, c(1,2), freqH)
    siteQ <- which( n_tH > 0, arr.ind=TRUE ) ; nQ <- dim(siteQ)[1]
    return( siteQ ) }
  thr   <- 0.25
  siteQ <- whereQ( x, z, thr )
  Q     <- matrix( 0, ncol=4, nrow=4 ) ; Q[siteQ] <- 1
  # plot( Q, key=NULL, asp=1, col=c("green","red") ) # Base R
  plot( rast(Q), main="Q", grid=T, col=c("green","red"),
        pax=list(at=1:4,ylabs=4:1,hadj=4) ) # R-package 'terra'

  PRA3 <- function( x=array(dim=c(nlon,nlat,n_t)),
                    z=array(dim=c(nlon,nlat,n_t)), thr.=thr ) {
    ns    <- prod( dim(x)[1:2] ) ; n_t <- dim(x)[3]
    freqH <- function(x,thr.=thr){ sum(x<thr.) }
    n_tH  <- apply(x, c(1,2), freqH)
    siteQ <- which( n_tH > 0, arr.ind=TRUE ) ; nQ <- dim(siteQ)[1]
    Q     <- nQ / ns ; s_Q <- sqrt( Q*(1-Q) / ns ) 
    xQ    <- sapply( 1:n_t, function(i){x[,,i][siteQ]} )
    zQ    <- sapply( 1:n_t, function(i){z[,,i][siteQ]} )
    PRAQ  <- PRA( c(xQ), c(zQ), thr. )
    pH    <- PRAQ["pH"]   ; V   <- PRAQ["V"]   ; R.Q   <- PRAQ["R"]
    s_pH  <- PRAQ["s_pH"] ; s_V <- PRAQ["s_V"] ; s_R.Q <- PRAQ["s_R"]
    R     <- Q * R.Q
    s_R   <- sqrt( s_Q^2*s_R.Q^2 + s_Q^2*R.Q^2 + Q^2*s_R.Q^2 )
    result           <- c(  Q ,  pH ,  V ,  R ,  s_Q ,  s_pH ,  s_V ,  s_R  )
    names ( result ) <- c( "Q", "pH", "V", "R", "s_Q", "s_pH", "s_V", "s_R" )
    return( result ) }

  PRAst <- PRA3( x, z, thr ) ; round(PRAst,3)

  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  curve( 1-exp(-x), xlim=c(0,3), ylab="", main="E[z|x]" )

  p <- function(x){ rep(1/3,length(x)) }
  v <- function(x){ exp(-x) }
  r <- function(x){ exp(-x) / 3 }
  n_x    <- 31       ; x.seq <- seq( 0, 3, length.out=n_x )
  p.seq  <- p(x.seq) ; v.seq <- v(x.seq) ; r.seq <- r(x.seq)
  vrange <- range(0,v.seq)
  par( mfrow=c(1,3), mar=c(5,3,2,0), bty="n", lwd=2 )
  plot( x.seq, p.seq, type="l", xlab="x", xlim=c(0,3), ylim=c(0,1), main="p[x]" )
  plot( x.seq, v.seq, type="l", xlab="x", xlim=c(0,3), ylim=c(0,1), main="v(x)" )
  plot( x.seq, r.seq, type="l", xlab="x", xlim=c(0,3), ylim=c(0,1), main="r(x)" )

  R.seq  <- ( exp(-x.seq[-n_x]) - exp(-x.seq[-1]) ) / 3
  pH.seq <- ( x.seq[-1] - x.seq[-n_x] ) / 3
  V.seq  <- R.seq / pH.seq

  par( mfrow=c(1,3), mar=c(5,3,2,0) )
  xtxt <- "Upper bound of interval" ; nms <- x.seq[-1]
  barplot( pH.seq, main="p[H]", xlab=xtxt, ylab="", ylim=c(0,0.1), names=nms )
  barplot( V.seq , main="V"   , xlab=xtxt, ylab="", ylim=c(0,1  ), names=nms )
  barplot( R.seq , main="R"   , xlab=xtxt, ylab="", ylim=c(0,0.1), names=nms )

  m     <- c(0,0) ; S  <- diag(1,2) ; S[1,2] <- S[2,1] <- 0.5
  Vx    <- S[1,1] ; Vz <- S[2,2]    ; rxz    <- S[1,2] / sqrt(Vx*Vz)
  mx    <- m[1]   ; mz <- m[2]
  px    <- function( x, m=mx, V=Vx ){ dnorm( x, m, sqrt(V) ) }
  x.seq <- seq( mx-2, mx+2, length.out=41 )
  Ez_x  <- function(x){ mz + (x-mx)*rxz }
  v     <- function(x,thr.=0){ EzVz_Gauss( m,S, thr.=0 )["Ez_xhi"] - Ez_x(x) }
  r     <- function(x,thr.=0){ px(x) * v(x,thr.) }
  thr   <- 1
  p.seq <- px(x.seq) ; v.seq <- v(x.seq,thr) ; r.seq <- r(x.seq,thr)

  par( mfrow=c(1,4), mar=c(5,3,2,0), bty="n", lwd=2 )
  plot( x.seq, p.seq, type="l", xlab="x", xlim=c(-2,1), main="p[x]", ylim=c(0,0.4) )
    abline( v=thr, col="red", lty=2 )
  plot( x.seq, v.seq, type="l", xlab="x", xlim=c(-2,1), main="v(x)" )
    abline( v=thr, col="red", lty=2 )
  plot( x.seq, r.seq, type="l", xlab="x", xlim=c(-2,1), main="r(x)" )
    abline( v=thr, col="red", lty=2 )
  plot( v.seq, p.seq, type="l", xlab="E[loss] = V", xlim=c(-0.5,2), main="p[ E[loss] ]", ylim=c(0,0.4)  )

# Example I: Sparse Data (n = 4)
## Sampling-based PRA
  PRA_Ss <- PRA( x_4, z_4, thr=0 )
## Distribution-based PRA
  set.seed(1)
  PRA_Ds  <- PRA_Gauss ( m_4, S_4, n_4, 0 )
  PRA_Ds5 <- PRA_Gauss5( m_4, S_4, n_4, 0 )
## Summary of various PRAs
  rbind( PRA_Ss, PRA_Ds, PRA_Ds5 )

# Example II: Linear Dataset
  xz  <- l_xz.L[[1]]   ; m  <- colMeans(xz) ; S  <- cov(xz)
  x   <- xz[,1]        ; z  <- xz[,2]       ; n  <- length(x)
  X   <- cbind( 1, x ) ; Vz <- 1            ; Sz <- diag(Vz,n)
  thr <- -1
## Sampling-based PRA
  PRA_Ss <- PRA( x, z, thr )
## Distribution-based PRA
  PRA_Ds <- PRA_Gauss( m, S, n, thr )
## Model-based PRA with LS72 for the linear data instead of Nimble
  # Prior:
  mb        <- c(0,0) ; Vb <- c(1.e4,1.e4) ; Sb <- diag(Vb)
  # Posterior:
  Sb_y_LS72 <- solve( solve(Sb) + t(X) %*% solve(Sz) %*% X )
  mb_y_LS72 <- Sb_y_LS72 %*% (solve(Sb) %*% mb + t(X) %*% solve(Sz) %*% z)

  par(mfrow=c(1,2))
  plot( x, z ) ; abline( mb_y_LS72, col="red" )
  nsmpl <- 10  ; smpl_b <- rmvnorm( nsmpl, mean=mb_y_LS72, sigma=Sb_y_LS72 )
  for( i in 1:nsmpl) { abline( smpl_b[i,], lty=2 ) }
### Model-based PRA using the posterior parameter distribution
  n_unc    <- 1e3
  theta    <- rmvnorm( n_unc, mean=mb_y_LS72, sigma=Sb_y_LS72 )
  lm.alpha <- theta[,1]             ; lm.beta <- theta[,2]
  lm.sigma <- rep( sqrt(Vz), n_unc )
  thr      <- seq( -2, 2, by=0.05 ) ; n_thr <- length(thr)
  EzHj     <- EzNotHj  <-   Rj  <-   Vj           <- numeric(n_unc)
  EzH      <- EzNotH   <-   R   <-   V   <-   pH  <- numeric(n_thr)
                          s_R   <- s_V   <- s_pH  <- numeric(n_thr)
                            R2  <-   V2  <-   pH2 <- numeric(n_thr) # remove later?
                          s_R2  <- s_V2  <- s_pH2 <- numeric(n_thr) # remove later?
  LCIR     <- UCIR     <- LCIV  <- UCIV           <- numeric(n_thr)
  LCIzNotH <- UCIzNotH <- LCIzH <- UCIzH          <- numeric(n_thr)
  qu       <- function( z, q=0.025 ){ quantile( z, q, na.rm=T ) }
  for(i in 1:n_thr) {
    i_H    <- which( x <  thr[i] ) ; n_H     <- length(i_H)
    i_NotH <- which( x >= thr[i] ) ; n_NotH  <- length(i_NotH)
    pH[i]  <- n_H / n              ; s_pH[i] <- sqrt( pH[i]*(1-pH[i]) / n )
    # Method 1: V, R from model expectations and UQ from law of total variance
    Ex_H    <- mean( x[i_H] )          ; Ex_NotH <- mean( x[i_NotH] )
    Ez_H    <- c(1,Ex_H) %*% mb_y_LS72 ; Ez_NotH <- c(1,Ex_NotH) %*% mb_y_LS72
    V[i]    <- Ez_NotH-Ez_H
    R[i]    <- pH[i] * V[i]
    Vzi     <- function(i){ t(c(1,x[i])) %*% Sb_y_LS72 %*% c(1,x[i]) + Vz }
    Vz_H    <- sum( sapply(i_H   ,Vzi) ) / n_H    + mb_y_LS72[2]^2 * var(x[i_H])
    Vz_NotH <- sum( sapply(i_NotH,Vzi) ) / n_NotH + mb_y_LS72[2]^2 * var(x[i_NotH])
    s_V[i]  <- sqrt( Vz_H / n_H + Vz_NotH / n_NotH )
    s_R[i]  <- sqrt( s_pH[i]^2 * s_V[i]^2 + s_pH[i]^2 * V[i]^2 + pH[i]^2 * s_V[i]^2 )
    # Method 2
    for(j in 1:n_unc) {
      zj      <- lm.alpha[j] + lm.beta[j] * x + rnorm( n, 0, lm.sigma[j] )
      EzHj[j] <- mean( zj[i_H] ) ; EzNotHj[j] <- mean( zj[i_NotH] )
      Vj[j]   <- EzNotHj[j] - EzHj[j]
      Rj[j]   <- pH[i] * Vj[j]
    }
    R2[i]   <- mean( Rj ) ; V2[i]   <- mean( Vj )
    s_R2[i] <- sd  ( Rj ) ; s_V2[i] <- sd  ( Vj )
    EzH[i]      <- mean( EzHj )  ; EzNotH[i]   <- mean( EzNotHj ) 
    LCIzNotH[i] <- qu( EzNotHj ) ; UCIzNotH[i] <- qu( EzNotHj, 0.975 )
    LCIzH[i]    <- qu( EzHj    ) ; UCIzH[i]    <- qu( EzHj   , 0.975 )
    LCIV[i]     <- qu( Vj      ) ; UCIV[i]     <- qu( Vj     , 0.975 )
    LCIR[i]     <- qu( Rj      ) ; UCIR[i]     <- qu( Rj     , 0.975 )
  }

  par    ( mfrow=c(1,2), mar=c(3,5,1,1) )
  plot   ( thr, EzH, type="l",
           ylim=range(c(LCIzNotH,LCIzH,UCIzNotH,UCIzH),na.rm=T),
           xlab="Threshold", ylab="(Conditional) Expectation of z")
  polygon( c(thr,thr[n_thr:1]), c(LCIzH,UCIzH[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIzNotH,UCIzNotH[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  abline ( h=m[2], col="blue", lty=3 )
  lines  ( thr, EzH, col="red") ; lines( thr, EzNotH, col="black" )
  legend ( "bottomright", legend=c("E[z|¬H]","E[z]","E[z|H]"),
           col=c("black","blue","red"), lty=c(1,3,1), cex=0.75 )

  plot   ( thr, V, col="black", type="l",
           ylim=range(c(LCIV,UCIV,LCIR,UCIR),na.rm=T),
           xlab="Threshold", ylab="V, R" )
  polygon( c(thr,thr[n_thr:1]), c(LCIV,UCIV[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIR,UCIR[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  lines  ( thr, V, col="black" ) ; lines( thr, R, col="red" )
  legend ( "bottomright", legend=c("V","R"), col=c("black","red"), lty=1, cex=0.75 )

  i       <- which(thr == -1)
  PRA_Ms  <- c( pH[i], V[i] , R[i] , s_pH[i], s_V[i] , s_R[i]  )
  PRA_Ms2 <- c( pH[i], V2[i], R2[i], s_pH[i], s_V2[i], s_R2[i] )

  PRA_LS72 <- function( x, z, thr=0, Vz=1 ) {
    n         <- length(x) ; X  <- cbind(1,x)
    mb        <- c(0,0)    ; Vb <- c(1.e4,1.e4)
    Sb        <- diag(Vb)  ; Sz <- diag(Vz,n)
    Sb_y_LS72 <- solve( solve(Sb) + t(X) %*% solve(Sz) %*% X )
    mb_y_LS72 <- Sb_y_LS72 %*% (solve(Sb) %*% mb + t(X) %*% solve(Sz) %*% z)
    i_H       <- which( x <  thr )       ; n_H     <- length(i_H)
    i_NotH    <- which( x >= thr )       ; n_NotH  <- length(i_NotH)
    pH        <- n_H / n                 ; s_pH    <- sqrt( pH*(1-pH) / n )
    Ex_H      <- mean( x[i_H] )          ; Ex_NotH <- mean( x[i_NotH] )
    Ez_H      <- c(1,Ex_H) %*% mb_y_LS72 ; Ez_NotH <- c(1,Ex_NotH) %*% mb_y_LS72
    V         <- Ez_NotH - Ez_H          ; R       <- pH * V
      Vzi     <- function(i){ t(c(1,x[i])) %*% Sb_y_LS72 %*% c(1,x[i]) + Vz }
    Vz_H      <- sum( sapply(i_H   ,Vzi) ) / n_H    + mb_y_LS72[2]^2 * var(x[i_H])
    Vz_NotH   <- sum( sapply(i_NotH,Vzi) ) / n_NotH + mb_y_LS72[2]^2 * var(x[i_NotH])
    s_V       <- sqrt( Vz_H / n_H + Vz_NotH / n_NotH )
    s_R       <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2 )
    return( list( mb  = mb_y_LS72, Sb = Sb_y_LS72,
                  PRA = c( pH=pH, V=V, R=R, s_pH=s_pH, s_V=s_V, s_R=s_R ) ) )
  }

  PRA_Ms_LS72 <- PRA_LS72( x, z, thr=-1, Vz=1 )
## Summary of various PRAs
  rbind( PRA_Ss, PRA_Ds, PRA_Ms_LS72$PRA, PRA_Ms, PRA_Ms2 )

# Example IV: A Sequence of Linear Datasets (n exponentially increasing)
  thr <- 0
  l_xz <- l_xz.L3 ; n_d <- length(l_xz)
  PRA.tbl <- t( sapply( 1:n_d, function(d){
    PRA ( l_xz[[d]][,1], l_xz[[d]][,2], thr=thr) } ) )
  i2 <- 2
  PRA.tbl2 <- t( sapply( i2:n_d, function(d){
    PRA_Gauss( m.=colMeans(l_xz[[d]]), cov(l_xz[[d]]), 2^d, thr) } ) )
  i3 <- 4 ; i3max <- 11
  PRA.tbl3 <- t( sapply( i3:i3max, function(d){
    PRA_LS72( l_xz[[d]][,1], l_xz[[d]][,2], thr=0, Vz=1 )$PRA } ) )

  par(mfrow=c(2,3), mar=c(5,2,2,1))
  plot( 1:n_d, PRA.tbl[,"pH"], main="pH", xlab="", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"pH"], pch=20, col="blue" )
  points( i3:i3max, PRA.tbl3[1:(i3max-i3+1),"pH"], pch= 3, col="green" )
  abline(h=0.5  ,col="red",lty=2)
  plot( 1:n_d, PRA.tbl[,"V"] , main="V" , xlab="", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"V" ], pch=20, col="blue" )
  points( i3:i3max, PRA.tbl3[1:(i3max-i3+1),"V" ], pch= 3, col="green" )
  abline(h=0.798,col="red",lty=2)
  plot( 1:n_d, PRA.tbl[,"R"] , main="R" , xlab="", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"R" ], pch=20, col="blue" )
  points( i3:i3max, PRA.tbl3[1:(i3max-i3+1),"R" ], pch= 3, col="green" )
  abline(h=0.399,col="red",lty=2)
  legend ( "topright", col=c("black","blue","green"), pch=c(1,20,3), cex=0.75,
           legend=c("Sampling-based PRA","Distribution-based PRA","Model-based PRA") )
  plot( 1:n_d, PRA.tbl[,"s_pH"], main="s_pH", xlab="log2(n_d)", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"s_pH"], pch=20, col="blue" )
  points( i3:i3max, PRA.tbl3[1:(i3max-i3+1),"s_pH"], pch= 3, col="green" )
  plot( 1:n_d, PRA.tbl[,"s_V"] , main="s_V" , xlab="log2(n_d)", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"s_V" ], pch=20, col="blue" )
  points( i3:i3max, PRA.tbl3[1:(i3max-i3+1),"s_V" ], pch= 3, col="green" )
  plot( 1:n_d, PRA.tbl[,"s_R"] , main="s_R" , xlab="log2(n_d)", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"s_R" ], pch=20, col="blue" )
  points( i3:i3max, PRA.tbl3[1:(i3max-i3+1),"s_R" ], pch= 3, col="green" )

# Example V: Nonlinear Dataset
  xz <- l_xz.NL[[1]] ; m <- colMeans(xz) ; S <- cov(xz)
  x  <- xz[,1]       ; z <- xz[,2]       ; n <- length(x)
## Sampling-based PRA
  thr <- 1
  PRA_Ss <- PRA( x, z, thr )
## Distribution-based PRA
  PRA_Ds <- PRA_Gauss( m, S, n, thr )
## Model-based PRA using Nimble
  Model1.Code <- nimbleCode({
    lm.alpha  ~ dnorm( 0, sd=100 )
    lm.beta   ~ dnorm( 0, sd=100 )
    lm.tau    ~ dgamma( 0.01, 0.01 )
    lm.sigma <- 1 / sqrt(lm.tau)
    for(i in 1:ndata){
      lm.mu[i] <- lm.alpha + lm.beta*exp(-x[i])
      z[i]      ~ dnorm( lm.mu[i], sd=lm.sigma )
    }
  } )
  Model1.Constants <- list( ndata=n, x=x )
  Model1.Data      <- list(z=z)
  Model1.Nimble    <- nimbleModel  ( Model1.Code, constants=Model1.Constants,
                                                  data=Model1.Data )
  Model1.Comp      <- compileNimble( Model1.Nimble )
  Model1.Conf      <- configureMCMC( Model1.Nimble, print=F )
  Model1.Conf$addMonitors( c("lm.sigma"), print=F )
  Model1.MCMC      <- buildMCMC    ( Model1.Conf )
  Model1.MCMC.Comp <- compileNimble( Model1.MCMC )
  ntheta           <- 1e4 ; nburnin <- 1e3 ; niter <- ntheta + nburnin
  set.seed(1)
  theta <- runMCMC( Model1.MCMC.Comp, nburnin=nburnin, niter=niter, prog=F )

  n_unc    <- 1e3
  itheta   <- sample( 1:ntheta, n_unc, replace=(ntheta<n_unc) )
  thr      <- seq( 0.1, 2.9, by=0.05 )  ;   n_thr <- length(thr)
  EzHj     <- EzNotHj  <-   Rj  <-   Vj          <- numeric(n_unc)
  EzH      <- EzNotH   <-   R   <-   V   <-   pH <- numeric(n_thr)
                          s_R   <- s_V   <- s_pH <- numeric(n_thr)
  LCIR     <- UCIR     <- LCIV  <- UCIV          <- numeric(n_thr)
  LCIzNotH <- UCIzNotH <- LCIzH <- UCIzH         <- numeric(n_thr)
  lm.alpha <- theta[,1] ; lm.beta <- theta[,2] ; lm.sigma <- theta[,3]
  for(i in 1:n_thr) {
    i_H   <- which( x < thr[i] ) ; n_H     <- length(i_H)
    pH[i] <- n_H / n             ; s_pH[i] <- sqrt( pH[i]*(1-pH[i]) / n )
    for(j in 1:n_unc) {
      zj      <- lm.alpha[itheta[j]] + lm.beta[itheta[j]] * exp(-x) +
                 rnorm( n, 0, lm.sigma[itheta[j]] )
      EzHj[j] <- mean( zj[i_H] ) ; EzNotHj[j] <- mean( zj[-i_H] )
      Vj[j]   <- EzNotHj[j] - EzHj[j]
      Rj[j]   <- pH[i] * Vj[j] }
    R[i]   <- mean( Rj ) ; V[i]   <- mean( Vj )
    s_R[i] <- sd  ( Rj ) ; s_V[i] <- sd  ( Vj )
    EzH[i]      <- mean( EzHj )  ; EzNotH[i]   <- mean( EzNotHj ) 
    LCIzNotH[i] <- qu( EzNotHj ) ; UCIzNotH[i] <- qu( EzNotHj, 0.975 )
    LCIzH[i]    <- qu( EzHj    ) ; UCIzH[i]    <- qu( EzHj   , 0.975 )
    LCIV[i]     <- qu( Vj      ) ; UCIV[i]     <- qu( Vj     , 0.975 )
    LCIR[i]     <- qu( Rj      ) ; UCIR[i]     <- qu( Rj     , 0.975 )
  }

  par    ( mfrow=c(1,2), mar=c(3,5,1,1) )
  plot   ( thr, EzH, type="l",
           ylim=range(c(LCIzNotH,LCIzH,UCIzNotH,UCIzH)),
           xlab="Threshold", ylab="(Conditional) Expectation of z")
  polygon( c(thr,thr[n_thr:1]), c(LCIzH,UCIzH[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIzNotH,UCIzNotH[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  abline ( h=m[2], col="blue", lty=3 )
  lines  ( thr, EzH, col="red") ; lines( thr, EzNotH, col="black" )
  legend ( "bottomright", legend=c("E[z|¬H]","E[z]","E[z|H]"),
           col=c("black","blue","red"), lty=c(1,3,1), cex=0.75 )
  plot   ( thr, V, col="black", type="l",
           ylim=range(c(LCIV,UCIV,LCIR,UCIR)),
           xlab="Threshold", ylab="V, R" )
  polygon( c(thr,thr[n_thr:1]), c(LCIV,UCIV[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIR,UCIR[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  lines  ( thr, V, col="black" ) ; lines( thr, R, col="red" )
  legend ( "bottomright", legend=c("V","R"), col=c("black","red"), lty=1, cex=0.75 )

  i      <- which(thr == 1)
  PRA_Ms <- c( pH[i], V[i], R[i], s_pH[i], s_V[i], s_R[i] )

## Summary of various PRAs
  rbind( PRA_Ss, PRA_Ds, PRA_Ms )

# Example VI: A Sequence of Nonlinear Datasets (n exponentially increasing)
  thr <- 1
  l_xz <- l_xz.NL3 ; n_d <- length(l_xz)
  PRA.tbl <- t( sapply( 1:n_d, function(d){
    PRA( l_xz[[d]][,1], l_xz[[d]][,2], thr=thr ) } ) )
  i2 <- 2
  PRA.tbl2 <- t( sapply( i2:n_d, function(d){
    PRA_Gauss( m.=colMeans(l_xz[[d]]), cov(l_xz[[d]]), 2^d, thr ) } ) )

  par(mfrow=c(2,3), mar=c(5,2,2,1))
  PRA_inf <- PRA( l_xz[[d]][,1], l_xz[[d]][,2], thr=thr )
  ylim_pH <- range(0,PRA.tbl[,"pH"],PRA.tbl2[,"pH"])
  ylim_V  <- range(0,PRA.tbl[,"V"] ,PRA.tbl2[,"V"] )
  ylim_R  <- range(0,PRA.tbl[,"R"] ,PRA.tbl2[,"R"] )
  plot( 1:n_d, PRA.tbl[,"pH"], main="pH", xlab="", ylab="", ylim=ylim_pH )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"pH"], pch=20, col="blue" )
    abline(h=PRA_inf["pH"],col="red",lty=2)
  plot( 1:n_d, PRA.tbl[,"V"] , main="V" , xlab="", ylab="", ylim=ylim_V )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"V" ], pch=20, col="blue" )
    abline(h=PRA_inf["V"] ,col="red",lty=2)
  plot( 1:n_d, PRA.tbl[,"R"] , main="R" , xlab="", ylab="", ylim=ylim_R )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"R" ], pch=20, col="blue" )
    abline(h=PRA_inf["R"] ,col="red",lty=2)
    legend ( "topright", col=c("black","blue"), pch=c(1,20), cex=0.75,
             legend=c("Sampling-based PRA","Distribution-based PRA") )
  plot( 1:n_d, PRA.tbl[,"s_pH"], main="s_pH", xlab="log2(n_d)", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"s_pH"], pch=20, col="blue" )
  plot( 1:n_d, PRA.tbl[,"s_V"] , main="s_V" , xlab="log2(n_d)", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"s_V" ], pch=20, col="blue" )
  plot( 1:n_d, PRA.tbl[,"s_R"] , main="s_R" , xlab="log2(n_d)", ylab="" )
  points( i2:n_d, PRA.tbl2[1:(n_d-i2+1),"s_R" ], pch=20, col="blue" )
  legend ( "topright", col=c("black","blue"), pch=c(1,20), cex=0.75,
           legend=c("Sampling-based PRA","Distribution-based PRA") )

# Example VII: German Forestry Data
## Sampling-based PRA
### Single-threshold PRA
  thr <- 250

  par ( mfrow=c(2,2), mar=c(4,4,1,2) )
  plot( x_r3, 100-z_Fs, main="Fs", xlab=""         , ylab="survival (%)",
        ylim=c(98,100))
    abline( v=thr, lty=2 )
  plot( x_r3, 100-z_Q , main="Q" , xlab=""         , ylab=""            ,
        ylim=c(98,100))
    abline( v=thr, lty=2 )
  plot( x_r3, 100-z_Pa, main="Pa", xlab="rain (mm)", ylab="survival (%)",
        ylim=c(95,100))
    abline( v=thr, lty=2 )
  plot( x_r3, 100-z_Ps, main="Ps", xlab="rain (mm)", ylab=""            ,
        ylim=c(98,100))
    abline( v=thr, lty=2 )

  PRA_Ss.Fs <- PRA( x_r3, 100-z_Fs, thr )
  PRA_Ss.Q  <- PRA( x_r3, 100-z_Q , thr )
  PRA_Ss.Pa <- PRA( x_r3, 100-z_Pa, thr )
  PRA_Ss.Ps <- PRA( x_r3, 100-z_Ps, thr )
  
  plotPRA( PRA( x_r3, 100-z_Pa, thr ) )
  
  PRA._ <- as.matrix( cbind(PRA_Ss.Fs, PRA_Ss.Q, PRA_Ss.Pa, PRA_Ss.Ps) )
  colnames(PRA._) <- c( "Fs", "Q", "Pa", "Ps" )
  m     <- PRA._[1:3,] ; s <- PRA._[4:6,]
  par(mfrow=c(1,1))
  col   <- c( "red", "cyan", "yellow" ) 
  bp    <- barplot( m, beside=T, ylim=c(-0.2,1), col=col,
                    legend.text=row.names(m),
                    args.legend=list(x="top",hor=T) )
  segments( bp, m - s, bp, m + s )
  ew    <- (bp[2,1]-bp[1,1])/4
  segments( bp - ew, m - s, bp + ew, m - s )
  segments( bp - ew, m + s, bp + ew, m + s )
### Multi-threshold PRA
  thr.seq <- c(220,250) ; n_thr <- length(thr.seq)

  par ( mfrow=c(2,2), mar=c(4,4,1,2) )
  plot( x_r3, 100-z_Fs, main="Fs", xlab=""         , ylab="survival (%)",
        ylim=c(98,100))
    abline( v=thr.seq, lty=2 )
  plot( x_r3, 100-z_Q , main="Q" , xlab=""         , ylab=""            ,
        ylim=c(98,100))
    abline( v=thr.seq, lty=2 )
  plot( x_r3, 100-z_Pa, main="Pa", xlab="rain (mm)", ylab="survival (%)",
        ylim=c(95,100))
    abline( v=thr.seq, lty=2 )
  plot( x_r3, 100-z_Ps, main="Ps", xlab="rain (mm)", ylab=""            ,
        ylim=c(98,100))
    abline( v=thr.seq, lty=2 )

# PRAm for survival (= 100 - mortality)
  PRA_Sm.Fs  <- PRAm(x_r3, 100-z_Fs, thr.seq)$seq
  PRA_Sm.Q   <- PRAm(x_r3, 100-z_Q , thr.seq)$seq
  PRA_Sm.Pa  <- PRAm(x_r3, 100-z_Pa, thr.seq)$seq
  PRA_Sm.Ps  <- PRAm(x_r3, 100-z_Ps, thr.seq)$seq
  PRA_Sm.tbl <- rbind( PRA_Sm.Fs, PRA_Sm.Q, PRA_Sm.Pa, PRA_Sm.Ps )
  pHm.tbl    <- matrix(PRA_Sm.tbl[,"pH"  ], nrow=n_thr)
  Vm.tbl     <- matrix(PRA_Sm.tbl[,"V"   ], nrow=n_thr)
  Rm.tbl     <- matrix(PRA_Sm.tbl[,"R"   ], nrow=n_thr)
  s_pHm.tbl  <- matrix(PRA_Sm.tbl[,"s_pH"], nrow=n_thr)
  s_Vm.tbl   <- matrix(PRA_Sm.tbl[,"s_V" ], nrow=n_thr)
  s_Rm.tbl   <- matrix(PRA_Sm.tbl[,"s_R" ], nrow=n_thr)
  m             <- rbind(   pHm.tbl,   Vm.tbl,   Rm.tbl )
  s             <- rbind( s_pHm.tbl, s_Vm.tbl, s_Rm.tbl )
  colnames( m ) <- c( "Fs" , "Q"  , "Pa", "Ps" )
  rownames( m ) <- c( "pH1", "pH2", "V1", "V2", "R1", "R2" )
  colnames( s ) <-                colnames( m )
  rownames( s ) <- paste0( "sd_", rownames( m ) )

  par(mfrow=c(1,2))
  col    <- c( "red", "cyan", "yellow" )
  list.i <- list( i1=c(1,3,5), i2=c(2,4,6) )
  sp     <- "Pa"
  for( j in 1:2 ) {
    i  <- list.i[[j]]
    bp <- barplot( m[i,sp], beside=T, ylim=c(-0.2,4.5), col=col,
                   legend.text=row.names(m)[i],
                   args.legend=list(x="topleft",hor=F) )
    segments( bp, (m-s)[i,sp], bp, (m+s)[i,sp] )
    ew <- (bp[2,1]-bp[1,1]) / 4
    segments( bp - ew, (m-s)[i,sp], bp + ew, (m-s)[i,sp] )
    segments( bp - ew, (m+s)[i,sp], bp + ew, (m+s)[i,sp] ) }

  par(mfrow=c(1,2))
  col <- c( "red", "cyan", "yellow" )
  list.i <- list( i1=c(1,3,5), i2=c(2,4,6) )
  for( j in 1:2 ) {
    i   <- list.i[[j]]
    bp  <- barplot( m[i,], beside=T, ylim=c(-0.2,4.5), col=col,
                    legend.text=row.names(m)[i],
                    args.legend=list(x="topleft",hor=F) )
    segments( bp, (m-s)[i,], bp, (m+s)[i,] )
    ew  <- (bp[2,1]-bp[1,1]) / 4
    segments( bp - ew, (m-s)[i,], bp + ew, (m-s)[i,] )
    segments( bp - ew, (m+s)[i,], bp + ew, (m+s)[i,] ) }
### Multi-threshold PRA for mortality rather than survival
  # PRAm for mortality
  PRA_Sm.Fs  <- PRAm(x_r3, z_Fs, thr.seq)$seq
  PRA_Sm.Q   <- PRAm(x_r3, z_Q , thr.seq)$seq
  PRA_Sm.Pa  <- PRAm(x_r3, z_Pa, thr.seq)$seq
  PRA_Sm.Ps  <- PRAm(x_r3, z_Ps, thr.seq)$seq
  PRA_Sm.tbl <- rbind( PRA_Sm.Fs, PRA_Sm.Q, PRA_Sm.Pa, PRA_Sm.Ps )
  PRA_Sm.tbl <- rbind( PRA_Sm.Fs, PRA_Sm.Q, PRA_Sm.Pa, PRA_Sm.Ps )
  pHm.tbl   <- matrix(PRA_Sm.tbl[,"pH"  ], nrow=n_thr)
  Vm.tbl    <- matrix(PRA_Sm.tbl[,"V"   ], nrow=n_thr)
  Rm.tbl    <- matrix(PRA_Sm.tbl[,"R"   ], nrow=n_thr)
  s_pHm.tbl <- matrix(PRA_Sm.tbl[,"s_pH"], nrow=n_thr)
  s_Vm.tbl  <- matrix(PRA_Sm.tbl[,"s_V" ], nrow=n_thr)
  s_Rm.tbl  <- matrix(PRA_Sm.tbl[,"s_R" ], nrow=n_thr)
  m             <- rbind(   pHm.tbl,   Vm.tbl,   Rm.tbl )
  s             <- rbind( s_pHm.tbl, s_Vm.tbl, s_Rm.tbl )
  colnames( m ) <- c( "Fs" , "Q"  , "Pa", "Ps" )
  rownames( m ) <- c( "pH1", "pH2", "V1", "V2", "R1", "R2" )
  colnames( s ) <-                colnames( m )
  rownames( s ) <- paste0( "sd_", rownames( m ) )

  par(mfrow=c(1,2))
  col <- c( "red", "cyan", "yellow" )
  list.i <- list( i1=c(1,3,5), i2=c(2,4,6) )
  for( j in 1:2 ) {
    i   <- list.i[[j]]
    bp  <- barplot( m[i,], beside=T, ylim=c(-4.5,1), col=col,
                    legend.text=row.names(m)[i],
                    args.legend=list(x="bottomleft",hor=F) )
    segments( bp, (m-s)[i,], bp, (m+s)[i,] )
    ew  <- (bp[2,1]-bp[1,1]) / 4
    segments( bp - ew, (m-s)[i,], bp + ew, (m-s)[i,] )
    segments( bp - ew, (m+s)[i,], bp + ew, (m+s)[i,] ) }
## Distribution-based PRA
  xz_Fs <- cbind( x_r3, 100-z_Fs )
  xz_Q  <- cbind( x_r3, 100-z_Q  )
  xz_Pa <- cbind( x_r3, 100-z_Pa )
  xz_Ps <- cbind( x_r3, 100-z_Ps )
  m_Fs <- colMeans(xz_Fs) ; S_Fs <- cov(xz_Fs) ; n <- length(x_r3)
  m_Q  <- colMeans(xz_Q ) ; S_Q  <- cov(xz_Q ) ; n <- length(x_r3)
  m_Pa <- colMeans(xz_Pa) ; S_Pa <- cov(xz_Pa) ; n <- length(x_r3)
  m_Ps <- colMeans(xz_Ps) ; S_Ps <- cov(xz_Ps) ; n <- length(x_r3)
  PRA_Ds.Fs <- PRA_Gauss( m_Fs, S_Fs, n, thr )
  PRA_Ds.Q  <- PRA_Gauss( m_Q , S_Q , n, thr )
  PRA_Ds.Pa <- PRA_Gauss( m_Pa, S_Pa, n, thr )
  PRA_Ds.Ps <- PRA_Gauss( m_Ps, S_Ps, n, thr )

  PRA._ <- as.matrix( cbind(PRA_Ds.Fs, PRA_Ds.Q, PRA_Ds.Pa, PRA_Ds.Ps) )
  colnames(PRA._) <- c( "Fs", "Q", "Pa", "Ps" )
  m     <- PRA._[1:3,] ; s <- PRA._[4:6,]
  par(mfrow=c(1,1))
  col   <- c( "red", "cyan", "yellow" ) 
  bp    <- barplot( m, beside=T, ylim=c(-0.2,1), col=col,
                    legend.text=row.names(m),
                    args.legend=list(x="topleft",hor=T) )
  segments( bp, m - s, bp, m + s )
  ew    <- (bp[2,1]-bp[1,1])/4
  segments( bp - ew, m - s, bp + ew, m - s )
  segments( bp - ew, m + s, bp + ew, m + s )
## Model-based conjugate PRA
### Survival
  set.seed(1)
  x  <- x_r3    ; z  <- 100-z_Pa     ; n  <- length(x)
  Vz <- 1       ; Sz <- diag(Vz,n)   ; X  <- cbind(1,x)
  mb <- c(99,0) ; Vb <- c(1.e4,1.e4) ; Sb <- diag(Vb)
  # Lindley & Smith (1972)
  Sb_y_LS72 <- solve( solve(Sb) + t(X) %*% solve(Sz) %*% X )
  mb_y_LS72 <- Sb_y_LS72 %*% (solve(Sb) %*% mb + t(X) %*% solve(Sz) %*% z)

  plot( x, z, main="Pa", ylim=c(96,101) ) ; abline( mb_y_LS72 )
  nsmpl  <- 10
  smpl_b <- rmvnorm( nsmpl, mean=mb_y_LS72, sigma=Sb_y_LS72 )
  for( i in 1:nsmpl) {
    abline( smpl_b[i,], lty=2 )
  }

  set.seed(1)
  thr <- 250
  i_H <- which( x<thr ) ; i_NotH <- which( x >= thr )
  n_H <- length( i_H )  ; n_NotH <- length( i_NotH )
  pH  <- n_H / n        ; s_pH   <- sqrt( pH*(1-pH) / n )
  # Method 1: V, R from model expectations and UQ from law of total variance
  Ex_H    <- mean( x[i_H ] )         ; Ex_NotH <- mean( x[i_NotH] )
  Ez_H    <- c(1,Ex_H) %*% mb_y_LS72 ; Ez_NotH <- c(1,Ex_NotH) %*% mb_y_LS72
  V       <- Ez_NotH-Ez_H
  R       <- pH * V
  Vzi     <- function(i){ t(c(1,x[i])) %*% Sb_y_LS72 %*% c(1,x[i]) + Vz }
  Vz_H    <- sum( sapply(i_H   ,Vzi) ) / n_H    + mb_y_LS72[2]^2 * var(x[i_H])
  Vz_NotH <- sum( sapply(i_NotH,Vzi) ) / n_NotH + mb_y_LS72[2]^2 * var(x[i_NotH])
  s_V     <- sqrt( Vz_H / n_H + Vz_NotH / n_NotH )
  s_R     <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2 )
  # Method 2: V, R and UQ from sampling
  n_unc    <- 1e3
  theta    <- rmvnorm( n_unc, mean=mb_y_LS72, sigma=Sb_y_LS72 )
  lm.alpha <- theta[,1] ; lm.beta <- theta[,2]
  lm.sigma <- rep( sqrt(Vz), n_unc )
  for(j in 1:n_unc) {
    zj      <- lm.alpha[j] + lm.beta[j] * x + rnorm( n, 0, lm.sigma[j] )
    EzHj[j] <- mean( zj[i_H] ) ; EzNotHj[j] <- mean( zj[i_NotH] )
    Vj[j]   <- EzNotHj[j] - EzHj[j]
    Rj[j]   <- pH * Vj[j]
  }
  R2   <- mean( Rj ) ; V2   <- mean( Vj )
  s_R2 <- sd  ( Rj ) ; s_V2 <- sd  ( Vj )
  PRA_Ms.Pa  <- c( pH, V , R , s_pH, s_V , s_R  )
  PRA_Ms2.Pa <- c( pH, V2, R2, s_pH, s_V2, s_R2 )

  PRA_LS72.Fs <- PRA_LS72( x_r3, 100-z_Fs, thr, Vz=1 )$PRA
  PRA_LS72.Q  <- PRA_LS72( x_r3, 100-z_Q , thr, Vz=1 )$PRA
  PRA_LS72.Pa <- PRA_LS72( x_r3, 100-z_Pa, thr, Vz=1 )$PRA
  PRA_LS72.Ps <- PRA_LS72( x_r3, 100-z_Ps, thr, Vz=1 )$PRA

  PRA._ <- as.matrix( cbind(PRA_LS72.Fs, PRA_LS72.Q, PRA_LS72.Pa, PRA_LS72.Ps) )
  colnames(PRA._) <- c( "Fs", "Q", "Pa", "Ps" )
  m     <- PRA._[1:3,] ; s <- PRA._[4:6,]
  par(mfrow=c(1,1))
  col   <- c( "red", "cyan", "yellow" ) 
  bp    <- barplot( m, beside=T, ylim=c(-0.2,1), col=col,
                    legend.text=row.names(m),
                    args.legend=list(x="topleft",hor=T) )
  segments( bp, m - s, bp, m + s )
  ew    <- (bp[2,1]-bp[1,1])/4
  segments( bp - ew, m - s, bp + ew, m - s )
  segments( bp - ew, m + s, bp + ew, m + s )
## Summary of various PRAs
  rbind( PRA_Ss.Pa, PRA_Ds.Pa, PRA_Ms.Pa, PRA_Ms2.Pa )

# Example VIII: Spatially Distributed Risk
## Chapter 14 of vO&B (2022)
### Data
  spdf_DEU <- readRDS( "data/GADM/gadm41_DEU_0_pk.rds" )
  spdf_DEU <- crop( spdf_DEU, ext(5,16,47,55) )
  ext_DEU  <- ext( spdf_DEU )
  dir_WC <- "/Users/marcelvanoijen/Documents/CAF2021_SEACAF/weather/WorldClim/"
  dir_WC_prec.2000_2009 <- paste0( dir_WC, "prec_2000_2009/" )
  dir_WC_prec.2010_2018 <- paste0( dir_WC, "prec_2010_2018/" )
  s_prec_DEU <- rast( "data/climate/s_prec_DEU.tif" )
  r_alt_DEU  <- rast( "data/elevation/r_alt_DEU.tif" )

  par( mfrow=c(1,2) )
  r_tree_DEU <- rast( "data/landuse/r_tree_DEU.tif" )
  plot( r_tree_DEU, col=terrain.colors(100,rev=T) )

  par( mfrow=c(1,3), mar=c(2,2,2,3) )
  plot( r_alt_DEU, col=terrain.colors(100), range=c(0,1500), main="Altitude (m)" )
  plot( s_prec_DEU[[ 1]], main="Prec. in 2000 (mm)",
        col=rainbow(100,end=0.7), range=c(0,1500), xaxt="n", yaxt="n" )
  plot( s_prec_DEU[[19]], main="Prec. in 2018 (mm)",
        col=rainbow(100,end=0.7), range=c(0,1500), xaxt="n", yaxt="n" )

### A simple model for forest yield class (YC)
  YC <- function( x1, x2 ) {
    alt <- x1  ; prec <- x2
    YC  <- max(0, 10 * (1-exp(-prec/1000)) * (1-alt/1000) )
    return( YC ) }

  r_prec_2000_DEU <- s_prec_DEU[[1]]
  r_YC_DEU <- xapp( r_alt_DEU, r_prec_2000_DEU, fun=YC )
  plot( r_YC_DEU     , main="YC: Model",
        col=terrain.colors(100,rev=T), range=c(0,9) )

## Chapter 15 of vO&B (2022)
  kz <- 30 ; r_U_DEU <- r_YC_DEU * kz
  IRRIG <- 500 ; ka <- 0.1
  r_YC.IRRIG_DEU    <- xapp( r_alt_DEU, r_prec_2000_DEU + IRRIG, fun=YC    )
  r_U.IRRIG_DEU     <- r_YC.IRRIG_DEU    * kz - IRRIG * ka

  par( mfrow=c(1,2), mar=c(3,3,3,1), cex=0.7 )
  plot( r_U_DEU                      , main="u: IRRIG=0" )
  plot( r_U.IRRIG_DEU    - r_U_DEU   , main="u-gain from IRRIG=500",
        zlim=c(-50,50) )
### Multiple action levels
  nI    <- 10 ; layers <- 1:nI ; IRRIG <- (layers-1)*50
  s_U.a_DEU <- NULL
  for(i in layers) {
    r_YC.a_DEU <- xapp( r_alt_DEU, r_prec_2000_DEU + IRRIG[i], fun=YC )
    r_U.a_DEU  <- r_YC.a_DEU * kz - IRRIG[i] * ka
    s_U.a_DEU  <- c( s_U.a_DEU, r_U.a_DEU ) }
  s_U.a_DEU <- rast(s_U.a_DEU)
  a.opt_DEU <- (which.max( s_U.a_DEU ) - 1) * 50

  par( mfrow=c(1,2), mar=c(2,2,2,3), cex=0.8 )
  plot( a.opt_DEU, main="Optimum irrigation (mm y-1)", xaxt="n", yaxt="n",
        col=terrain.colors(100,rev=T))

## Chapter 17 of vO&B (2022)
  gadm5  <- world(resolution=5, path="data/gadm")
  IRRIG <- 500
  s_YC_DEU       <- rast( sapply(1:19, function(i){
    xapp( r_alt_DEU, s_prec_DEU[[i]]      , fun=YC )}) )
  s_YC.IRRIG_DEU <- rast( sapply(1:19, function(i){
    xapp( r_alt_DEU, s_prec_DEU[[i]]+IRRIG, fun=YC )}) )

  kz <- 30 ; ka <- 0.1
  s_U_DEU       <- s_YC_DEU       * kz
  s_U.IRRIG_DEU <- s_YC.IRRIG_DEU * kz - IRRIG * ka

  r_U_DEU       <- mean( s_U_DEU )
  r_U.IRRIG_DEU <- mean( s_U.IRRIG_DEU )

  PRA.Ez_notH <- function( x, z, thr=1000, rel=F ) {
    n       <- length(x) ; H  <- which(x < thr) ; notH   <- which(x >= thr)
    n_H     <- length(H) ; pH <- n_H / n        ; n_notH <- length(notH)
    Ez_H    <- mean( z[H]    ) ; s_Ez_H    <- sqrt( var(z[H   ]) / n_H    )
    Ez_notH <- mean( z[notH] ) ; s_Ez_notH <- sqrt( var(z[notH]) / n_notH )
    if(n_H==0){ Ez_H <- min(z) } ; if(n_notH==0){ Ez_notH <- max(z) }
    V       <- Ez_notH - Ez_H  ; R         <- pH * V
    s_pH    <- sqrt( pH*(1-pH) / n )
    s_V     <- sqrt( s_Ez_H^2 + s_Ez_notH^2 )
    s_R     <- sqrt( s_pH^2*s_V^2 + s_pH^2*V^2 + pH^2*s_V^2 )
    if(rel){V <- V/Ez_notH ; R=R/Ez_notH; s_V=s_V/Ez_notH; s_R=s_R/Ez_notH}
    return( c(Ez_notH=Ez_notH, pH=pH, V=V, R=R, s_pH=s_pH, s_V=s_V, s_R=s_R) )
  }

  s_PRA_DEU           <- xapp( s_prec_DEU      , s_U_DEU      , fun=PRA.Ez_notH )
  s_PRA.IRRIG_DEU     <- xapp( s_prec_DEU+IRRIG, s_U.IRRIG_DEU, fun=PRA.Ez_notH )
  r_Ez_notH_DEU       <- s_PRA_DEU$Ez_notH
  r_Ez_notH.IRRIG_DEU <- s_PRA.IRRIG_DEU$Ez_notH
  r_pH_DEU            <- s_PRA_DEU$pH
  r_pH.IRRIG_DEU      <- s_PRA.IRRIG_DEU$pH
  r_pH_DEU            <- mask(r_pH_DEU,r_alt_DEU)
  r_pH.IRRIG_DEU      <- mask(r_pH.IRRIG_DEU,r_alt_DEU)
  r_V_DEU <- s_PRA_DEU$V ; r_V.IRRIG_DEU <- s_PRA.IRRIG_DEU$V
  r_R_DEU <- s_PRA_DEU$R ; r_R.IRRIG_DEU <- s_PRA.IRRIG_DEU$R

  r_Rc_DEU       <- r_R_DEU       - r_Ez_notH_DEU
  r_Rc.IRRIG_DEU <- r_R.IRRIG_DEU - r_Ez_notH.IRRIG_DEU

  par( mfrow=c(2,3), mar=c(1,1,1,3) )
  plot( r_pH_DEU      , main="pH", xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0, 1) ) ; plot( gadm5, add=T)
  plot( r_V_DEU       , main="V" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,50) ) ; plot( gadm5, add=T)
  plot( r_R_DEU       , main="R" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,50) ) ; plot( gadm5, add=T)
  plot( r_pH.IRRIG_DEU, main="pH.IRRIG", xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0, 1) ) ; plot( gadm5, add=T)
  plot( r_V.IRRIG_DEU , main="V.IRRIG" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,50) ) ; plot( gadm5, add=T)
  plot( r_R.IRRIG_DEU , main="R.IRRIG" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,50) ) ; plot( gadm5, add=T)

## Chapter 18 of vO&B (2022)
  exposure <- function( pH, ... ){
    pH <- pH[!is.na(pH)]        ; n <- length(pH)
    nE <- length( pH[pH>1/21] ) ; Q <- nE/n
    return( Q ) }
  r_Q4_DEU    <- aggregate( r_pH_DEU, fact=4, fun=exposure )
  s_prec4_DEU <- aggregate( s_prec_DEU, fact=4, fun="c" )
  s_U4_DEU    <- aggregate( s_U_DEU   , fact=4, fun="c" )
    msk_DEU   <- ifel( r_Q4_DEU>0 & !is.na(s_prec4_DEU[[1]]), 1, NA )
  s_PRA4_DEU  <- xapp( s_prec4_DEU, s_U4_DEU, fun=PRA.Ez_notH )
  r_pH4_DEU   <- s_PRA4_DEU$pH ; r_pH4_DEU <- mask( r_pH4_DEU, msk_DEU )
  r_V4_DEU    <- s_PRA4_DEU$V  ; r_V4_DEU  <- mask( r_V4_DEU , msk_DEU )
  r_R4_DEU    <- s_PRA4_DEU$R  ; r_R4_DEU  <- mask( r_R4_DEU , msk_DEU )

  par( mfrow=c(2,3), mar=c(1,1,1,3) )
  plot( mean(s_prec4_DEU), main="mean(prec)", xaxt="n", yaxt="n",
        col=map.pal("ryb"), range=c(0,1700))   ; plot( gadm5, add=T)
  plot( mean(s_U4_DEU   ), main="mean(u)"  , xaxt="n", yaxt="n",
        col=map.pal("ryg"), range=c(0, 200))   ; plot( gadm5, add=T)
  plot( r_Q4_DEU , main="Q" , xaxt="n", yaxt="n",
        col=map.pal("gyr"), range=c(   0, 1) ) ; plot( gadm5, add=T)
  plot( r_pH4_DEU, main="pH", xaxt='n',yaxt='n',
        col=map.pal("gyr"), range=c(   0, 1) ) ; plot( gadm5, add=T)
  plot( r_V4_DEU , main="V" , xaxt='n',yaxt='n',
        col=map.pal("gyr"), range=c( -30,50) ) ; plot( gadm5, add=T)
  plot( r_R4_DEU , main="R" , xaxt='n',yaxt='n',
        col=map.pal("gyr"), range=c( -30,50) ) ; plot( gadm5, add=T)

# BDT
  pA <- 0.9 ; pB_A <- 0.9 ; pB_a <- 0.6 ; pC_B <- 0.9 ; pC_b <- 0.2
  # Analytical solutions using Law of Total Probability (LTP) and Bayes' Theorem (BT):
  pB   <- pA*pB_A + (1-pA)*pB_a     # = 0.9 *0.9  + 0.1 *0.6 = 0.87        # LTP
  pC   <- pB*pC_B + (1-pB)*pC_b     # = 0.87*0.9  + 0.13*0.2 = 0.809       # LTP
  # pB_a <- 0.6                                                            # Prior
  pC_a <- pB_a*pC_B + (1-pB_a)*pC_b # = 0.6 *0.9  + 0.4 *0.2 = 0.62        # LTP
  pA_b <- pA*(1-pB_A) / (1-pB)      # = 0.9 *0.1  / 0.13     = 9/13        # BT
  # pC_b <- 0.2                                                            # Prior
    pc_A <- pB_A*(1-pC_B) + (1-pB_A)*(1-pC_b) # = 0.9*0.1 + 0.1*0.8 = 0.17 # LTP
  pA_c <- pA*pc_A     / (1-pC)      # = 0.9 *0.17 / 0.191    = 153/191     # BT
  pB_c <- pB*(1-pC_B) / (1-pC)      # = 0.87*0.1  / 0.191    = 87/191      # BT
  solA <- rbind( c(pB,pC), c(pB_a,pC_a), c(pA_b,pC_b), c(pA_c,pB_c) )
  # Numerical solutions:
  set.seed(1)
  pA   <- 0.9 ; pB_A <- 0.9 ; pB_a <- 0.6 ; pC_B <- 0.9 ; pC_b <- 0.2
  n    <- 1e6
  x    <- rbinom( n, 1, pA )
  y    <- rbinom( n, 1, x*pB_A + (1-x)*pB_a )
  z    <- rbinom( n, 1, y*pC_B + (1-y)*pC_b )
  pB   <- sum(y==1) / n                ; pC   <- sum(z==1) / n
  pB_a <- sum(x==0 & y==1) / sum(x==0) ; pC_a <- sum(x==0 & z==1) / sum(x==0)
  pA_b <- sum(x==1 & y==0) / sum(y==0) ; pC_b <- sum(y==0 & z==1) / sum(y==0)
  pA_c <- sum(x==1 & z==0) / sum(z==0) ; pB_c <- sum(y==1 & z==0) / sum(z==0) 
  solN <- rbind( c(pB,pC), c(pB_a,pC_a), c(pA_b,pC_b), c(pA_c,pB_c) )
  # Comparison:
  round( cbind( solA, solN ), 3 )

  u <- function( a, x=1, t=1, e=0, ka=0.2, kz=1 ) {
    z    <- t*(1-exp(-a-x)) + e
    cost <- ka*a ; benefit <- kz*z
    return( benefit - cost ) }

  par( mfrow=c(1,1), mar=c(4,2,3,0) )
  na    <- 31
  a.seq <- seq( 0, 3, length.out=na )
  u.seq <- u(a.seq)
  imaxu <- which( u.seq == max(u.seq) ) ; amaxu <- a.seq[imaxu]
  plot( a.seq, u.seq, xlab="a", ylab="", main="Utility\n(no uncertainty)",
        type="l" )
  abline( v=amaxu, col="red", lty=2 )
  text( amaxu, mean(u.seq), labels=paste("a_opt =",amaxu ), pos=4 )

  set.seed(1)
  np     <- 5e3
  x.smp  <- rnorm( np, 1  , 1   ) ; t.smp  <- rnorm( np, 1  , 0.5 )
  e.smp  <- rnorm( np, 0  , 1   )
  ka.smp <- runif( np, 0.1, 0.3 ) ; kz.smp <- runif( np, 0.5, 1.5 )
  Qlo    <- 0.25 ; Qhi <- 0.75

  kunc    <- c( 1/1.5, 1, 1.5 )
  x.smp1  <- 1   + (x.smp  - 1   ) * kunc[1]
  t.smp1  <- 1   + (t.smp  - 1   ) * kunc[1]
  e.smp1  <-        e.smp          * kunc[1]
  ka.smp1 <- 0.2 + (ka.smp - 0.2 ) * kunc[1]
  kz.smp1 <- 1   + (kz.smp - 1   ) * kunc[1]
  x.smp2  <- 1   + (x.smp  - 1   ) * kunc[2]
  t.smp2  <- 1   + (t.smp  - 1   ) * kunc[2]
  e.smp2  <-        e.smp          * kunc[2]
  ka.smp2 <- 0.2 + (ka.smp - 0.2 ) * kunc[2]
  kz.smp2 <- 1   + (kz.smp - 1   ) * kunc[2]
  x.smp3  <- 1   + (x.smp  - 1   ) * kunc[3]
  t.smp3  <- 1   + (t.smp  - 1   ) * kunc[3]
  e.smp3  <-        e.smp          * kunc[3]
  ka.smp3 <- 0.2 + (ka.smp - 0.2 ) * kunc[3]
  kz.smp3 <- 1   + (kz.smp - 1   ) * kunc[3]
  u.tbl1  <- u.tbl2 <- u.tbl3 <- NULL
  for(i in 1:np) {
    u.tbl1 <- rbind( u.tbl1,
      u( a.seq, x.smp1[i], t.smp1[i], e.smp1[i], ka.smp1[i], kz.smp1[i] ) )
    u.tbl2 <- rbind( u.tbl2,
      u( a.seq, x.smp2[i], t.smp2[i], e.smp2[i], ka.smp2[i], kz.smp2[i] ) )
    u.tbl3 <- rbind( u.tbl3,
      u( a.seq, x.smp3[i], t.smp3[i], e.smp3[i], ka.smp3[i], kz.smp3[i] ) )
  }
  uQlo.seq1 <- sapply( 1:na, function(i) { quantile( u.tbl1[,i], probs=Qlo ) } )
  uQlo.seq2 <- sapply( 1:na, function(i) { quantile( u.tbl2[,i], probs=Qlo ) } )
  uQlo.seq3 <- sapply( 1:na, function(i) { quantile( u.tbl3[,i], probs=Qlo ) } )
  umn.seq1  <- colMeans(u.tbl1)
  umn.seq2  <- colMeans(u.tbl2)
  umn.seq3  <- colMeans(u.tbl3)
  uQhi.seq1 <- sapply( 1:na, function(i) { quantile( u.tbl1[,i], probs=Qhi ) } )
  uQhi.seq2 <- sapply( 1:na, function(i) { quantile( u.tbl2[,i], probs=Qhi ) } )
  uQhi.seq3 <- sapply( 1:na, function(i) { quantile( u.tbl3[,i], probs=Qhi ) } )
  imaxmnu.1 <- which( umn.seq1 == max(umn.seq1) ) ; amaxmnu.1 <- a.seq[imaxmnu.1]
  imaxmnu.2 <- which( umn.seq2 == max(umn.seq2) ) ; amaxmnu.2 <- a.seq[imaxmnu.2]
  imaxmnu.3 <- which( umn.seq3 == max(umn.seq3) ) ; amaxmnu.3 <- a.seq[imaxmnu.3]

  par( mfrow=c(1,3), mar=c(4,2,3,0) )
  u.range <- range( uQlo.seq1, umn.seq1, uQhi.seq1,
                    uQlo.seq2, umn.seq2, uQhi.seq2,
                    uQlo.seq3, umn.seq3, uQhi.seq3 )
  plot  ( a.seq, umn.seq1 , xlab="a", ylab="", type="l", 
          main="Utility\n(low uncertainty)", ylim=u.range )
  points( a.seq, uQlo.seq1, type="l", lty=2 )
  points( a.seq, uQhi.seq1, type="l", lty=2 )
  abline( v=amaxmnu.1, col="red", lty=2 )
  text( amaxmnu.1, min(umn.seq1), labels=paste("a_opt =",amaxmnu.1 ), pos=2 )
  plot  ( a.seq, umn.seq2 , xlab="a", ylab="", type="l", 
          main="Utility", ylim=u.range )
  points( a.seq, uQlo.seq2, type="l", lty=2 )
  points( a.seq, uQhi.seq2, type="l", lty=2 )
  abline( v=amaxmnu.2, col="red", lty=2 )
  text( amaxmnu.2, min(umn.seq2), labels=paste("a_opt =",amaxmnu.2 ), pos=4 )
  plot  ( a.seq, umn.seq3 , xlab="a", ylab="", type="l", 
          main="Utility\n(high uncertainty)", ylim=u.range )
  points( a.seq, uQlo.seq3, type="l", lty=2 )
  points( a.seq, uQhi.seq3, type="l", lty=2 )
  abline( v=amaxmnu.3, col="red", lty=2 )
  text( amaxmnu.3, min(umn.seq3), labels=paste("a_opt =",amaxmnu.3 ), pos=2 )

  