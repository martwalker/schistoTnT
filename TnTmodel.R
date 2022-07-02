################################################################################
## TnT model
################################################################################

################################################################################
## model equations
################################################################################
W.dyn <- function(t,y,par){
  with(as.list(c(y,par)),{
    
    # is frequency of treatment divisible by dt
    if((parameters["freq.tx"]/parameters["dt"])%%1!=0) {
      stop("frequency of treatment must be divisible bt dt - change dt")
    }
    
     R0<-parameters["R0"]
     R0sh <- R0/R0hs
    
    # y and dy/dt empty matrices
    ymat<- matrix(y, ncol=1, nrow=nw)
    dymat <- matrix(0, ncol=1, nrow=nw)
    
    # indicator variables
    W.id <- seq(1:nw)
    Wj <- vector("numeric", 1)
    
    ## total contribution to transmission
    Wj = sum(ymat)  #sum up all nw compartments
   
    # dynamic overdispersion parameter for the distribution of worms (NBD)
    kdyn <- as.vector(e$kdyn)
    kpost <- as.vector(e$kpost)
    kinf <- as.vector(e$kinf)
    Wpost <- as.vector(e$Wpost)
    Winf <- as.vector(e$Winf)
   
    # all worm prevalence
    prevold<- 1-(1+(Wj)/kdyn)^-kdyn 
    
    # mating probability
    mmat<- matingf(Wj,kdyn)
    
    # force of infection
    foi= foi_(mmat, R0sh, R0hs, N1N2,mu0,
              mu1, mu2, Wj)
    
    # effective reproduction number
    Re = as.vector(Re_(mmat, R0sh, R0hs, N1N2, mu0, mu1, mu2, Wj))
    
    for (i in 1:nw) { #for worm compartment (i) in worm compartments W1:Wnw
      if (i==1) { #i.e. W1
        dymat[i,1] = foi - (nw*(mu1+mu0))*ymat[i,1]
      } else if (i>1) { #i.e. W1..Wnw
        dymat[i,1] = nw*(mu0+mu1)*ymat[i-1,1] - nw*(mu0+mu1)*ymat[i,1]
      }
    }
    
    if (parameters["dotx"]==1) {
      ## is it a treatment time?
      if (t>=parameters["start.tx"] & t<(parameters["start.tx"]+parameters["n.tx"]*parameters["freq.tx"]) & 
          t %% parameters["freq.tx"]==0 ) {
        ## update k dynamic, Wjpost and kpost
        kdyn <- kinf*Wj/( (1+kinf)*Winf - kinf*Wj)
        Wpost <- Wj
        kpost <- kdyn
      } 
      else if(t>parameters["start.tx"]) {
       kdyn <- Wj^2*(Winf-Wpost)^2/
         ((Winf^2/kinf)*(Wj-Wpost)^2 + (Wpost^2/kpost)*(Wj-Winf)^2)
      }
   
    }
    
    #update prevalence
    prevnew<- 1-(1+(Wj*2)/kdyn)^-kdyn  
    prevold <- prevnew
  
    ## update all dynamic parameters
    e$kdyn <<- kdyn
    e$kinf <<- kinf
    e$kpost <<- kpost
    e$Wpost <<- Wpost
    e$Winf <<- Winf 
    
    ## PPV & NPV
    PPV <- as.numeric( parameters["se"]*prevold/( (parameters["se"]*prevold)+(1-parameters["sp"])*(1-prevold) ) )
    NPV <- as.numeric( parameters["sp"]*(1-prevold)/( (1-parameters["se"])*prevold+(parameters["sp"])*(1-prevold) ))
    
    return(list(rbind(dymat),
                Wj=Wj, prev=prevnew, PPV=PPV, NPV=NPV, k=kdyn))
  })
}


################################################################################
## mating probability, FoI and Re
################################################################################

## mating probability
matingf <- function(M,k){
  
  alpha= M/(k+M)
  
  integrand <- function(theta)  
  {(1-cos(theta))/(((1+alpha*cos(theta))^(1+k)))}
  intsol <- integrate(integrand, lower=0, upper=2*pi)
  
  int <- intsol$value
  I= (((1-alpha)^(1+k))/(2*pi))
  I1= 1-I*int
  
  return(I1)
}

## force of infection
foi_= function(matingf, R0sh, R0hs, N1N2,mu0,
               mu1, mu2, W){
  
  
  FOIN<- R0hs*R0sh*(mu0+mu1)*mu2*matingf*W
  FOID<- (R0hs*(mu0+mu1)*N1N2*W*matingf+mu2)
  FOI<- FOIN/FOID
  
  return(FOI)
}

## Re
Re_= function(matingf, R0sh, R0hs, N1N2,mu0,
              mu1, mu2, W){
  R0sh*R0hs*mu2*matingf/(R0hs*(mu0+mu1)*N1N2*(W)*matingf+mu2)
}

####################################################
## treatment event function
#####################################################

eventfun <- function(t, y, parameters) {
  with(as.list(y), {
    
    kdyn <- e$kdyn
    
    ymat <- matrix(y, ncol=1, nrow=parameters["nw"])
    
    Wj = sum(ymat)  #Sum up all nw compartments
    
    prev<- 1-(1+(Wj*2)/kdyn)^-kdyn
    
    if (parameters["dotx"]==1) {
    ## is it a treatment time?
    if (t>=parameters["start.tx"] & t<(parameters["start.tx"]+parameters["n.tx"]*parameters["freq.tx"]) & 
        t %% parameters["freq.tx"]==0 ) {
      ## probability of a positive diagnostic result
      prob.detect <- 1-(1-parameters["se"]*prev-(1-parameters["sp"])*(1-prev))^parameters["n"]
      treatment <- prob.detect
        ymat[,1] <- ymat[1,] - ymat[,1]*parameters["coverage"]*parameters["efficacy"]*treatment
    } else {
      ymat[,1] <- ymat[,1]
    }
    }
    return(c(ymat))
  })
  
}

####################################################
## function to stop simulation if Re<1
#####################################################
rootfun <- function(t, y, parameters) {
  with(as.list(y), {
    ymat<- matrix(y, ncol=1, nrow=parameters["nw"])
    kdyn <- as.vector(e$kdyn)
    dd <- 1  
    Wj = sum(ymat)
    mmat<- matingf(Wj*2,kdyn)
    Re = as.vector(Re_(mmat, parameters["R0sh"], parameters["R0hs"], parameters["N1N2"],
                       parameters["mu0"], parameters["mu1"], parameters["mu2"], Wj))
    return(Re-0.99)
       
  })
}
  
################################################################################
## main function that runs model
################################################################################
runmod <- function(parameters, inits=NULL) {
  
  if (parameters["dotx"]==1 & is.null(inits)) {
    stop("Simulating treatment requires initial values to be equilibrium worm burden")
  }
  
  if (is.null(inits)) {
    ## approximate initial values
    Wj0 <- (parameters["mu2"]*(parameters["R0"]^2-1)/
               ( (parameters["mu0"]+parameters["mu1"])*parameters["R0hs"]*parameters["N1N2"] ))

    inits <- as.vector(rep( Wj0/parameters["nw"], 
      parameters["nw"]))
  }
  
  ## initialize dynamic variable environment
  e <<- new.env()
  e$kdyn <- as.vector(parameters["k"])
  e$kinf <- as.vector(parameters["k"])
  e$kpost <- as.vector(parameters["k"])
  e$Winf <- sum(inits)
  e$Wpost <- sum(inits)
  
  if (parameters["dotx"]==1) {
    if(parameters["stop.t"] <= parameters["start.tx"] + parameters["n.tx"]*parameters["freq.tx"]) {
      stop("treatments run longer than simulation time")
    }
    
    events = list(func= eventfun, time= seq(parameters["start.tx"],
                                            parameters["start.tx"]+parameters["n.tx"]*parameters["freq.tx"],
                                            parameters["freq.tx"]), parameters= parameters, root=T)
    
  } 
  
  if (is.vector(inits)) {
    y1 <- rep(as.numeric(inits), (parameters["nw"]))
  } else {
    y1<- rep(as.numeric(1/parameters["nw"]), (parameters["nw"]))
  }
  
  t <- seq(0,parameters["stop.t"],by=parameters["dt"])
  
  if(parameters["dotx"]==1) {
    out<-lsoda(y=y1,times=t, func=W.dyn, par=parameters,  events = events, rootfunc = rootfun)
  } else {
     out<-lsoda(y=y1,times=t, func=W.dyn, par=parameters)
  }
  
  out <- as.data.frame(out)
  
  every <- floor(parameters["dtout"]/parameters["dt"])
  
  out[seq(1,nrow(out), every),]
}


