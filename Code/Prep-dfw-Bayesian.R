# using NIMBLE ------------------------------------------------------------

## NIMBLE model
## Leroux model
Leroux_model <- nimbleCode({
  for(i in 1:Nzip){
    O[i] ~ dbern(p[i])
    logit(p[i]) <- alpha +delta[i]+beta[1]*dep[i]+beta[2]*gini[i]+beta[3]*black[i]+beta[4]*hisp[i]+beta[5]*asian[i]+beta[6]*hiv[i]+beta[7]*popdens[i]
    delta[i] <- s[i] - mean.s
  }
  s[1:Nzip] ~ dmnorm(mean=mu.s[1:Nzip],prec=prec.s[1:Nzip,1:Nzip])
  prec.s[1:Nzip,1:Nzip] <- tau.s*(rho*Q[1:Nzip,1:Nzip]+(1-rho)*diag(Nzip))
  mean.s <- sum(s[1:Nzip])/Nzip
  
  alpha ~ dflat()
  for(k in 1:Nbeta){
    beta[k] ~ dnorm(0,sd=100)
  }
  
  # sd.s ~ T(dnorm(0,sd=1),0,)
  # tau.s <- 1/((sd.s)^2)
  tau.s ~ dgamma(1,0.01)
  sd.s <- sqrt(1/tau.s)
  
  rho ~ dunif(0,1)
})

## NIMBLE data
O <- sapply(dfw.zcta.output2$transit.dis,function(x) ifelse(x<=30,1,0))
# O <- sapply(dfw.zcta.output2$transit.dis,function(x) ifelse(x<=60,1,0))
Nzip <- nrow(dfw.zcta.output2)
Nbeta <- 7
dep <- scale(dfw.zcta.output2$deprive,T,T)[,1]
gini <- scale(dfw.zcta.output2$ACS_GINI_INDEX_ZC,T,T)[,1]
uninsured <- scale(dfw.zcta.output2$ACS_PCT_UNINSURED_ZC,T,T)[,1]
black <- scale(dfw.zcta.output2$ACS_PCT_BLACK_ZC,T,T)[,1]
hisp <- scale(dfw.zcta.output2$ACS_PCT_HISPANIC_ZC,T,T)[,1]
asian <- scale(dfw.zcta.output2$ACS_PCT_ASIAN_ZC,T,T)[,1]
hiv <- scale(dfw.zcta.output2$risk,T,T)[,1]
popdens <- scale(dfw.zcta.output2$CEN_POPDENSITY_ZC,T,T)[,1]

neighs <- poly2nb(dfw.zcta.output2,queen = T)
adj <- unlist(neighs)
sumNumNeigh <- length(adj)
num <- lengths(neighs)

W <- nb2mat(neighs, zero.policy = TRUE, style = "B")
Q <- matrix(0,nrow = nrow(W),ncol = ncol(W))

diag(Q) <- apply(W,1,sum)
Q <- Q-W

constant.nimble <- list(dep=dep,gini=gini,popdens=popdens,black=black,hisp=hisp,asian=asian,hiv=hiv,Nbeta=Nbeta,Nzip=Nzip,Q=Q,mu.s=rep(0,Nzip)) 

data.nimble <- list(O=O)

## NIMBLE initials
initials.leroux <- list(alpha=rnorm(1),beta=rnorm(Nbeta),s=rnorm(Nzip),tau.s=rgamma(1,1,0.01),rho=runif(1,0,1))

nimbleOptions('showCompilerOutput' = TRUE)

## NIMBLE model running
set.seed(1)
##Build and compile the model
bym.model <- nimbleModel(code=Leroux_model,constants=constant.nimble,data=data.nimble,inits = initials.leroux)
bym.model$calculate() ##calculate the log probability
bym.cModel <- compileNimble(bym.model)
##Configure, build, and compile MCMC
bym.confMCMC <- configureMCMC(bym.model,enableWAIC = T)
bym.confMCMC$addMonitors(c("alpha","beta","s","tau.s","rho"))
bym.MCMC <- buildMCMC(bym.confMCMC)
bym.cMCMC <- compileNimble(bym.MCMC,project = bym.cModel)
##Run MCMC via "runMCMC"
bym.samples <- runMCMC(bym.cMCMC,nburnin = 500000,niter = 1000000,nchains = 2,thin = 100,setSeed = T,samplesAsCodaMCMC = T)

## Model diagnosis
ggs_traceplot(ggs(bym.samples,family = "alpha"))
ggs_traceplot(ggs(bym.samples,family = "beta"))
ggs_traceplot(ggs(bym.samples,family = "tau.s"))
ggs_traceplot(ggs(bym.samples,family = "rho"))
ggs_traceplot(ggs(bym.samples,family = "s\\[1\\]"))
ggs_traceplot(ggs(bym.samples,family = "u\\[1\\]"))
ggs_traceplot(ggs(bym.samples,family = "s\\[21\\]"))
ggs_traceplot(ggs(bym.samples,family = "u\\[12\\]"))

gelman.diag(subset(bym.samples,pars = "alpha"))
gelman.diag(subset(bym.samples,pars = "beta"))
gelman.diag(subset(bym.samples,pars = "tau.s"))
gelman.diag(subset(bym.samples,pars = "rho"))

GR.diag <- gelman.diag(bym.samples, multivariate = FALSE)
all(GR.diag$psrf[,"Point est."] < 1.1) 

which(GR.diag$psrf[,"Point est."]>1.1)

##WAIC
calculateWAIC(bym.cMCMC)

for(i in 1:Nbeta){
  beta.family <- paste0("beta\\[",i)
  beta.samples <- ggs(bym.samples,family = beta.family)
  print(mean(beta.samples$value))
  print(quantile(beta.samples$value,probs = c(0.025,0.975)))
}

sd.s.samples <- ggs(bym.samples,family = "tau.s")
mean(sd.s.samples$value)

rho.samples <- ggs(bym.samples,family="rho")
mean(rho.samples$value)
