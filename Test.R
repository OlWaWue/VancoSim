library(mrgsolve)
library(ggplot2)
library(dplyr)
library(plyr)

vanc_modl <- mread("Goti_et_al")


dosing <- ev(id=1, 
             amt=1000,
             rate=1000/1,
             cmt=2,
             addl=0,
             ii=12,
             CrCL=120,
             WT=70,
             DIAL=0)

sim_res <- vanc_modl %>%
  ev(dosing) %>%
  zero_re %>%
  mrgsim(end=72)

data <- as.data.frame(sim_res)


ggplot(data) + geom_line(aes(x=time, y=IPRED)) + geom_point(data=obs, aes(x=time, y=y))

library('rjags')

doses_prior_to_t <- list()

dosing_time = c(0,12,24)

dosing_time < 15

times=seq(0,72,0.1)

tdm_times = c(10, 23)

temp_times <- NULL
for(i in 1:length(tdm_times)){
  temp_times <- c(temp_times, which(times == tdm_times[i]))
}

a <- estimate_etas_empiricalBayes()



estimate_etas_empiricalBayes <- function(mod = vanc_modl,
                                         init = c(ETA1=0.1, ETA2=0.1, ETA3=0.1),
                                         data = dosing,
                                         obs = data.frame(time=5, y=17)) {
  
  ### This is the objective function to be minimized
  mapbayes <- function(eta, y, d, m){
    
    eta <- as.list(eta)
    names(eta) <- names(init)
    eta_m <- eta %>% unlist %>% matrix(nrow=1)
    m %<>% param(eta)
    
    out <- m %>%  
      param(eta) %>%
      obsonly %>%               #Simulate only for timepoint where tdm was performed
      zero_re() %>%             #Omega and Sigma Matrices are dropped from the model
      data_set(d)  %>%          #Data set used to simulate the Data 
      mrgsim
    
    sig2j <- (out$IPRED*(sigma[1,1])+(sigma[2,2]))^2 ## Residual variance
    sqwres <- log(sig2j) + (1/sig2j)*(y-out$IPRED)^2 ## square sum for -2LL estimate Bauer et al. (2007) AAPSJ E60
    nOn <- diag(eta_m %*% omega.inv %*% t(eta_m)) ## 
    return(sum(sqwres)+nOn)                       ## -2LL empirical bayes estimate objective function value, see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
  }
  
  
  #Retrieve omega matrix from model for MAP Bayesian estimation
  omega <- (omat(mod, make=T))
  
  
  omega <- diag(c((omega[1,1]),(omega[2,2]), (omega[3,3])))
  #Invert Omega Matrix for Standard Error estimation
  omega.inv <- solve(omega)
  
  sigma <- diag(0.227, 3.4)
  

  
  t.obs=obs$time
  y.obs=obs$y
  #simulate concentration at time of concentration measurement 
  #to get a data frame containing information for the fitting process
  sim <- mod %>% 
    omat(omega) %>%       #Omega matrix is reset, if it has been deleted from the model earlier
    data_set(data) %>%
    carry.out(amt, evid, cmt) %>% #important to carry out EVERY information on dosing and covariates
    mrgsim(end=-1, add=t.obs) #only simulate for time(s) with an observed concentration
  
  #Dataframe used for simulation process
  sim <- as.data.frame(sim)
  
  #Dataframe used for fitting IMPORTANT: drop parameters to be predicted
  fit_data <- sim %>% select(ID, time, evid, amt, cmt,  IPRED)
  
  fit_data
  #minimize the mapbayes function and save ETAs to fit$par
  #y= measured concentration(s)
  #m= model object to be used
  #d= data.frame containing information on dosing and patient
  fit <- optim(init, fn=mapbayes, y=y.obs, m=mod, d=fit_data)
  # save individual etas for the report
  fit$par
  est_eta <- fit$par
  names(est_eta) <- c("ETA1", "ETA2", "ETA3")
  
  ret = list(post_etas = est_eta)
}


sim_res <- vanc_modl %>%
  ev(dosing) %>%
  param(a$post_etas) %>%
  zero_re %>%
  mrgsim(end=72)

data <- as.data.frame(sim_res)


ggplot(data) + geom_line(aes(x=time, y=IPRED)) + geom_point(data=obs, aes(x=time, y=y)) + ylim(c(0,69)) + xlim(c(0,24))


tdm_data = data.frame(time=tdm_times, conc=c(10, 20))

jags.mod <- jags.model('Goti_et_al.bug',
                  data = list('c' = c(10, 20),
                              'amt' = c(1000), 
                              'dosing_time' = dosing_time,
                              't_inf' = c(0.5),
                              'CLCR'=time_dep_cov$ClCr, 
                              'WT'=time_dep_cov$WT, 
                              'DIAL'=time_dep_cov$DIAL,
                              'times'= times,
                              'tdm_times' = temp_times,
                              "theta"=thetas,
                              "omega"=omegas,
                              'sigma'=sigma ),
                  n.chains = 4,
                  n.adapt = 1000)


d <- coda.samples(jags.mod,
                  c('eta1', 'eta2', 'eta3'),
                  1000, thin=1)

df <- as.data.frame(d[[4]])

mcmc_se <- list()

for (i in 1:nrow(df)) {

  temp_dat <- pk_2cmt_infusion(theta = thetas,
                        CLCR=time_dep_cov$ClCr, WT=time_dep_cov$WT, DIAL=time_dep_cov$DIAL,
                        eta = c(df$eta1[i], df$eta2[i], df$eta3[i]),
                        dosing_events = dosing_events,
                        times=seq(0,72,0.1))
  
  mcmc_se[[i]] <- temp_dat$IPRED
}



df_temp <- NULL

for (k in 1:nrow(df)) {
  df_temp <- rbind(df_temp, mcmc_se[[k]])
  
}

s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.05, 0.10, 0.15, 0.20, 0.8, 0.85, 0.9, 0.95, 0.5)))

## Combine individual PK data and quantils in a data.frame
pk_data <- data.frame(time=seq(0,72,0.1),
                      s1=s[1,],s2=s[8,], 
                      s3=s[2,],s4=s[7,],
                      s5=s[3,],s6=s[6,],
                      s7=s[4,],s8=s[5,],
                      max=s[9,]) # median 


## Build raw individual PK plot
p <- ggplot(pk_data) + 
  geom_ribbon(aes(ymin=s1, ymax=s2, x=time), fill="red", alpha=0.15) + 
  geom_ribbon(aes(ymin=s3, ymax=s4, x=time), fill="red", alpha=0.15) + 
  geom_ribbon(aes(ymin=s5, ymax=s6, x=time), fill="red", alpha=0.15) + 
  geom_ribbon(aes(ymin=s7, ymax=s8, x=time), fill="red", alpha=0.15) + 
  geom_line(aes(y=max, x=time))+
  geom_point(data=tdm_data, aes(x=time, y=conc)) +
  geom_line(data=res, aes(x=times, y=IPRED), linetype=2)

plot(p)

