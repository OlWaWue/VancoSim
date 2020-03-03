library(mrgsolve)
library(ggplot2)
library(dplyr)
library(plyr)
library(dmutate)
library(magrittr)


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
  mrgsim(end=12)

data <- as.data.frame(sim_res)

obs = data.frame(time=5, y=17)

ggplot(data) + geom_line(aes(x=time, y=IPRED), size=1) + 
  geom_point(data=obs, aes(x=time, y=y), size=2, shape=1, colour="red", stroke=2.5) + theme_bw() +
  xlab("Time since last dose [h]") + ylab("Vancomycin Plasma concentration [mg/L]") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + ylim(0,60)




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
    
    sig2j <- (out$IPRED^2)*(0.227^2)+(3.4)^2 ## Residual variance
    sqwres <- log(sig2j) + (1/sig2j)*(y-out$IPRED)^2 ## square sum for -2LL estimate Bauer et al. (2007) AAPSJ E60
    nOn <- diag(eta_m %*% omega.inv %*% t(eta_m)) ## 
    return(sum(sqwres)+nOn)                       ## -2LL empirical bayes estimate objective function value, see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
  }
  
  
  #Retrieve omega matrix from model for MAP Bayesian estimation
  omega <- (omat(mod, make=T))
  
  
  omega <- diag(c((omega[1,1]),(omega[2,2]), (omega[3,3])))
  #Invert Omega Matrix for Standard Error estimation
  omega.inv <- solve(omega)
  
  sigma <- c(0.227, 3.4)
  

  
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


a <- estimate_etas_empiricalBayes()

sim_res <- vanc_modl %>%
  ev(dosing) %>%
  param(a$post_etas) %>%
  zero_re %>%
  mrgsim(end=12)

new_data <- as.data.frame(sim_res)




ggplot(data) + geom_line(aes(x=time, y=IPRED), linetype=2, size=1) + 
  geom_line(data=new_data, aes(x=time, y=IPRED), colour="red", size=1) + 
  geom_point(data=obs, aes(x=time, y=y), size=2, shape=1, colour="red", stroke=2.5) + theme_bw() +
  xlab("Time since last dose [h]") + ylab("Vancomycin Plasma concentration [mg/L]") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + ylim(0,100)


mcmc_data <- read.csv("temp.csv")

head(mcmc_data)

ggplot(data) + geom_line(aes(x=time, y=IPRED), linetype=2, size=1) +  
  geom_point(data=obs, aes(x=time, y=y), size=2, shape=1, colour="red", stroke=2.5) + theme_bw() +
  xlab("Time since last dose [h]") + ylab("Vancomycin Plasma concentration [mg/L]") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  geom_line(data = mcmc_data, aes(x=time, y=max), colour="red", size=1) +
  geom_ribbon(data= mcmc_data, aes(x=time, ymin=s1, ymax=s2), fill="red", alpha=0.1) +
  geom_ribbon(data= mcmc_data, aes(x=time, ymin=s3, ymax=s4), fill="red", alpha=0.1) +
  geom_ribbon(data= mcmc_data, aes(x=time, ymin=s5, ymax=s6), fill="red", alpha=0.1) +
  geom_ribbon(data= mcmc_data, aes(x=time, ymin=s7, ymax=s8), fill="red", alpha=0.1) +
  xlim(0,12) + ylim(0,100)


data <- data.frame(times=seq(0,12, by=0.1))


jags.mod <- jags.model('Goti_et_al.bug',
                  data = list('c' = c(17),
                              'amt' = c(1000), 
                              'dosing_time' = c(0),
                              't_inf' = c(1),
                              'CLCR'=c(120), 
                              'WT'=c(70), 
                              'DIAL'=c(0),
                              'times'= seq(0,12,by=0.1),
                              'tdm_times' = c(5),
                              "theta"=c(4.5,  ## THETA1 -  clearance (L/h)
                                        58.4, ## THETA2 - volume of central compartment (L)
                                        38.4, ## THETA3 - volume of peripheral compartment (L)
                                        6.5,  ## THETA4 - intercompartmental clearance (L/h)
                                        0.8,  ## THETA5 - CrCl on CL
                                        0.7,  ## THETA6 - DIAL on CL
                                        0.5),  ## THETA7 - DIAL on Vc
                              "omega"= c(0.398,
                                         0.816,
                                        0.571),
                              'sigma'=c(0.227 , 3.4) ),
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

