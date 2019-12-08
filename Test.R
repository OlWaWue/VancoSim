library(mrgsolve)
library(ggplot2)

vanc_modl <- mread("Goti_et_al")


dosing <- ev(id=1, 
             amt=1000,
             rate=1000/0.5,
             cmt=2,
             addl=2,
             ii=12,
             CrCL=120,
             WT=70,
             DIAL=0)

sim_res <- vanc_modl %>%
  ev(dosing) %>%
  zero_re %>%
  mrgsim(end=72)

data <- as.data.frame(sim_res)


ggplot(data) + geom_line(aes(x=time, y=IPRED))

library('rjags')

doses_prior_to_t <- list()

dosing_time = c(0,12,24)

dosing_time < 15

times = c(10)



tdm_data = data.frame(time=times, conc=c(10))

ags <- jags.model('Goti_et_al.bug',
                  data = list('c' = c(10),
                              'amt' = c(1000,1000,1000), 
                              'dosing_time' = dosing_time,
                              't_inf' = c(0.5,0.5,0.5),
                              'times'= times,
                              'params'=params,
                              "theta"=thetas,
                              "omega"=omegas,
                              'sigma'=sigma ),
                  n.chains = 4,
                  n.adapt = 1000)

d <- coda.samples(ags,
                  c('eta1', 'eta2', 'eta3'),
                  1000, thin=1)

df <- as.data.frame(d[[4]])

mcmc_se <- list()

for (i in 1:nrow(df)) {

  temp_dat <- pk_2cmt_infusion(theta = thetas,
                        params = params,
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
