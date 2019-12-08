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

times = c(10, 18)

for(j in 1:length(times)){
  temp = NULL
    for(i in 1:length(dosing_time)){
      
      if (dosing_time[i]<times[j]){
        temp <- c(temp, i)
      }
      
      
    }
  doses_prior_to_t[[j]] <- temp
}



ags <- jags.model('Goti_et_al.bug',
                  data = list('c' = c(10,15),
                              'amt' = c(1000,1000,1000), 
                              'dosing_time' = dosing_time,
                              't_inf' = c(0.5,0.5,0.5),
                              'times'= times,
                              'params'=params,
                              "theta"=thetas,
                              "omega"=omegas,
                              'sigma'=sigma ),
                  n.chains = 4,
                  n.adapt = n.iter)

d <- coda.samples(jags,
                  c('eta1', 'eta2', 'eta3', 'eta4'),
                  n.iter, thin=1)