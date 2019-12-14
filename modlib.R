thetas = c(4.5,  ## THETA1 -  clearance (L/h)
            58.4, ## THETA2 - volume of central compartment (L)
            38.4, ## THETA3 - volume of peripheral compartment (L)
            6.5,  ## THETA4 - intercompartmental clearance (L/h)
            0.8,  ## THETA5 - CrCl on CL
            0.7,  ## THETA6 - DIAL on CL
            0.5)  ## THETA7 - DIAL on Vc

params = c(70,    ## Body weight in kg
           120,   ## CrCl in mL/min
           0)     ## Dialysis yes or no

omegas = c(0.398,
           0.816,
           0.571)

sigma = c(3.4)

dosing_events = data.frame(time=c(0, 12, 24),
                           amt=c(1000, 1000, 1000),
                           t_inf=c(0.5, 0.5, 0.5))

pat_CLCR = data.frame(time=c(0,19),
                        value=c(120,120))

pat_WT = data.frame(time=c(0,24),
                     value=c(70,75))

pat_DIAL = data.frame(time=c(0,8, 12),
                     value=c(0, 1, 0))

times=seq(0,72,0.1)





etas = c(0,0,0)




time_dep_cov <- data.frame(time=times)
time_dep_cov$ClCr <- 0
time_dep_cov$WT <- 0
time_dep_cov$DIAL <- 0



for(i in 1:nrow(pat_CLCR)){
  time_dep_cov$ClCr <- ifelse(pat_CLCR$time[i] <= time_dep_cov$time, pat_CLCR$value[i], time_dep_cov$ClCr)
}

if(nrow(pat_CLCR>1)){
  for(i in 1:(nrow(pat_CLCR)-1) ){
    temp_data <- data.frame(time=c(pat_CLCR$time[i], pat_CLCR$time[i+1]), value=c(pat_CLCR$value[i], pat_CLCR$value[i+1]))
    
    temp_lin_mod <- lm(value~time, data = temp_data)
    
    predicted_CrCL <- (predict(temp_lin_mod, newdata=data.frame(time=time_dep_cov$time)))
    
    time_dep_cov$ClCr <- ifelse(pat_CLCR$time[i+1] >= time_dep_cov$time, predicted_CrCL, time_dep_cov$ClCr)
    
    
  }
}

for(i in 1:nrow(pat_WT)){
  time_dep_cov$WT <- ifelse(pat_WT$time[i] <= time_dep_cov$time, pat_WT$value[i], time_dep_cov$WT)
}


if(nrow(pat_WT>1)){
  for(i in 1:(nrow(pat_WT)-1) ){
    temp_data <- data.frame(time=c(pat_WT$time[i], pat_WT$time[i+1]), value=c(pat_WT$value[i], pat_WT$value[i+1]))
    
    
    print(temp_data)
    temp_lin_mod <- lm(value~time, data = temp_data)
    
    predicted_WT<- (predict(temp_lin_mod, newdata=data.frame(time=time_dep_cov$time)))
    
    time_dep_cov$WT <- ifelse(pat_WT$time[i+1] >= time_dep_cov$time, predicted_WT, time_dep_cov$WT)
    
    
  }
}



for(i in 1:nrow(pat_DIAL)){
  time_dep_cov$DIAL <- ifelse(pat_DIAL$time[i] <= time_dep_cov$time, pat_DIAL$value[i], time_dep_cov$DIAL)
}


print(time_dep_cov)

res <- pk_2cmt_infusion(theta = thetas,
                        CLCR=time_dep_cov$ClCr, WT=time_dep_cov$WT, DIAL=time_dep_cov$DIAL,
                 eta = etas,
                 dosing_events = dosing_events,
                 times=seq(0,72,0.1))

res_alt <- pk_2cmt_infusion_alt(theta = thetas,
                        params = params,
                        eta = etas,
                        dosing_events = dosing_events,
                        times=seq(0,72,0.1))


ggplot(res) + geom_line(aes(x=times, y=IPRED)) + geom_line(aes(x=times, y=CLCr)) + geom_line(aes(x=times, y=DIAL)) + geom_line(aes(x=times, y=WT))
ggplot(res_alt) + geom_line(aes(x=times, y=IPRED))

## analytical solution of 2cmt model multiple doses 
pk_2cmt_infusion <- function(theta, CLCR, WT, DIAL, eta, dosing_events, times){
  
  dosing_time <- dosing_events[,1]
  amt <- dosing_events[,2]
  t_inf <- dosing_events[,3]
  

  
  

  
  
  IPRED <- vector()
  
  for(x in 1:length(times)){
    
    Cl <- theta[1]*( (CLCR[x]/120)^theta[5]) * (theta[6]^DIAL[x])*exp(eta[1])
    V1 <- theta[2]*(WT[x]/70) * (theta[7]^DIAL[x])*exp(eta[2])
    V2 <- theta[3] * exp(eta[3])
    Q <- theta[4]
    
    k=Cl/V1
    k12 = Q/V1
    k21 = Q/V2
    
    
    beta = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - 4 * k21 * k))
    
    alpha = (k21 * k) / beta
    
    A = 1/V1 *(alpha-k21)/(alpha-beta)
    B = 1/V1 *(beta-k21)/(beta-alpha)

      t=times[x]
        if(t==0) {
          temp_ipred =0
          } else{     
              

             
            doses_prior_to_t <- NULL
            
            for(i in 1:length(dosing_time)){
              if (dosing_time[i]<t){
                doses_prior_to_t <- c(doses_prior_to_t, i)
              }
                          
            }

            
            temp_delta_td <- dosing_time[length(doses_prior_to_t)]

            temp_delta_tinf <- t_inf[length(doses_prior_to_t)]

            temp_amt <- amt[length(doses_prior_to_t)]
          
            
            Di <- amt[doses_prior_to_t] 
            Di <- Di[-length(amt[doses_prior_to_t])]
            
            # Di = amt[dosing_time<t][-length(amt[dosing_time<t])]
            
            temp_delta_tdi <- dosing_time[doses_prior_to_t]
            temp_delta_tdi <- temp_delta_tdi[-length(temp_delta_tdi)]
            
            
            tinf_i <- t_inf[doses_prior_to_t]
            tinf_i <- tinf_i[-length(amt[doses_prior_to_t])]

            
            temp_ipred <- ifelse((t-temp_delta_td)<=temp_delta_tinf,
                                 sum( ( Di/tinf_i *(A/alpha * (1-exp(-alpha*tinf_i )) * exp(-alpha*(t-temp_delta_tdi-tinf_i ) ) + B/beta  * (1-exp(-beta*tinf_i ))  * exp(-beta*(t-temp_delta_tdi-tinf_i ) )) ) ) + ( temp_amt /temp_delta_tinf * (A/alpha * (1-exp(-alpha*(t-temp_delta_td))) + B/beta * (1-exp(-beta*(t-temp_delta_td))))),
                                 sum( amt[dosing_time<t]/t_inf[dosing_time<t] * (A/alpha * (1-exp(-alpha*t_inf[dosing_time<t])) * exp(-alpha*(t-dosing_time[dosing_time<t]-t_inf[dosing_time<t])) + B/beta * (1-exp(-beta*t_inf[dosing_time<t]))  * exp(-beta*(t-dosing_time[dosing_time<t]-t_inf[dosing_time<t]))) )
                                 )  
        }
        IPRED[x] <- temp_ipred
  }
  
  dat <- data.frame(times=times, IPRED=IPRED, CL_i=Cl, Vc_i=V1, Vp_i=V2, Q_i=Q, CLCr=time_dep_cov$ClCr, WT=time_dep_cov$WT, DIAL=time_dep_cov$DIAL)
  
  return(dat)
}


pk_2cmt_infusion_alt <- function(theta, params, eta, dosing_events, times){
  
  dosing_time <- dosing_events[,1]
  amt <- dosing_events[,2]
  t_inf <- dosing_events[,3]
  
  WT <- params[1]
  CrCL <- params[2]
  DIAL <- params[3]
  
  
  Cl <- theta[1]*( (CrCL/120)^theta[5]) * (theta[6]^DIAL)*exp(eta[1])
  V1 <- theta[2]*(WT/70) * (theta[7]^DIAL)*exp(eta[2])
  V2 <- theta[3] * exp(eta[3])
  Q <- theta[4]
  
  k=Cl/V1
  k12 = Q/V1
  k21 = Q/V2
  
  
  beta = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - 4 * k21 * k))
  
  alpha = (k21 * k) / beta
  
  A = 1/V1 *(alpha-k21)/(alpha-beta)
  B = 1/V1 *(beta-k21)/(beta-alpha)
  
  IPRED <- vector()
  
  
  for(x in 1:length(times)){
    
    t=times[x]
   
      
      
      
      temp_conc <- c(nrow = length(amt))
      for(i in 1:length(amt)) {
        
        temp_conc[i] <- ifelse((times[x]-dosing_time[i]) <= 0,
                               0 , 
                                ifelse(times[x]-dosing_time[i] <= t_inf[i], 
                                       amt[i]/t_inf[i] * ( (A/alpha)*(1-exp(-alpha*(times[x]-dosing_time[i]))) + (B/beta) * (1-exp(-beta*(times[x]-dosing_time[i])))), 
                                       amt[i]/t_inf[i] * (  (A/alpha)*(1-exp(-alpha*t_inf[i]) )*exp(-alpha*(times[x]-dosing_time[i]-t_inf[i])) + (B/beta)*(1-exp(-beta*t_inf[i]) )*exp(-beta*(times[x]-dosing_time[i]-t_inf[i])) )
                                       )
                               ) 
        
      }
    
      IPRED[x] <- sum(temp_conc)
  }
  
  
  dat <- data.frame(times=times, IPRED=IPRED, CL_i=Cl, Vc_i=V1, Vp_i=V2, Q_i=Q)
  
  return(dat)
}

