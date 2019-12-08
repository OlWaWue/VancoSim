library('rjags')
library('ggplot2')
library('readxl')
library('shiny')
library('gridExtra')

THETAS = c(4.5,  ## THETA1 -  clearance (L/h)
           58.4, ## THETA2 - volume of central compartment (L)
           38.4, ## THETA3 - volume of peripheral compartment (L)
           6.5,  ## THETA4 - intercompartmental clearance (L/h)
           0.8,  ## THETA5 - CrCl on CL
           0.7,  ## THETA6 - DIAL on CL
           0.5)  ## THETA7 - DIAL on Vc

PARAMS = c(70,    ## Body weight in kg
           120,   ## CrCl in mL/min
           0)     ## Dialysis yes or no

OMEGAS = c(0.398,
           0.816,
           0.571)

SIGMA = c(0.227)

dosing_events = data.frame(time=c(0, 12, 24),
                           amt=c(1000, 1000, 1000),
                           t_inf=c(0.5, 0.5, 0.5))

ETAS = c(0,0,0)





## analytical solution of 2cmt model multiple doses 
pk_2cmt_infusion <- function(theta, params, eta, dosing_events, times){
  
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
  
  dat <- data.frame(times=times, IPRED=IPRED, CL_i=Cl, Vc_i=V1, Vp_i=V2, Q_i=Q)
  
  return(dat)
}

## Process a NONMEM like dataset
process_data_set <- function(pk_data = data.frame(time=c(0,4,6,12,30,50),
                                                  amt=c(100,".",".",100,".","."),
                                                  conc=c(".", 2, 3, ".", 1.5, 0.72),
                                                  evid=c(1, 0, 0, 1, 0, 0)),
                             n.burn=200, n.iter=1000, 
                             thetas  = c(4.5,  ## THETA1 -  clearance (L/h)
                                         58.4, ## THETA2 - volume of central compartment (L)
                                         38.4, ## THETA3 - volume of peripheral compartment (L)
                                         6.5,  ## THETA4 - intercompartmental clearance (L/h)
                                         0.8,  ## THETA5 - CrCl on CL
                                         0.7,  ## THETA6 - DIAL on CL
                                         0.5), ## THETA7 - DIAL on Vc
                             omegas = c(0.398,
                                        0.816,
                                        0.571),
                             params,
                             TIME = seq(0, 72, by=1), SIGMAS=c(0.227)) {
  
  
  
  # inner function to generate individual PK plot
  do_plot <- function(jags_result, nburn=n.burn){
    
    ## Use the fourth chain for the simulation
    df <- as.data.frame(jags_result[[4]])
    
    ## remove burnin iterations
    df <- df[-(1:nburn),]
    
    mcmc_se <- list()
    
    mcmc_Cl_ind <- list()
    mcmc_V1_ind <- list()
    mcmc_V2_ind <- list()
    
    ## Simulate with the etas obtained by sampling from the posterior distribution
    withProgress(message = "Processing MCMC results", max = nrow(df), {
      for (i in 1:nrow(df)) {
        
          
          temp_dat <- pk_2cmt_infusion(theta = thetas,
                                       params = params,
                                       eta = c(df$eta1[i], df$eta2[i], df$eta3[i]),
                                       dosing_events = dosing_events,
                                       times=TIME)
          
          mcmc_se[[i]] <- temp_dat$IPRED
          mcmc_Cl_ind[[i]] <- temp_dat$CL_i
          mcmc_V1_ind[[i]] <- temp_dat$Vc_i
          mcmc_V2_ind[[i]] <- temp_dat$Vp_i
        
        incProgress(1)
      }
    })
    
    
    df_temp <- NULL
    
    withProgress(message = "Calculating Prediction interval", max = nrow(df), {
      ## bind simulations in a data.frame
      for (k in 1:nrow(df)) {
        df_temp <- rbind(df_temp, mcmc_se[[k]])
        incProgress(1)
      }
    })
    
    ## Generate Quantils
    s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.05, 0.10, 0.15, 0.20, 0.8, 0.85, 0.9, 0.95, 0.5)))
    
    ## Combine individual PK data and quantils in a data.frame
    pk_data <- data.frame(time=TIME,
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
      geom_point(data=tdm_data, aes(x=time, y=conc))
    
    ## Get the last simulated concentration for the boxplot
    c_at_tlast <- df_temp[,ncol(df_temp)]
    
    ## carry on max value in PK plot for adjusting y-axis
    ind_y_max <- max(pk_data$s2)
    
    ind_y_min <- min(pk_data$s1[pk_data$s1>0])
    
    return(list(p, c_at_tlast, ind_y_max, ind_y_min, mcmc_Cl_ind, mcmc_V1_ind, mcmc_V2_ind ))
  }
  
  ### Generate MCMC trace plots and distributions for every eta
  mcmc_diagnosticplots <- function(chain=1, jags_result, nburn = n.burn, omega, colour="red") {
    
    df <- as.data.frame(jags_result[[chain]]) ## Dataframe with etas
    
    nb.etas <- ncol(jags_result[[chain]]) 
    
    niter <- nrow(df)
    
    df$iteration <- (1:niter) 
    
    df <- df[-(1:nburn),]
    
    
    ## I will never again understand this...
    
    
    ## Problems that COULD happen: never refer aesthetics to columns via df[,n.et]
    ## Do not pass data to ggplot, but to the geom layer => problems with additional
    ## layers that are supposed to use fixed values like geom_vline
    plot_list <- list()
    
    for(n.et in 1:nb.etas) {
      
      dens_post <- density(df[,n.et], adjust=2)
      max_eta <- dens_post$x[which.max(dens_post$y)]
      
      dens_post <- data.frame(ETA=dens_post$x, freq=dens_post$y, max_eta=max_eta)
      x_this_eta <- seq(min(dens_post$ETA), max(dens_post$ETA), length.out = 512)
      
      dens_prior <- dnorm(x=x_this_eta, mean=0, sd=sqrt(omega[n.et]))
      dens_prior <- data.frame(ETA=x_this_eta, freq=dens_prior)
      
      current_data <- data.frame(iteration=(nburn+1):niter,eta=df[,n.et])
      
      y_name <- paste("ETA", n.et)
      
      
      ## trace plot for this eta in this chain
      
      plot_list[[paste("p_iter_ETA", n.et, sep="")]] <- ggplot(data=current_data) + 
        geom_line(aes(x=iteration,y=eta), colour=colour)+ 
        ylim(min(x_this_eta),
             max(x_this_eta)) +
        theme_bw() + ylab(y_name) + xlab("Iteration")
      
      
      ## 90° flipped density plot for this eta in this chain
      plot_list[[paste("p_dens_ETA", n.et, sep="")]] <- ggplot(data = dens_post) +
        geom_vline(aes(xintercept=max_eta), colour="red", linetype=2) +
        geom_line(aes(x=ETA, y=freq), colour="red") + 
        geom_line(data=dens_prior, aes(x=ETA,y=freq), colour = "blue") + 
        coord_flip() + 
        xlim(min(x_this_eta),
             max(x_this_eta))+
        annotate(geom="text", 
                 y=max(dens_post$freq)+0.2, 
                 x=max_eta, 
                 label="posterior", colour="red") +
        annotate(geom="text", y=max(dens_prior$freq)+0.2, x=0, label="prior", colour="blue") +
        ylim(0,2)+ ylab(paste("Density of", y_name)) +
        theme_bw() + theme(axis.title.y = element_blank()) 
      
      
    }
    
    df <- df[,-ncol(df)]
    
    plot_list[["chain_data"]] <- df
    
    
    
    return(plot_list)
  }
  
  
  
  ### separate dosing events ...
  dosing_events <- data.frame(time=as.numeric(as.character(pk_data[pk_data$evid==1,]$time)),
                              amt=as.numeric(as.character(pk_data[pk_data$evid==1,]$amt)),
                              dur=as.numeric(as.character(pk_data[pk_data$evid==1,]$dur)))
  
  
  ### ... and TDM events
  tdm_data <- data.frame(conc=as.numeric(as.character(pk_data[pk_data$evid==0,]$conc)),
                         time=as.numeric(as.character(pk_data[pk_data$evid==0,]$time)))
  
  
  ### Use JAGS model to sample from posterior distribution
  ### According to the selected model
  print(SIGMAS)
  #### ---- For some reason, JAGS runs more stable when SD is submitted to the model and is ----#
  #### Calculated back to variance in the model file --- ?? 
  jags <- jags.model('Goti_et_al.bug',
                    data = list('c' = tdm_data$conc,
                                'amt' = dosing_events$amt, 
                                'dosing_time' = dosing_events$time,
                                't_inf' = dosing_events$dur,
                                'times'= tdm_data$time,
                                'params'=params,
                                "theta"=thetas,
                                "omega"=omegas,
                                'sigma'=SIGMAS ),
                    n.chains = 4,
                    n.adapt = n.iter)
  
  d <- coda.samples(jags,
                    c('eta1', 'eta2', 'eta3'),
                    n.iter, thin=1)
  
  
  
  pk_profile <- (do_plot(d))
  
  
  mcmc_plots_1 <- mcmc_diagnosticplots(1, d, nburn=n.burn, omega=omegas, "red")
  mcmc_plots_2 <- mcmc_diagnosticplots(2, d, nburn=n.burn, omega=omegas, "orange")
  mcmc_plots_3 <- mcmc_diagnosticplots(3, d, nburn=n.burn, omega=omegas, "yellow")
  
  ## Use the fourth chain for the simulation
  df <- as.data.frame(d[[4]])
  
  ## remove burnin iterations
  df <- df[-(1:n.burn),]
  
  
  result = list(mcmc_plots_1, 
                mcmc_plots_2, 
                mcmc_plots_3, 
                pk_profile=pk_profile[[1]], 
                c_at_tlast=pk_profile[[2]], 
                ind_y_max=pk_profile[[3]], 
                ind_y_min=pk_profile[[4]], 
                mcmc_etas = df,
                mcmc_Cl_ind=pk_profile[[5]],
                mcmc_V1_ind=pk_profile[[6]],
                mcmc_V2_ind=pk_profile[[7]])
  
  return(result)
}


perform_mc_simulation <- function(n.mc, omegas, thetas, app_data, t_from, t_to) {
  

    mc_eta2 <- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[2])))
    mc_eta3<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[3])))
    mc_eta1<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[1])))
    
    all_etas <- data.frame(ETA1=mc_eta1, ETA2=mc_eta2, ETA3=mc_eta3)
  
  
  dat_mc <- NULL
  
  
  mc_se <- list()
  
  mc_Cl_ind <- list()
  mc_V1_ind <- list()
  mc_V2_ind <- list()
  

  
  ## Do with progress
  withProgress(message = "Performing Monte Carlo simulation", max = n.mc, {
    for(i in 1:n.mc){
      ### Simulate population PK profiles for every monte carlo sample according to the PK model currently
      ### Selected
      temp_dat <- pk_2cmt_infusion(theta = thetas,
                                   params = app_data$params,
                                   eta = c(all_etas$ETA1[i], all_etas$ETA2[i], all_etas$ETA3[i]),
                                   dosing_events = dosing_events,
                                   times=seq(t_from,t_to,0.2))
      
  
      
      mc_se[[i]] <- temp_dat$IPRED
      mc_Cl_ind[[i]] <- temp_dat$CL_i
      mc_V1_ind[[i]] <- temp_dat$Vc_i
      mc_V2_ind[[i]] <- temp_dat$Vp_i
        

      dat_mc <- cbind(dat_mc, mc_se[[i]])
      incProgress(1)
    }
  })
  
  ## transpose data -> rows into columns
  dat_mc <- t(dat_mc)
  
  ## Get quantile from monte carlo result
  s <- apply(dat_mc,2,function(x) quantile(x,probs=c(0.05,0.5,0.95)))
  
  ## Data for population PK Plot
  plot_dat <- data.frame(TIME=seq(t_from, t_to, by=0.2),CP_min=s[1,],CP=s[2,],CP_max=s[3,], DELTA=(s[3,]-s[1,]))
  return(list(plot_dat, dat_mc, all_etas,mc_Cl_ind, mc_V1_ind, mc_V2_ind))
}
