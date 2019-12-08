
shinyServer(function(input, output, session) {
  
  ## Load an Excel file with data
  observeEvent(input$loadPT,{
    inFile <- input$loadPT
    
    if (is.null(inFile))
      return(NULL)
    
    tryCatch({
      app_data$data_set <- read_xlsx(inFile$datapath)
    },
    error = function(e){
      showModal(modalDialog(
        title = label_error,
        HTML(paste("File not recognized!<br>Details:<br><br>",e)),
        easyClose = TRUE,
        footer = NULL
      ))
    })
    
  })
  
  app_data <- reactiveValues(
    
    ## AppData used in simulation
    
    mcmc_result = NULL,
    mc_result= NULL,
    data_set = data.frame(time=c(0,4,6,12,24,30,36,48),
                          amt=c(1000,".",".",1000,1000,".",1000,"."),
                          conc=c(".", 10, 15, ".",".", 16,".", 5),
                          evid=c(1, 0, 0, 1,1, 0, 1, 0),
                          dur=c("0.5",".",".","0.5","0.5",".","0.5",".")),
    params= c(70,    ## Body weight in kg
      120,   ## CrCl in mL/min
      0),     ## Dialysis yes or no
    
    demo_loaded = FALSE, # Flag shows whether demo simulation has been loaded
    pk_plots = NULL
  )
  
  ## This function takes care of all the simulations
  ## MCMC is outsourced to global.R
  updatePKPlot <- function(){
    app_data$mc_result <- perform_mc_simulation(100, ## number of simulations
                                                OMEGAS, ## omegas
                                                THETAS, ## thetas
                                                app_data, ## App Data for Dosing / TDM Data
                                                0, 72) ## Time to simulate
    
    app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = 200, n.burn = 20,
                                            thetas = THETAS,
                                            omegas = OMEGAS,
                                            params = PARAMS,
                                            TIME =seq(0, 72, by=0.2), SIGMAS=3.4) 
    
    
    plot_dat <- app_data$mc_result [[1]] ## get Plot data ...
    dat_mc <- app_data$mc_result [[2]]  ### .. and raw results from the mc simulation to complete the plots
    
    ## Get lowest value (non-zero <- log-scale) and max value in PK plot to adjust y-axis
    pop_y_max <- max(plot_dat$CP_max)
    pop_y_min <- min(plot_dat$CP_min[plot_dat$CP_min >0])
    ind_y_max <- app_data$mcmc_result[[6]]
    ind_y_min <- app_data$mcmc_result[[7]]
    
    
    ## get tdm data from the table
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    
    ## Check if TDM concentration is above upper prediction interval to readjust the plot limits on y-axis
    pop_y_max <- ifelse(max(tdm_data$conc) > pop_y_max, max(tdm_data$conc), pop_y_max)
    ind_y_max <- ifelse(max(tdm_data$conc) > ind_y_max, max(tdm_data$conc), ind_y_max)
    
    ## prepare individual boxplot
    ind_boxplot <- ggplot(data=data.frame(conc=app_data$mcmc_result[[5]], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank())+ 
      ggtitle(paste("C last at ", input$TIME[2], " h"), "Individual") + xlab("") + ylim(c(0,ind_y_max*1.2))
    
    ## prepare individual PK plot
    
    ind_plot <- app_data$mcmc_result[[4]] + theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") +
      ggtitle("MCMC Result including data (posterior)", "80/85/90/95% PI") + 
      geom_line(data=plot_dat, aes(x=TIME, y=CP), colour="blue", linetype=2) + ylim(c(0,ind_y_max*1.2))
    
    ## prepare population boxplot
    pop_boxplot <- ggplot(data=data.frame(conc=dat_mc[,ncol(dat_mc)], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank())+ ggtitle(paste("C last at ", input$TIME[2], " h") , "Population") + 
      xlab("") + ylim(c(0,pop_y_max*1.2))
    
    ## Prepare population PK plot
    pop_plot <- ggplot(data=plot_dat)  + geom_line(aes(x=TIME, y=CP), colour="blue") +
      geom_ribbon(aes(x=TIME, ymax=CP_max, ymin=CP_min), alpha=0.15, fill="blue") +
      theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") + ggtitle("MC Result without data (prior)", "95% PI") +
      geom_point(data=tdm_data, aes(x=time, y=conc)) + ylim(c(0,pop_y_max*1.2))
    
    
    plots <- list(ind_plot, pop_plot, ind_boxplot, pop_boxplot)
  }
  
  output$pkPlot <- renderPlot({
    
    
    grid.arrange(app_data$pk_plots[[2]], app_data$pk_plots[[4]], 
                 app_data$pk_plots[[1]], app_data$pk_plots[[3]], nrow=2, ncol=2,widths=c(4,1))
  })
  
  output$data_set <- renderTable({
    
    display_data <-app_data$data_set
    
    display_data
    
  })
  
  output$cov_plot <- renderPlot({
    
    
    
  })
  
  output$pop_cov_plot <- renderPlot({
    
    
  })
  
  output$par_dist <- renderPlot({
    
   # print(app_data$mcmc_result[[9]])
    
#    dist_plots[[1]] <- ggplot() +theme_bw()+ 
#      geom_density(data=pop_pars, aes(x=pop_ka, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
#      geom_density(data=ind_pars, aes(x=ind_ka, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1) +
#      xlab("Absorption rate constant (ka) [1/h]") + ylab("Frequency")
#    dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
#      geom_density(data=pop_pars, aes(x=pop_Vc, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
#      geom_density(data=ind_pars,aes(x=ind_Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
#      xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")
#    dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
#      geom_density(data=pop_pars, aes(x=pop_Cl, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
#      geom_density(data=ind_pars,aes(x=ind_Cl, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
#      xlab("Clearance (Cl) [L/h]") + ylab("Frequency")
#    dist_plots[[4]] <- ggplot(  ) +theme_bw()+ 
#      geom_density(data=pop_pars, aes(x=pop_Vp, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
#      geom_density(data=ind_pars,aes(x=ind_Vp, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
#      xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")
#    dist_plots[[5]] <- ggplot(  ) +theme_bw()+ 
#      geom_density(data=pop_pars, aes(x=pop_Q, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
#      geom_density(data=ind_pars,aes(x=ind_Q, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
#      xlab("Intercompartmental Clearance (Q) [L/h]") + ylab("Frequency")
#    dist_plots[[6]] <- ggplot(  ) +theme_bw()+ 
#      geom_density(data=pop_pars, aes(x=pop_F_oral, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
#      geom_density(data=ind_pars,aes(x=ind_F_oral, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
#      xlab("Systemically available fraction") + ylab("Frequency")
    
#    grid.arrange(dist_plots[[1]],dist_plots[[2]],dist_plots[[3]],dist_plots[[4]],dist_plots[[5]],dist_plots[[6]], nrow=3, ncol=2)
  
  })
  
  output$traceplot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    ## Has to be generalized 
    
    
    ## 
      gridExtra::grid.arrange(app_data$mcmc_result[[chain]]$p_iter_ETA1, app_data$mcmc_result[[chain]]$p_dens_ETA1,    
                              app_data$mcmc_result[[chain]]$p_iter_ETA2, app_data$mcmc_result[[chain]]$p_dens_ETA2, 
                              app_data$mcmc_result[[chain]]$p_iter_ETA3, app_data$mcmc_result[[chain]]$p_dens_ETA3, nrow=4, ncol=2, widths=c(3,1))  
    
    
    
  })
  
  
  output$modelfile <- renderText({
    
    ## Show the correct JAGS model file
    

      includeText("Goti_et_al.bug")
    
    
  })
  
  output$info <- renderText({
    
    ## Show a short info about the currently selected PK Model

    paste("Goti et al. Ther Drug Monit (2018) 40:212–221 ")

  })
  
  observeEvent(input$submit,
               ## Submit changes
               app_data$pk_plots <- updatePKPlot()
  )
  
  
  
})