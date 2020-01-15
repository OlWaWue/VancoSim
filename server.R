
shinyServer(function(input, output, session) {
  
  ## Load an Excel file with data
  observeEvent(input$loadPT,{
    inFile <- input$loadPT
    
    if (is.null(inFile))
      return(NULL)
    
    tryCatch({
      
      path_info <- read_excel(inFile$datapath, sheet = "PATHOGEN")
      
      pat_info <- read_excel(inFile$datapath, sheet = "PATIENT")
      
      updateTextInput(session, inputId = "pat_ID", value = pat_info$ID[1])
      updateNumericInput(session, inputId = "WT", value = as.numeric(pat_info$WT[1]))
      updateNumericInput(session, inputId = "CRCL", value = as.numeric(pat_info$CRCL[1]))
      updateCheckboxInput(session, inputId = "has_dialysis", value=pat_info$DIAL[1])
      
      updateNumericInput(session, inputId = "MIC", value = as.numeric(path_info$MIC[1]))
      updateSelectInput(session, inputId = "choose_pathogen", selected = as.numeric(path_info$PATHOGEN[1]))
      
      temp <- convert_xlsx_to_NMTRAN(read_excel(inFile$datapath, sheet = "DATA_SET") )
      
      
      
      app_data$data_set <- temp$conv_data
      app_data$user_data_set <- temp$original_data
      app_data$time_reference <- temp$time_reference
      
      
      ii = 12
      doses = temp$conv_data[temp$conv_data$evid==1,]
      if(nrow(doses)>1){
        ii = doses$time[nrow(doses)] - doses$time[nrow(doses)-1]
      }
      
      if(nrow(temp$conv_data[temp$conv_data$evid==0,])>=1){
        app_data$tdm_samples_available <- TRUE
      } else {
        app_data$tdm_samples_available <- FALSE
      }
      
      app_data$last_known_dose <- tail(temp$conv_data[temp$conv_data$evid==1,],1)
      app_data$last_known_dose_orig <- tail(temp$original_data[temp$original_data$EVID==1,],1)
      
      app_data$last_known_dose$ii <- ii
      app_data$last_known_dose_orig$II <- ii
      
      
      
    },
    error = function(e){
      showModal(modalDialog(
        title = "ERROR",
        HTML(paste("File not recognized!<br>Details:<br><br>",e)),
        easyClose = TRUE,
        footer = NULL
      ))
    })
    
  })
  
  
  get_cov_from_dataset <- function(dataset, times) {
    
    
    
    WT <- dataset[dataset$evid == 2 & dataset$wt !=".",]
    CLCR <- dataset[dataset$evid == 2 & dataset$crcl !=".",]
    DIAL <- dataset[dataset$evid == 2 & dataset$dial !=".",]
    
    if(nrow(CLCR)==0){
      pat_CLCR = data.frame(time=c(0),
                            value=c(input$CRCL))
    
    } else {
      pat_CLCR = data.frame(time=c(0,CLCR$time),
                            value=c(input$CRCL, as.numeric(as.character(CLCR$crcl))))
    }
    
    if(nrow(WT)==0) {
      pat_WT = data.frame(time=c(0),
                          value=c(input$WT))
    } else {
      pat_WT = data.frame(time=c(0,WT$time),
                          value=c(input$WT,as.numeric(as.character(WT$wt))))
    }
    if(nrow(DIAL)==0){
      pat_DIAL = data.frame(time=c(0),
                            value=c(as.numeric(input$has_dialysis)))
    } else {
      pat_DIAL = data.frame(time=c(0,DIAL$time),
                            value=c(as.numeric(input$has_dialysis), as.numeric(as.character(DIAL$dial))))
    }
    
    
    
    time_dep_cov <- data.frame(time=times)
    time_dep_cov$CLCR <- 0
    time_dep_cov$WT <- 0
    time_dep_cov$DIAL <- 0
    
    
    
    for(i in 1:nrow(pat_CLCR)){
      time_dep_cov$CLCR <- ifelse(pat_CLCR$time[i] <= time_dep_cov$time, pat_CLCR$value[i], time_dep_cov$CLCR)
    }
    
    
    
    if(nrow(pat_CLCR)>1){
      for(i in 1:(nrow(pat_CLCR)-1) ){
        temp_data <- data.frame(time=c(pat_CLCR$time[i], pat_CLCR$time[i+1]), value=c(pat_CLCR$value[i], pat_CLCR$value[i+1]))
        
        temp_lin_mod <- lm(value~time, data = temp_data)
        
        predicted_CrCL <- (predict(temp_lin_mod, newdata=data.frame(time=time_dep_cov$time)))
        
        time_dep_cov$CLCR <- ifelse(pat_CLCR$time[i+1] >= time_dep_cov$time, predicted_CrCL, time_dep_cov$CLCR)

      }
    }
    
    for(i in 1:nrow(pat_WT)){
      time_dep_cov$WT <- ifelse(pat_WT$time[i] <= time_dep_cov$time, pat_WT$value[i], time_dep_cov$WT)
    }
    
    
    if(nrow(pat_WT)>1){
      for(i in 1:(nrow(pat_WT)-1) ){
        temp_data <- data.frame(time=c(pat_WT$time[i], pat_WT$time[i+1]), value=c(pat_WT$value[i], pat_WT$value[i+1]))
        
        
        temp_lin_mod <- lm(value~time, data = temp_data)
        
        predicted_WT<- (predict(temp_lin_mod, newdata=data.frame(time=time_dep_cov$time)))
        
        time_dep_cov$WT <- ifelse(pat_WT$time[i+1] >= time_dep_cov$time, predicted_WT, time_dep_cov$WT)
        
        
      }
    }
    
    
    
    for(i in 1:nrow(pat_DIAL)){
      time_dep_cov$DIAL <- ifelse(pat_DIAL$time[i] <= time_dep_cov$time, pat_DIAL$value[i], time_dep_cov$DIAL)
    }
   
    return(time_dep_cov)
  }
  
  observeEvent(input$but.resetApp, {

    resetApp()
    
  })
  
  resetApp <- function(){
    app_data$times_to_calculate = NULL
    app_data$mcmc_result = NULL
    app_data$mc_result= NULL
    app_data$disc_shown = F
    app_data$user_data_set = NULL
    app_data$data_set = NULL
    app_data$time_reference = NULL
    app_data$params= c(WT=NA,    
                       CRCL=NA,   
                       DIAL=NA)     
    
    app_data$covariates = NULL
    
    app_data$forecast_data_set = NULL
    app_data$user_forecast_data_set = NULL
    app_data$forecast_advanced = F
    
    app_data$step_one_completed = FALSE
    app_data$step_two_completed = FALSE 
    app_data$step_three_completed = FALSE 
    
    app_data$TDM_recommendation = NULL

    app_data$adapted_pk_plot_withTDM = NULL
    
    app_data$pk_plots = NULL
    app_data$adapted_pk_plot = NULL
    app_data$dist_plots = NULL
    app_data$tdm_samples_available =T
    app_data$last_known_dose = NULL
    app_data$last_known_dose_orig = NULL
    app_data$standard_y_zoom = NULL
    app_data$standard_x_zoom = NULL
    
    app_data$user_y_zoom = NULL
    app_data$user_x_zoom = NULL
    
    app_data$standard_y_zoom_pop = NULL
    app_data$standard_x_zoom_pop = NULL
    
    app_data$user_y_zoom_pop = NULL
    app_data$user_x_zoom_pop = NULL
    
    app_data$standard_y_zoom_adapt = NULL
    app_data$standard_x_zoom_adapt = NULL
    
    app_data$user_y_zoom_adapt = NULL
    app_data$user_x_zoom_adapt = NULL
    
    app_data$covariate_plot = NULL
    
    
    
    updateTextInput(session, inputId = "pat_ID", value = "PAT0001")
    updateNumericInput(session, inputId = "WT", value = 70)
    updateNumericInput(session, inputId = "CRCL", value = 120)
    updateNumericInput(session, inputId = "MIC", value = 5)
    
    updateNumericInput(session, inputId = "adapt.dose", value = 1000)
    updateNumericInput(session, inputId = "adapt.ii", value = 12)
    updateNumericInput(session, inputId = "adapt.dur", value = 60)
    updateNumericInput(session, inputId = "adapt.n", value = 5)
    
    updateSelectInput(session, inputId = "adapt.for", selected = 1)
    updateSelectInput(session, inputId = "adapt.what", selected = 1)
    updateSelectInput(session, inputId = "choose_pathogen", selected = 1)
    updateSelectInput(session, inputId = "choose_recommendation", selected = 1)
    updateSelectInput(session, inputId = "select_chain", selected = 1)
    
    updateCheckboxInput(session, inputId = "additional_tdm", value=F)
    updateCheckboxInput(session, inputId = "add_dur_info", value=F)
    updateCheckboxInput(session, inputId = "has_dialysis", value=F)
    
    updateTextAreaInput(session, inputId = "report_comment", value="")
    
    updateTabsetPanel(session, inputId = "mainpage", selected = "Enter Patient data")
  }
  
  app_data <- reactiveValues(
    
    ## AppData used in simulation
    
    mcmc_result = NULL,
    times_to_calculate = NULL,
    mc_result= NULL,
    covariates = NULL,
    disc_shown = F,
    user_data_set = NULL,
    data_set = NULL,
    time_reference = NULL,
    params= c(WT=NA,    ## Body weight in kg
              CRCL=NA,   ## CrCl in mL/min
              DIAL=NA),     ## Dialysis yes or no
    
    
    # Flags showing which stages were completed
    step_one_completed = FALSE, 
    step_two_completed = FALSE, 
    step_three_completed = FALSE, 
    
    forecast_data_set = NULL,
    user_forecast_data_set = NULL,
    forecast_advanced = F,
    
    TDM_recommendation = NULL,
    
    pk_plots = NULL,
    adapted_pk_plot = NULL,
    adapted_pk_plot_withTDM = NULL,
    dist_plots = NULL,
    tdm_samples_available =T,
    last_known_dose = NULL,
    last_known_dose_orig = NULL,
    standard_y_zoom = NULL,
    standard_x_zoom = NULL,
    
    user_y_zoom = NULL,
    user_x_zoom = NULL,
    
    standard_y_zoom_pop = NULL,
    standard_x_zoom_pop = NULL,
    
    user_y_zoom_pop = NULL,
    user_x_zoom_pop = NULL,
    
    standard_y_zoom_adapt = NULL,
    standard_x_zoom_adapt = NULL,
    
    user_y_zoom_adapt = NULL,
    user_x_zoom_adapt = NULL,
    
    covariate_plot = NULL
    
  )
  
  reset_adapt_to_last_known_dose <- function(){

    updateNumericInput(session, inputId = "adapt.ii", value = round(as.numeric(as.character(app_data$last_known_dose_orig$II)),2) )
    updateNumericInput(session, inputId = "adapt.dose", value = round(as.numeric(as.character(app_data$last_known_dose_orig$AMT)),2) )
    updateNumericInput(session, inputId = "adapt.dur", value = round(as.numeric(as.character(app_data$last_known_dose_orig$DUR)),2) )
    
    ### Change recommendation to stay on current therapy
    updateSelectizeInput(session, inputId = "choose_recommendation", selected = 1 )
    
  }
  
  convert_xlsx_to_NMTRAN <- function(data){
    date <- (as.character(data$DATE))
    
    hour <- (as.character(data$TIME))
    

    hour_time <- as.POSIXct(paste(date, hour, sep= " "))
    
    
    
    orig_data <- data.frame(DATE = date,
                            TIME = hour,
                            AMT = data$AMT,
                            CONC = data$CONC,
                            EVID = data$EVID,
                            DUR = data$DUR,
                            WT = data$WT,
                            CRCL = data$CRCL,
                            DIAL = data$DIAL) 


    
    t_ref <- hour_time[1]
    
    hour_time <- as.numeric(hour_time)
    

    
    
    hour_time <- (hour_time-min(hour_time))/3600

    conv_data <- data.frame(time=hour_time,
                            amt = as.numeric(data$AMT),
                            conc = as.numeric(data$CONC),
                            evid = data$EVID,
                            dur = as.numeric(data$DUR)/60,
                            wt = data$WT,
                            crcl = data$CRCL,
                            dial = data$DIAL)
    

    
    return(list(original_data=orig_data,
                conv_data=conv_data,
                time_reference=t_ref))
  }
  
  ## This function takes care of all the simulations
  ## MCMC is outsourced to global.R
  updatePKPlot <- function(){
    
    times <- app_data$data_set$time
    
    if(length(times) == 1){
      times <- c(times, times+12)
    }
    
    x_min <- as.POSIXct.numeric(min(times)*3600,origin=app_data$time_reference)
    x_max <- as.POSIXct.numeric((max(times)+input$simulate.t)*3600,origin=app_data$time_reference)
    
    ## get tdm data from the table
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    

    
    # --- Update parameters from GUI
    
    app_data$params <- c(input$WT, input$CRCL, input$has_dialysis)
    
    
    ## --- get time variant covariates
    

    temp_time <- seq(min(times), max(times)+input$simulate.t, by=input$delta.t)
    
    temp_time <- c(temp_time, tdm_data$time)
    
    temp_time <- unique(temp_time)
    
    temp_time <- sort(temp_time)
    
    app_data$times_to_calculate <- temp_time
    
    cov <- get_cov_from_dataset(app_data$data_set, temp_time)
    
    app_data$covariates <- cov
    
    app_data$mc_result <- perform_mc_simulation(input$mc.iter, ## number of simulations
                                                OMEGAS, ## omegas
                                                THETAS, ## thetas
                                                app_data$data_set, ## App Data for Dosing / TDM Data
                                                covariates=cov,
                                                times=temp_time) ## Time to simulate
    
    
    
    plot_dat <- app_data$mc_result [[1]] ## get Plot data ...
    dat_mc <- app_data$mc_result [[2]]  ### .. and raw results from the mc simulation to complete the plots
    
    ## Get lowest value (non-zero <- log-scale) and max value in PK plot to adjust y-axis
    pop_y_max <- max(plot_dat$CP_max)
    pop_y_min <- min(plot_dat$CP_min[plot_dat$CP_min >0])
    
    tdm_data$time <- as.POSIXct.numeric(tdm_data$time*3600,origin=app_data$time_reference)
    
    plot_dat$TIME <- as.POSIXct.numeric(plot_dat$TIME*3600, origin=app_data$time_reference)
    
    if(app_data$tdm_samples_available) {
    
        app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = input$mcmc.iter, n.burn = input$mcmc.burn,
                                                thetas = THETAS,
                                                omegas = OMEGAS,
                                                covariates=cov,
                                                TIME =temp_time, 
                                                SIGMAS=3.4, time_reference=app_data$time_reference) 
        ind_y_max <- app_data$mcmc_result[[7]]
        ind_y_min <- app_data$mcmc_result[[8]]
        
        ind_y_max <- ifelse(max(tdm_data$conc) > ind_y_max, max(tdm_data$conc), ind_y_max)
        
        
        
        app_data$standard_y_zoom <- c(0,ind_y_max*1.1)
        app_data$standard_x_zoom <- c(x_min, x_max)
        
        app_data$user_y_zoom <- c(0,ind_y_max*1.1)
        app_data$user_x_zoom <- c(x_min, x_max)
        
        ## prepare individual boxplot ### Currently not used
        ind_boxplot <- ggplot(data=data.frame(conc=app_data$mcmc_result[[6]], time="")) + geom_boxplot(aes(x=time, y=conc)) + plot_theme  +
          theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
                axis.ticks.y = element_blank())+ 
          ggtitle("C last [mg/L]", "Individual") + xlab("\n") + coord_cartesian(ylim=c(0,ind_y_max*1.1), expand = F) 
        
        ## prepare individual PK plot
        
        
        
        ind_plot_data <- app_data$mcmc_result[[5]]
        
        ind_plot_data$time <- as.POSIXct.numeric(ind_plot_data$time*3600,origin=app_data$time_reference)
        
        ind_plot_data$low_target <- input$low.target
        
        ind_plot_data$high_target <-input$high.target
        
        
        ### These data can be used to present the simulated covariates
        temp_cov <- app_data$covariates
        temp_cov$time <- (as.POSIXct.numeric(app_data$covariates$time*3600,origin=app_data$time_reference))
       
        
        ## Build raw individual PK plot
        ind_plot <- ggplot(ind_plot_data)  +
          geom_ribbon(aes(ymin=low_target, ymax=high_target, x=time, fill="target"), alpha=0.3) + 
          geom_ribbon(aes(ymin=s1, ymax=s2, x=time,fill="s1"), alpha=0.15, show.legend = T) + 
          geom_ribbon(aes(ymin=s3, ymax=s4, x=time,fill="s2"),  alpha=0.15) + 
          geom_ribbon(aes(ymin=s5, ymax=s6, x=time,fill="s3"), alpha=0.15) + 
          geom_ribbon(aes(ymin=s7, ymax=s8, x=time,fill="s4"),  alpha=0.15) + 
          geom_line(aes(y=max, x=time, colour="ind"), show.legend = T)  + 
          geom_point(data=tdm_data, aes(x=time, y=conc, colour="tdm"), size=3, shape=1, stroke=2) + plot_theme + xlab("") + ylab("Vancomycin Plasma Concentration [mg/L]") +
          ggtitle("Individual Prediction Using TDM Data and Covariates", "Select area and double click for zoom and unzoom") + 
          coord_cartesian(ylim=app_data$standard_y_zoom, xlim=app_data$standard_x_zoom, expand = F) +
          geom_line(data=plot_dat, aes(x=TIME, y=CP, colour="pop"), linetype=2) +  
          geom_hline(aes(yintercept=input$MIC, colour="mic"), linetype=3, size=1) +
          scale_x_datetime(labels = date_format("%a %d.%m.%Y\n%H:%M", tz = "CET")) + 
          scale_colour_manual(values=c("ind"="blue",
                                        "tdm"="firebrick",
                                       "pop"="black",
                                       "mic"="black"),
                                        guide = guide_legend(override.aes = list(
                                          linetype =  c(1,2,3,0),
                                          shape = c(NA, NA,NA,1),
                                          fill = c("white", "white","white", "white")
                                          ),
                                          title=""),
                                        labels=c("Individual Prediction", 
                                                 "Population Prediction", "MIC","TDM Data")
                              ) + 
          scale_fill_manual(values=c("s1"="blue","s2"="blue","s3"="blue","s4"="blue","target"="darkgreen"
                                     ),
                                     guide = guide_legend(override.aes = list(
                                       alpha=c(0.2,0.175,0.15,0.1,0.1),
                                        linetype =  c(0,0,0,0,0),
                                        shape = c(NA, NA,NA,NA,NA)),
                                        title=""),
                                        labels=c("80 % Interval", 
                                                 "85 % Interval",
                                                 "90 % Interval",
                                                 "95 % Interval","Target for Cmin"
                                                 )) + 
          theme(legend.position = "left") 
        
        
        ind_pars <- (app_data$mcmc_result[[10]])
    
    }
    
    
    
    
    
    ## Check if TDM concentration is above upper prediction interval to readjust the plot limits on y-axis
    pop_y_max <- ifelse(max(tdm_data$conc) > pop_y_max, max(tdm_data$conc), pop_y_max)
    
    app_data$standard_y_zoom_pop <- c(0,pop_y_max*1.1)
    app_data$standard_x_zoom_pop <- c(x_min, x_max)
       
    ## prepare population boxplot
    pop_boxplot <- ggplot(data=data.frame(conc=dat_mc[,ncol(dat_mc)], time="")) + geom_boxplot(aes(x=time, y=conc)) + plot_theme +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank())+ ggtitle("C last [mg/L]" , "Population") + 
      xlab("\n") + coord_cartesian(ylim=c(0,pop_y_max*1.1), expand = F)
    
    ## Prepare population PK plot
    pop_plot <- ggplot(data=plot_dat)  + geom_line(aes(x=TIME, y=CP, colour="pop")) +
      geom_ribbon(aes(x=TIME, ymax=CP_max, ymin=CP_min, fill="s1"), alpha=0.15) +
      plot_theme + xlab("") + ylab("Vancomycin Plasma Concentration [mg/L]") + ggtitle("Population Prediction Using Patient Covariates", "Select area and double click for zoom and unzoom") +
      geom_hline(aes(yintercept=input$MIC, colour="mic"), linetype=3, size=1) +
      geom_ribbon(aes(ymin=input$low.target, ymax=input$high.target, x=TIME,fill="target"), alpha=0.3) +  scale_x_datetime( labels = date_format("%a %d.%m.%Y\n%H:%M", tz = "CET")) + 
      coord_cartesian(ylim=app_data$standard_y_zoom_pop, xlim=app_data$standard_x_zoom_pop, expand = F) +
      scale_colour_manual(values=c("mic"="black", "pop"="blue"),
                          guide = guide_legend(override.aes = list(
                            linetype =  c(3,1),
                            shape = c(NA,NA),
                            fill = c("white","white")
                          ),
                          title=""),
                          labels=c("MIC", "Population Prediction")
      ) + 
      scale_fill_manual(values=c("s1"="blue", "target"="darkgreen"
      ),
      guide = guide_legend(override.aes = list(
        alpha=c(0.1,0.1),
        linetype =  c(NULL,NULL),
        shape = c(NA,NA)),
        title=""),
        labels=c("95 % Interval","Target for Cmin")
      ) + 
      theme(legend.position = "left")
      
    
    
    ### Prepare Distribution plots
    
    
    pop_pars <- (app_data$mc_result[[4]])
    
    dist_plots <- list()
    
    if(app_data$tdm_samples_available) {
    
        dist_plots[[1]] <- ggplot() +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=CL, y=..density.., fill="pop"), colour="blue", size=.5,alpha=0.25, linetype=1) +
          geom_density(data=ind_pars, aes(x=CL, y=..density.., fill="ind"), colour="red", size=.5,alpha=0.25, linetype=1) +
          xlab("Clearance (Cl) [L/h]") + ylab("Frequency") +
          scale_fill_manual(values=c("ind"="red", "pop"="blue"
          ),
          guide = guide_legend(override.aes = list(
            linetype =  c(NULL,NULL),
            colour = c("red", "blue"),
            alpha = c(0.25, 0.25),
            shape = c(NA,NA)),
            title=""),
          labels=c("Individual","Population")
          ) + 
          theme(legend.position = "bottom")
        
        dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=V1, y=..density.., fill="pop"), colour="blue", size=.5,alpha=0.25, linetype=1) +
          geom_density(data=ind_pars,aes(x=V1, y=..density.., fill="ind"), colour="red", size=.5,alpha=0.25, linetype=1)+
          xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")+
          scale_fill_manual(values=c("ind"="red", "pop"="blue"
          ),
          guide = guide_legend(override.aes = list(
            linetype =  c(NULL,NULL),
            colour = c("red", "blue"),
            alpha = c(0.25, 0.25),
            shape = c(NA,NA)),
            title=""),
          labels=c("Individual","Population")
          ) + 
          theme(legend.position = "bottom")
        
        dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=V2, y=..density.., fill="pop"), colour="blue", size=.5,alpha=0.25, linetype=1)+
          geom_density(data=ind_pars,aes(x=V2, y=..density.., fill="ind"), colour="red", size=.5,alpha=0.25, linetype=1)+
          xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")+
          scale_fill_manual(values=c("ind"="red", "pop"="blue"
          ),
          guide = guide_legend(override.aes = list(
            linetype =  c(NULL,NULL),
            colour = c("red", "blue"),
            alpha = c(0.25, 0.25),
            shape = c(NA,NA)),
            title=""),
          labels=c("Individual","Population")
          ) + 
          theme(legend.position = "bottom")
        
    } else {
      dist_plots[[1]] <- ggplot() +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=CL, y=..density.., fill="pop"), colour="blue", size=.5,alpha=0.25, linetype=1) +
        xlab("Clearance (Cl) [L/h]") + ylab("Frequency")+
        scale_fill_manual(values=c("pop"="blue"
        ),
        guide = guide_legend(override.aes = list(
          linetype =  c(NULL),
          colour = c("blue"),
          alpha = c(0.25),
          shape = c(NA)),
          title=""),
        labels=c("Population")
        ) + 
        theme(legend.position = "bottom")
      
      dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=V1, y=..density.., fill="pop"), colour="blue", size=.5,alpha=0.25, linetype=1) +
        xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")+
        scale_fill_manual(values=c("pop"="blue"
        ),
        guide = guide_legend(override.aes = list(
          linetype =  c(NULL),
          colour = c("blue"),
          alpha = c(0.25),
          shape = c(NA)),
          title=""),
        labels=c("Population")
        ) + 
        theme(legend.position = "bottom")
      
      dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=V2, y=..density.., fill="pop"), colour="blue", size=.5,alpha=0.25, linetype=1)+
        xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")+
        scale_fill_manual(values=c("pop"="blue"
        ),
        guide = guide_legend(override.aes = list(
          linetype =  c(NULL),
          colour = c("blue"),
          alpha = c(0.25),
          shape = c(NA)),
          title=""),
        labels=c("Population")
        ) + 
        theme(legend.position = "bottom")
    }
    
    app_data$dist_plots <- dist_plots
    
    if(app_data$tdm_samples_available){
    
      return(plots <- list(ind_plot, pop_plot, ind_boxplot, pop_boxplot))
    } else {
      return(plots <- list(NULL, pop_plot, NULL, pop_boxplot))
    }

  }
  
  output$adapted_pkPlot <- renderPlot({
    
    if(is.null(app_data$adapted_pk_plot)){
      return()
    }
    

    app_data$adapted_pk_plot
    
  })
  
  output$adapted_pkPlot_withTDM <- renderPlot({
    
    if(is.null(app_data$adapted_pk_plot_withTDM)){
      return()
    }

    app_data$adapted_pk_plot_withTDM
    
  })
  
  output$pkPlot <- renderPlot({
    
    if(is.null(app_data$pk_plots)){
      return()
    } 
    
    if(app_data$tdm_samples_available) {
      app_data$pk_plots[[1]]
    } else {
      app_data$pk_plots[[2]]
    }
    
    
  })
  
  
  
  output$data_set <- DT::renderDataTable({
    
    if(!app_data$disc_shown){
      app_data$disc_shown = T
      showModal(modalDialog(
        title = "Disclaimer",
        HTML(paste("<B>daGama - VancoSim has been created for personal use only! </B><BR>", 
                   "There is no guarantee for the correctness of results obtained via VancoSim Using data generated by the application is at the risk of the daGama - VancoSim user alone!", 
                   "Any data entered is only saved temporarly and deleted upon closing the application. ", 
                   "<BR>",
                   "The user accepts <a href=\"http://www.osc-lab.de/impressum.html\"> our legal notices</a> automatically when using the application. ")),
        easyClose = TRUE,
        footer = NULL
      ))
    }
    
    if(is.null(app_data$data_set)){
    
      app_data$user_data_set <- read_excel("./blank.xlsx")

    }
    
    display_data <-app_data$user_data_set
    
    display_data
    
    
  })
  
  output$DT_FORE <- DT::renderDataTable({
    

    
    if(is.null(app_data$user_forecast_data_set)){
      
      app_data$user_forecast_data_set <- read_excel("./blank_forecast.xlsx")
      
    }
    
    display_data <-app_data$user_forecast_data_set
    
    display_data
    
    
  })
  
  output$cov_plot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    if(is.null(app_data$mcmc_result)){
      return()
    }
    
    
    ## Show correlation matrix - the lazy way
    if(app_data$tdm_samples_available){
      chart.Correlation(app_data$mcmc_result[[chain]]$chain_data, histogram=TRUE)
    }
    
  })
  
  output$pop_cov_plot <- renderPlot({
    ## Show correlation matrix - the lazy way
    
    if(is.null(app_data$mc_result)){
      return()
    }
    
    chart.Correlation(app_data$mc_result[[3]], histogram=TRUE)
    
  })
  
  output$par_dist <- renderPlot({
    

    if(is.null(app_data$mc_result)){
      return()
    }
    
    
    
    grid.arrange(app_data$dist_plots[[1]],app_data$dist_plots[[2]],app_data$dist_plots[[3]], nrow=1, ncol=3)
  
  })
  
  output$traceplot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    
    if(is.null(app_data$mcmc_result)){
      return()
    }
    
    ## Has to be generalized 
    
    if(app_data$tdm_samples_available){
    ## 
      gridExtra::grid.arrange(app_data$mcmc_result[[chain]]$p_iter_ETA1, app_data$mcmc_result[[chain]]$p_dens_ETA1,    
                              app_data$mcmc_result[[chain]]$p_iter_ETA2, app_data$mcmc_result[[chain]]$p_dens_ETA2, 
                              app_data$mcmc_result[[chain]]$p_iter_ETA3, app_data$mcmc_result[[chain]]$p_dens_ETA3, nrow=3, ncol=2, widths=c(3,1))  
    }
    
    
    
  })
  
  ## Please comment ASAP
  
  updateAdapted_pk_plot <- function(){
    
    
    if (app_data$tdm_samples_available) {
        current_etas <- app_data$mcmc_result$mcmc_etas
    } else {
      current_etas <- app_data$mc_result[[3]]
    }
    
    
    
    
    adapted_dosing_events <- app_data$data_set
    

    adapted_dosing_events$time <- as.numeric(as.character(adapted_dosing_events$time))
    adapted_dosing_events$amt <- as.numeric(as.character(adapted_dosing_events$amt))
    adapted_dosing_events$dur <- as.numeric(as.character(adapted_dosing_events$dur))
    
    
    t <- max(adapted_dosing_events$time)    
    
    
    time <- NULL
    amt <- NULL
    dur <- NULL
    
    
    
    for(i in 1:input$adapt.n){
      
      
      time <- c(time, t+(input$adapt.ii)*i)
      amt <- c(amt,input$adapt.dose)
      dur <- c(dur, input$adapt.dur/60)
      
      
    }
    
    x_min <- as.POSIXct.numeric(min(time)*3600,origin=app_data$time_reference)
    x_max <- as.POSIXct.numeric(max(time+input$adapt.ii)*3600,origin=app_data$time_reference)
    
    adapted_dosing_events <- adapted_dosing_events[adapted_dosing_events$evid==1,]
    
    
    
    adapted_dosing_events <- data.frame(time = adapted_dosing_events$time,
                                        amt = adapted_dosing_events$amt,
                                        dur = adapted_dosing_events$dur)
    
 
    
    if(app_data$forecast_advanced){
      
      if(nrow(app_data$user_forecast_data_set)==0){
        showModal(modalDialog(
          title = "ERROR",
          HTML("No Entries to simulate not found!"),
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
      new_dosing_events <- app_data$forecast_data_set
      
      x_min <- as.POSIXct.numeric(min(app_data$forecast_data_set$time)*3600,origin=app_data$time_reference)
      x_max <- as.POSIXct.numeric(max(app_data$forecast_data_set$time+12)*3600,origin=app_data$time_reference)
      
    } else {
      new_dosing_events <- data.frame(time = time, amt = amt, dur = dur)
    }
    
    
    
    adapted_dosing_events <- rbind(adapted_dosing_events, new_dosing_events)
    
    
    
    ## Remove Factors
    adapted_dosing_events$time <- as.numeric(as.character(adapted_dosing_events$time))
    adapted_dosing_events$amt <- as.numeric(as.character(adapted_dosing_events$amt))
    adapted_dosing_events$dur <- as.numeric(as.character(adapted_dosing_events$dur))
    
    
    
    times <- adapted_dosing_events$time
    
    TIME <- seq(min(times), max(times)+input$adapt.ii, by=input$delta.t)
    
    
    
    cov <- get_cov_from_dataset(app_data$data_set, TIME)
     
    mcmc_se <- list()
   
    ## get tdm data from the table
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    
    withProgress(message = "Performing Simulation...", max = nrow(current_etas), {
      for (i in 1:nrow(current_etas)) {
        
        
        temp_dat <- pk_2cmt_infusion(theta = THETAS,
                                     CLCR=cov$CLCR, WT=cov$WT, DIAL=cov$DIAL,
                                     eta = c(current_etas$eta1[i], current_etas$eta2[i], current_etas$eta3[i]),
                                     dosing_events = adapted_dosing_events,
                                     times=TIME)
        
        mcmc_se[[i]] <- temp_dat$IPRED
        
        
        incProgress(1)
      }
    })
    
    
    df_temp <- NULL
    
    withProgress(message = "Calculating Prediction Interval...", max = nrow(current_etas), {
      ## bind simulations in a data.frame
      for (k in 1:nrow(current_etas)) {
        df_temp <- rbind(df_temp, mcmc_se[[k]])
        incProgress(1)
      }
    })
    
    ## Generate Quantils
    s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.025, 0.05, 0.075, 0.10, 0.9, 0.925, 0.95, 0.975, 0.5)))
    
    ## Combine individual PK data and quantils in a data.frame
    pk_data <- data.frame(time=TIME,
                          s1=s[1,],s2=s[8,], 
                          s3=s[2,],s4=s[7,],
                          s5=s[3,],s6=s[6,],
                          s7=s[4,],s8=s[5,],
                          max=s[9,]) # median 
    
    
    auc_data_24 <- tail(pk_data, ceiling(24/input$delta.t))
    
    auc_data_24$time <- auc_data_24$time - min(auc_data_24$time)
    
    auc_temp_low <- NULL
    auc_temp_median <- NULL
    auc_temp_high <- NULL
    
    for(i in 2:nrow(auc_data_24)){
      sum_c_low <- auc_data_24$s1[i]+auc_data_24$s1[i-1]
      sum_c_median <- auc_data_24$max[i]+auc_data_24$max[i-1]
      sum_c_high <- auc_data_24$s2[i]+auc_data_24$s2[i-1]
      
      delta_t <- auc_data_24$time[i]-auc_data_24$time[i-1]
      
      auc_temp_low <- c(auc_temp_low, delta_t*sum_c_low/2)
      auc_temp_median <- c(auc_temp_median, delta_t*sum_c_median/2)
      auc_temp_high <- c(auc_temp_high, delta_t*sum_c_high/2)
    }
    
    AUC_24 <- c(median=sum(auc_temp_median), lower=sum(auc_temp_low), higher=sum(auc_temp_high))

    
    middle_of_the_plot <- as.POSIXct.numeric( (as.numeric(x_max)-as.numeric(x_min))/2, origin = x_min)
    y_max <- max(pk_data$s2)
    
    pk_data$time <- as.POSIXct.numeric(pk_data$time*3600,origin=app_data$time_reference)
    
    app_data$standard_y_zoom_adapt = c(0, max(pk_data$s2)*1.1)
    app_data$standard_x_zoom_adapt = c(x_min, x_max)
    
    app_data$user_y_zoom_adapt = c(0, max(pk_data$s2)*1.1)
    app_data$user_x_zoom_adapt = c(x_min, x_max)
    

    ## Build raw individual PK plot
    p <- ggplot(pk_data) + 
      geom_hline(aes(yintercept=input$MIC, colour="mic"), linetype=3, size=1) +
      annotate("text", x=middle_of_the_plot, y=y_max, label=paste("AUC 24h: ", round(AUC_24[1],1) , " (95% PI: ", round(AUC_24[2],1), "-", round(AUC_24[3],1), ") mg x h /L", sep="" ), size=6 ) +
      annotate("text", x=middle_of_the_plot, y=y_max-(y_max*0.1), label=paste("AUC 24h/MIC: ", round(AUC_24[1]/input$MIC,1) , " (95% PI: ", round(AUC_24[2]/input$MIC,1), "-", round(AUC_24[3]/input$MIC,1), ") (mg x h /L)/(mg/L)", sep="" ), size=6 ) +
      geom_ribbon(aes(ymin=s1, ymax=s2, x=time, fill="s1" ), alpha=0.15) + 
      geom_ribbon(aes(ymin=s3, ymax=s4, x=time, fill="s2" ), alpha=0.15) + 
      geom_ribbon(aes(ymin=s5, ymax=s6, x=time, fill="s3" ), alpha=0.15) + 
      geom_ribbon(aes(ymin=s7, ymax=s8, x=time, fill="s4" ), alpha=0.15) + 
      geom_line(aes(y=max, x=time, colour="pred"))+
      plot_theme + xlab("") + ylab("Vancomycin Plasma Concentration [mg/L]") + ggtitle("Prediction of new Dosing Scheme beginning at last event (Dose or TDM)", "Select area and double click for zoom and unzoom") + 
      geom_ribbon(aes(ymin=input$low.target, ymax=input$high.target, x=time,fill="target"), alpha=0.3) +  scale_x_datetime(labels = date_format("%a %d.%m.%Y\n%H:%M", tz = "CET")) +
      coord_cartesian(xlim=app_data$standard_x_zoom_adapt, ylim=app_data$standard_y_zoom_adapt, expand = F)
    
    if (app_data$tdm_samples_available){
         p <- p+ scale_colour_manual(values=c("mic"="black",
                                       "pred"="red"),
                              guide = guide_legend(override.aes = list(
                                linetype =  c(3,1),
                                shape = c(NA,NA),
                                fill = c("white","white")
                              ),
                              title=""),
                              labels=c("MIC","Individual Prediction")
          ) + 
          scale_fill_manual(values=c("s1"="red","s2"="red","s3"="red","s4"="red","target"="darkgreen"
          ),
          guide = guide_legend(override.aes = list(
            alpha=c(0.2,0.175,0.15,0.1,0.1),
            linetype =  c(0,0,0,0,0),
            shape = c(NA, NA,NA,NA,NA)),
            title=""),
          labels=c("80 % Interval", 
                   "85 % Interval",
                   "90 % Interval",
                   "95 % Interval","Target for Cmin"
          )) + 
          theme(legend.position = "left")
    } else {
      p <- p+ scale_colour_manual(values=c("mic"="black",
                                           "pred"="red"),
                                  guide = guide_legend(override.aes = list(
                                    linetype =  c(3,1),
                                    shape = c(NA,NA),
                                    fill = c("white","white")
                                  ),
                                  title=""),
                                  labels=c("MIC","Population Prediction")
      ) + 
        scale_fill_manual(values=c("s1"="red","s2"="red","s3"="red","s4"="red","target"="darkgreen"
        ),
        guide = guide_legend(override.aes = list(
          alpha=c(0.2,0.175,0.15,0.1,0.1),
          linetype =  c(0,0,0,0,0),
          shape = c(NA, NA,NA,NA,NA)),
          title=""),
        labels=c("80 % Interval", 
                 "85 % Interval",
                 "90 % Interval",
                 "95 % Interval","Target for Cmin"
        )) + 
        theme(legend.position = "left")
    }
      
    
    
    app_data$adapted_pk_plot <- p
  }
  
  
  output$modelfile <- renderText({
    
    ## Show the correct JAGS model file
    

      includeText("Goti_et_al.bug")
    
    
  })
  
  output$info <- renderText({
    
    
    paste("<h4>1. Please enter below:</h4><h5><BR>Patient covariates, dosing and if available TDM data. Then, click the <B>\"Analyze Data\"</B> Button.</h5>")

  })
  
  output$info.pk <- renderText({
    
    
    paste("<h4>2. Time-Plasmaconcentration-Curves:</h4><h5><BR>Click the <B>\"Automatic Dose Adaptation\"</B> Button for automatic dosing proposal or choose <B>\"Manual Dose Adaptation\"</B>.</h5>")
    
  })
  
  output$info.adapt <- renderText({
    
    
    paste("<h4>3. Dose Adaptation</h4><h5><BR>Explore and modify dose adaptation. click the <B>\"Refresh Simulation\"</B> Button to inspect the result. Then, click the <B>\"Continue\"</B> Button to produce a PDF report.</h5>")
    
  })
  
  output$info.report <- renderText({
    
    
    paste("<h4>4. Create and Export PDF Report:</h4><h5><BR>Click the <B>\"Download Report\"</B> Button. After downloading, <B>\"Reset\"</B> the App to analyze additional data.</h5>")
    
  })
  
  output$info.mcmc <- renderText({
    
    
    paste("<h4>MCMC Diagnostic Plots</h4><h5><BR>Explore the MCMC sampling and the resulting ETA if TDM data was available. These plots show, <B>how</B> different this patient is compared to the typical patient</h5>")
    
  })
  
  output$info.mc <- renderText({
    
    
    paste("<h4>MC Diagnostic Plots</h4><h5><BR>These plots show the population predictions for this patient based solely on the popPK model</h5>")
    
  })
  
  output$info.dist <- renderText({
    
    
    paste("<h4>Individual Parameter Distribution</h4><h5><BR>These plots show the distributions of population and if available individual pharmacokinetic parameters at baseline.</h5>")
    
  })
  
  output$info.model <- renderText({
    
    
    paste("<h4>Model File</h4><h5><BR>This is a JAGS compatible model file of the popPK Model <b>Goti et al. Ther Drug Monit (2018) 40:212–221</B></h5>")
    
  })
  
  output$info.settings <- renderText({
    
    
    paste("<h4>Application Settings</h4><h5><BR>Expert user <B>only!</B></h5>")
    
  })
  
  output$info.about <- renderText({
    
    
    paste("<h4>About</h4><h5><BR>daGama - VancoSim is free and open source: at <a href=\"https://github.com/OlWaWue/VancoSim\">GitHub</a></h5>")
    
  })
  
  observeEvent(input$ADV_FORECAST, {
    if(input$ADV_FORECAST == 2){
      app_data$forecast_advanced = T
    } else {
      app_data$forecast_advanced = F
    }
  })
  
  observeEvent(input$submit, {
    
    if(is.null(app_data$data_set)){
      showModal(modalDialog(
        title = "ERROR",
        HTML("No data entry present!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return()
    }
               
    app_data$step_one_completed = T
    
    app_data$pk_plots <- updatePKPlot()
               
    reset_adapt_to_last_known_dose()
               
    updateTabsetPanel(session, inputId = "mainpage", selected = "PK Plots")

  })
  
  
  observeEvent(input$but.adapt, {
    ## Use dosage adaptation algorithm
    if (app_data$tdm_samples_available) {
      current_etas <- app_data$mcmc_result$mcmc_etas
    } else {
      current_etas <- app_data$mc_result[[3]]
    }
    
    
    current_etas <- tail(current_etas, (trunc(nrow(current_etas)/10)) )
    
    obj_fun_all <- function(par, data_set, N) {
      
      AMT=par[1]
      II=par[2]
      DUR=par[3]
      
      dos_ev <- data_set[data_set$evid==1,]
      
      dos_ev$time <- as.numeric(as.character(dos_ev$time))
      dos_ev$amt <- as.numeric(as.character(dos_ev$amt))
      dos_ev$dur <- as.numeric(as.character(dos_ev$dur))
      
      t <- max(dos_ev$time)  
      
      time <- NULL
      amt <- NULL
      dur <- NULL
      
      
      for(i in 1:N){
        
        
        time <- c(time, t+(II)*i)
        amt <- c(amt, AMT)
        dur <- c(dur, DUR/60)
        
        
      }
      
      
      
      temp <- data.frame(time=time, amt=amt, dur=dur)
      
      
      dos_ev <- data.frame(time=dos_ev$time,
                           amt=dos_ev$amt,
                           dur=dos_ev$dur)
      
      tested_dosing_events <- rbind(dos_ev, temp)
      
      
      times = seq(min(tested_dosing_events$time), max(tested_dosing_events$time)+II, by=1)
      
      cov <- get_cov_from_dataset(app_data$data_set, times)
      
      cmins_adapt <- NULL
      
      withProgress(message = "Autoadaptation...", max = nrow(current_etas), {
        for (i in 1:nrow(current_etas)) {
          
          
          temp_dat <- pk_2cmt_infusion(theta = THETAS,
                                       CLCR=cov$CLCR, WT=cov$WT, DIAL=cov$DIAL,
                                       eta = c(current_etas$eta1[i], current_etas$eta2[i], current_etas$eta3[i]),
                                       dosing_events = tested_dosing_events,
                                       times=times)
          
          cmins_adapt <- c(cmins_adapt, temp_dat$IPRED[nrow(temp_dat)])
          
          
          incProgress(1)
        }
      })
      

      ipre <- mean(cmins_adapt)
      
      ther_tar <- (input$high.target+input$low.target)/2
      
      return((ther_tar-ipre)^2)
    }
    
    obj_fun_dose <- function(par, II, DUR, data_set, N) {
      
      AMT=par[1]
      
      dos_ev <- data_set[data_set$evid==1,]
      
      dos_ev$time <- as.numeric(as.character(dos_ev$time))
      dos_ev$amt <- as.numeric(as.character(dos_ev$amt))
      dos_ev$dur <- as.numeric(as.character(dos_ev$dur))
      
      t <- max(dos_ev$time)  
      
      time <- NULL
      amt <- NULL
      dur <- NULL
      
      
      for(i in 1:N){
        
        
        time <- c(time, t+(II)*i)
        amt <- c(amt, AMT)
        dur <- c(dur, DUR/60)
        
        
      }
      
      
      
      temp <- data.frame(time=time, amt=amt, dur=dur)
      
      
      
      dos_ev <- data.frame(time=dos_ev$time,
                           amt=dos_ev$amt,
                           dur=dos_ev$dur)
      
      tested_dosing_events <- rbind(dos_ev, temp)
      
      
      times = seq(min(tested_dosing_events$time), max(tested_dosing_events$time)+II, by=1)
      
      cov <- get_cov_from_dataset(app_data$data_set, times)
      
      cmins_adapt <- NULL
      
      withProgress(message = "Autoadaptation...", max = nrow(current_etas), {
        for (i in 1:nrow(current_etas)) {
          
          
          temp_dat <- pk_2cmt_infusion(theta = THETAS,
                                       CLCR=cov$CLCR, WT=cov$WT, DIAL=cov$DIAL,
                                       eta = c(current_etas$eta1[i], current_etas$eta2[i], current_etas$eta3[i]),
                                       dosing_events = tested_dosing_events,
                                       times=times)
          
          cmins_adapt <- c(cmins_adapt, temp_dat$IPRED[nrow(temp_dat)])
          
          
          incProgress(1)
        }
      })
      
      ipre <- mean(cmins_adapt)
      
      ther_tar <- (input$high.target+input$low.target)/2
      
      return((ther_tar-ipre)^2)
    }
    obj_fun_ii <- function(par, AMT, DUR, data_set, N) {
      
      II=par[1]
      
      dos_ev <- data_set[data_set$evid==1,]
      
      dos_ev$time <- as.numeric(as.character(dos_ev$time))
      dos_ev$amt <- as.numeric(as.character(dos_ev$amt))
      dos_ev$dur <- as.numeric(as.character(dos_ev$dur))
      
      t <- max(dos_ev$time)  
      
      time <- NULL
      amt <- NULL
      dur <- NULL
      
      
      for(i in 1:N){
        
        
        time <- c(time, t+(II)*i)
        amt <- c(amt, AMT)
        dur <- c(dur, DUR/60)
        
        
      }
      
      
      
      temp <- data.frame(time=time, amt=amt, dur=dur)
      
      
      dos_ev <- data.frame(time=dos_ev$time,
                           amt=dos_ev$amt,
                           dur=dos_ev$dur)
      
      tested_dosing_events <- rbind(dos_ev, temp)
      
      
      times = seq(min(tested_dosing_events$time), max(tested_dosing_events$time)+II, by=1)
      
      cov <- get_cov_from_dataset(app_data$data_set, times)
      
      cmins_adapt <- NULL
      
      withProgress(message = "Autoadaptation...", max = nrow(current_etas), {
        for (i in 1:nrow(current_etas)) {
          
          
          temp_dat <- pk_2cmt_infusion(theta = THETAS,
                                       CLCR=cov$CLCR, WT=cov$WT, DIAL=cov$DIAL,
                                       eta = c(current_etas$eta1[i], current_etas$eta2[i], current_etas$eta3[i]),
                                       dosing_events = tested_dosing_events,
                                       times=times)
          
          cmins_adapt <- c(cmins_adapt, temp_dat$IPRED[nrow(temp_dat)])
          
          
          incProgress(1)
        }
      })
      
      ipre <- mean(cmins_adapt)
      
      ther_tar <- (input$high.target+input$low.target)/2
      
      return((ther_tar-ipre)^2)
    }
    
    obj_fun_dur <- function(par, AMT, II, data_set, N) {
      
      DUR=par[1]
      
      dos_ev <- data_set[data_set$evid==1,]
      
      dos_ev$time <- as.numeric(as.character(dos_ev$time))
      dos_ev$amt <- as.numeric(as.character(dos_ev$amt))
      dos_ev$dur <- as.numeric(as.character(dos_ev$dur))
      
      t <- max(dos_ev$time)  
      
      time <- NULL
      amt <- NULL
      dur <- NULL
      
      
      for(i in 1:N){
        
        
        time <- c(time, t+(II)*i)
        amt <- c(amt, AMT)
        dur <- c(dur, DUR/60)
        
        
      }
      
      
      
      temp <- data.frame(time=time, amt=amt, dur=dur)
      
      
      dos_ev <- data.frame(time=dos_ev$time,
                           amt=dos_ev$amt,
                           dur=dos_ev$dur)
      
      tested_dosing_events <- rbind(dos_ev, temp)
      
      
      times = seq(min(tested_dosing_events$time), max(tested_dosing_events$time)+II, by=1)
      
      cov <- get_cov_from_dataset(app_data$data_set, times)
      
      cmins_adapt <- NULL
      
      withProgress(message = "Autoadaptation...", max = nrow(current_etas), {
        for (i in 1:nrow(current_etas)) {
          
          
          temp_dat <- pk_2cmt_infusion(theta = THETAS,
                                       CLCR=cov$CLCR, WT=cov$WT, DIAL=cov$DIAL,
                                       eta = c(current_etas$eta1[i], current_etas$eta2[i], current_etas$eta3[i]),
                                       dosing_events = tested_dosing_events,
                                       times=times)
          
          cmins_adapt <- c(cmins_adapt, temp_dat$IPRED[nrow(temp_dat)])
          
          
          incProgress(1)
        }
      })
      
      ipre <- mean(cmins_adapt)
      
      ther_tar <- (input$high.target+input$low.target)/2
      
      return((ther_tar-ipre)^2)
    }
    
    if(input$adapt.what==4) {
    
      opt_res <- optim(par=c(AMT=1000, II=12, DUR=60), fn=obj_fun_all, data_set=app_data$data_set, N=input$adapt.n)
      
      updateNumericInput(session, inputId = "adapt.ii", value = round(as.numeric(as.character(opt_res$par[2])),2) )
      updateNumericInput(session, inputId = "adapt.dose", value = round(as.numeric(as.character(opt_res$par[1])),2) )
      updateNumericInput(session, inputId = "adapt.dur", value = round(as.numeric(as.character(opt_res$par[3])),2) )
    } else if(input$adapt.what==1) {
      opt_res <- optim(par=c(AMT=input$adapt.dose), fn=obj_fun_dose, data_set=app_data$data_set, N=input$adapt.n, II=input$adapt.ii, DUR=input$adapt.dur, method="Brent", lower=0,upper=5000,
                       control=list(maxit=500))
      
      updateNumericInput(session, inputId = "adapt.dose", value = round(as.numeric(as.character(opt_res$par[1])),2) )
    } else if(input$adapt.what==2) {
      opt_res <- optim(par=c(II=input$adapt.ii), fn=obj_fun_ii, data_set=app_data$data_set, N=input$adapt.n, AMT=input$adapt.dose, DUR=input$adapt.dur, method="Brent", lower=6,upper=24,
                       control=list(maxit=500))
      
      updateNumericInput(session, inputId = "adapt.ii", value = round(as.numeric(as.character(opt_res$par[1])),2) )
    } else if(input$adapt.what==3) {
      opt_res <- optim(par=c(DUR=input$adapt.dur), fn=obj_fun_dur, data_set=app_data$data_set,II=input$adapt.ii, N=input$adapt.n, AMT=input$adapt.dose, method="Brent", lower=60,upper=120)
      
      updateNumericInput(session, inputId = "adapt.dur", value = round(as.numeric(as.character(opt_res$par[1])),2) )
    }
    
    #
    # ---- Perform consequences in GUI and app_data
    #
    
    app_data$step_two_completed = T
    
    updateSelectizeInput(session, inputId = "choose_recommendation", selected = 2 )
    
    updateTabsetPanel(session, inputId = "mainpage", selected = "Dose Adaptation")
    delay(1000,
      updateAdapted_pk_plot()
    )
  })
  
  observeEvent(input$but.refresh, {
    ## Use dosage adaptation algorithm
    if(input$adapt.dur < 60){
      showModal(modalDialog(
        title = "CAVE",
        HTML(paste("<B>Infusion duration should exceed 60 minutes! </B><BR>", 
                   "Risk of Red-Man-Syndrome!", 
                   "<BR>",
                   "See: <a href=\"https://www.ncbi.nlm.nih.gov/books/NBK482506/\" target=\"_blank\"> Red Man Syndrome</a>")),
        easyClose = TRUE,
        footer = NULL
      ))
    }
    
    updateAdapted_pk_plot()
    
  })
  
  observeEvent(input$but.reset, {
    ## Use dosage adaptation algorithm
    reset_adapt_to_last_known_dose()
    
    delay(1000,
          updateAdapted_pk_plot()
          )
    
  })
  
  observeEvent(input$but.man_adapt, {
    
    ## Use manual adaptation adaptation algorithm
    
    app_data$step_two_completed = T
    
    reset_adapt_to_last_known_dose()
    
    delay(1000,{
        #  updateAdapted_pk_plot()
          updateTabsetPanel(session, inputId = "mainpage", selected = "Dose Adaptation")
      }
    )
    
    
    
  })
  
  observeEvent(input$pk_doubleclick, {
    
    if(!is.null(input$pk_brush)){
        app_data$user_y_zoom <- c(input$pk_brush$ymin, input$pk_brush$ymax)
        app_data$user_x_zoom <- c(input$pk_brush$xmin, input$pk_brush$xmax)
        
        app_data$user_y_zoom_pop <- c(input$pk_brush$ymin, input$pk_brush$ymax)
        app_data$user_x_zoom_pop <- c(input$pk_brush$xmin, input$pk_brush$xmax)

        app_data$user_x_zoom<- as.POSIXct.numeric(app_data$user_x_zoom,origin=as.POSIXct(strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")))
        app_data$user_x_zoom_pop<- as.POSIXct.numeric(app_data$user_x_zoom_pop,origin=as.POSIXct(strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")))

        session$resetBrush("pk_brush")
        
        if(!is.null(app_data$pk_plots[[1]])){
          app_data$pk_plots[[1]] <- app_data$pk_plots[[1]]+ coord_cartesian(ylim=app_data$user_y_zoom, xlim=app_data$user_x_zoom, expand = F)
        } else {
          app_data$pk_plots[[2]] <- app_data$pk_plots[[2]]+ coord_cartesian(ylim=app_data$user_y_zoom_pop, xlim=app_data$user_x_zoom_pop, expand = F)
        }
    } else {
      if(!is.null(app_data$pk_plots[[1]])){
        app_data$pk_plots[[1]] <- app_data$pk_plots[[1]]+ coord_cartesian(ylim=app_data$standard_y_zoom, xlim=app_data$standard_x_zoom, expand = F)
      } else {
        app_data$pk_plots[[2]] <- app_data$pk_plots[[2]]+ coord_cartesian(ylim=app_data$standard_y_zoom_pop, xlim=app_data$standard_x_zoom_pop, expand = F)
      }
    }
     
  })
  
  observeEvent(input$pk_adapt_doubleclick, {
    
    if(!is.null(input$pk_adapt_brush)){
      app_data$user_y_zoom_adapt <- c(input$pk_adapt_brush$ymin, input$pk_adapt_brush$ymax)
      app_data$user_x_zoom_adapt <- c(input$pk_adapt_brush$xmin, input$pk_adapt_brush$xmax)

      app_data$user_x_zoom_adapt<- as.POSIXct.numeric(app_data$user_x_zoom_adapt,origin=as.POSIXct(strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")))

      
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + coord_cartesian(ylim=app_data$user_y_zoom_adapt, xlim=app_data$user_x_zoom_adapt, expand = F)
      app_data$adapted_pk_plot <- app_data$adapted_pk_plot + coord_cartesian(ylim=app_data$user_y_zoom_adapt, xlim=app_data$user_x_zoom_adapt, expand = F)
      
      session$resetBrush("pk_adapt_brush")
      
    } else {
      app_data$user_y_zoom_adapt <- NULL
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + coord_cartesian(ylim=app_data$standard_y_zoom_adapt, xlim=app_data$standard_x_zoom_adapt, expand = F)
      app_data$adapted_pk_plot <- app_data$adapted_pk_plot + coord_cartesian(ylim=app_data$standard_y_zoom_adapt, xlim=app_data$standard_x_zoom_adapt, expand = F)
    }
    
    if(!is.null(app_data$TDM_recommendation)) {
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot
        if(is.null(app_data$user_y_zoom_adapt)) {
          app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + geom_vline(aes(xintercept=app_data$TDM_recommendation), size=2, colour="blue", alpha=0.25)+
            annotate("text", x=app_data$TDM_recommendation, y=(max(app_data$standard_y_zoom_adapt)/2), label=paste("Take TDM sample here:", "\n", as.character(app_data$TDM_recommendation, format = "%a %d.%m.%Y\n%H:%M", tz = "CET") ), size=6)
        } else {
          app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + geom_vline(aes(xintercept=app_data$TDM_recommendation), size=2, colour="blue", alpha=0.25)+
            annotate("text", x=app_data$TDM_recommendation, y=(max(app_data$user_y_zoom_adapt)/2), label=paste("Take TDM sample here:", "\n", as.character(app_data$TDM_recommendation, format = "%a %d.%m.%Y\n%H:%M", tz = "CET") ), size=6)
        }
    }
    
  })
  
  observeEvent(input$pk_adapt_doubleclick_withTDM, {
    if(!is.null(input$pk_adapt_brush_withTDM)){
      app_data$user_y_zoom_adapt <- c(input$pk_adapt_brush_withTDM$ymin, input$pk_adapt_brush_withTDM$ymax)
      app_data$user_x_zoom_adapt <- c(input$pk_adapt_brush_withTDM$xmin, input$pk_adapt_brush_withTDM$xmax)
      
      app_data$user_x_zoom_adapt<- as.POSIXct.numeric(app_data$user_x_zoom_adapt,origin=as.POSIXct(strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")))
      
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + coord_cartesian(ylim=app_data$user_y_zoom_adapt, xlim=app_data$user_x_zoom_adapt, expand = F)
      app_data$adapted_pk_plot <- app_data$adapted_pk_plot + coord_cartesian(ylim=app_data$user_y_zoom_adapt, xlim=app_data$user_x_zoom_adapt, expand = F)
      
      session$resetBrush("pk_adapt_brush_withTDM")
      
    } else {
      app_data$user_y_zoom_adapt <- NULL
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + coord_cartesian(ylim=app_data$standard_y_zoom_adapt, xlim=app_data$standard_x_zoom_adapt, expand = F)
      app_data$adapted_pk_plot <- app_data$adapted_pk_plot + coord_cartesian(ylim=app_data$standard_y_zoom_adapt, xlim=app_data$standard_x_zoom_adapt, expand = F)
    }
    if(!is.null(app_data$TDM_recommendation)) {
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot
      if(is.null(app_data$user_y_zoom_adapt)) {
        app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + geom_vline(aes(xintercept=app_data$TDM_recommendation), size=2, colour="blue", alpha=0.25)+
          annotate("text", x=app_data$TDM_recommendation, y=(max(app_data$standard_y_zoom_adapt)/2), label=paste("Take TDM sample here:", "\n", as.character(app_data$TDM_recommendation, format = "%a %d.%m.%Y\n%H:%M", tz = "CET") ), size=6)
      } else {
        app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + geom_vline(aes(xintercept=app_data$TDM_recommendation), size=2, colour="blue", alpha=0.25)+
          annotate("text", x=app_data$TDM_recommendation, y=(max(app_data$user_y_zoom_adapt)/2), label=paste("Take TDM sample here:", "\n", as.character(app_data$TDM_recommendation, format = "%a %d.%m.%Y\n%H:%M", tz = "CET") ), size=6)
      }
    }
  })
  
  observeEvent(input$but.report, {
    ## Use dosage adaptation algorithm
    
    app_data$step_three_completed = T
    
    updateAdapted_pk_plot()
    
    app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot
    
    updateTabsetPanel(session, inputId = "mainpage", selected = "Clinical Report")
    
  })

  output$vers <- renderText({
    includeText("version.txt")
  })
  
  
  
  observe({
    if (req(input$mainpage) == "Clinical Report"){
      if(!(app_data$step_one_completed & app_data$step_two_completed & app_data$step_three_completed)){
        showModal(modalDialog(
          title = "Error",
          HTML(paste("<B>Not all three stages have been completed! </B><BR>", 
                     "Please complete steps 1 to 3 prior to report creation.")), 
          
          easyClose = TRUE,
          footer = NULL
        ))
        updateTabsetPanel(session, inputId = "mainpage", selected = "Enter Patient data")
      }
    }

  })
  
  
  observe({
    if (req(input$mainpage) == "PK Plots"){
      if(!(app_data$step_one_completed)){
        showModal(modalDialog(
          title = "Error",
          HTML(paste("<B>No data submitted! </B><BR>", 
                     "Please submit data")), 
          
          easyClose = TRUE,
          footer = NULL
        ))
        updateTabsetPanel(session, inputId = "mainpage", selected = "Enter Patient data")
      }
    }
    
  })
  
  observe({
    if (req(input$mainpage) == "Dose Adaptation"){
      if(!(app_data$step_one_completed & app_data$step_two_completed)){
        showModal(modalDialog(
          title = "Error",
          HTML(paste("<B>Not all three stages have been completed! </B><BR>", 
                     "Please submit and analyze data prior to dose adaptation")), 
          
          easyClose = TRUE,
          footer = NULL
        ))
        updateTabsetPanel(session, inputId = "mainpage", selected = "Enter Patient data")
      }
    }
    
  })
  
  output$downPT <- downloadHandler(
    filename = function() {
      paste(input$pat_ID, ".xlsx", sep = "")
    },
    content = function(file) {
      
      
      tryCatch({
        
        new_xlsx_data <- list("DATA_SET"=app_data$user_data_set,
                              "PATIENT"=data.frame('ID'=input$pat_ID,
                                                   'WT'=input$WT,
                                                   'CRCL'=input$CRCL,
                                                   'DIAL'=as.numeric(input$has_dialysis)),
                              "PATHOGEN"=data.frame('PATHOGEN'=input$choose_pathogen,
                                                    'MIC'=input$MIC))
        
        writexl::write_xlsx(new_xlsx_data, file)
        
        
      },
      error = function(e){
        showModal(modalDialog(
          title = "ERROR",
          HTML(paste("Error while saving Patient!<br>Details:<br><br>",e)),
          easyClose = TRUE,
          footer = NULL
        ))
        stop(safeError(e))
      })
    }
  )
  
  output$but.download = downloadHandler(
    
    
    
    ## Filename includes patient ID
    filename = paste("report_", input$pat_ID, ".pdf", sep=""),
    content = function(file) {
      withProgress(message = "Compiling report ...", style="notification", value =0 ,{
        
        
        incProgress(0.33)
        out <- render('report.Rmd', output_format=pdf_document(latex_engine = "xelatex"))
        incProgress(0.33)
        file.rename(out, file)
        incProgress(0.34)
      })
    }
  )
  
  observeEvent(input$BUT_ADD_TDM, {
    new_date = as.character(input$ADD_TDM_DATE)
    new_time = (strsplit((as.character(input$ADD_TDM_TIME)), " ")[[1]][2])
    
    app_data$TDM_recommendation <-as.POSIXct(as.POSIXlt.character(paste(new_date, new_time, "")), format = "%a %d.%m.%Y\n%H:%M", tz = "CET")
    
    
    
    if(is.null(app_data$user_y_zoom_adapt)) {
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + geom_vline(aes(xintercept=app_data$TDM_recommendation), size=2, colour="blue", alpha=0.25)+
        annotate("text", x=app_data$TDM_recommendation, y=(max(app_data$standard_y_zoom_adapt)/2), label=paste("Take TDM sample here:", "\n", as.character(app_data$TDM_recommendation, format = "%a %d.%m.%Y\n%H:%M", tz = "CET") ), size=6)
    } else {
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot_withTDM + geom_vline(aes(xintercept=app_data$TDM_recommendation), size=2, colour="blue", alpha=0.25)+
        annotate("text", x=app_data$TDM_recommendation, y=(max(app_data$user_y_zoom_adapt)/2), label=paste("Take TDM sample here:", "\n", as.character(app_data$TDM_recommendation, format = "%a %d.%m.%Y\n%H:%M", tz = "CET") ), size=6)
    }
    

  })
  
  observeEvent(input$additional_tdm,{
    if(!input$additional_tdm){
      app_data$TDM_recommendation <- NULL
      app_data$adapted_pk_plot_withTDM <- app_data$adapted_pk_plot
    }
  })
  
  observeEvent(input$REM_FORE_OK,{
    
    if(input$REM_FORE>nrow(app_data$user_forecast_data_set)){
      showModal(modalDialog(
        title = "ERROR",
        HTML("Entry not found!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return()
    }
    
    if(nrow(app_data$user_forecast_data_set)>1){  
      
      app_data$user_forecast_data_set <- app_data$user_forecast_data_set[-as.numeric(input$REM_FORE),]
      app_data$user_forecast_data_set <- app_data$user_forecast_data_set[with(app_data$user_forecast_data_set, order(DATE, TIME)),]
      
      write_xlsx(app_data$user_forecast_data_set, "./temp_fore.xlsx")
      app_data$user_forecast_data_set <-(read_xlsx("./temp_fore.xlsx"))

    } else if(nrow(app_data$user_forecast_data_set)<=1){
      app_data$user_forecast_data_set <- read_excel("./blank_forecast.xlsx")
    }
    
    
    
    
    if(nrow(app_data$user_forecast_data_set)>=1){
      fore_data_set <- data.frame()
      for(i in 1:nrow(app_data$user_forecast_data_set)){
        cur_date <- app_data$user_forecast_data_set$DATE[i]
        cur_time <- app_data$user_forecast_data_set$TIME[i]
        cur_date_time <- as.POSIXlt.character(paste(cur_date, cur_time, ""))
        
        cur_time <- (as.numeric(cur_date_time)-as.numeric(app_data$time_reference))/3600
        
       
        
        fore_data_set <- rbind(fore_data_set, c('time'=cur_time,
                                                'amt'=as.numeric(app_data$user_forecast_data_set$AMT[i]),
                                                'dur'=as.numeric(app_data$user_forecast_data_set$DUR[i])/60))
        
      }
      colnames(fore_data_set) <- c('time', 'amt', 'dur')

      app_data$forecast_data_set <- fore_data_set
    }
    
    
    toggleModal(session, "REM_FORE_MODAL", toggle = "close")
  })
  
  observeEvent(input$ADD_FORE_OK,{
    
    ## Get Entry type
    new_date = as.character(input$ADD_FORE_DATE)
    new_time = (strsplit((as.character(input$ADD_FORE_TIME)), " ")[[1]][2])

    new_entry = data.frame('DATE'=new_date,
                           'TIME'=new_time,
                           'AMT'=input$ADD_FORE_AMT,
                           'DUR'=input$ADD_FORE_DUR)

    
    
    if(nrow(app_data$user_forecast_data_set)>0){  
      app_data$user_forecast_data_set$DATE <- as.character(app_data$user_forecast_data_set$DATE )
      app_data$user_forecast_data_set$TIME <- as.character(app_data$user_forecast_data_set$TIME )
      app_data$user_forecast_data_set$DUR <- as.character(app_data$user_forecast_data_set$DUR )

      app_data$user_forecast_data_set <- rbind(app_data$user_forecast_data_set, new_entry)
      
      app_data$user_forecast_data_set <- app_data$user_forecast_data_set[with(app_data$user_forecast_data_set, order(DATE, TIME)),]
      
    } else {
      app_data$user_forecast_data_set <- new_entry
    }
    
    write_xlsx(app_data$user_forecast_data_set, "./temp_fore.xlsx")
    app_data$user_forecast_data_set <-(read_xlsx("./temp_fore.xlsx"))

    if(nrow(app_data$user_forecast_data_set)>=1){
      fore_data_set <- data.frame()
      for(i in 1:nrow(app_data$user_forecast_data_set)){
        cur_date <- app_data$user_forecast_data_set$DATE[i]
        cur_time <- app_data$user_forecast_data_set$TIME[i]
        cur_date_time <- as.POSIXlt.character(paste(cur_date, cur_time, ""))
        
 
        cur_time <- (as.numeric(cur_date_time)-as.numeric(app_data$time_reference))/3600
        
        fore_data_set <- rbind(fore_data_set, c('time'=cur_time,
                                                'amt'=as.numeric(app_data$user_forecast_data_set$AMT[i]),
                                                'dur'=as.numeric(app_data$user_forecast_data_set$DUR[i])/60))
        
      }
      colnames(fore_data_set) <- c('time', 'amt', 'dur')
      
      
      app_data$forecast_data_set <- fore_data_set
    }
    
    toggleModal(session, "ADD_FORE_MODAL", toggle = "close")
  })
  
  observeEvent(input$REM_OK,{
    
    if(input$REM_ENT>nrow(app_data$user_data_set)){
        showModal(modalDialog(
          title = "ERROR",
          HTML("Entry not found!"),
          easyClose = TRUE,
          footer = NULL
        ))
        return()
    }
    
    if(nrow(app_data$user_data_set)>1){  
      
      app_data$user_data_set <- app_data$user_data_set[-as.numeric(input$REM_ENT),]
      
      new_xlsx_data <- list("DATA_SET"=app_data$user_data_set,
                            "PATIENT"=data.frame('ID'=input$pat_ID,
                                                 'WT'=input$WT,
                                                 'CRCL'=input$CRCL,
                                                 'DIAL'=as.numeric(input$has_dialysis)),
                            "PATHOGEN"=data.frame('PATHOGEN'=input$choose_pathogen,
                                                  'MIC'=input$MIC))
      
      writexl::write_xlsx(new_xlsx_data, "./temp.xlsx")
      temp <- convert_xlsx_to_NMTRAN(read_xlsx("./temp.xlsx"))

      
      app_data$data_set <- temp$conv_data
      app_data$user_data_set <- temp$original_data
      app_data$time_reference <- temp$time_reference
      
      
      
      
      if(nrow(temp$conv_data[temp$conv_data$evid==0,])>=1){
        app_data$tdm_samples_available <- TRUE
      } else {
        app_data$tdm_samples_available <- FALSE
      }
      
      ii = 12
      doses = temp$conv_data[temp$conv_data$evid==1,]
      if(nrow(doses)>1){
        ii = doses$time[nrow(doses)] - doses$time[nrow(doses)-1]
      } 
      
      
      app_data$last_known_dose <- tail(temp$conv_data[temp$conv_data$evid==1,],1)
      app_data$last_known_dose_orig <- tail(temp$original_data[temp$original_data$EVID==1,],1)
      
      app_data$last_known_dose$ii <- ii
      app_data$last_known_dose_orig$II <- ii
      
      
    } else if(nrow(app_data$user_data_set)<=1){
      resetApp()
    }
    toggleModal(session, "REM_MODAL", toggle = "close")
  })
  
  observeEvent(input$ADD_OK,{
    
    ## Get Entry type
    new_date = as.character(input$ADD_DATE)
    new_time = (strsplit((as.character(input$ADD_TIME)), " ")[[1]][2])

    
    if(input$SELECT_TYPE==1) { ## New Dose

      new_entry = data.frame('DATE'=new_date,
                    'TIME'=new_time,
                    'AMT'=input$ADD_AMT,
                    'CONC'=".",
                    'EVID'=1,
                    'DUR'=input$ADD_DUR,
                    'WT'=".",
                    'CRCL'=".",
                    'DIAL'=".")
      
    } else if(input$SELECT_TYPE==2) { ## TDM Measurement
      new_entry = data.frame('DATE'=new_date,
                    'TIME'=new_time,
                    'AMT'=".",
                    'CONC'=input$ADD_TDM,
                    'EVID'=0,
                    'DUR'=".",
                    'WT'=".",
                    'CRCL'=".",
                    'DIAL'=".")
    } else if(input$SELECT_TYPE==3) { ## CLCR Entry
      new_entry = data.frame('DATE'=new_date,
                    'TIME'=new_time,
                    'AMT'=".",
                    'CONC'=".",
                    'EVID'=2,
                    'DUR'=".",
                    'WT'=".",
                    'CRCL'=input$ADD_CLCR,
                    'DIAL'=".")
    } else if(input$SELECT_TYPE==4) { ## WT Measurement
      new_entry = data.frame('DATE'=new_date,
                    'TIME'=new_time,
                    'AMT'=".",
                    'CONC'=".",
                    'EVID'=2,
                    'DUR'=".",
                    'WT'=input$ADD_WT,
                    'CRCL'=".",
                    'DIAL'=".")
    } else if(input$SELECT_TYPE==5){ ## DIAL Entry
      new_entry = data.frame('DATE'=new_date,
                    'TIME'=new_time,
                    'AMT'=".",
                    'CONC'=".",
                    'EVID'=2,
                    'DUR'=".",
                    'WT'=".",
                    'CRCL'=".",
                    'DIAL'=as.numeric(input$ADD_DIAL))
    }
     
     
    if(nrow(app_data$user_data_set)>0){  
      app_data$user_data_set$DATE <- as.character(app_data$user_data_set$DATE )
      app_data$user_data_set$TIME <- as.character(app_data$user_data_set$TIME )
      app_data$user_data_set$CONC <- as.character(app_data$user_data_set$CONC )
      app_data$user_data_set$DUR <- as.character(app_data$user_data_set$DUR )
      app_data$user_data_set$WT <- as.character(app_data$user_data_set$WT )
      app_data$user_data_set$CRCL <- as.character(app_data$user_data_set$CRCL )
      app_data$user_data_set$DIAL <- as.character(app_data$user_data_set$DIAL )
      app_data$user_data_set <- rbind(app_data$user_data_set, new_entry)
      
      app_data$user_data_set <- app_data$user_data_set[with(app_data$user_data_set, order(DATE, TIME)),]
      
    } else {
      app_data$user_data_set <- new_entry
    }
    
    
    new_xlsx_data <- list("DATA_SET"=app_data$user_data_set,
                          "PATIENT"=data.frame('ID'=input$pat_ID,
                                               'WT'=input$WT,
                                               'CRCL'=input$CRCL,
                                               'DIAL'=as.numeric(input$has_dialysis)),
                          "PATHOGEN"=data.frame('PATHOGEN'=input$choose_pathogen,
                                                'MIC'=input$MIC))

    writexl::write_xlsx(new_xlsx_data, "./temp.xlsx")
    temp <- convert_xlsx_to_NMTRAN(read_xlsx("./temp.xlsx", sheet="DATA_SET"))
    
    app_data$data_set <- temp$conv_data
    app_data$user_data_set <- temp$original_data
    app_data$time_reference <- temp$time_reference
    
    
    
    
    if(nrow(temp$conv_data[temp$conv_data$evid==0,])>=1){
      app_data$tdm_samples_available <- TRUE
    } else {
      app_data$tdm_samples_available <- FALSE
    }
    
    ii = 12
    doses = temp$conv_data[temp$conv_data$evid==1,]
    if(nrow(doses)>1){
      ii = doses$time[nrow(doses)] - doses$time[nrow(doses)-1]
    }
    
    
    app_data$last_known_dose <- tail(temp$conv_data[temp$conv_data$evid==1,],1)
    app_data$last_known_dose_orig <- tail(temp$original_data[temp$original_data$EVID==1,],1)
    
    app_data$last_known_dose$ii <- ii
    app_data$last_known_dose_orig$II <- ii
    
    toggleModal(session, "ADD_MODAL", toggle = "close")
    })

  
})