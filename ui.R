

shinyUI(navbarPage("VancoSim - by Oliver Scherf-Clavel (c) 2019 - JMU Wuerzburg", id="mainpage",theme = shinytheme("united"),
                   
                   tabPanel("Enter Patient data",
                            htmlOutput("info"),hr(),
                            sidebarLayout(
                              sidebarPanel(
                                useShinyjs(), 
                                
                                actionButton("submit", label = "Analyze Data", icon = icon("chart-bar"), width = NULL),
                                
                                br(),br(),
                                wellPanel("Patient covariates at baseline",br(), br(),
                                          textInput(inputId="pat_ID", label="Patient ID", value="PAT0001"),
                                          numericInput(inputId="WT", label="Body weight [kg]", value =70),
                                          numericInput(inputId="CRCL", label="Creatinine Clearance [mL/min]", value =120),
                                          checkboxInput("has_dialysis", "Dialysis?", value = F)
                                ),
                                wellPanel("Pathogen Information", br(), br(),
                                          selectInput("choose_pathogen", "Pathogen:", selected=1, GLOB_PATHOGENS),
                                          numericInput(inputId="MIC", label="MIC [mg/L]", value =5)
                                ),
                                bsModal("ADD_MODAL", title="Add new Data Entry", trigger="BUT_ADD_ENTRY", size = "medium", 
                                        
                                        dateInput(inputId="ADD_DATE", label="Date:", value = today()),
                                        timeInput(inputId = "ADD_TIME", label = "Time [HH:MM]:", value = Sys.time(), seconds = F),
                                        
                                        selectInput("SELECT_TYPE", label="Entry type:", 
                                                    choices = list("Dosing Event" = 1,
                                                                   "TDM Measurement" = 2,
                                                                   "Creatinine Clearance Measurement" = 3,
                                                                   "Weight Measurement" = 4,
                                                                   "Dialysis" = 5)
                                        ),
                                        conditionalPanel(condition="input.SELECT_TYPE==1",
                                                         numericInput("ADD_AMT", "Dose [mg]:", min = 0, max = 5000, value = 1000),
                                                         numericInput("ADD_DUR", "Infusion duration [min]:", min = 0, max = 1440, value = 60)
                                        ),
                                        conditionalPanel(condition="input.SELECT_TYPE==2",
                                                         numericInput("ADD_TDM", "Measured Vancomycin Concentration [mg/L]:", min = 0, max = 200, value = 0)
                                        ),
                                        conditionalPanel(condition="input.SELECT_TYPE==3",
                                                         numericInput("ADD_CLCR", "Creatinine Clearance [mL/min]:", min = 0, max = 180, value = 100)
                                        ),
                                        conditionalPanel(condition="input.SELECT_TYPE==4",
                                                         numericInput("ADD_WT", "Body Weight [kg]:", min = 10, max = 180, value = 70)
                                        ),
                                        conditionalPanel(condition="input.SELECT_TYPE==5",
                                                         selectInput("ADD_DIAL", label="Dialysis:", 
                                                                     choices = list("End" = 0, "Start" = 1))
                                        ),
                                        br(),
                                        actionButton("ADD_OK", "OK")
                                ),
                                bsModal("REM_MODAL", title="Remove Data Entry", trigger="BUT_REM_ENTRY", size = "small", 
                                        
                                        
                                        numericInput("REM_ENT", "Remove Entry No.:", min = 0, max = 5000, value = 1),
                                        
                                        br(),
                                        actionButton("REM_OK", "OK")
                                )
                              ),
                              
                              
                              
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                
                                fileInput("loadPT", "Load data", multiple = FALSE, 
                                          accept = c(".xlsx", ".xls"),                                     
                                          width = NULL,buttonLabel = "Browse...", 
                                          placeholder = "No file selected"),
                                actionButton("BUT_ADD_ENTRY", "Add Data Entry", icon = icon("plus-square")),
                                actionButton("BUT_REM_ENTRY", "Remove Data Entry", icon = icon("minus-square")),
                                downloadButton(outputId="downPT", "Download Data"), br(),br(),
                                wellPanel("Dataset used in the computation",
                                          DT::dataTableOutput("data_set")
                                )
                              )
                            )
                   ),
                   tabPanel("PK Plots", htmlOutput("info.pk"),hr(),
                            sidebarLayout(
                              sidebarPanel(
                                actionButton("but.man_adapt", label = "Manual Adaptation", icon = icon("hand-paper"), width = NULL),br(),br(),
                                actionButton("but.adapt", label = "Automatic Adaptation", icon = icon("calculator"), width = NULL),br(),br(),
                                selectInput("adapt.for", "Adapt for:", selected=1, GLOB_ADAPT_FOR),
                                selectInput("adapt.what", "Adapt ...", selected=1, GLOB_ADAPT_WHAT)
                                
                              ),
                              mainPanel(
                                plotOutput("pkPlot", height = 600,
                                           brush=brushOpts(id="pk_brush"),
                                           dblclick="pk_doubleclick")
                              )
                            )
                    ),
                   tabPanel("Dose Adaptation", htmlOutput("info.adapt"),hr(),
                            sidebarLayout(
                              sidebarPanel(
                                actionButton("but.report", label = "Continue", icon = icon("file-alt"), width = NULL),br(),br(),
                                selectInput(inputId="ADV_FORECAST", label = "Level:", choices = list("Simple"=1,
                                                                                                     "Advanced"=2)),
                                conditionalPanel(condition="input.ADV_FORECAST==1",
                                              numericInput(inputId="adapt.dose", label="New Dose [mg]", value =1000, step=100, min=0),
                                              numericInput(inputId="adapt.ii", label="New Interdose Interval [h]", value = 12, step=1, min = 6),
                                              numericInput(inputId="adapt.dur", label="New Duration of Infusion [min]", value = 60, step=10, min=10),
                                              numericInput(inputId="adapt.n", label="Number of Dosing Events to simulate", value = 5, step=1, min=3)
                                  ),
                                conditionalPanel(condition="input.ADV_FORECAST==2",
                                                 actionButton("BUT_ADD_FORE", "Add Simulation Entry", icon = icon("plus-square")),br(),br(),
                                                 actionButton("BUT_REM_FORE", "Remove Simulation Entry", icon = icon("minus-square")),br(),br(),
                                                 dataTableOutput(outputId="DT_FORE")
                                ),
                                actionButton("but.reset", label = "Reset to last known dose", icon = icon("undo"), width = NULL),br(),br(),
                                actionButton("but.refresh", label = "Refresh Simulation", icon = icon("refresh"), width = NULL),
                                bsModal("ADD_FORE_MODAL", title="Add new Simulation Entry", trigger="BUT_ADD_FORE", size = "small", 
                                        dateInput(inputId="ADD_FORE_DATE", label="Date:", value = today()),
                                        timeInput(inputId = "ADD_FORE_TIME", label = "Time [HH:MM]:", value = Sys.time(), seconds = F),
                                        numericInput("ADD_FORE_AMT", "Dose [mg]:", min = 0, max = 5000, value = 1000),
                                        numericInput("ADD_FORE_DUR", "Infusion duration [min]:", min = 0, max = 1440, value = 60),
                                        
                                        br(),
                                        actionButton("ADD_FORE_OK", "OK")
                                        
                                        
                                        
                                ),
                                bsModal("REM_FORE_MODAL", title="Remove Simulation Entry", trigger="BUT_REM_FORE", size = "small", 
                                        
                                        
                                        numericInput("REM_FORE", "Remove Entry No.:", min = 0, max = 5000, value = 1),
                                        
                                        br(),
                                        actionButton("REM_FORE_OK", "OK")
                                )
                              ),
                              mainPanel(
                                plotOutput("adapted_pkPlot", height = 600,
                                           brush=brushOpts(id="pk_adapt_brush"),
                                           dblclick="pk_adapt_doubleclick")
                              )
                              
                          
                        )
                   ),
                   tabPanel("Clinical Report",htmlOutput("info.report"),hr(),
                            
                            
                            sidebarLayout(
                              sidebarPanel(
                                selectInput("choose_recommendation", "Recommendation:", selected=1, GLOB_RECOMMENDATIONS),
                                textAreaInput("report_comment", "Additional comment: ", width="100%"),
                                checkboxInput("additional_tdm", "Recommend additional TDM?", value = F),
                                conditionalPanel(condition="input.additional_tdm",
                                                 dateInput(inputId="ADD_TDM_DATE", label="Date:", value = today()),
                                                 timeInput(inputId = "ADD_TDM_TIME", label = "Time [HH:MM]:", value = Sys.time(), seconds = F),
                                                 actionButton("BUT_ADD_TDM", "Add TDM Recommendation", icon = icon("plus-square"))
                                ),
                                checkboxInput("add_dur_info", "Add note about Infusion duration?", value = F),
                                downloadButton("but.download", "Download Report"), br(), br(),
                                actionButton("but.resetApp", label = "Reset Application", icon = icon("undo"), width = NULL)
                              ),
                              mainPanel(
                                plotOutput("adapted_pkPlot_withTDM", height = 600,
                                           brush=brushOpts(id="pk_adapt_brush_withTDM"),
                                           dblclick="pk_adapt_doubleclick_withTDM")
                              )
                            )  
                   ),
                   tabPanel("MCMC plots", htmlOutput("info.mcmc"),hr(),
                            selectInput("select_chain", "MCMC Chain:", choices = c(1,2,3,4)),br(),
                            wellPanel("Trace plot of MCMC",br(),
                                      
                                      plotOutput("traceplot", height = 600)
                            ),
                            wellPanel("Correlation between random effects",br(),
                                      plotOutput("cov_plot", height = 600)
                            )### TODO: More diagnostic plots
                   ),
                   tabPanel("MC plots",htmlOutput("info.mc"),hr(),
                            
                            wellPanel("Correlation between random effects",br(),
                                      plotOutput("pop_cov_plot", height = 600)
                            )### TODO: More diagnostic plots
                   ),
                   tabPanel("Parameter distributions",htmlOutput("info.dist"),hr(),
                            
                            plotOutput("par_dist", height = 600)
                            
                   ),
                   tabPanel("Model File", htmlOutput("info.model"),hr(),
                            verbatimTextOutput("modelfile")),
                   tabPanel("Settings", htmlOutput("info.settings"),hr(),
                             numericInput(inputId="mcmc.iter",   label="Iterations MCMC", value =1000),
                             numericInput(inputId="mc.iter",     label="Iterations MC", value =1000),
                             numericInput(inputId="mcmc.burn",   label="Burn-in Iterations MCMC", value =200),
                             numericInput(inputId="delta.t",     label="Delta time [h]", value =0.25),
                             numericInput(inputId="simulate.t",  label="Simulate time [h]", value =12),
                             numericInput(inputId="low.target",  label="Lower Limit Cthrough [mg/L]", value =10),
                             numericInput(inputId="high.target", label="Upper Limit Cthrough [mg/L]", value =20)
                    ),
                    tabPanel("About", htmlOutput("info.about"),hr(),
                                                     img(src="daGama_logo.png", height = 90, width = 300),
                                                     withTags({
                                                       div(class="header", checked=NA, 
                                                           h4("Information"),
                                                           p("Version: 0.0.1 Alpha"),
                                                           p("This is an early version for demonstration only!"),
                                                           p("Entered data is NOT saved!"),
                                                           p("All rights reserved: 09. November 2019"),hr(),
                                                           h4("Resources"),
                                                           a(href="http://go.uniwue.de/osc-group", "Official Homepage of OSC-Group"),br(),
                                                           a(href="https://doi.org/10.1097/FTD.0000000000000490", "Goti et al. Ther Drug Monit (2018) 40:212–221"),br(),
                                                           a(href="http://www.osc-lab.de/example_dat.xlsx", "Get Example Excel File with multiple dosing with TDM samples"),br(),
                                                           a(href="http://www.osc-lab.de/example_dat_2.xlsx", "Get Example Excel File with multiple dosing WITHOUT TDM samples"),br(),
                                                           a(href="http://www.osc-lab.de/example_dat_3.xlsx", "Get Example Excel File with single dosing with TDM samples"),br(),
                                                           a(href="http://www.osc-lab.de/example_dat_4.xlsx", "Get Example Excel File with single dosing WITHOUT TDM samples"),br(),
                                                           a(href="http://www.osc-lab.de/example_dat_5.xlsx", "Get Example Excel File with time variable covariate"),br(),
                                                           hr(),
                                                           h4("Version history"),
                                                           verbatimTextOutput("vers"))
                                                     })
                          
                    )
            )
)