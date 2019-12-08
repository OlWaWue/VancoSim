shinyUI(navbarPage("VancoSim - by Oliver Scherf-Clavel (c) 2019 - JMU Wuerzburg",
                   
                   tabPanel("Enter Patient data",
                            htmlOutput("info"),br(),hr(),
                            sidebarLayout(
                              sidebarPanel(
                                
                                h6("Use an Excel File with similar structure (see Table)"),
                                fileInput("loadPT", "Load data", multiple = FALSE, 
                                          accept = c(".xlsx", ".xls"),                                     
                                          width = NULL,buttonLabel = "Browse...", 
                                          placeholder = "No file selected"),
                                
                                wellPanel("Patient Information",
                                          textInput(inputId="pat_ID", label="Patient ID", value="PAT0001"),
                                          numericInput(inputId="WT", label="Body weight [kg]", value =70),
                                          numericInput(inputId="CRCL", label="Creatinine Clearance [mL/min]", value =120),
                                          checkboxInput("has_dialysis", "Dialysis?", value = F)
                                ),
                                actionButton("submit", label = "Apply Changes", icon = NULL, width = NULL)
                              ),
                              
                              
                              
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                
                                
                                wellPanel("Dataset used in the computation",
                                          tableOutput("data_set")
                                )
                              )
                            )
                   ),
                   tabPanel("PK Plots",
                            plotOutput("pkPlot", height = 800)
                            ),
                   tabPanel("MCMC plots",
                            selectInput("select_chain", "MCMC Chain:", choices = c(1,2,3,4)),br(),
                            wellPanel("Trace plot of MCMC",br(),
                                      
                                      plotOutput("traceplot", height = 600)
                            ),
                            wellPanel("Correlation between random effects",br(),
                                      plotOutput("cov_plot", height = 600)
                            )### TODO: More diagnostic plots
                   ),
                   tabPanel("MC plots",
                            
                            wellPanel("Correlation between random effects",br(),
                                      plotOutput("pop_cov_plot", height = 600)
                            )### TODO: More diagnostic plots
                   ),
                   tabPanel("Parameter distributions",
                            
                            plotOutput("par_dist", height = 600)
                            
                   ),
                   tabPanel("Model File",
                            verbatimTextOutput("modelfile")),
                    tabPanel("Settings",
                             numericInput(inputId="mcmc.iter", label="Iterations MCMC", value =200),
                             numericInput(inputId="mc.iter", label="Iterations MC", value =200),
                             numericInput(inputId="mcmc.burn", label="Burn-in Iterations MCMC", value =10),
                             numericInput(inputId="delta.t", label="Delta time [h]", value =0.5),
                             numericInput(inputId="simulate.t", label="Simulate time [h]", value =0)
                             )
            )
)