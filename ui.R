

shinyUI(navbarPage("VancoSim - by Oliver Scherf-Clavel (c) 2019 - JMU Wuerzburg", id="mainpage",
                   
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
                                wellPanel("Pathogen Information",
                                          selectInput("choose_pathogen", "Pathogen:", selected=1, list("MRSA"=1, 
                                                                                                       "Unidentified"=2, 
                                                                                                       "Mixed infection"=3)),
                                          numericInput(inputId="MIC", label="MIC [mg/L]", value =5)
                                ),
                                actionButton("submit", label = "Apply Changes", icon = NULL, width = NULL)
                              ),
                              
                              
                              
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                
                                disabled(actionButton("but.pkplots", label="Next ...")),
                                
                                wellPanel("Dataset used in the computation",
                                          tableOutput("data_set")
                                )
                              )
                            )
                   ),
                   tabPanel("PK Plots", 
                            plotOutput("pkPlot", height = 800)
                            ),
                   tabPanel("Dose Adaptation"
                            
                   ),
                   tabPanel("Clinical Report",
                            selectInput("choose_recommendation", "Recommendation:", selected=1, list("Continue on this dose"=1, 
                                                                                      "Increase the dose"=2, 
                                                                                      "Decrease the dose"=3,
                                                                                      "Increase the dosing interval"=4,
                                                                                      "Decrease the dosing interval"=5)),
                            textAreaInput("report_comment", "Additional comment: "),
                            checkboxInput("additional_tdm", "Recommend additional TDM?", value = F),
                            checkboxInput("add_dur_info", "Add note about Infusion duration?", value = F)
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
                             numericInput(inputId="simulate.t", label="Simulate time [h]", value =0),
                             numericInput(inputId="low.target", label="Target Throughconcentration [mg/L]", value =10),
                             numericInput(inputId="high.target", label="Limit Cmax [mg/L]", value =20)
                             )
            )
)