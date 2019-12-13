

shinyUI(navbarPage("VancoSim - by Oliver Scherf-Clavel (c) 2019 - JMU Wuerzburg", id="mainpage",theme = shinytheme("united"),
                   
                   tabPanel("Enter Patient data",
                            htmlOutput("info"),hr(),
                            sidebarLayout(
                              sidebarPanel(
                                useShinyjs(), 
                                
                                actionButton("submit", label = "Analyze Data", icon = icon("chart-bar"), width = NULL),
                                
                                br(),br(),
                                wellPanel("Patient Information",br(), br(),
                                          textInput(inputId="pat_ID", label="Patient ID", value="PAT0001"),
                                          numericInput(inputId="WT", label="Body weight [kg]", value =70),
                                          numericInput(inputId="CRCL", label="Creatinine Clearance [mL/min]", value =120),
                                          checkboxInput("has_dialysis", "Dialysis?", value = F)
                                ),
                                wellPanel("Pathogen Information", br(), br(),
                                          selectInput("choose_pathogen", "Pathogen:", selected=1, GLOB_PATHOGENS),
                                          numericInput(inputId="MIC", label="MIC [mg/L]", value =5)
                                )
                              ),
                              
                              
                              
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                
                                fileInput("loadPT", "Load data", multiple = FALSE, 
                                          accept = c(".xlsx", ".xls"),                                     
                                          width = NULL,buttonLabel = "Browse...", 
                                          placeholder = "No file selected"),
                                h4("Use an Excel File with similar structure (see Table)"),br(),br(),
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
                                numericInput(inputId="adapt.dose", label="New Dose [mg]", value =1000, step=100, min=0),
                                numericInput(inputId="adapt.ii", label="New Interdose Interval [h]", value = 12, step=1, min = 6),
                                numericInput(inputId="adapt.dur", label="New Duration of Infusion [min]", value = 60, step=10, min=10),
                                numericInput(inputId="adapt.n", label="Number of Dosing Events to simulate", value = 5, step=1, min=3),
                                actionButton("but.reset", label = "Reset to last known dose", icon = icon("undo"), width = NULL),br(),br(),
                                actionButton("but.refresh", label = "Refresh Simulation", icon = icon("refresh"), width = NULL)
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
                                
                                checkboxInput("additional_tdm", "Recommend additional TDM?", value = F),
                                checkboxInput("add_dur_info", "Add note about Infusion duration?", value = F),
                                downloadButton("but.download", "Download Report"), br(), br(),
                                actionButton("but.resetApp", label = "Reset Application", icon = icon("undo"), width = NULL)
                              ),
                              mainPanel(
                                textAreaInput("report_comment", "Additional comment: ", width="100%")
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
                             numericInput(inputId="simulate.t",  label="Simulate time [h]", value =0),
                             numericInput(inputId="low.target",  label="Target Throughconcentration [mg/L]", value =15),
                             numericInput(inputId="high.target", label="Limit Cmax [mg/L]", value =20)
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
                                                           hr(),
                                                           h4("Version history"),
                                                           verbatimTextOutput("vers"))
                                                     })
                          
                    )
            )
)