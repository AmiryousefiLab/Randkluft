build_app_ui <- function(html_code) shinyUI(fluidPage(
  shinyjs::useShinyjs(),
  tags$header(HTML(html_code)),


  tags$style(HTML("
  .pagination-button {
    float: right; /* This will float the button to the right */
    margin-right: 10px; /* You can adjust the margin to control spacing */
  }
")),

tags$style(HTML("
  .pagination-button-back {
    float: left; /* This will float the button to the right */
    margin-left: 10px; /* You can adjust the margin to control spacing */
  }
")),

tags$head(tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css")),


  tags$div(
    id = "title",
    tags$h1(id = "tool_name", "Randk\\uft", style = "color:navy;"),
    tags$h4(
      id = "tool_exp",
      tags$em("Unitary gating of the CyCIF markers"),
      style = "color:black;"
    )
  ),

  tags$div(
    id = "tool",
    tabsetPanel(
      tabPanel(
        strong("Home"),
        br(),
        tags$div(includeMarkdown("./documents/home.md"), style = "max-width:800px;")
      ),

      tabPanel(
        strong("Randkluft"),
        br(),
        sidebarLayout(
          # Sidebar with a slider input
          sidebarPanel(width = 4, br(),
                       tabsetPanel(

                         tabPanel(
                          value=1,
                           strong("Upload file"),
                           br(),

                          #  p("Please provide your .csv cell-marker table."),

                          p("Please upload a .csv file with cell-marker data following the format described in",
                                                   strong("Help,"), "
                                        or download the example data set found in",
                                                   downloadLink('download_example_data', strong(' here.') )
                                                 ),


                          #  fluidRow(column(
                          #    width = 8,  fileInput("cell_file", label = "Cell-Marker File")
                          #  ),

                           fluidRow(
                              column(width = 8,
                                    div(id = "file_input_div",
                                        fileInput("cell_file", label = "Upload a CSV file", accept = ".csv", buttonLabel = "Browse"))),
                                        tags$div(id = "message_loading_data", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),

                              column(width = 4, actionButton("remove_file", "Remove", icon = icon("trash")))
                            ),


                         ),


                         tabPanel(strong("Randkluft"),
                          value=3,
                          tabsetPanel(id='sidebar2',
                            tabPanel(
                          value=1,
                           strong("Essential"),
                           br(),

                          #  column(
                          #    width = 4, actionButton("remove_file", "Remove", icon = icon("trash"))
                          #  )),
                           actionButton(
                             inputId = "run_gate",
                             "Randkluft",
                           ),

                          br(),
                          br(),

                           checkboxGroupInput("selected_columns", "Select Markers", choices = NULL, selected = NULL),
                           br(),


                           materialSwitch(inputId = "gen_hist_plots_on_off",
                                                                       label = "Show gates",
                                                                       status = "danger",
                                                                       right=TRUE,
                                                                       value = TRUE ),


                          # Add the progress bar div element here
                           tags$div(id = "message", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),
                           tags$div(id = "message_gen_hist", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),


                          p('Click button to download a .csv  file with the gate estimations.'),
                          downloadButton(outputId = "downloadEstimations", label = "Download Gate Estimates"),

                          p('Click button to download current set of plots that are displayed.'),
                          downloadButton(outputId = "downloadcurrent", label = "Download Current Plot"),

                          p('Click button to download pdf file of plots for all selected markers.'),
                          downloadButton(outputId = "downloadall", label = "Download All Plots")


                             ),


                              tabPanel(
                          value=2,

                           strong("Bivariate"),
                           br(),

                             varSelectInput("xvar", "X variable", NULL, selected = NULL),
                             varSelectInput("yvar", "Y variable", NULL, selected = NULL),

                             p('Click button to download the current bivariate plot.'),
                             downloadButton(outputId = "download_bivariate_current", label = "Download Current Plot"),

                             p('Click button to download bivariate plots for all selected marker pairs.'),
                             downloadButton(outputId = "download_bivariate_all", label = "Download All Plots"),

                             hr(),
                         ),

                            tabPanel(
                          value=3,

                           strong("Trivariate"),
                           br(),

                             varSelectInput("xvarTri", "X variable", NULL, selected = NULL),
                             varSelectInput("yvarTri", "Y variable", NULL, selected = NULL),
                            varSelectInput("zvar", "Z variable", NULL, selected = NULL),


                             hr(),
                         )
                             )

                         ),
                         tabPanel(
                         value=4,
                          strong("Extra"),
                          tabsetPanel(id='sidebar_post',
                            tabPanel(
                          value=1,
                           strong("Phenotyping"),
                            br(),

                          #  p("Please provide your .csv cell-marker table."),

                          p("Please upload your phenotyping workflow following the format described in",
                                                   strong("Help,"), "
                                        or download the example workflow found in",
                                                   downloadLink('download_example_data2', strong(' here.') )
                                                 ),

                           fluidRow(
                              column(width = 8,
                                    div(id = "file_input_div2",
                                        fileInput("phen_wfl", label = "Upload a CSV file", accept = ".csv", buttonLabel = "Browse"))),
                                        tags$div(id = "message_loading_data2", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),
                                        br(),
                              column(width = 4, actionButton("remove_file2", "Remove", icon = icon("trash")))
                            ),


                          actionButton("define_phenotype_AUTO", "Phenotype my data"),
                          br(),
                          br(),

                          uiOutput("marker_checkboxes"),
                          materialSwitch(inputId = "any_indicator",
                                                            label = "Any positive",
                                                            status = "info",
                                                            right= TRUE,
                                                            value = FALSE),
                          # Marker selection input
                          textInput("phenotype_name", "Phenotype Name"),
                          actionButton("define_phenotype", "Add phenotype definition"),
                          br(),
                          br(),
                          p('You can download the your original CSV file with phenotypes added as a new column, as well as the phenotype workflow you defined.'),
                          downloadButton(outputId = "downloadPhenotypes", label = "Download Phenotyped Data"),
                          downloadButton(outputId = "downloadPhenotypeTable", label = "Download Workflow"),
                          br(),


                        )


                        )


                        ),
                    id = "sidebartab"
                       )),


                     mainPanel(width=8,

                     tabsetPanel(
            conditionalPanel(
              condition = "input.sidebar2 == 1 && input.sidebartab == 3",
              br(),
              numericInput('intercept', 'Type Gate Value', value = ""),
              actionButton(inputId = "updateGates", label = "Update Gate"),


              conditionalPanel(
                condition = "input.intercept == ''",
                div(
                  class = "alert alert-danger",
                  "Please enter a valid gate value."
                )
              ),

        # plotOutput("gated_histogram_on_page"),
        plotOutput("gated_histogram_on_page", height = "1000px", width = "1000px"),

         # Pagination buttons
	div(
	  class = "pagination-button",
	  actionButton(inputId = "nextMarker", label = "Next Marker ", icon("arrow-right"))
	),

div(
  class = "pagination-button-back",
  actionButton(inputId = "prevMarker", label = "Previous Marker ", icon("arrow-left"))
)


      ),

      conditionalPanel(
        condition = "input.sidebartab == 3 && input.sidebar2 == 2",
        # Add numeric input fields and a button
        numericInput("gate_xvar_update", "Enter Gate for X Variable:", value = ""),
        numericInput("gate_yvar_update", "Enter Gate for Y Variable:", value = ""),
        actionButton("update_gates_bivariate", "Update Gates"),
        plotOutput("plot2", height = "1000px", width = "1000px"),
        verbatimTextOutput("prop_summary")


      ),
       conditionalPanel(
        condition = "input.sidebartab == 3 && input.sidebar2 == 3",
        # Add numeric input fields and a button
        numericInput("gate_xvar_updateTri", "Enter Gate for X Variable:", value = ""),
        numericInput("gate_yvar_updateTri", "Enter Gate for Y Variable:", value = ""),
        numericInput("gate_zvar_update", "Enter Gate for Y Variable:", value = ""),
        actionButton("update_gates_trivariate", "Update Gates"),
        plotlyOutput("plot_trivariate", height = "1000px", width = "1000px"), # learn this to be min and max, put 2 percent offsset multipled by like 2 percent plus minus


      ),
      conditionalPanel(
        condition = "input.sidebartab == 2 && input.sidebar1 == 1",
        # plotOutput("image_garage_output", brush="plot_brush", height = "1000px", width = "1000px"),
        plotlyOutput("image_garage_output", height = "1000px", width = "1000px"),

        verbatimTextOutput("subset_summary")
      ),
       conditionalPanel(
        condition = "input.sidebartab == 2 && input.sidebar1 == 2",
        plotOutput("tissue_score_out"),
        plotOutput("modes_plot_output"),
        # plotOutput("histogram_plots_cyles"),
        # plotOutput("CRUDEINDEX"),
        # plotOutput("suggested_removal"),
        plotOutput("klPLOToutput"),
        plotlyOutput("quality_gauge"),
        # plotlyOutput("quality_gauge_kl")

      ),


       conditionalPanel(
        condition = "input.sidebartab == 4 && input.sidebar_post == 1",
        # textOutput("phenotype_output"),
              # Create two sections: Left and Right

                column(width = 6, tableOutput("phenotypeTable")),
                column(width = 6, plotOutput("pheno_bar", height = "500px", width = "500px")),

              verbatimTextOutput("post_statistics"),


      # tableOutput("phenotypeTable")# Add your table output here

      ),


      conditionalPanel(
        condition = "input.sidebartab == 4 && input.sidebar_post == 2",
        plotOutput("icaAnalysis2"),
      ),


       conditionalPanel(
        condition = "input.sidebartab == 4 && input.sidebar_post == 0",
        # plotOutput("icaAnalysis"),
      ),


    ),


        )


        )


      ),
      tabPanel(
        strong("Help"),
        br(),
        tags$div(includeMarkdown("./documents/help.md"), style = "max-width:800px;")
      ),
      tabPanel(
        strong("FAQ"),
        br(),
        column(width = 1, ""),
        br(),
        column(
          width = 6,
          h4(strong("Q:"), tags$em(strong(
            "Why are some proportions or total sample numbers zero?"
          ))),
          p(strong("A:"),
            "Randkluft searches for a positively skewed signal emerging from background noise.
  In cases where the marker distribution is already negatively skewed or lacks a discernible positive tail,
  the algorithm terminates early without estimating a gate.
  In these situations, we recommend visual inspection of the distribution and manual gating."
          ),
          br(),

          h4(strong("Q:"), tags$em(strong(
            "Why do I get 'Disconnected from the server' after uploading my data?"
          ))),
          p(strong("A:"),
            "This typically indicates that the uploaded file does not conform to the expected input format.
  Please consult the Help section and ensure that column names, data types, and required fields
  strictly follow the documented input structure before re-uploading."
          ),
          br(),

          h4(strong("Q:"), tags$em(strong(
            "Why are upload and analysis slow?"
          ))),
          p(strong("A:"),
            "Randkluft treats all numeric columns—including spatial coordinates and DNA/Hoechst channels—as potential gating targets.
  If your input file contains many columns that are not required for analysis, removing them before upload
  can significantly improve performance and reduce processing overhead."
          ),
          br(),

          h4(strong("Q:"), tags$em(strong(
            "Which data are used to detect the gates?"
          ))),
          p(strong("A:"),
            "At each step of the workflow, Randkluft internally stores the active dataset associated with the selected panel.
  All subsequent analyses use this updated data.
  Both the modified datasets and the resulting gate estimates can be downloaded at each stage of the analysis."
          ),
          br()
        )
      ),
      tabPanel(
        strong("Contact"),
        br(),
        tags$div(includeMarkdown("./documents/contact.md"), style = "max-width:800px;")
      )
    )
  )
))
