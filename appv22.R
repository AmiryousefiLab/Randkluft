library(umap)
library(NMF)
library(tsne)
library(plotly)
library(readxl)
library(ComplexHeatmap)
library(shiny)
library(shinyjs)
library(shinyvalidate)
library(shinyWidgets)
library(shinyalert)
library(purrr)
library(seqinr)
library(colourpicker)
library(markdown)
library(rmarkdown)
library(httpuv)
library(shinydashboard)
library(bslib)
library(mclust)
library(moments)
library(multimode)
library(data.table)
library(viridis)
library(dplyr)
library(scales)
library(philentropy)
library(HistogramTools)
#library(DescToolo)
library(fastICA)
library(tiff)
library(datasets)
data(iris)
library(umap)
library(stringr)
library(tidyr)
library(PEkit)
library(sqldf)
library(ggplot2)
library(gridExtra)

# source codes
source("utils.R")

df_example = read.csv('exemplar-001--unmicst_cell.csv')

ground_truth_gates_loaded = read.csv('tuulia_data_GT.csv')

df_example_PHENOTYPE = read.csv("phenotype_table_help.csv")

options(shiny.trace = TRUE,
        shiny.maxRequestSize = 30 * 1024 ^ 3) # change max file limit
        # options(shiny.maxRequestSize = 10 * 1024 ^ 3)  # Set to 10 GB


ui <- shinyUI(fluidPage(
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
                            # tags$style(".shiny-file-input-progress {display: none}"),
                              column(width = 8, 
                                    div(id = "file_input_div", 
                                        fileInput("cell_file", label = "Upload a CSV file", accept = ".csv", buttonLabel = "Browse"))),
                                        tags$div(id = "message_loading_data", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),

                              column(width = 4, actionButton("remove_file", "Remove", icon = icon("trash")))
                            ),
                            

                         ),
                         
                         ##COMMENTED OUT THE PRE_Randkluft COMPLETELY

                        #  tabPanel(strong("Pre-Randkluft"),
                        #  value=2,
                        #  tabsetPanel(id="sidebar1",
                        #  tabPanel(value=1,
                        #  strong("Remove Cells"),
                        #  br(),
                        #     radioButtons(
                        #       "patient_number_im_gargage",
                        #       "Select Patients",
                        #       choices = character(0),
                        #       selected = character(0)
                        #     ),
                        #     hr(),
                        #      radioButtons(
                        #       "precrev_marker",
                        #       "Select Marker",
                        #       choices = character(0),
                        #       selected = character(0)
                        #     ),
                        #     sliderInput("opacity_slider", "Marker Opacity", min = 0.25, max = 1, value = 0.5, step = 0.25, ticks=FALSE),
                        #     # Add a button to capture the subset
                        #     # actionButton("RetainSubset", "Retain"),
                        #     # actionButton("RemoveSubset", "Remove"),
                        #     actionButton("RetrySubset", "Retry"),
                        #     # actionButton("SaveSubset", "Save"),
                        #     actionButton("UseSubset", "Use")
                        #  ),
                        # 
                        # 
                        # tabPanel(value=2,
                        #  strong("Tissue Scoring"),
                        #  br(),
                        #     radioButtons(
                        #       "patients_ts",
                        #       "Select Patients",
                        #       choices = character(0),
                        #       selected = character(0)
                        #     ),
                        # radioButtons(
                        #       "cycle_ts",
                        #       "Select Cycle",
                        #       choices = character(0),
                        #       selected = character(0)
                        #     ),
                        #     hr(),
                        #   actionButton("acceptChangesQC", "Accept Changes")
                        # 
                        #  )
                        # 
                        # 
                        #  )
                        # ),

                          

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

                           checkboxGroupInput("unique_patients_gating", "Select Patients", choices = NULL, selected = NULL),



                                                                             br(),


                           materialSwitch(inputId = "gen_hist_plots_on_off", 
                                                                       label = "Show gates",
                                                                       status = "danger",
                                                                       right=TRUE,
                                                                       value = FALSE ),

                            # materialSwitch(inputId = "GMM_gate_on_off", 
                            #                                            label = "GMM",
                            #                                            status = "success",
                            #                                            right=TRUE,
                            #                                            value = FALSE ),


                           
                          # Add the progress bar div element here
                           tags$div(id = "message", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),
                           tags$div(id = "message_gen_hist", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),


                          p('Click button to download a .csv  file with the gate estimations.'),
                          downloadButton(outputId = "downloadEstimations", label = "Download Gate Estimates"),

                          p('Click button to download current set of plots that are displayed.'),
                          downloadButton(outputId = "downloadcurrent", label = "Download Current Plot"),

                          p('Click button to download pdf file of plots for all selected markers and patients.'),
                          downloadButton(outputId = "downloadall", label = "Download All Plots")

                          # p('Click button to download pdf file of plots for all selected markers and patients.'),
                          # downloadButton(outputId = "downloadPDF_4panel", label = "Download 4 Panel")


                             ), 
                             
                             
                              tabPanel(
                          value=2,

                           strong("Bivariate"),
                           br(),
                           
                             varSelectInput("xvar", "X variable", NULL, selected = NULL),
                             varSelectInput("yvar", "Y variable", NULL, selected = NULL),
                            #  checkboxGroupInput(
                            #   "patient_number_plot2",
                            #   "Select unique patients",
                            #   choices = NULL,
                            #   selected = NULL
                            # ),
                            radioButtons(
                              "patient_number_plot2",
                              "Select Patients",
                              choices = character(0),
                              selected = character(0)
                            ),
                             hr(),
                             # Add a horizontal rule
                            #  checkboxInput("by_species", "Show individuals", TRUE),
                            #  checkboxInput("show_margins", "Show marginal plots", TRUE),
                            #  checkboxInput("smooth", "Add smoother"),
                           
                         ),

                            tabPanel(
                          value=3,

                           strong("Trivariate"),
                           br(),
                           
                             varSelectInput("xvarTri", "X variable", NULL, selected = NULL),
                             varSelectInput("yvarTri", "Y variable", NULL, selected = NULL),
                            varSelectInput("zvar", "Z variable", NULL, selected = NULL),

                            #  checkboxGroupInput(
                            #   "patient_number_plot2",
                            #   "Select unique patients",
                            #   choices = NULL,
                            #   selected = NULL
                            # ),
                            radioButtons(
                              "trivar_patient",
                              "Select Patients",
                              choices = character(0),
                              selected = character(0)
                            ),
                             hr(),
                             # Add a horizontal rule
                            #  checkboxInput("by_species", "Show individuals", TRUE),
                            #  checkboxInput("show_margins", "Show marginal plots", TRUE),
                            #  checkboxInput("smooth", "Add smoother"),
                           
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
                            # tags$style(".shiny-file-input-progress {display: none}"),
                              column(width = 8, 
                                    div(id = "file_input_div2", 
                                        fileInput("phen_wfl", label = "Upload a CSV file", accept = ".csv", buttonLabel = "Browse"))),
                                        tags$div(id = "message_loading_data2", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),
                                        br(),  
                              column(width = 4, actionButton("remove_file2", "Remove", icon = icon("trash")))
                            ),
                            

                        #  # Marker selection input
                        #   textInput("phenotype_name", "Phenotype Name"),
                          # checkboxGroupInput("selected_columns_phenotyping", "Select Markers", choices = NULL, selected = NULL),
                          
                          # checkboxInput("marker_positive", "Is it a Positive Marker?"),
                          # actionButton("add_marker", "Add Marker"),
                          # actionButton("define_phenotype", "Define Phenotype"),
                          actionButton("define_phenotype_AUTO", "Phenotype my data"),
                          br(),
                          br(),
                          # p('Click button to download a .csv  file with the phenotypes'),
                          # downloadButton(outputId = "downloadPhenotypes", label = "Download Phenotype"),
                          # downloadButton(outputId = "downloadPhenotypeTable", label = "Download Workflow"),
                          # br(),
                          # checkboxGroupInput("selected_columns_phenotyping", "Select Markers", choices = NULL, selected = NULL)
                          # checkboxGroupInput("selected_columns_phenotyping", "Select Markers"),
                          # tags$div(id = "updateProgressPhenotypePost", style = "font-size: 20px; position: fixed; bottom: 0; right: calc(5% + 10px);"),
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




                          #  # Checkbox group for marker selection with individual material switches
                          #         fluidRow(
                          #           # tags$style(".shiny-file-input-progress {display: none}"),
                          #           column(width = 8,
                          #             checkboxGroupInput("selected_columns_phenotyping", "Select Markers", choices = NULL, selected = NULL)
                          #           ),
                          #           column(width = 4,
                          #             # Generate three material switches dynamically
                          #             lapply(reactive({ seq_along(input$selected_columns_phenotyping) }), function(i) {
                          #               materialSwitch(
                          #                 inputId = paste0("pos_neg_switch_", i),
                          #                 label = "- / +",
                          #                 status = "danger",
                          #                 right = TRUE,
                          #                 value = FALSE
                          #               )
                          #             })
                          #           )
                          #         )
                        )

                        
                  

                        




                        )
                        
                        
                        
                        ),                       
                     ,
                        
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
        # actionButton(inputId = "nextPatient", label = "Next Patient"),
        div(
  class = "pagination-button",
  actionButton(inputId = "nextPatient", label = "Next Patient ", icon("arrow-right"))
),

div(
  class = "pagination-button-back",
  actionButton(inputId = "prevPatient", label = "Previous Patient ", icon("arrow-left"))
),

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


      # conditionalPanel(
      #   condition = "input.sidebartab == 1",
      #   verbatimTextOutput("sum")
      # ), 

  
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


# options(shiny.maxRequestSize=100*1024^2)
options(shiny.maxRequestSize=100*1024^3)
server <- shinyServer(function(input, output, session) {

  remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
  low_threshold <- quantile(target, low_percentile / 100)
  high_threshold <- quantile(target, high_percentile / 100)
  target <- target[target >= low_threshold & target <= high_threshold]
  return(target) }

    #  # Define a function to find the mode of a vector
    #   find_mode <- function(x) {
    #     ux <- unique(x)
    #     ux[which.max(tabulate(match(x, ux)))]
    #   }

    
     # Define a function to find the mode of a vector
      find_mode2 <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }


  outputinterceptreactive <- reactiveVal(NULL)

  plot_gating_individual <- function(target, output, marker, GroundTruthInput, do_ggplot = TRUE) {
    
    if (!is.null(outputinterceptreactive())) {
      output$cutoff <- outputinterceptreactive()
      output$N_removed <- sum(target > output$cutoff)
      output$percentage_removed <- round(output$N_removed / length(target), 3)
    } 
    
    print("mode is here")
    print(find_mode2(target))
    
    if (do_ggplot == TRUE) {
      plot <- ggplot(data.frame(x = target), aes(x = x)) +
        geom_histogram(aes(y = after_stat(density)), bins = 100,
                       fill = "lightblue",
                       color = "black") +
        geom_vline(xintercept = output$cutoff, color = "red") +
        geom_vline(xintercept = GroundTruthInput, color = "blue") +
        geom_density() +
        geom_text(x = output$cutoff, size = 5, # Size of the text
                  y = find_mode2(target), # Set y to the mode of the target
                  label = paste0("Gate:", round(output$cutoff, 2)),
                  hjust = 0,
                  vjust = -0.5, # Adjust vertical alignment
                  fill = "white",
                  color = "black") +
        labs(
          x = "",
          y = "",
          title = paste(
            "N+=", output$N_removed, 
            " ", 
            "+R=", round(output$percentage_removed, 3),
            " ", 
            "Gate=", round(output$cutoff, 2)
          )
        ) +
        xlab(marker) +
        theme(
          plot.title = element_text(size = 18, face = "bold", family = "Arial"),
          axis.title = element_text(size = 16, family = "Arial"),
          axis.text = element_text(size = 14, family = "Arial")
        )
      return(plot)
    } 
  }
  
  plot_gating_grid <- function(df, unique_ids, column_indices) {
    plot_list <- list()
    
    names_columns <- names(df)
    n_plots_per_row <- length(column_indices)  # Adjust this as needed
    n_plots <- length(unique_ids) * length(column_indices)
    n_rows <- ceiling(n_plots / n_plots_per_row)
    
    print(n_rows)
    print(n_plots)
    print(n_plots_per_row)
    
    for (unique_id in unique_ids) {
      par(mfrow = c(5, 4), mar = c(5.1, 3, 4.1, 2))
      
      for (j in column_indices) {
        target <- df[df$imageid == unique_id, j]
        target <- remove_outliers2(target)
        target <- target[!is.infinite(target)]
        
        marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
        output <- skew_gate(target, 0.01)
        gtGate <- ground_truth_gates_loaded$Gate[ground_truth_gates_loaded$Patient == unique_id & ground_truth_gates_loaded$Marker == names_columns[j]]
        p_temp <- plot_gating_individual(target, output, marker, gtGate)
        plot_list <- c(plot_list, list(p_temp))
      }
    }
    
    p <- gridExtra::grid.arrange(grobs = plot_list, ncol = n_plots_per_row, nrow = n_rows)
    return(p)
  }
  

observeEvent(input$updateGates, {
  # Get the new gate value from the numeric input
  new_gate_value <- input$intercept

   if (!is.na(new_gate_value) && !is.null(new_gate_value) && new_gate_value != "") {
     # Update the reactive value with the new gate value
      outputinterceptreactive(new_gate_value)
  } else {
     # Update the reactive value with the new gate value
      outputinterceptreactive(NULL)
  }

  # # Update the reactive value with the new gate value
  # outputinterceptreactive(new_gate_value)

  # Re-run the code for generating plots
  generatePlot()

  selected_columns_man <- input$selected_columns
  selected_patients_man <- input$unique_patients_gating

  #outputinterceptreactive(NULL)
  chosen_patient <- selected_patients_man[current_patient()]
  chosen_marker <- selected_columns_man[current_marker()]

 
  updated_df <- csv_save_file_path()
    updated_df$Gate[
      updated_df$Marker == chosen_marker &
      updated_df$Patient == chosen_patient
    ] <- new_gate_value
    csv_save_file_path(updated_df)

})






plot_gating_grid_pdf_all = function(df, unique_ids, column_names) {
  plot_list = list()

  names_columns = names(df)
  # Calculate the number of rows and columns for mfrow
  n_plots_per_row = length(column_indices)  # Adjust this as needed
  n_plots = length(unique_ids) * length(column_indices)
  n_rows = ceiling(n_plots / n_plots_per_row)

  print(n_rows)
  print(n_plots)
  print(n_plots_per_row)

  for (unique_id in unique_ids) {
    par(mfrow=c(5,4), mar = c(5.1, 3, 4.1, 2))

  # for (a in c(1:2)) {
    # par(mfrow = c(n_rows, n_plots_per_row), mar = c(5.1, 4.1, 4.1, 2.1))
    # par(mfrow = c(1, 1))

    for (j in column_indices) {
    target <- df[df$imageid == unique_id, j]




    target <- remove_outliers2(target)


  

    target <- target[!is.infinite(target)]

    # print(target)

    marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
    output <- skew_gate(target, 0.01)
    
    p_temp = plot_gating_individual(target, output, marker)
    plot_list = c(plot_list, list(p_temp))
  }
  }

  # Limit the number of plots to a maximum of 180
  # plot_list <- plot_list[1:min(length(plot_list), 1)]
  
  # p = gridExtra::grid.arrange(grobs = plot_list, ncol = n_plots_per_row, nrow= n_rows)
  p = gridExtra::grid.arrange(grobs = plot_list)

  return(p)
}

  skew_gate <- function(x, alpha=0.01) {
    sk <- moments::skewness(x)
    
    n <- length(x)
    a <- locmodes(x)$locations[1]
    
    b <- max(x)
    
    if (sk < 0) {
      message("The skewness is negative!")


      outputGMM <- Mclust(x, G = 2)
      gmm_gate <- mean(outputGMM$parameters$mean)

       n_removed <- sum(x > gmm_gate)
       perc_removed <- round(n_removed / n, 3)

      return(list(
      skewness = sk,
      cutoff = gmm_gate,
      N = n,
      N_removed = n_removed,
      percentage_removed = perc_removed,
      returnvalplot = x[which(x> a+(b-a)/2)]
    ))
    }

    iteration <- 0
    
    while (abs(sk) > alpha & iteration <= 100) {
      # print(iteration)
      if (sk >= 0) {
        b <- a + (b - a) / 2
      } else {
        a <- a + (b - a) / 2
      }
      a <- min(a, b)
      b <- max(a, b)
      
      sk <- skewness(x[which(x < b)])
      if (is.nan(sk)) {
        message("Warning: skewness is NaN")
        break
      }
      iteration <- iteration + 1
    }
    
    n_removed <- sum(x > b)
    perc_removed <- round(n_removed / n, 3)

    print(n_removed)
    print(perc_removed)
    list(
      skewness = sk,
      cutoff = b,
      N = n,
      N_removed = n_removed,
      percentage_removed = perc_removed,
      returnvalplot = x[which(x> a+(b-a)/2)]
    )
  }

  




  get_gates_csv <- function(dataframe, csv_save_file_path) {


  data <- dataframe # Replace with your actual CSV file path


  # Get unique values from the imageid column
  unique_imageids <- unique(data$imageid)


  alpha <- 0.01  # Set your desired alpha value

  # Create an empty data frame to store the results
  results_df <- data.frame(Patient = character(),
                          Marker = character(),
                          Gate = numeric(),
                          stringsAsFactors = FALSE)


# Before the loop, set the total number of iterations
total_iterations <- length(unique_imageids) * ncol(data)

progress <- Progress$new(session, min = 1, max = total_iterations)
progress$set(message = "Randkluft in action...", value = 0)

# Loop through each unique imageid
for (imageid in unique_imageids) {
  sub_data <- data[data$imageid == imageid, -1]

  # Loop through each column in sub_data
  for (col_idx in 1:ncol(sub_data)) {
    col_name <- colnames(sub_data)[col_idx]
    values <- sub_data[, col_idx]

    values <- remove_outliers2(values)
    values <- values[!is.infinite(values)]

    result <- skew_gate(values, alpha)$cutoff

    # Append row to results_df
    results_df <- rbind(results_df, data.frame(Patient = imageid, Marker = col_name, Gate = result))

    cat("Image ID:", imageid, "Marker:", col_name, "Result:", result, "\n")

    # Calculate the current iteration
    current_iteration <- (match(imageid, unique_imageids) - 1) * ncol(sub_data) + col_idx

    # Update the progress bar
    progress$set(message = "Randkluft in Action...", value = current_iteration)
  }
}

progress$close()

return(results_df)




  }

  get_gates_csv_single <- function(dataframe, csv_save_file_path) {


  remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
  low_threshold <- quantile(target, low_percentile / 100)
  high_threshold <- quantile(target, high_percentile / 100)
  target <- target[target >= low_threshold & target <= high_threshold]
  return(target)
}


  data <- dataframe # Replace with your actual CSV file path




  # List of markers
  markers <- colnames(data)[-1]

  # Create a list to store ggplot objects
  histogram_list <- list()

   # Define skew_gate function above

  alpha <- 0.01  # Set your desired alpha value

  # Create an empty data frame to store the results
  results_df <- data.frame(Patient = character(),
                          Marker = character(),
                          Gate = numeric(),
                          stringsAsFactors = FALSE)

  # Loop over unique imageids
  unique_imageids <- unique(data$imageid)
  for (imageid in unique_imageids) {
    for (marker in markers) {
      # Subset data for the current imageid and marker
      subset_data <- data[data$imageid == imageid, c("imageid", marker)]

    

  
      # Calculate gate value using skew_gate function
      result_to_plot <- skew_gate(subset_data[[marker]])  # Replace with your desired alpha
      skewness_value <- result_to_plot$skewness
      cutoff_value <- result_to_plot$cutoff

       # Append row to results_df
      results_df <- rbind(results_df, data.frame(Patient = imageid, Marker = marker, Gate = cutoff_value))
      
      cat("Image ID:", imageid, "Marker:", marker, "Result:", cutoff_value, "\n")
    }}




  return(results_df)


  }
  
  

  #shinyjs::hide("GMM_gate_on_off")



  shinyjs::hide("patient_number")
  
  observe({

    req(uploaded_df())

    runjs("
      function animateDots() {
        var messageDiv = document.getElementById('message_loading_data');
        var dots = '';
        setInterval(function() {
          if (dots.length === 3) dots = '';
          else dots += '.';
          messageDiv.innerText = 'Loading Data' + dots;
        }, 500);
      }
      animateDots();
      document.getElementById('message_loading_data').style.display = 'block';
    ")

    # Show the checkboxGroupInput
    shinyjs::show("selected_columns")

    # Show the checkboxGroupInput
    shinyjs::show("unique_patients_gating")
    
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", "DNA2", "Hoechst1", "Hoechst2", "Hoechst3", 
                            "Hoechst4", "Hoechst5", "Hoechst6", "Hoechst7", "Hoechst8", "Hoechst9", "Hoechst10", "Hoechst_1", "Hoechst_2", "Hoechst_3", 
                            "Hoechst_4", "Hoechst_5", "Hoechst_6", "Hoechst_7", "Hoechst_8", "Hoechst_9", "Hoechst_10", "DAPI1", "DAPI2", "DAPI3", "DAPI4", "DAPI5", "DAPI6", "DAPI7", "DAPI8", "DAPI9", "DAPI10",
                            "DAPI_1", "DAPI_2", "DAPI_3", "DAPI_4", "DAPI_5", "DAPI_6", "DAPI_7", "DAPI_8", "DAPI_9", "DAPI_10", "DNA1", "DNA3", "DNA2", "DNA4", "DNA5", "DNA6","DNA7","DNA8", "DNA9", "DNA10", "DNA11",
                            "DNA12", "DNA13", "DNA_1", "DNA_3", "DNA_2", "DNA_4", "DNA_5", "DNA_6","DNA_7","DNA_8", "DNA_9", "DNA_10", "DNA_11",
                            "DNA_12", "DNA_13", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
                            "MinorAxisLength", "Extent", "Solidity", "Orientation", "") # List the columns to exclude

      # Define the regular expression pattern to identify columns to be excluded
  exclude_pattern <- "DNA|DAPI|Hoechst"

  # Use grep to get the column names matching the pattern
  exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

  # Add the excluded columns to the original 'columns_to_exclude' vector
  columns_to_exclude <- c(columns_to_exclude, exclude_columns)

     # Define the regular expression pattern to identify columns to be excluded
  exclude_pattern2 <- "_positivity"

  # Use grep to get the column names matching the pattern
  exclude_columns2 <- grep(exclude_pattern2, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

  columns_to_exclude <- c(columns_to_exclude, exclude_columns2)



    # Select all column names except the ones to exclude
    column_names <- setdiff(names(uploaded_df()), columns_to_exclude)
    reactive_markers_selected(column_names)


    updateCheckboxGroupInput(session, "selected_columns", choices = column_names, selected = column_names)
    updateCheckboxGroupInput(session, "selected_columns_phenotyping", choices = column_names, selected = column_names)
    updateCheckboxGroupInput(session, "unique_patients_gating", choices = unique_patients(), selected = unique_patients())






    

    runjs("document.getElementById('message_loading_data').style.display = 'none';")

  })

# Define a reactiveVal for plot height
#   plotHeight <- reactiveVal("2060px")


  reactive_markers_selected <- reactiveVal(NULL)

  uploaded_df <- reactiveVal(NULL)
# Define a reactive value for the GMM gate switch
  gmm_gate_switch <- reactiveVal(TRUE)  # Initialize with the desired default value
  # You can set the initial value based on the user's preference or the default value

  histplot_gate_switch <- reactiveVal(FALSE)  # Initialize with the desired default value
  # You can set the initial value based on the user's preference or the default value

  # Observe the change in the GMM gate switch input and update the reactive value
  observeEvent(input$GMM_gate_on_off, {
    gmm_gate_switch(input$GMM_gate_on_off)
  })

  # Observe the change in the GMM gate switch input and update the reactive value
  observeEvent(input$gen_hist_plots_on_off, {

    histplot_gate_switch(input$gen_hist_plots_on_off)
  })

 # This button will reset the inFileLoad
  observeEvent(input$remove_file, {
    reset("cell_file")  # reset is a shinyjs function
    # reset("selected_columns")
    # shinyjs::hide("gated_histogram_on_page")

     # Hide the checkboxGroupInput
    shinyjs::hide("selected_columns")
    shinyjs::hide("unique_patients_gating")

    # shinyjs::hide("marker")
    # shinyjs::hide("xvar")
    # shinyjs::hide("yvar")
    shinyjs::hide("gated_histogram_on_page")
    shinyjs::hide("patient_number")
    shinyjs::hide("patient_number_plot2")
    shinyjs::hide("histogram_plot")
    # shinyjs::hide("plot2")
    shinyjs::reset("patient_number_plot2")

    shinyjs::hide("summary_output")
    shinyjs::hide("image_garage_output")
    shinyjs::hide("selected_columns_phenotyping")

  

    

    # shinyjs::toggleState(id = "gen_hist_plots_on_off")  # Toggle the switch

    # output$plot2 <- NULL 

    # Clear the select inputs and update their choices


    updateSelectInput(session, "yvar", choices = character(0), selected = character(0))
    updateSelectInput(session, "xvar", choices = character(0), selected = character(0))
    updateSelectInput(session, "marker", choices = character(0), selected = character(0))



  })

  #  # This button will reset the inFileLoad
  # observeEvent(input$image_garage, {
    
  #   req(uploaded_df())

  #   print(colnames(uploaded_df()))



  # })

  # Define a function to check if any numerical column has values greater than 20
# check_if_logged <- function(df) {
#   numeric_cols <- sapply(df, is.numeric)
#   max_value <- max(df[numeric_cols])
#   return(max_value > 20)
# }

check_if_logged <- function(data) {
  max_vals <- sapply(data, max, na.rm = TRUE)
  any(max_vals > 20)
}

sub_data_logged <- reactiveVal(NULL)

 read_batch_with_progress = function(file_path,nrows,no_batches){
    progress = Progress$new(session, min = 1,max = no_batches)
    progress$set(message = "Processing File...")
    seq_length = ceiling(seq.int(from = 2, to = nrows-2,length.out = no_batches+1))
    seq_length = seq_length[-length(seq_length)]
  
    
    #read the first line
    
    
    
    # Generate a sample of 10 numbers from 1 to 20 without replacement
    #sample_numbers <- sample(1:nrows, 10000, replace = FALSE)
    # Create the SQL query string
    #query <- sprintf("SELECT * FROM file WHERE CellID IN (%s)", paste(sample_numbers, collapse = ", ")) #the first row should be called CellID
    #df <- read.csv.sql(file_path, sql = query, dbname = tempfile())
    
    #col_names = colnames(df)

    for(i in seq_along(seq_length)){
      progress$set(value = i)
      if(i == no_batches) chunk_size = -1 else chunk_size = seq_length[i+1] - seq_length[i]

      # df_temp = read.csv(file_path, skip = seq_length[i], nrows = chunk_size, header = FALSE, stringsAsFactors = FALSE)
      df_temp = read.csv.sql(file_path, sql = query, dbname = tempfile())
       # read.csv(file_path, skip = seq_length[i], nrows = chunk_size, header = TRUE, check.names=FALSE)

      colnames(df_temp) = col_names
      df = rbind(df,df_temp)
    }

    progress$close()
    return(df)
  }


  # observe({
  #   output$sum <- renderPrint({
  #     print(head(df(), 10))
  #   }) 
  # }) 

  

  # Event handler for the button click to read the uploaded CSV file
  observeEvent(input$cell_file, {
    # runjs("
    #   function animateDots() {
    #     var messageDiv = document.getElementById('message_loading_data');
    #     var dots = '';
    #     setInterval(function() {
    #       if (dots.length === 3) dots = '';
    #       else dots += '.';
    #       messageDiv.innerText = 'Loading Data' + dots;
    #     }, 500);
    #   }
    #   animateDots();
    #   document.getElementById('message_loading_data').style.display = 'block';
    # ")

    n_rows = length(count.fields(input$cell_file$datapath))
    

    # df_out = read_batch_with_progress(input$file1$datapath,n_rows,10)
    # Generate a sample of 10 numbers from 1 to 20 without replacement
    if (n_rows > 70000){
    sample_numbers <- sample(1:n_rows, 70000, replace = FALSE)}
    else{sample_numbers <- seq(1:n_rows)}
    # Create the SQL query string
    query <- sprintf("SELECT * FROM file WHERE CellID IN (%s)", paste(sample_numbers, collapse = ", ")) #the first row should be called CellID
    uploaded_df(read.csv.sql(input$cell_file$datapath, sql = query, dbname = tempfile()))
    
    

    # uploaded_df(read.csv(input$cell_file$datapath, header = TRUE, check.names=FALSE))

    # uploaded_df(read(input$cell_file$datapath))

    print("DONE LOADING CSV INITIAL")


      new_column_name <- "imageid"
      old_column_name <- "imageID"

      # Get current value of the reactiveVal
      current_df <- uploaded_df()

      if (old_column_name %in% colnames(current_df)) {
        colnames(current_df)[colnames(current_df) == old_column_name] <- new_column_name
        uploaded_df(current_df)  # Update the reactiveVal with the modified dataframe
      } else {
        print("Column not found")
      }

      print(colnames(current_df))

      # Get the indices of columns with empty string as name
      empty_string_cols <- which(colnames(current_df) == "")

      # Remove columns with empty string as name
      if (length(empty_string_cols) > 0) {
        current_df <- current_df[,-empty_string_cols, drop = FALSE]
        uploaded_df(current_df)
      }

    print('colnames after drop blank')

    print(colnames(uploaded_df()))
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
    "MinorAxisLength", "Extent", "Solidity", "Orientation", "", "DNA6a")  # List the columns to exclude
    # Select all column names except the ones to exclude

    # # Define the regular expression pattern to identify columns to be excluded
    # exclude_pattern <- "DNA|DAPI"

    # # Use grep to get the column names matching the pattern
    # exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

    # Add the excluded columns to the original 'columns_to_exclude' vector
    # columns_to_exclude <- c(columns_to_exclude, exclude_columns)


    column_names <- setdiff(names(uploaded_df()), columns_to_exclude)

    print(column_names)



    print("pre")
    max_values <- sapply(uploaded_df()[column_names], max)
    print("post")

  

    uploaded_intermediate <- uploaded_df()
    uploaded_intermediate[column_names] <- sapply(uploaded_intermediate[column_names], as.numeric)
    uploaded_intermediate[column_names] <- data.frame(lapply(uploaded_intermediate[column_names], function(x) ifelse(x >= 0 & x < 10, 10, x)))

    if (any(max_values > 20)) {
        cat("Working with raw data\n")
        uploaded_intermediate[column_names] <- sapply(uploaded_intermediate[column_names], log)
        uploaded_df(uploaded_intermediate)
      } else {
        cat("Working with logged values\n")
      }

      print("crash")
      

        #uploaded_intermediate <- uploaded_df()
      # Assuming uploaded_df() is your uploaded dataframe

      if (!"imageid" %in% colnames(uploaded_df())) {
        uploaded_intermediate$imageid <- rep("image", nrow(uploaded_df()))
        uploaded_df(uploaded_intermediate)
      }
      # #removing random rows based on the size of the file more than 10 megab
      #  if (object.size(uploaded_intermediate)/1024/1024 > 10) {
      #    #sample_size<- round(dim(uploaded_df(uploaded_intermediate))[1] * (10/object.size(uploaded_df(uploaded_intermediate))/1024/1024))
      #    uploaded_df(uploaded_intermediate[sample(dim(uploaded_intermediate)[1], 400000), ])
      #  }
      if (any(c("DNA1", "DNA_1", "DAPI1", "DAPI_1", "Hoechst1", "Hoechst_1") %in% colnames(uploaded_df()))) {
        uploaded_intermediate$DNA1 <- rep(1, nrow(uploaded_df()))
        uploaded_df(uploaded_intermediate)
      }

      if ("X" %in% colnames(uploaded_df())) {
        colnames(uploaded_intermediate)[which(names(uploaded_intermediate) == "X")] <- "X_centroid"
        uploaded_df(uploaded_intermediate)
      }


        if ("Y" %in% colnames(uploaded_df())) {
        colnames(uploaded_intermediate)[which(names(uploaded_intermediate) == "Y")] <- "Y_centroid"
        uploaded_df(uploaded_intermediate)
        }

     
      

          print("DONE LOADING CSV FINAL after POST")


          print(colnames(uploaded_df()))

          #print(uploaded_df()$ELANE)


    showNotification("File Ready for Use", duration = 10, id = "message")


     

    # Show the checkboxGroupInput
    shinyjs::show("selected_columns")
    shinyjs::show("unique_patients_gating")
        shinyjs::show("selected_columns_phenotyping")




  })
  




  pdf_file_path <- reactiveVal(NULL)
  csv_save_file_path <- reactiveVal(NULL)
  resultdf_reactive <- reactiveVal(NULL)



  observeEvent(input$run_gate, {
    req(uploaded_df())

    selected_columns <- input$selected_columns

    print(colnames(uploaded_df()))


    print(selected_columns)

    selcol_len = length(selected_columns)

    unique_patients <- input$unique_patients_gating  # Assuming this function returns unique image IDs

    print(unique_patients)


  # # Subsetting the data based on selected columns
  #   sub_data <- uploaded_df()[, c("imageid", selected_columns)]

    # Subsetting the data based on selected columns and unique image IDs
    sub_data <- uploaded_df() %>%
        filter(imageid %in% unique_patients) %>%
        select(imageid, all_of(selected_columns))


     print(head(sub_data))

    pdf_file_name <- "temp_histograms.pdf"  # Get user-provided file name
    csv_save_file_name <- "temp_gates.csv"  # Get user-provided file name

    # Show the progress bar
    
    
    #  runjs("
    #   function animateDots() {
    #     var messageDiv = document.getElementById('message');
    #     var dots = '';
    #     setInterval(function() {
    #       if (dots.length === 3) dots = '';
    #       else dots += '.';
    #       messageDiv.innerText = 'Randkluft in Action' + dots;
    #     }, 500);
    #   }
    #   animateDots();
    #   document.getElementById('message').style.display = 'block';
    # ")


  
   
    if (length(selected_columns)==1){
      print("hi")
      print(colnames(sub_data))
      resultdf_to_save = get_gates_csv_single(sub_data, csv_save_file_name)
    }
    else {
      resultdf_to_save = get_gates_csv(sub_data, csv_save_file_name)
      resultdf_reactive(resultdf_to_save)
    }


    print(resultdf_to_save)
    
    csv_save_file_path(resultdf_to_save)

    # showNotification("Randkluft Found", duration = 1, id = "message")

    # Clear the message
    # runjs("document.getElementById('message').style.display = 'none';")
    # runjs("document.getElementById('message').style.display = 'none';")

      # Hide the progress bar when the task is complete

      dataframe_pos <- uploaded_df()

    # Get all marker column names (excluding non-marker columns)
    marker_columns <- selected_columns

    # Iterate over each marker column and add positivity/negativity column
    for (chosen_patient in unique_patients) {
    for (marker_col in selected_columns) {
      gate_value <- resultdf_reactive()$Gate[
        resultdf_reactive()$Marker == marker_col &
        resultdf_reactive()$Patient == chosen_patient
      ]

      
      # Add marker_positivity column based on the marker and its gate value
      dataframe_pos[[paste0(marker_col, "_positivity")]] <- sapply(dataframe_pos[[marker_col]], determine_positivity, gate_value)
    }
    }

    # Print first few rows to check the changes
    print(head(dataframe_pos))

    print(unique(dataframe_pos$ELANE_positivity))

    uploaded_df(dataframe_pos)
      

    shinyjs::show("gated_histogram_on_page")
    

    # Hide the progress bar when the task is complete


    # Display a completion message to the user
    shinyalert::shinyalert(
      title = "Randkluft Found",
      text = "Gates Saved in CSV file and set on histograms. You may download these files in the Download section.",
      type = "success"
    )

  })

  # Function to determine marker positivity/negativity based on gate value and marker intensity
determine_positivity <- function(marker_intensity, gate_value) {
  if (marker_intensity > gate_value) {
    return("+")  # Marker intensity above gate value is positive
  } else {
    return("-")  # Marker intensity below gate value is negative
  }
}

  # Define reactive values
current_marker <- reactiveVal(1)
current_patient <- reactiveVal(1)


current_4panel <- reactiveValues(p1=NULL,p2=NULL,p3=NULL, p4=NULL)

observeEvent(input$nextPatient, {

  outputinterceptreactive(NULL)

  
  observe({
  updateNumericInput(session, "intercept", value = '')
})
  current_patient((current_patient() %% length(unique_patients_gating())) + 1)
})

observeEvent(input$prevPatient, {
 
   outputinterceptreactive(NULL)


    observe({
  updateNumericInput(session, "intercept", value = '')
})



  current_patient(ifelse(current_patient() == 1, length(unique_patients_gating()), current_patient() - 1))
})

# Define a reactiveVal to store the user-defined intercept value
user_defined_intercept <- reactiveVal(NULL)

observeEvent(input$nextMarker, {

  outputinterceptreactive(NULL)

  updateNumericInput(session, "intercept", value = '') 

  current_marker((current_marker() %% length(selected_columns())) + 1)
})


observeEvent(input$prevMarker, {

    outputinterceptreactive(NULL)

 
      observe({
  updateNumericInput(session, "intercept", value = '')
})



  current_marker(ifelse(current_marker() == 1, length(selected_columns()), current_marker() - 1))
})

   get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


          observeEvent(input$gen_hist_plots_on_off, {
  req(uploaded_df())
  
  if (histplot_gate_switch()) {
    # Run the code for generating histograms
    output$gated_histogram_on_page <- renderPlot({
      generatePlot()
    })
    shinyjs::show("gated_histogram_on_page")

  }

  else {
        shinyjs::hide("gated_histogram_on_page")
  }
})


  generatePlot <- function() {
    req(uploaded_df())

      runjs("
      function animateDots() {
        var messageDiv = document.getElementById('message_gen_hist');
        var dots = '';
        setInterval(function() {
          if (dots.length === 3) dots = '';
          else dots += '.';
          messageDiv.innerText = 'Generating plots' + dots;
        }, 500);
      }
      animateDots();
      document.getElementById('message_gen_hist').style.display = 'block';
    ")

    df <- uploaded_df()
  
 # Get the column indices for the selected columns
    column_indices <- sapply(selected_columns(), function(col_name) {
      which(names(df) == col_name)
    })

    unique_imageids <- unique(df$imageid)
    
    selected_patients_reactively <- unique_patients_gating()

    plot_list <- list()  # Create a list to store plots

    unique_id <- selected_patients_reactively[current_patient()]
    column_index <- column_indices[current_marker()]

    print(selected_patients_reactively)
    print(unique_patients_gating())
    print(column_indices)
    print(unique_id)
    print(column_index)

    print(current_patient())
    print(current_marker())


    histogram_plot <- plot_gating_grid(df, selected_patients_reactively[current_patient()], column_indices[current_marker()])

      selected_columns_man <- input$selected_columns
      selected_patients_man <- input$unique_patients_gating

          filtered_data <- uploaded_df() %>%
    filter(imageid %in% selected_patients_man[current_patient()]) %>%
    select(imageid, all_of(selected_columns_man[current_marker()]))
          
    filtered_data_xy <- uploaded_df() %>%
    filter(imageid %in% selected_patients_man[current_patient()]) 
    
 

chosen_patient <- selected_patients_man[current_patient()]
chosen_marker <- selected_columns_man[current_marker()]
digrepresentation <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid, color = filtered_data_xy[[chosen_marker]])) +
  geom_point(shape = 20, size = 0.3, alpha = 0.5) +
  scale_color_gradient(low = "grey30", high = "white") +
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
        axis.title = element_text(size = 16, family = "Arial"),
        axis.text = element_text(size = 14, family = "Arial")) +
  xlab("X Centroid") +
  ylab("Y Centroid") +
  labs(title = 'Digital Representation') 

gate_value <- resultdf_reactive()$Gate[
  resultdf_reactive()$Marker == chosen_marker &
    resultdf_reactive()$Patient == chosen_patient
]

if (!is.null(outputinterceptreactive())) {
  gate_value <- outputinterceptreactive()
}

temp_df <- filtered_data_xy
temp_df[[chosen_marker]][temp_df[[chosen_marker]] < gate_value] <- 0

density_values <- get_density(filtered_data_xy$X_centroid, filtered_data_xy$Y_centroid, n = 100)

contour_plot <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
  geom_point(
    aes(color = ifelse(filtered_data_xy[[chosen_marker]] > gate_value,
                       density_values, 0)),
    size = 0.3,
    alpha = 0.35
  ) +
  geom_density_2d(
    data = subset(filtered_data_xy, filtered_data_xy[[chosen_marker]] > gate_value),
    color = "black"
  ) +  
  scale_color_viridis(
    option = "turbo",
    name = "Intensity"
  ) +
  theme(
    legend.position = c(0.98, 0.98),        # top-right inside plot
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      color = "grey70"
    ),
    legend.key.height = unit(12, "pt"),
    legend.key.width  = unit(6, "pt"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
    axis.title = element_text(size = 16, family = "Arial"),
    axis.text  = element_text(size = 14, family = "Arial"),
    panel.background = element_rect(fill = "white"),  # very pale grey for the plotting area
    plot.background = element_rect(fill = "white")
  ) +
  labs(
    title = "Positive Density",
    x = "X Centroid",
    y = "Y Centroid"
  )
overlay_plot2 <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
  geom_point(aes(color = ifelse(filtered_data_xy[[chosen_marker]] > gate_value, "Positive", "Negative")), size = 0.3, alpha = 0.35) +
  scale_color_manual(guide = guide_legend(title = ""), values = c("Positive" = "red", "Negative" = "grey")) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
        axis.title = element_text(size = 16, family = "Arial"),
        axis.text = element_text(size = 14, family = "Arial"),
        text = element_text(family = "Arial"),
        panel.background = element_rect(fill = "white"),  # very pale grey for the plotting area
        plot.background = element_rect(fill = "white")) +  # catches any remaining text elements
  labs(title = 'Positive Cells', x = "X Centroid", y = "Y Centroid")

arranged_plots <- grid.arrange(histogram_plot, digrepresentation, overlay_plot2, contour_plot, ncol = 2)

  # Print the arranged plots
  print(arranged_plots)

  current_4panel$p1 <- histogram_plot
  current_4panel$p2 <- overlay_plot2
  current_4panel$p3 <- digrepresentation
  current_4panel$p4 <- contour_plot

    output$downloadcurrent = downloadHandler(
    filename = function() {"plots.pdf"},
    content = function(file) {
     pdf(file, width = 11, height = 11, onefile = TRUE)
     grid.arrange(current_4panel$p1, current_4panel$p2, current_4panel$p3, current_4panel$p4, ncol = 2) 
     dev.off()
    })

  runjs("document.getElementById('message_gen_hist').style.display = 'none';")
   shinyjs::show("gated_histogram_on_page")

}


  output$downloadall <- downloadHandler(
   filename = "all_plots.pdf", 
    content = function(file) {
      selected_columns_save <- input$selected_columns
      selected_patients_save <- input$unique_patients_gating
    # Subsetting the data based on selected columns and unique image IDs
    sub_data <- uploaded_df() %>%
      filter(imageid %in% selected_patients_save) %>%
      select(imageid, all_of(selected_columns_save))


    print(colnames(sub_data))

    pdf_file_name <- "all_plots.pdf"  # Get user-provided file name

     runjs("
    function animateDots() {
      var messageDiv = document.getElementById('message_gen_hist');
      var dots = '';
      setInterval(function() {
        if (dots.length === 3) dots = '';
        else dots += '.';
        messageDiv.innerText = 'Generating Histograms' + dots;
      }, 500);
    }
    animateDots();
    document.getElementById('message_gen_hist').style.display = 'block';
  ")

    generate_histogram_pdf(sub_data, pdf_file_name)

    runjs("document.getElementById('message_gen_hist').style.display = 'none';")

    
        file.copy(pdf_file_path(), file)
      }
    
)


generate_four_panel_plot <- function(data, imageid, marker) {
  # Your code to generate the four-panel plot here for a specific imageid and marker
  # Replace this code with your existing code that generates the four-panel plot

  # # Get the column indices for the selected columns
  #       column_indices <- sapply(selected_columns, function(col_name) {
  #         which(names(filtered_data) == col_name)
  #       })

  
  selected_columns_save <- input$selected_columns

  print(selected_columns_save)

  # Get the column indices for the selected columns
  column_indices <- sapply(selected_columns_save, function(col_name) {
    which(names(data) == col_name)
  })

  
  histogram_plot <- plot_gating_grid(data, imageid, column_indices[marker])

  filtered_data <- uploaded_df() %>%
    filter(imageid %in% imageid) %>%
    select(imageid, all_of(marker))
  
  filtered_data_xy <- uploaded_df() %>%
    filter(imageid %in% imageid) 

  chosen_patient <- imageid
  chosen_marker <- marker

  digrepresentation <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid, color = filtered_data_xy[[chosen_marker]])) +
    geom_point(shape = 20, size = 0.3) +
    scale_color_gradient(low = "grey30", high = "white") +
    theme(panel.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("X Centroid") +
    ylab("Y Centroid") +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14)) +
    labs(title = 'Digital Representation') 
  
  gate_value <- resultdf_reactive()$Gate[
    resultdf_reactive()$Marker == chosen_marker &
      resultdf_reactive()$Patient == chosen_patient
  ]
  
  if (!is.null(outputinterceptreactive())) {
    gate_value <- outputinterceptreactive()
  }
  
  temp_df <- filtered_data_xy
  temp_df[[chosen_marker]][temp_df[[chosen_marker]] < gate_value] <- 0
  
  density_values <- get_density(filtered_data_xy$X_centroid, filtered_data_xy$Y_centroid, n = 100)
  
  contour_plot <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
    geom_point(aes(color = ifelse(filtered_data_xy[[chosen_marker]] > gate_value, density_values, 0)), size = 0.3) +
    geom_density_2d(data = subset(filtered_data_xy, filtered_data_xy[[chosen_marker]] > gate_value), color = 'black') +  
    scale_color_viridis(option = "turbo") + 
    theme(legend.position = "Bottom",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14)) +
    labs(title = 'Positive Density', x = "X Centroid", y = "Y Centroid")
  
  overlay_plot2 <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
    geom_point(aes(color = ifelse(filtered_data_xy[[chosen_marker]] > gate_value, "Positive", "Negative")), size = 0.3) +
    scale_color_manual(guide = guide_legend(title = ""), values = c("Positive" = "red", "Negative" = "grey")) +
    theme(legend.position = "Bottom",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14)) +
    labs(title = 'Positive Cells', x = "X Centroid", y = "Y Centroid")
  
  arranged_plots <- grid.arrange(histogram_plot, overlay_plot2, digrepresentation, contour_plot, ncol = 2)
  
  return(arranged_plots)
}




# # Create a reactive value to store histograms and gate values

# # Function to generate the PDF with histograms
# generate_histogram_pdf <- function(subdata_to_plot, pdf_file_name) {

#   data <- subdata_to_plot

#   print('colnames of subdata hist here')

#   print(colnames(data))
#   # List of markers
#   markers <- colnames(data)[-1]
#   print(markers)

#   # Create a list to store ggplot objects
#   histogram_list <- list()

#   # Loop over unique imageids
#   unique_imageids <- unique(data$imageid)
#   print(unique_imageids)
#   for (imageid in unique_imageids) {
#     print(imageid)
#     for (marker in markers) {
#       # Subset data for the current imageid and marker
#       # subset_data <- data[data$imageid == imageid, c("imageid", marker)]
#       print(marker)


#     # Get the index of the column
#     col_index <- which(names(data) == marker)

      
#     #     unique_id <- unique(df$imageid)[current_patient()]
#     #     column_index <- column_indices[current_marker()]
#     # histogram_plot <- plot_gating_grid(df, unique_id, column_index)

#       histogram <- plot_gating_grid(data, imageid, col_index)
#       # Store the histogram plot in the list
#       histogram_list[[paste(imageid, marker)]] <- histogram

#     }

      
#       # # Store the histogram plot in the list
#       # histogram_list[[paste(imageid, marker)]] <- histogram
#     }
  

#   # Arrange and save histograms across multiple pages
#   pdf_file <- pdf_file_name
#   pdf(file = pdf_file, width = 14, height = 10)
#   histogram_count <- length(histogram_list)
#   histograms_per_page <- 6
#   pages <- ceiling(histogram_count / histograms_per_page)

#   for (page in 1:pages) {
#     start_idx <- (page - 1) * histograms_per_page + 1
#     end_idx <- min(page * histograms_per_page, histogram_count)
#     page_histograms <- histogram_list[start_idx:end_idx]
    
#     do.call(grid.arrange, c(page_histograms, ncol = 3))
#     if (page < pages) {
#       cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")  # Add page breaks
#     }
#   }
#   dev.off()
#   pdf_file_path(pdf_file)

#   }

# generate_histogram_pdf <- function(subdata_to_plot, pdf_file_name) {
#   data <- subdata_to_plot
  
#   # Create a list to store ggplot objects (four-panel plots)
#   four_panel_plot_list <- list()
  
#   # Loop over unique imageids and markers
#   unique_imageids <- unique(data$imageid)
#   markers <- colnames(data)[-1]
  
#   for (imageid in unique_imageids) {
#     for (marker in markers) {
#       # Generate the four-panel plot for each combination of imageid and marker
#       four_panel_plot <- generate_four_panel_plot(data, imageid, marker)
      
#       # Store the four-panel plot in the list
#       four_panel_plot_list[[paste(imageid, marker)]] <- four_panel_plot
#     }
#   }
  
#   # Arrange and save four-panel plots across multiple pages in the PDF
#   pdf_file <- pdf_file_name
#   pdf(file = pdf_file, width = 14, height = 10)
  
#   # Number of four-panel plots and plots per page
#   plot_count <- length(four_panel_plot_list)
#   plots_per_page <- 6
#   pages <- ceiling(plot_count / plots_per_page)
  
#   for (page in 1:pages) {
#     start_idx <- (page - 1) * plots_per_page + 1
#     end_idx <- min(page * plots_per_page, plot_count)
#     page_plots <- four_panel_plot_list[start_idx:end_idx]
    
#     # Modify the arrangement logic here to display the plots properly in the PDF
#     # Example: do.call(grid.arrange, c(page_plots, ncol = 3))
    
#     if (page < pages) {
#       cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")  # Add page breaks
#     }
#   }
  
#   dev.off()
#   pdf_file_path(pdf_file)
# }

# generate_histogram_pdf <- function(subdata_to_plot, pdf_file_name) {
#   data <- subdata_to_plot
  
#   # Create a list to store ggplot objects (four-panel plots)
#   four_panel_plot_list <- list()
  
#   # Loop over unique imageids and markers
#   unique_imageids <- unique(data$imageid)
#   markers <- colnames(data)[-1]
  
#   for (imageid in unique_imageids) {
#     for (marker in markers) {
#       # Generate the four-panel plot for each combination of imageid and marker
#       four_panel_plot <- generate_four_panel_plot(data, imageid, marker)
      
#       # Store the four-panel plot in the list
#       four_panel_plot_list[[paste(imageid, marker)]] <- four_panel_plot
#     }
#   }
  
#   # Arrange and save four-panel plots across multiple pages in the PDF
#   pdf_file <- pdf_file_name
#   pdf(file = pdf_file, width = 14, height = 10)
  
#   # Number of four-panel plots and plots per page
#   plot_count <- length(four_panel_plot_list)
#   plots_per_page <- 4  # Set to 4 for four plots per page
#   pages <- ceiling(plot_count / plots_per_page)
  
#   for (page in 1:pages) {
#     start_idx <- (page - 1) * plots_per_page + 1
#     end_idx <- min(page * plots_per_page, plot_count)
#     page_plots <- four_panel_plot_list[start_idx:end_idx]
    
#     # Arrange plots in a 2x2 grid on each page
#     do.call(grid.arrange, c(page_plots))
#     # do.call(grid.arrange, c(page_plots, nrow = 2, ncol = 2))

    
#     if (page < pages) {
#       cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")  # Add page breaks
#     }
#   }
  
#   dev.off()
#   pdf_file_path(pdf_file)
# }

# generate_histogram_pdf <- function(subdata_to_plot, pdf_file_name) {
#   data <- subdata_to_plot
  
#   # Create a list to store ggplot objects (four-panel plots)
#   four_panel_plot_list <- list()
  
#   # Loop over unique imageids and markers
#   unique_imageids <- unique(data$imageid)
#   markers <- colnames(data)[-1]
  
#   for (imageid in unique_imageids) {
#     for (marker in markers) {
#       # Generate the four-panel plot for each combination of imageid and marker
#       four_panel_plot <- generate_four_panel_plot(data, imageid, marker)
      
#       # Store the four-panel plot in the list
#       four_panel_plot_list[[paste(imageid, marker)]] <- four_panel_plot
#     }
#   }
  
#   # Arrange and save four-panel plots across multiple pages in the PDF
#   pdf_file <- pdf_file_name
#   pdf(file = pdf_file, width = 14, height = 10)
  
#   # Number of four-panel plots
#   plot_count <- length(four_panel_plot_list)
  
#   for (page in 1:plot_count) {
#     page_plot <- four_panel_plot_list[[page]]
    
#     # Display one plot per page
#     print(page_plot)
    
#     if (page < plot_count) {
#       cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")  # Add page breaks
#     }
#   }
  
#   dev.off()
#   pdf_file_path(pdf_file)
# }

generate_histogram_pdf <- function(subdata_to_plot, pdf_file_name) {
  data <- subdata_to_plot
  
  # Create a list to store ggplot objects (four-panel plots)
  four_panel_plot_list <- list()
  
  # Loop over unique imageids and markers
  unique_imageids <- unique(data$imageid)
  markers <- colnames(data)[-1]
  
  for (imageid in unique_imageids) {
    for (marker in markers) {
      # Generate the four-panel plot for each combination of imageid and marker
      four_panel_plot <- generate_four_panel_plot(data, imageid, marker)
      
      # Store the four-panel plot in the list
      four_panel_plot_list[[paste(imageid, marker)]] <- four_panel_plot
    }
  }
  
  # Arrange and save four-panel plots across multiple pages in the PDF
  pdf_file <- pdf_file_name
  pdf(file = pdf_file, width = 14, height = 10)
  
  # Number of four-panel plots
  plot_count <- length(four_panel_plot_list)
  
  for (page in 1:plot_count) {
    page_plot <- four_panel_plot_list[[page]]
    
    # Use do.call to place one plot per page
    do.call(grid.arrange, list(page_plot, ncol = 1))
    
    if (page < plot_count) {
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")  # Add page breaks
    }
  }
  
  dev.off()
  pdf_file_path(pdf_file)
}




# library(gridExtra)

# generatePlot <- function() {
#   # Add your previous code here to setup and prepare the required data
#   # ...

#   # Initialize a list to store the plots for all combinations of columns and patients
#   all_plots <- list()

#   for (i in seq_along(input$selected_columns)) {
#     for (j in seq_along(input$unique_patients_gating)) {
#       # Your existing logic for generating individual plots for each selected column and patient here...
#       # ... (Your existing logic for generating plots here)
      
#       # Store each plot in the all_plots list
#       all_plots[[length(all_plots) + 1]] <- grid.arrange(
#         histogram_plot, overlay_plot2, digrepresentation, contour_plot,
#         ncol = 2, main = paste("Plots for", input$selected_columns[i], "and", input$unique_patients_gating[j])
#       )
#     }
#   }

#   # Combine all plots into a single grid
#   combined_plots <- do.call(grid.arrange, c(all_plots, ncol = 2))

#   # Save the combined plots into a PDF file
#   output$downloadcurrent = downloadHandler(
#     filename = function() { "combined_plots.pdf" },
#     content = function(file) { 
#       pdf(file)
#       print(combined_plots)
#       dev.off()
#     }
#   )
# }


# output$icaAnalysis <- renderPlot({
#   # Assuming 'uploaded_df()' is your dataframe

#   # Select numeric columns for ICA
#   numeric_data <- uploaded_df()[sapply(uploaded_df(), is.numeric)]

#   # Apply FastICA
#   ica_result <- fastICA(numeric_data, n.comp = 5)  # Adjust the number of components as needed

#   # Extract the independent components
#   independent_components <- ica_result$S

# #   # Generate a plot for the independent components
# #   plot_ica <- plot(independent_components, type = "l")  # Adjust the plot type as per your preference

# #   plot_ica
# # })


# output$icaAnalysis <- renderPlot({
#   # Assuming 'uploaded_df()' is your dataframe

#   # Select numeric columns for ICA
#   numeric_data <- uploaded_df()[sapply(uploaded_df(), is.numeric)]

#   # # Replace infinite values with NA
#   # numeric_data[is.infinite(numeric_data)] <- NA

#   # Replace missing values with column means
#   # numeric_data[is.na(numeric_data)] <- colMeans(numeric_data, na.rm = TRUE)

#   numeric_data <- numeric_data[!is.infinite(numeric_data)]


#   # Apply FastICA
#   ica_result <- fastICA(numeric_data, n.comp = 5)  # Adjust the number of components as needed

#   # Extract the independent components
#   independent_components <- ica_result$S

#   # Generate a plot for the independent components
#   plot_ica <- plot(independent_components, type = "l")  # Adjust the plot type as per your preference

#   plot_ica
# })

# output$icaAnalysis <- renderPlot({
#   # Assuming 'uploaded_df()' is your dataframe

#   # Select numeric columns for ICA
#   numeric_data <- uploaded_df()[sapply(uploaded_df(), is.numeric)]

#   # Replace infinite values with NA
#   numeric_data[is.infinite(numeric_data)] <- NA

#   # Replace missing values with column means
#   numeric_data[is.na(numeric_data)] <- colMeans(numeric_data, na.rm = TRUE)

#   # Convert numeric data to a matrix
#   numeric_matrix <- as.matrix(numeric_data)

#   # Apply FastICA
#   ica_result <- fastICA(numeric_matrix, n.comp = 5)  # Adjust the number of components as needed

#   # Extract the independent components
#   independent_components <- ica_result$S

#   # Generate a plot for the independent components
#   plot_ica <- plot(independent_components, type = "l")  # Adjust the plot type as per your preference

#   plot_ica
# })

# Load TIFF file (replace "your_file.tiff" with your TIFF file path)
#   tiff_data <- readTIFF("wsi_t.tif")

  # # Display TIFF image in the plotOutput
  # output$icaAnalysis <- renderPlot({
  #     # Perform UMAP dimensionality reduction
  #     umap_result <- umap(iris[, -5])  # Exclude the species column for UMAP
  #     print(umap_result)
  #     print(head(iris))

  #     # Add UMAP results back to the dataset
  #     iris_with_umap <- cbind(iris, umap_result$layout)

  #     print(iris_with_umap)

  #     # Plotting UMAP results
  #     ggplot(iris_with_umap, aes(x = iris_with_umap[['1']], y = iris_with_umap[['2']], color = Species)) +
  #       geom_point() +
  #       labs(title = "UMAP Plot of Iris Dataset")
  # })


  # # # Display TIFF image in the plotOutput
  # output$icaAnalysis <- renderPlot({

  #   req(uploaded_df())
  #   colnames(uploaded_df())
  #   # Define the column names you want to exclude
  #   columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", "DNA1", "DNA3", "DNA2", "Hoechst"
  #   , "DAPI1", "DAPI2", "DNA13", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
  #   "MinorAxisLength", "Extent", "Solidity", "Orientation", "DNA_6", "DNA_7", "DNA_8", "DNA_9", "")  # List the columns to exclude
  #   # Select all column names except the ones to exclude
  #   # Define the regular expression pattern to identify columns to be excluded
  #   exclude_pattern <- "DNA|DAPI"

  #   # Use grep to get the column names matching the pattern
  #   exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

  #   # Add the excluded columns to the original 'columns_to_exclude' vector
  #   columns_to_exclude <- c(columns_to_exclude, exclude_columns)

  #   column_names <- setdiff(names(uploaded_df()), columns_to_exclude)
  #   print(column_names)
  #     # Perform UMAP dimensionality reduction
  #     # Perform UMAP dimensionality reduction
  #     # Drop rows with any NA values in the selected columns
  #     your_matrix <- data.matrix(uploaded_df()[, 'ELANE'])

  #       # Perform UMAP dimensionality reduction on the cleaned dataset
  #       umap_result <- umap(your_matrix)


  #     print('hi')
      
  #     # Combine UMAP results with original dataframe
  #     umap_df <- cbind(uploaded_df(), umap_result$layout)

  #     # Plotting UMAP results
  #     ggplot(umap_df, aes(x = umap_df[['1']], y = umap_df[['2']], color = factor(1:nrow(umap_df)))) +
  #       geom_point() +
  #       scale_color_discrete(guide = FALSE) +
  #       labs(title = "UMAP Plot of Selected Columns")
  # })


  # Display TIFF image in the plotOutput
  output$icaAnalysis <- renderPlot({

    req(uploaded_df())
    colnames(uploaded_df())
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", "DNA2", "Hoechst1", "Hoechst2", "Hoechst3", 
                            "Hoechst4", "Hoechst5", "Hoechst6", "Hoechst7", "Hoechst8", "Hoechst9", "Hoechst10", "Hoechst_1", "Hoechst_2", "Hoechst_3", 
                            "Hoechst_4", "Hoechst_5", "Hoechst_6", "Hoechst_7", "Hoechst_8", "Hoechst_9", "Hoechst_10", "DAPI1", "DAPI2", "DAPI3", "DAPI4", "DAPI5", "DAPI6", "DAPI7", "DAPI8", "DAPI9", "DAPI10",
                            "DAPI_1", "DAPI_2", "DAPI_3", "DAPI_4", "DAPI_5", "DAPI_6", "DAPI_7", "DAPI_8", "DAPI_9", "DAPI_10", "DNA1", "DNA3", "DNA2", "DNA4", "DNA5", "DNA6","DNA7","DNA8", "DNA9", "DNA10", "DNA11",
                            "DNA12", "DNA13", "DNA_1", "DNA_3", "DNA_2", "DNA_4", "DNA_5", "DNA_6","DNA_7","DNA_8", "DNA_9", "DNA_10", "DNA_11",
                            "DNA_12", "DNA_13", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
                            "MinorAxisLength", "Extent", "Solidity", "Orientation", "") # List the columns to exclude
    # Select all column names except the ones to exclude
    # Define the regular expression pattern to identify columns to be excluded
    exclude_pattern <- "DNA|DAPI|Hoechst"

    # Use grep to get the column names matching the pattern
    exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

    # Add the excluded columns to the original 'columns_to_exclude' vector
    columns_to_exclude <- c(columns_to_exclude, exclude_columns)

    column_names <- setdiff(names(uploaded_df()), columns_to_exclude)
    print(column_names)
      # Perform UMAP dimensionality reduction
      # Perform UMAP dimensionality reduction
      # Drop rows with any NA values in the selected columns

      df_to_umapify <- uploaded_df()
      df_to_umapify <- df_to_umapify[, column_names]
      # Remove rows with any infinite values across all columns
      # Loop through the list and replace infinite values in each vector with 0
        for (i in seq_along(df_to_umapify)) {
          df_to_umapify[[i]][is.infinite(df_to_umapify[[i]])] <- 0
        }

      # df_to_umapify[!is.infinite(df_to_umapify)] <- 0
     

      print(colnames(df_to_umapify))

      your_matrix <- data.matrix(df_to_umapify[, column_names])

        # Perform UMAP dimensionality reduction on the cleaned dataset
        umap_result <- umap(your_matrix)


      print('hi')
      
      # Combine UMAP results with original dataframe
      umap_df <- cbind(uploaded_df(), umap_result$layout)

      # Plotting UMAP results
      ggplot(umap_df, aes(x = umap_df[['1']], y = umap_df[['2']], color = factor(1:nrow(umap_df)))) +
        geom_point() +
        scale_color_discrete(guide = FALSE) +
        labs(title = "UMAP Plot of Selected Columns")
  })

  #   # Display TIFF image in the plotOutput
  # output$icaAnalysis <- renderPlot({
  #   req(uploaded_df())
  #   colnames(uploaded_df())
  #   # Define the column names you want to exclude
  #   columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", "DNA1", "DNA3", "DNA2", "Hoechst"
  #   , "DAPI1", "DAPI2", "DNA13", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
  #   "MinorAxisLength", "Extent", "Solidity", "Orientation", "DNA_6", "DNA_7", "DNA_8", "DNA_9", "")  # List the columns to exclude
  #   # Select all column names except the ones to exclude
  #   # Define the regular expression pattern to identify columns to be excluded
  #   exclude_pattern <- "DNA|DAPI"

  #   # Use grep to get the column names matching the pattern
  #   exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

  #   # Add the excluded columns to the original 'columns_to_exclude' vector
  #   columns_to_exclude <- c(columns_to_exclude, exclude_columns)

  #   column_names <- setdiff(names(uploaded_df()), columns_to_exclude)
  #   print(column_names)
  #     # Perform UMAP dimensionality reduction
  #     # Perform UMAP dimensionality reduction
  #     # Drop rows with any NA values in the selected columns

  #     df_to_umapify <- uploaded_df()
  #     df_to_umapify <- df_to_umapify[, column_names]
  #     # Remove rows with any infinite values across all columns
  #     # Loop through the list and replace infinite values in each vector with 0
  #       for (i in seq_along(df_to_umapify)) {
  #         df_to_umapify[[i]][is.infinite(df_to_umapify[[i]])] <- 0
  #       }

  #     # df_to_umapify[!is.infinite(df_to_umapify)] <- 0
     

  #     print(colnames(df_to_umapify))

  #     your_matrix <- data.matrix(df_to_umapify[, column_names])

  #       # # Perform UMAP dimensionality reduction on the cleaned dataset
  #       # umap_result <- umap(your_matrix)

  #     print('hi')
      
      

  #        # Perform K-means clustering with 3 clusters
  #   kmeans_fit <- kmeans(your_matrix, centers = 3)

  #   # Add cluster labels to the UMAP dataframe
  #   umap_df$cluster <- kmeans_fit$cluster

  #   # Combine UMAP results with original dataframe
  #         umap_df <- cbind(uploaded_df(), umap_df$cluster)


  #             ggplot(umap_df, aes(x = umap_df[['1']], y = umap_df[['2']], color = factor(cluster))) +
  #         geom_point() +
  #         scale_color_discrete(palette = c("red", "green", "blue"), # Customize color palette
  #                             labels = c("Cluster 1", "Cluster 2", "Cluster 3"), # Specify custom labels
  #                             guide = guide_legend(title = "Pseudo Clusters")) +
  #         labs(title = "UMAP Plot of Selected Columns")
  # })






  unique_patients <- reactive({
  req(uploaded_df())
  unique(uploaded_df()$imageid)
  })

#    observe({
#   updateNumericInput(session, "intercept", value = '')
# })

  observe({
  shinyjs::show("patient_number")
  updateRadioButtons(session, "patient_number", choices = unique_patients())
})

 observe({
  #shinyjs::show("patient_number")
  updateRadioButtons(session, "patient_number_im_gargage", choices = unique_patients())
})

 observe({
  shinyjs::show("patients_ts")
  updateRadioButtons(session, "patients_ts", choices = unique_patients())
})


 observe({
  shinyjs::show("cycle_ts")
  updateRadioButtons(session, "cycle_ts", choices = cycle_detected())
})

 observe({
  shinyjs::show("patient_number_plot2")
  updateRadioButtons(session, "patient_number_plot2", choices = unique_patients())
})

observe({
  shinyjs::show("trivar_patient")
  updateRadioButtons(session, "trivar_patient", choices = unique_patients())
})


 observe({
  shinyjs::show("unique_patients_gating")
  updateCheckboxGroupInput(session, "unique_patients_gating", choices = unique_patients(), selected = unique_patients())
})

 observe({
  shinyjs::show("precrev_marker")
  updateRadioButtons(session, "precrev_marker", choices = unique_markers_w_DNA())
})

  unique_markers_w_DNA <- reactive({
    req(uploaded_df())
    colnames(uploaded_df())
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", 
     "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
    "MinorAxisLength", "Extent", "Solidity", "Orientation", "")  # List the columns to exclude
    # Select all column names except the ones to exclude
    # Define the regular expression pattern to identify columns to be excluded
    # exclude_pattern <- "DNA|DAPI"

    # # Use grep to get the column names matching the pattern
    # exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

    # Add the excluded columns to the original 'columns_to_exclude' vector
    # columns_to_exclude <- c(columns_to_exclude, exclude_columns)
      # Define the regular expression pattern to identify columns to be excluded
      exclude_pattern2 <- "_positivity"

      # Use grep to get the column names matching the pattern
      exclude_columns2 <- grep(exclude_pattern2, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

      columns_to_exclude <- c(columns_to_exclude, exclude_columns2)

    column_names <- setdiff(names(uploaded_df()), columns_to_exclude)
    
  })


  unique_markers <- reactive({
    req(uploaded_df())
    colnames(uploaded_df())
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y", "DNA2", "Hoechst1", "Hoechst2", "Hoechst3", 
                            "Hoechst4", "Hoechst5", "Hoechst6", "Hoechst7", "Hoechst8", "Hoechst9", "Hoechst10", "Hoechst_1", "Hoechst_2", "Hoechst_3", 
                            "Hoechst_4", "Hoechst_5", "Hoechst_6", "Hoechst_7", "Hoechst_8", "Hoechst_9", "Hoechst_10", "DAPI1", "DAPI2", "DAPI3", "DAPI4", "DAPI5", "DAPI6", "DAPI7", "DAPI8", "DAPI9", "DAPI10",
                            "DAPI_1", "DAPI_2", "DAPI_3", "DAPI_4", "DAPI_5", "DAPI_6", "DAPI_7", "DAPI_8", "DAPI_9", "DAPI_10", "DNA1", "DNA3", "DNA2", "DNA4", "DNA5", "DNA6","DNA7","DNA8", "DNA9", "DNA10", "DNA11",
                            "DNA12", "DNA13", "DNA_1", "DNA_3", "DNA_2", "DNA_4", "DNA_5", "DNA_6","DNA_7","DNA_8", "DNA_9", "DNA_10", "DNA_11",
                            "DNA_12", "DNA_13", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
    "MinorAxisLength", "Extent", "Solidity", "Orientation", "")  # List the columns to exclude
    # Select all column names except the ones to exclude
    # Define the regular expression pattern to identify columns to be excluded
    exclude_pattern <- "DNA|DAPI|Hoechst"

    # Use grep to get the column names matching the pattern
    exclude_columns <- grep(exclude_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

    # Add the excluded columns to the original 'columns_to_exclude' vector
    columns_to_exclude <- c(columns_to_exclude, exclude_columns)

    column_names <- setdiff(names(uploaded_df()), columns_to_exclude)
    
  })

    cycle_detected <- reactive({
    req(uploaded_df())
    colnames(uploaded_df())
    cycle_pattern <- "DNA|DAPI|Hoechst"  # Define your pattern


    column_names <- grep(cycle_pattern, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)

    column_names
    
  })


 selected_columns <- reactive({
    input$selected_columns
  })

  selected_columns_phenotyping <- reactive({
    input$selected_columns_phenotyping
  })


  unique_patients_gating <- reactive({
    input$unique_patients_gating
  })

  inputreactiveintercept_USERDEF <- reactive({
    input$intercept
  })


 selected_opacity <- reactive({
    input$opacity_slider
  })


   

  # inputreactiveintercept <- reactive({
  #   input
  # })

  observe({
    updateSelectInput(session, "marker", choices = unique_markers())
})

 observe({
    updateSelectInput(session, "marker_im_garage", choices = unique_markers())
})



   observe({
    updateSelectInput(session, "xvar", choices = unique_markers())
})


 observe({
    updateSelectInput(session, "yvar", choices = unique_markers())
})


   observe({
    updateSelectInput(session, "xvarTri", choices = unique_markers())
})


 observe({
    updateSelectInput(session, "yvarTri", choices = unique_markers())
})

 observe({
    updateSelectInput(session, "zvar", choices = unique_markers())
})



  
  # newdf <- reactiveVal(NULL)

  # newdf <- reactiveVal(NULL)

  brushed_data <- reactiveVal(NULL)
  filtered_data_original <- reactiveVal(NULL)
  brushed_data_removal <- reactiveVal(NULL)

    filtered_data_reactive <- reactiveVal(NULL)


    observeEvent(c(input$precrev_marker, input$patient_number_im_gargage), {
  output$image_garage_output <- renderPlotly({
    req(uploaded_df())
    req(input$patient_number_im_gargage)
    req(input$precrev_marker)
    req(input$opacity_slider)

      selected_patient <- input$patient_number_im_gargage
      selected_marker <- input$precrev_marker
      
      # Filter data based on selected patient
      filtered_data <- uploaded_df()[uploaded_df()$imageid %in% selected_patient, ]


      imageid <- filtered_data$imageid  # Extract the desired column
      filtered_data <- filtered_data[, -which(names(filtered_data) == "imageid")]  # Remove the desired column from its original position
      filtered_data <- cbind(imageid, filtered_data)  # Add the desired column as the first column

      filtered_data_original(filtered_data)  # Store the original data


      filtered_data_reactive(filtered_data)

    ggplot(filtered_data_reactive(), aes(x = X_centroid, y = Y_centroid, color = filtered_data_reactive()[[selected_marker]])) +
      geom_point(alpha = selected_opacity()) +
      xlab('X Centroid') + ylab('Y Centroid') + 
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18)) +
      labs(title = "Digital Marker Overlay") +
      scale_color_continuous(name = "Intensity (logged)")
  })
})




# Function to filter columns based on a pattern
filter_columns_by_pattern <- function(data, pattern) {
  selected_columns <- grep(pattern, colnames(data), value = TRUE)
  return(data[, selected_columns])
}

# crudeway_react <- reactiveVal(NULL)
regression_mode_react <- reactiveVal(NULL)
sugremReact <- reactiveVal(NULL)
klPlotReact <- reactiveVal(NULL)
filtered_data_updated <- reactiveVal(NULL)

    output$tissue_score_out <- renderPlot({
     req(uploaded_df())
     req(input$patients_ts)
     req(input$cycle_ts)


      # Get selected marker and patient
      selected_patient <- input$patients_ts

      selected_cycle <- input$cycle_ts


      filtered_data <- subsetted_precrev_ts()

      # Assuming 'df' is your dataframe and 'column_name' is the name of the column to drop
      #filtered_data <- filtered_data[, !colnames(filtered_data) %in% c("DNA", "DAPI", "Hoechst")]


      print(selected_patient)

  
      cycle_pattern <- "DNA|DAPI|Hoechst"  # Define your pattern

      print(colnames(filtered_data))


        selected_columns <- grep(cycle_pattern, colnames(filtered_data), value = TRUE, ignore.case = TRUE)

         # Get the column indices for the selected columns
        column_indices <- sapply(selected_columns, function(col_name) {
          which(names(filtered_data) == col_name)
        })

      print(selected_columns)

  
     print(column_indices[selected_columns])
     print(colnames(uploaded_df()))
      histogram_plots <- plot_gating_grid(filtered_data, selected_patient, column_indices[selected_columns])



      # output$histogram_plots_cyles 

      #  output$histogram_plots_cyles <- renderPlot({
      #       histogram_plots
      #     })

          # Define a function to find the mode of a vector
      find_mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }

        # Extract histogram data and find modes
        modes <- list()
        kldivs <- list()
        # crudeway <- list()
        # cut_threshold_list <- list()

        for (column_name in selected_columns) {
          histogram_data <- filtered_data[[column_name]]  # Assuming the histogram is the first grob

          histogram_data <- remove_outliers2(histogram_data)
          histogram_data <- histogram_data[!is.infinite(histogram_data)]


          # Find mode of the histogram data
          mode_value <- find_mode(histogram_data)

        

          print(mode_value)

          # print(kldiv)
          
          # Store mode values
          modes[[column_name]] <- mode_value
        }


        # Combine modes into one data frame
        modes_df <- data.frame(Index = seq_along(modes), Value = unlist(modes))

        print(modes_df)

        print("DNA 1 value here")
        print(modes_df[1,2])


         # Plot the modes separately
        plot_modes <- ggplot(modes_df, aes(x = Index, y = Value)) +
        geom_point() + 
          geom_line() +
          labs(x = "Cycle Number", y = "Mode") +
          ggtitle("Regression of Quality") +
                theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 18), legend.position = "none"
          ) +  scale_y_continuous(label = scales::comma)

          regression_mode_react(plot_modes)


        print(colnames(filtered_data_updated()))

        print(unique(filtered_data_updated()$Removed_FOXP3_DNA_8))

            # Calculate the number of columns
      num_columns <- length(selected_columns)

      # Initialize a list to store KL divergences
      kl_divergences <- list()

      # Define the range for iteration based on odd or even number of columns
      iter_range <- ifelse(num_columns %% 2 == 0, num_columns - 2, num_columns - 1)

      # Iterate through consecutive pairs of columns
      for (i in seq(1, iter_range, by = 1)) {
        # Extract consecutive pairs of columns
        intensity_vector1 <- filtered_data[[selected_columns[i]]]
        intensity_vector2 <- filtered_data[[selected_columns[i + 1]]]
        
        # Exponentiate log values to obtain non-logarithmic values
        non_log_values1 <- exp(intensity_vector1)
        non_log_values2 <- exp(intensity_vector2)
        
        # Normalize the resulting values to ensure they represent probability distributions
        normalize <- function(x) x / sum(x)
        prob_distribution1 <- normalize(non_log_values1)
        prob_distribution2 <- normalize(non_log_values2)
        
        # Bind the vectors
        vectors <- rbind(prob_distribution1, prob_distribution2)
        
        # Calculate KL divergence between histograms
        kl_divergence <- KL(vectors, unit = 'log')

        print(kl_divergence)
        
        # Store the KL divergence in the list
        kl_divergences[[paste(selected_columns[i], selected_columns[i + 1], sep = "_")]] <- kl_divergence
      }

      # Convert kl_divergences to a data frame
      kl_df <- data.frame(
        Index = seq_along(kl_divergences),
        KL_Divergence = unlist(kl_divergences)
      )

      # Plot the trend of KL divergences
      plot_kl_divergences <- ggplot(kl_df, aes(x = Index, y = KL_Divergence)) +
        geom_point() + 
        geom_line() +
        labs(x = "Pair Number", y = "KL Divergence") +
        ggtitle("Trend of KL Divergence") +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          legend.position = "none"
        ) +
        scale_y_continuous(label = scales::comma)

      klPlotReact(plot_kl_divergences)       


            # Example: Calculate mode differences
        mode_diffs <- diff(modes_df$Value)

        print("mode diffs here")

        print(mode_diffs)


        # # Normalize mode differences
        # ideal_mode_diff <- 0  # Ideal mode difference
        # worst_mode_diff <- max(abs(mode_diffs))  # Worst-case mode difference

        # # Assign scores for mode differences
        # ideal_score_mode_diff <- 100
        # worst_score_mode_diff <- 0

        # # Normalize mode differences based on ideal and worst-case scenarios
        # normalize_mode_diff <- function(value, ideal, worst) {
        #   (value - worst) / (ideal - worst) * 100
        # # }

        # # Normalize mode differences
        # # normalized_mode_diffs <- normalize_mode_diff(mode_diffs, ideal_mode_diff, worst_mode_diff)

        # # Calculate overall quality score (average of KL divergence and mode difference scores)
        # overall_quality_score <- mean(normalized_mode_diffs)

        # # Determine gauge chart bounds based on overall quality score
        # lower_bound <- min(overall_quality_score, 0)  # Lower bound
        # upper_bound <- max(overall_quality_score, 100)  # Upper bound


            # Calculate an overall score using KL divergence and mode values
      # Calculate an overall score using KL divergence and mode values
# # You might want to standardize these values to have comparable ranges
#           kl_divergence_score <- mean(kl_df$KL_Divergence)  # Mean of KL divergences
#           mode_score <- max(modes_df$Value) - min(modes_df$Value)  # Range of mode values

#           # Combine scores to create an overall quality metric
#           quality_score <- (kl_divergence_score + mode_score) / 2  # You can use a weighted sum based on importance


          # # Create a gauge plot using Plotly
          # gauge_plot <- plot_ly(
          #   type = "indicator",
          #   mode = "gauge+number",
          #   value = overall_quality_score,
          #   title = "Quality Gauge",
          #   gauge = list(
          #     axis = list(range = list(0, upper_bound)),
          #     # steps = list(
          #     #   list(range = c(0, 2), color = "red"),
          #     #   list(range = c(2, 4), color = "orange"),
          #     #   list(range = c(4, 6), color = "yellow"),
          #     #   list(range = c(6, 8), color = "green"),
          #     #   list(range = c(8, 10), color = "darkgreen")
          #     # ),
          #     steps = list(
          #       list(range = c(0, 0.2*upper_bound), color = "red"),
          #       list(range = c(0.2*upper_bound, 0.4*upper_bound), color = "orange"),
          #       list(range = c(0.4*upper_bound, 0.6*upper_bound), color = "yellow"),
          #       list(range = c(0.6*upper_bound, 0.8*upper_bound), color = "green"),
          #       list(range = c(0.8*upper_bound, upper_bound), color = "darkgreen")
          #     )
          #     # threshold = list(
          #     #   line = list(color = "black", width = 4),
          #     #   thickness = 0.75,
          #     #   value = 5
          #     # )
          #   )
          # )

          # library(plotly)

          # Replace 0.5 with the actual mean of your KL divergences
          mean_KL_divergences <- mean(kl_df$KL_Divergence)


          # library(plotly)

          # Replace 0.5 with the actual mean of your KL divergences
          # mean_KL_divergences <- 0.5

          # Subtracting the value from 1 to invert the scale and normalize between 0 and 1
          normalized_value <- 1 - mean_KL_divergences

          first_cycle <- modes_df[1,2]

          # Get the minimum value from the second column of modes_df
          min_value_cycle <- min(modes_df[, 2], na.rm = TRUE)

          print(min_value_cycle)

          # Get the index of the minimum value from the second column of modes_df
          min_index <- which.min(modes_df[, 2])

          # Get the name of the row corresponding to the minimum value
          min_row_name <- rownames(modes_df)[min_index]

  #         gauge_plot <- gauge_plot %>%
  # layout(
  #   annotations = list(
  #     text = "Your Text Here",
  #     x = 0.5,  # X-coordinate position (0 to 1)
  #     y = 0.5,  # Y-coordinate position (0 to 1)
  #     showarrow = FALSE  # Set to TRUE if you want an arrow
  #   )
  # )

          # Print the name of the row corresponding to the minimum value
          print('min name')
          print(min_row_name)

          print("here is kl divs")
          print(kl_df)

          gauge_plot <- plot_ly(
            domain = list(x = c(0, 1), y = c(0, 1)),
            value = min_value_cycle,  # Set the value to the mean of KL divergences
            title = list(text = "Data Quality - Mode Deviation"),
            type = "indicator",
            mode = "gauge+number",
            gauge = list(
              steps = list(
                list(range = c(0, 2*first_cycle), color = "red"),
                list(range = c(0.2*first_cycle, 0.4*first_cycle), color = "orange"),
                list(range = c(0.4*first_cycle, 0.6*first_cycle), color = "yellow"),
                list(range = c(0.6*first_cycle, 0.8*first_cycle), color = "green"),
                list(range = c(0.8*first_cycle, first_cycle), color = "darkgreen")
              ),
              axis = list(range = list(0, modes_df[1,2])),  # Set the range from 0 to 1
              bar = list(color = "black"),  # Customize the color of the gauge
              threshold = list(
              line = list(color = "black", width = 4),
              thickness = 0.75,
              value = min_value_cycle)
              
              ))


            gauge_plot <- gauge_plot %>%
              layout(
                annotations = list(
                  text = paste0("lowest cycle intensity: ", min_row_name),
                  x = 0.5,  # X-coordinate position (0 to 1)
                  y = 0.5,  # Y-coordinate position (0 to 1)
                  showarrow = FALSE  # Set to TRUE if you want an arrow
                )
              )


              gauge_plot <- gauge_plot %>%
              layout(
                annotations = list(
                  text = paste0("First Cycle Intensity"),
                                    # text = paste0("First Cycle Intensity: ", modes_df[1,2]),

                  x = 1,  # X-coordinate position (0 to 1)
                  y = 0,  # Y-coordinate position (0 to 1)
                  showarrow = FALSE  # Set to TRUE if you want an arrow
                )
              )


            # axis =list(range = list(NULL, 500)),

          # modes_df[1,2]
    
            


          # gauge_plot <- gauge_plot %>%
          #   layout(
          #     margin = list(l=20,r=30),
          #     gauge = list(
          #       axis = list(range = list(1, 0.5)),  # Set the range from 0 to 1
          #       bar = list(color = "darkblue"),  # Customize the color of the gauge
          #       steps = list(
          #         list(range = c(0, 1), color = "lightgray")  # Customize steps if needed
          #       )
          #     )
          #   )

          # Render the gauge plot
          output$quality_gauge <- renderPlotly({
            gauge_plot
          })


          # library(plotly)

            # Define upper bound as 0 and lower bound as the maximum value from kl_df
            upper_bound <- 0
            # lower_bound <- kl_df[max_index, 2]


            # Get the index of the maximum value from the second column of kl_df
            max_index <- which.max(kl_df[, 2])

            lower_bound <- kl_df[max_index, 2]


            # Get the name of the row corresponding to the maximum value
            max_row_name <- rownames(kl_df)[max_index]

            # Print the name of the row corresponding to the maximum value
            print(max_row_name)

#             # Calculate percentiles for the KL divergences
# percentiles <- quantile(kl_df[, 2], probs = seq(0, 1, 0.2))

# # Extract the values for each percentile
# percentile_values <- percentiles[order(percentiles, decreasing = TRUE)]

# library(plotly)

# # Define upper bound as 0 and lower bound as the highest KL divergence value
# upper_bound <- 0
# lower_bound <- percentile_values[1]  # The highest KL divergence value

# # Create the gauge plot
# gauge_plot_kl <- plot_ly(
#   domain = list(x = c(0, 1), y = c(0, 1)),
#   value = lower_bound,  # Set the value to the highest KL divergence value
#   title = list(text = "Data Quality - KL Divergence"),
#   type = "indicator",
#   mode = "gauge+number",
#   gauge = list(
#     steps = list(
#       list(range = c(0, 2 * lower_bound), color = "red"),
#       list(range = c(0.2 * lower_bound, 0.4 * lower_bound), color = "orange"),
#       list(range = c(0.4 * lower_bound, 0.6 * lower_bound), color = "yellow"),
#       list(range = c(0.6 * lower_bound, 0.8 * lower_bound), color = "green"),
#       list(range = c(0.8 * lower_bound, lower_bound), color = "darkgreen")
#     ),
#     axis = list(range = list(upper_bound, lower_bound)),  # Set the range from upper_bound to lower_bound
#     bar = list(color = "black"),  # Customize the color of the gauge
#     threshold = list(
#       line = list(color = "black", width = 4),
#       thickness = 0.75,
#       value = lower_bound
#     )
#   )
# )

# # Print the gauge plot
# gauge_plot_kl



            # Create the gauge plot
            gauge_plot_kl <- plot_ly(
              domain = list(x = c(0, 1), y = c(0, 1)),
              value = lower_bound,  # Set the value to the maximum value from kl_df
              title = list(text = "Data Quality - KL Divergence"),
              type = "indicator",
              mode = "gauge+number",
              gauge = list(
                steps = list(
                  list(range = c(0, 2 * lower_bound), color = "red"),
                  list(range = c(0.2 * lower_bound, 0.4 * lower_bound), color = "orange"),
                  list(range = c(0.4 * lower_bound, 0.6 * lower_bound), color = "yellow"),
                  list(range = c(0.6 * lower_bound, 0.8 * lower_bound), color = "green"),
                  list(range = c(0.8 * lower_bound, 0), color = "darkgreen")
                ),
                axis = list(range = list(NULL, 0)),  # Set the range from upper_bound to lower_bound
                bar = list(color = "black"),  # Customize the color of the gauge
                threshold = list(
                  line = list(color = "black", width = 4),
                  thickness = 0.75,
                  value = lower_bound
                )
              )
            )

          #     # Render the gauge plot
          # output$quality_gauge_kl <- renderPlotly({
          #   gauge_plot_kl 
          # })


            
                })

 # Render the plot
        output$modes_plot_output <- renderPlot({
          regression_mode_react()
        })
        
         output$klPLOToutput <- renderPlot({
          klPlotReact()
        })

    
#        observeEvent(input$RetainSubset, {
#   brushed_data(brushedPoints(filtered_data_reactive()[, c("CellID", "X_centroid", "Y_centroid")], input$plot_brush, xvar = "X_centroid", yvar = "Y_centroid"))
#   filtered_data_reactive(brushed_data())
#   output$image_garage_output <- renderPlot({
#     ggplot(filtered_data_reactive(), aes(x = X_centroid, y = Y_centroid)) +
#       geom_point() +
#       xlab('X Centroid') + ylab('Y Centroid')
#   })
# })

    observeEvent(input$RetrySubset, {
      selected_marker <- input$precrev_marker
  brushed_data(NULL)
  filtered_data_reactive(filtered_data_original())
  output$image_garage_output <- renderPlotly({
    ggplot(filtered_data_reactive(), aes(x = X_centroid, y = Y_centroid, color = filtered_data_reactive()[[selected_marker]])) +
      geom_point(alpha = selected_opacity()) +
      xlab('X Centroid') + ylab('Y Centroid') + 
      labs(title = "Digital Marker Overlay") +
        scale_color_continuous(name = "Intensity (logged)") +
           theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ) 
  })
})


observeEvent(input$RemoveSubset, {
  selected_marker <- input$precrev_marker
  brushed_data(brushedPoints(filtered_data_reactive()[, c("CellID", "X_centroid", "Y_centroid")], input$plot_brush, xvar = "X_centroid", yvar = "Y_centroid"))
  filtered_data_new <- filtered_data_reactive()[!(filtered_data_reactive()$CellID %in% brushed_data()$CellID), ]
  filtered_data_reactive(filtered_data_new)
  output$image_garage_output <- renderPlotly({
    ggplot(filtered_data_reactive(), aes(x = X_centroid, y = Y_centroid, color = filtered_data_reactive()[[selected_marker]])) +
      geom_point(alpha = selected_opacity()) +
      xlab('X Centroid') + ylab('Y Centroid') +
      labs(title = "Digital Marker Overlay") +
        scale_color_continuous(name = "Intensity (logged)") +
           theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ) 
  })
  
})


# observeEvent(event_data("plotly_selected"), {
#     selected_marker <- input$precrev_marker

#   # event_data("plotly_selected") contains the selected points' information
#   selected_points <- event_data("plotly_selected")
#   print("selected points printed below!!")
#   print(selected_points)


#   selected_points$
  
#   if (!is.null(selected_points)) {
#     filtered_data_reactive(filtered_data_reactive() %>% 
#                              filter(!CellID %in% selected_points$CellID))
    
#     output$image_garage_output <- renderPlotly({
#       ggplot(filtered_data_reactive(), aes(x = X_centroid, y = Y_centroid, color = filtered_data_reactive()[[selected_marker]])) +
#         geom_point(alpha = selected_opacity()) +
#         xlab('X Centroid') + ylab('Y Centroid')
#     })
#   }
# })


observeEvent(event_data("plotly_selected"), {
  selected_marker <- input$precrev_marker

  # event_data("plotly_selected") contains the selected points' information
  selected_points <- event_data("plotly_selected")
  print("selected points printed below!!")
  print(selected_points)

  if (!is.null(selected_points)) {
    # Filter the reactive dataframe based on X_centroid and Y_centroid
    filtered_data_reactive(filtered_data_reactive() %>%
                             filter(!X_centroid %in% selected_points$x & !Y_centroid %in% selected_points$y))

    output$image_garage_output <- renderPlotly({
      ggplot(filtered_data_reactive(), aes(x = X_centroid, y = Y_centroid, color = filtered_data_reactive()[[selected_marker]])) +
        geom_point(alpha = selected_opacity()) +
      xlab('X Centroid') + ylab('Y Centroid') +
      labs(title = "Digital Marker Overlay") +
        scale_color_continuous(name = "Intensity (logged)") +
           theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ) 
    })
  }
})




    observeEvent(input$UseSubset, {
    # #Get the unique imageids from filtered_data_reactive
    # selected_imageid <- unique(filtered_data_reactive()$imageid)

    # # Update uploaded_df by replacing rows with matching imageid
    # updated_df <- uploaded_df()
    # updated_df[updated_df$imageid %in% selected_imageid, ] <- filtered_data_reactive()

    # # Update uploaded_df with the updated data
    # uploaded_df(updated_df)

    #   # Get the common identifier
    # common_ids <- filtered_data_reactive()$imageid
    
    # print(common_ids)
    # print(uploaded_df()$imageid)
    
    # # Subset the rows in uploaded_df() based on the common identifier
    # uploaded_df()[uploaded_df()$imageid %in% common_ids, ] <- filtered_data_reactive()

    uploaded_df(filtered_data_reactive())
    

  # showNotification("Data subset loaded for patient!", duration = NULL, id = "loadedsubset")
    # Display a completion message to the user
    shinyalert::shinyalert(
      title = "Subset Loaded",
      text = "You may use this subsetted data in the Randkluft Tab for analysis.",
      type = "success"
    )

})

# observeEvent(input$UseSubset, {
#   # Get the common identifier
#   common_ids <- filtered_data_reactive()$imageid
  
#   # Subset the rows in uploaded_df() based on the common identifier
#   uploaded_df()[uploaded_df()$imageid %in% common_ids, ] <- filtered_data_reactive()
  
#   print(colnames(uploaded_df()))
# })



# observeEvent(input$UseSubset, {
#   # Get the unique imageids from filtered_data_reactive
#   selected_imageid <- unique(filtered_data_reactive()$imageid)

#   # Update uploaded_df by replacing rows with matching imageid
#   updated_df <- uploaded_df()
#   updated_df[updated_df$imageid %in% selected_imageid, ] <- filtered_data_reactive()

#   # Update uploaded_df with the updated data
#   uploaded_df(updated_df)
# })


 

# bursher <- reactive({
#       input$plotBrush
#     })
# ... Your previous code ...

# Inside the server function


# output$image_garage_output <- renderPlotly({
#   # Get selected marker and patient
#   selected_patient <- input$patient_number_im_gargage
  
#   # Filter data based on selected patient
#   filtered_data <- uploaded_df()[uploaded_df()$imageid %in% selected_patient, ]
#   filtered_data_react(filtered_data)
  
#   p <- plot_ly(x = ~filtered_data$X_centroid, y = ~filtered_data$Y_centroid, type = 'scatter', mode = 'markers')
  
#   p <- layout(p, xaxis = list(title = 'X Centroid'), yaxis = list(title = 'Y Centroid'))
  
#   p <- p %>% highlight(on = 'plotly_selected', selectize = TRUE)
  
#   p
# })

# observeEvent(event_data("plotly_selected"), {
#   selected_points <- event_data("plotly_selected")
#   if (!is.null(selected_points)) {
#     selected_points <- filtered_data_react[selected_points$pointNumber + 1, ] # +1 because pointNumber is 0-based
#     newdf(selected_points)
#     print(newdf())
#   }
# })

  
  # ... Your previous code ...
  


# ... Your remaining code ...

# Download button handler to save the Phenotype Table as a CSV file
  output$downloadPhenotypeTable <- downloadHandler(
    filename = function() {
      "phenotype_table.csv"
    },
    content = function(file) {
      write.csv(phenotype_df(), file, row.names = FALSE)
    }
  )



      # Functions that return statistics from the models in list, save into csv and download
  output$downloadPhenotypes <- downloadHandler(
    filename = function() {
      if (!is.null(csv_save_file_path())) {
        if (!is.null(input$csv_name)) {
          return(paste0(input$csv_name, ".csv"))
        } else {

          return('Phenotypes_included.csv')
        }

      }
      showNotification("Run Randkluft first!", duration = NULL, id = "warngate")

    },
    content = function(file) {
      mydftosave = uploaded_df()
      write.csv(mydftosave, file, row.names = FALSE) 
    }
  )






      # Functions that return statistics from the models in list, save into csv and download
  output$downloadEstimations <- downloadHandler(
    filename = function() {
      if (!is.null(csv_save_file_path())) {
        if (!is.null(input$csv_name)) {
          return(paste0(input$csv_name, ".csv"))
        } else {

          return('Gates.csv')
        }

      }
      showNotification("Run Randkluft first!", duration = NULL, id = "warngate")

    },
    content = function(file) {
      mydftosave = csv_save_file_path()
      write.csv(mydftosave, file, row.names = FALSE) 
    }
  )

  output$download_example_data <- downloadHandler(
    filename = function() {
      'ExampleCellMarkerData.csv'
    },
    content = function(con) {
      write.csv(df_example, con, row.names=FALSE)
    }
  )


     output$download_example_data2 <- downloadHandler(
     filename = function() {
       'Phenotype_table_example.csv'
     },
     content = function(con) {
       write.csv(df_example_PHENOTYPE, con, row.names=FALSE)
     }
   )



observeEvent(c(input$xvar, input$yvar), {
  # Clear the input fields
  updateTextInput(session, "gate_xvar_update", value = "")
  updateTextInput(session, "gate_yvar_update", value = "")

  # Set react_xgate and react_ygate to NULL
  react_xgate(NULL)
  react_ygate(NULL)
})


 
  subsetted <- reactive({
    # shinyjs::show("plot2")

    req(input$patient_number_plot2)
    uploaded_df() |> filter(imageid %in% input$patient_number_plot2)  
   
  })

  
  subsetted_precrev_ts <- reactive({
    # shinyjs::show("plot2")

    req(input$patients_ts)
    uploaded_df() |> filter(imageid %in% input$patients_ts)  
   
  })

    subsetted_tri <- reactive({
    # shinyjs::show("plot2")

    req(input$trivar_patient)
    uploaded_df() |> filter(imageid %in% input$trivar_patient)  
   
  })

  react_xgate <- reactiveVal(NULL)
  react_ygate <- reactiveVal(NULL)

  output$plot2 <- renderPlot({
      generatePlotBivar()
    })


      output$plot_trivariate <- renderPlotly({
      generatePlotTrivar()
    })


      generatePlotBivar <- function() {

    output$plot2 <- renderPlot({
        dfplot2 <- subsetted()

        print(colnames(dfplot2))

        patient_selected <- input$patient_number_plot2

        xvar <- input$xvar
        yvar <- input$yvar

        print(patient_selected)
        print(xvar)
        print(yvar)

        print(resultdf_reactive())

        print(resultdf_reactive()$Marker == xvar)
        print(resultdf_reactive()$Patient == patient_selected)

        gate_xvar <- resultdf_reactive()$Gate[
            resultdf_reactive()$Marker == xvar &
            resultdf_reactive()$Patient == patient_selected
        ]

        gate_yvar <- resultdf_reactive()$Gate[
            resultdf_reactive()$Marker == yvar &
            resultdf_reactive()$Patient == patient_selected
        ]

        if (!is.null(react_xgate())) {
            gate_xvar <- react_xgate()
            if (xvar == yvar) {
                gate_yvar <- gate_xvar
            }
        }

        if (!is.null(react_ygate())) {
            gate_yvar <- react_ygate()
            if (xvar == yvar) {
                gate_xvar <- gate_yvar
            }
        }

        print(gate_xvar)
        print(gate_yvar)

        num_points <- nrow(dfplot2)
        num_pp <- sum(dfplot2[[xvar]] > gate_xvar & dfplot2[[yvar]] > gate_yvar)
        num_pm <- sum(dfplot2[[xvar]] > gate_xvar & dfplot2[[yvar]] <= gate_yvar)
        num_mp <- sum(dfplot2[[xvar]] <= gate_xvar & dfplot2[[yvar]] > gate_yvar)
        num_mm <- sum(dfplot2[[xvar]] <= gate_xvar & dfplot2[[yvar]] <= gate_yvar)

        prop_pp <- num_pp / num_points
        prop_pm <- num_pm / num_points
        prop_mp <- num_mp / num_points
        prop_mm <- num_mm / num_points

        complete_rows <- complete.cases(dfplot2[[xvar]], dfplot2[[yvar]])
        dfplot2 <- dfplot2[complete_rows, ]

        finite_rows <- is.finite(dfplot2[[xvar]]) & is.finite(dfplot2[[yvar]])
        dfplot2 <- dfplot2[finite_rows, ]

        density_values_bi <- get_density(dfplot2[[xvar]], dfplot2[[yvar]], n = 100)

        p <- ggplot(dfplot2, aes(!!input$xvar,!!input$yvar)) + 
            geom_vline(xintercept = gate_xvar, color = "orange") +
            geom_hline(yintercept = gate_yvar, color = "orange") +
            scale_color_viridis(option="turbo") + 
            geom_point(aes(color = density_values_bi), size = 0.3, alpha = 0.35) +
            theme(
                legend.position = "none",
                plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
                axis.title = element_text(size = 16, family = "Arial"),
                axis.text = element_text(size = 14, family = "Arial"),
                panel.background = element_rect(fill = "#f8f8f8"),  # very pale grey for the plotting area
                plot.background = element_rect(fill = "#f8f8f8")
            ) +
            labs(title = "Bivariate Density") 

        p <- p +
            geom_text(x = max(dfplot2[[xvar]])-max(dfplot2[[xvar]])/40, y = max(dfplot2[[yvar]])-max(dfplot2[[yvar]])/40, label = paste(round(prop_pp,3)), color = "red", size = 6) +
            geom_text(x = max(dfplot2[[xvar]])-max(dfplot2[[xvar]])/40, y = min(dfplot2[[yvar]])+min(dfplot2[[yvar]])/40, label = paste(round(prop_pm,3)), color = "green", size = 6) +
            geom_text(x = min(dfplot2[[xvar]])+min(dfplot2[[xvar]])/40, y = max(dfplot2[[yvar]])-max(dfplot2[[yvar]])/40, label = paste(round(prop_mp,3)), color = "blue", size = 6) +
            geom_text(x = min(dfplot2[[xvar]])+min(dfplot2[[xvar]])/40, y = min(dfplot2[[yvar]])+min(dfplot2[[yvar]])/40, label = paste(round(prop_mm,3)), color = "black", size = 6)

        overlay_plot2 <- ggplot(dfplot2, aes(x = X_centroid, y = Y_centroid)) +
            geom_point(aes(
                color = case_when(
                    dfplot2[[xvar]] > gate_xvar & dfplot2[[yvar]] > gate_yvar ~ "+/+",
                    dfplot2[[xvar]] <= gate_xvar & dfplot2[[yvar]] <= gate_yvar ~ "-/-",
                    dfplot2[[xvar]] > gate_xvar & dfplot2[[yvar]] <= gate_yvar ~ "+/-",
                    dfplot2[[xvar]] <= gate_xvar & dfplot2[[yvar]] > gate_yvar ~ "-/+",
                    TRUE ~ "Other"
                )
            ), size = 0.3, alpha = 0.7) +
            scale_color_manual(
                guide = guide_legend(title = "", override.aes = list(size = 5)),
                values = c(
                    "+/+" = "red",
                    "-/-" = "grey",
                    "+/-" = "green",
                    "-/+" = "blue",
                    "Other" = "black"
                )
            ) +
            labs(title = "Bivariate Gating", x = "X Centroid", y = "Y Centroid") +
            theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
                axis.title = element_text(size = 16, family = "Arial"),
                axis.text = element_text(size = 14, family = "Arial"),
                legend.position = "none",
                panel.background = element_rect(fill = "white"),  # very pale grey for the plotting area
                plot.background = element_rect(fill = "white")
            ) 

        prop_above_cutoff_xvar <- sum(dfplot2[[xvar]] > gate_xvar) / nrow(dfplot2)

        contour_plot_xvar <- ggplot(dfplot2, aes(x = X_centroid, y = Y_centroid)) +
            geom_point(aes(color = ifelse(dfplot2[[xvar]] > gate_xvar, "Positive", "Negative")), size = 0.3, alpha = 0.35) +
            scale_color_manual(guide = guide_legend(title = ""), values = c("Positive" = "green", "Negative" = "grey")) +
            geom_density_2d(data = subset(dfplot2, dfplot2[[xvar]] > gate_xvar), color = 'black') +  
            theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
                axis.title = element_text(size = 16, family = "Arial"),
                axis.text = element_text(size = 14, family = "Arial"),
                legend.position = "bottom",
                panel.background = element_rect(fill = "white"),  # very pale grey for the plotting area
                plot.background = element_rect(fill = "white")
            ) +
            labs(title = paste(xvar, "+ cell=", round(prop_above_cutoff_xvar, 3), sep = ""), x = "X Centroid", y = "Y Centroid")

        prop_above_cutoff_yvar <- sum(dfplot2[[yvar]] > gate_yvar) / nrow(dfplot2)

        contour_plot_yvar <- ggplot(dfplot2, aes(x = X_centroid, y = Y_centroid)) +
            geom_point(aes(color = ifelse(dfplot2[[yvar]] > gate_yvar, "Positive", "Negative")), size = 0.3, alpha = 0.35) +
            scale_color_manual(guide = guide_legend(title = ""), values = c("Positive" = "blue", "Negative" = "grey")) +
            geom_density_2d(data = subset(dfplot2, dfplot2[[yvar]] > gate_yvar), color = 'black') +  
            theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
                axis.title = element_text(size = 16, family = "Arial"),
                axis.text = element_text(size = 14, family = "Arial"),
                legend.position = "bottom",
                panel.background = element_rect(fill = "white"),  # very pale grey for the plotting area
                plot.background = element_rect(fill = "white")
            ) +
            labs(title = paste(yvar, "+ cell=", round(prop_above_cutoff_yvar, 3), sep = ""), x = "X Centroid", y = "Y Centroid")

        arranged_plots <- grid.arrange(p, overlay_plot2, contour_plot_xvar, contour_plot_yvar, ncol = 2)

        # Save the arranged plots with higher resolution
        ggsave("arranged_plots.png", arranged_plots, dpi = 1000, width = 12, height = 8)

        print(arranged_plots)
    })
}

      

generatePlotTrivar <- function() {

    output$plot_trivariate <- renderPlotly({

  
    dfplot2 <- subsetted_tri()


  patient_selected <- input$trivar_patient

    xvar <- input$xvarTri
    yvar <- input$yvarTri
    zvar <- input$zvar


    # print(patient_selected)
    # print(xvar)
    # print(yvar)

    # print(resultdf_reactive())

    # print(resultdf_reactive()$Marker == xvar)
    # print(resultdf_reactive()$Patient == patient_selected)

    gate_xvar <- resultdf_reactive()$Gate[
    resultdf_reactive()$Marker == xvar &
    resultdf_reactive()$Patient == patient_selected
  ]

  gate_yvar <- resultdf_reactive()$Gate[
    resultdf_reactive()$Marker == yvar &
    resultdf_reactive()$Patient == patient_selected
  ]

  gate_zvar <- resultdf_reactive()$Gate[
    resultdf_reactive()$Marker == zvar &
    resultdf_reactive()$Patient == patient_selected
  ]

    # Check if react_xgate or react_ygate is defined
    # if (!is.null(react_xgate())) {
    #   # Code to execute when react_xgate is defined
    #   gate_xvar <- react_xgate()
    #   if (xvar == yvar) {
    #   gate_yvar <- gate_xvar
    #   }
    # }


    # if (!is.null(react_ygate())) {
    #   # Code to execute when react_ygate is defined
    #   gate_yvar <- react_ygate()
    #   if (xvar == yvar) {
    #   gate_xvar <- gate_yvar
    #   }
    # }


  

    print(gate_xvar)
    print(gate_yvar)

    print(xvar)
    print(yvar)
    print(zvar)

        axx <- list(
      title = paste(xvar)
    )

    axy <- list(
      title = paste(yvar)
    )

    axz <- list(
      title = paste(zvar)
    )



    # p <- plot_ly(dfplot2, x = ~dfplot2[[xvar]], y = ~dfplot2[[yvar]], z = ~dfplot2[[zvar]], type = "scatter3d") %>%
    #   layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz), 
    #   # Add vertical lines (intercepts)
    #     add_trace(data = data.frame(x = gate_xvar, y = dfplot2[[yvar]], z = dfplot2[[zvar]]), type = "scatter3d", mode = "lines", line = list(color = "red")),
    #     add_trace(data = data.frame(x = dfplot2[[xvar]], y = gate_yvar, z = dfplot2[[zvar]]), type = "scatter3d", mode = "lines", line = list(color = "blue")),
    #     add_trace(data = data.frame(x = dfplot2[[xvar]], y = dfplot2[[yvar]], z = gate_zvar), type = "scatter3d", mode = "lines", line = list(color = "green")))

#     # gate_zvar_matrix = data.frame(x = range(dfplot2[[xvar]]), y = range(dfplot2[[yvar]]), z = gate_zvar)
# # Suppose gate_zvar is the z-value for the plane
# gate_zvar_plane <- expand.grid(
#   x = range(dfplot2[[xvar]]),
#   y = range(dfplot2[[yvar]])
# )

# # Create a matrix to represent the plane with constant z-values
# # Define the dimensions (rows and columns) based on the length of x and y coordinates
# rows <- length(unique(gate_zvar_plane$y))  # Number of rows
# cols <- length(unique(gate_zvar_plane$x))  # Number of columns

# # Create a matrix with z-values filled with gate_zvar
# gate_z_matrix <- matrix(rep(gate_zvar, rows * cols), nrow = rows, ncol = cols, byrow = TRUE)

# # Assign the z-values to gate_zvar_plane
# gate_zvar_plane$z <- as.vector(t(gate_z_matrix))

# print(dim)

# # Get the desired z-value
# z_value <- gate_zvar  # Replace this with the desired z-value

# # Create a matrix with the same dimensions as your data with the specified z-value
# z_values <- matrix(rep(z_value, nrow(dfplot2) * ncol(dfplot2)), nrow = nrow(dfplot2[[xvar]]), ncol = ncol(dfplot2[[yvar]]))

# Create a Plotly surface plot
# fig <- plot_ly(z = ~z_values, type = "surface", showscale = FALSE)



     p <- plot_ly(dfplot2, x = ~dfplot2[[xvar]], y = ~dfplot2[[yvar]], z = ~dfplot2[[zvar]], type = "scatter3d") %>%
    layout(scene = list(
      xaxis = axx,
      yaxis = axy,
      zaxis = axz)
      # Add planes to separate regions
      # add_trace(data = data.frame(x = c(gate_xvar, max(dfplot2[[xvar]])), y = dfplot2[[yvar]], z = dfplot2[[zvar]]), type = "mesh3d", color = I("rgba(255, 0, 0, 0.3)")) %>%
      #   add_trace(data = data.frame(x = dfplot2[[xvar]], y = c(gate_yvar, max(dfplot2[[yvar]])), z = dfplot2[[zvar]]), type = "mesh3d", color = I("rgba(0, 255, 0, 0.3)")) %>%
      #   add_trace(data = data.frame(x = dfplot2[[xvar]], y = dfplot2[[yvar]], z = c(gate_zvar, max(dfplot2[[zvar]]))), type = "mesh3d", color = I("rgba(0, 0, 255, 0.3)"))
      )
      

  # p <- p %>% add_surface(z = ~z_values, type = "surface")

  # p <- plot_ly(data = xy_grid, x = ~x, y = ~y, z = ~z, type = "surface")



  


p



#     z <- c(
#   c(8.83,8.89,8.81,8.87,8.9,8.87),
#   c(8.89,8.94,8.85,8.94,8.96,8.92),
#   c(8.84,8.9,8.82,8.92,8.93,8.91),
#   c(8.79,8.85,8.79,8.9,8.94,8.92),
#   c(8.79,8.88,8.81,8.9,8.95,8.92),
#   c(8.8,8.82,8.78,8.91,8.94,8.92),
#   c(8.75,8.78,8.77,8.91,8.95,8.92),
#   c(8.8,8.8,8.77,8.91,8.95,8.94),
#   c(8.74,8.81,8.76,8.93,8.98,8.99),
#   c(8.89,8.99,8.92,9.1,9.13,9.11),
#   c(8.97,8.97,8.91,9.09,9.11,9.11),
#   c(9.04,9.08,9.05,9.25,9.28,9.27),
#   c(9,9.01,9,9.2,9.23,9.2),
#   c(8.99,8.99,8.98,9.18,9.2,9.19),
#   c(8.93,8.97,8.97,9.18,9.2,9.18)
# )
# dim(z) <- c(15,6)
# z2 <- z + 1
# z3 <- z - 1

# fig <- plot_ly(showscale = FALSE)
# fig <- fig %>% add_surface(z = ~z)
# fig <- fig %>% add_surface(z = ~z2, opacity = 0.98)
# fig <- fig %>% add_surface(z = ~z3, opacity = 0.98)

# fig

    
# p





# p <- plot_ly(dfplot2, x = ~dfplot2[[xvar]], y = ~dfplot2[[yvar]], z = ~dfplot2[[zvar]], type = "scatter3d") %>%
#   layout(scene = list(
#     xaxis = axx,
#     yaxis = axy,
#     zaxis = axz
#   )) %>%
#   add_trace(data = data.frame(x = c(gate_xvar, gate_xvar), y = range(dfplot2[[yvar]]), z = range(dfplot2[[zvar]])), type = "scatter3d", mode = "lines", line = list(color = "red", width = 2)) %>%
#   add_trace(data = data.frame(x = range(dfplot2[[xvar]]), y = c(gate_yvar, gate_yvar), z = range(dfplot2[[zvar]])), type = "scatter3d", mode = "lines", line = list(color = "green", width = 2)) %>%
#   add_trace(data = data.frame(x = range(dfplot2[[xvar]]), y = range(dfplot2[[yvar]]), z = c(gate_zvar, gate_zvar)), type = "scatter3d", mode = "lines", line = list(color = "blue", width = 2))



# p <- plot_ly(dfplot2, x = ~dfplot2[[xvar]], y = ~dfplot2[[yvar]], z = ~dfplot2[[zvar]], type = "scatter3d") %>%
#   layout(scene = list(
#     xaxis = axx,
#     yaxis = axy,
#     zaxis = axz
#   )) %>%
#   add_trace(data = data.frame(x = c(min(dfplot2[[xvar]]), max(dfplot2[[xvar]])), y = c(gate_yvar, gate_yvar), z = c(gate_zvar, gate_zvar), type = "mesh3d", color = I("rgba(0, 255, 0, 0.3)"))) %>%
#   add_trace(data = data.frame(x = c(gate_xvar, gate_xvar), y = c(min(dfplot2[[yvar]]), max(dfplot2[[yvar]])), z = c(gate_zvar, gate_zvar), type = "mesh3d", color = I("rgba(255, 0, 0, 0.3)"))) %>%
#   add_trace(data = data.frame(x = c(gate_xvar, gate_xvar), y = c(gate_yvar, gate_yvar), z = c(min(dfplot2[[zvar]]), max(dfplot2[[zvar]])), type = "mesh3d", color = I("rgba(0, 0, 255, 0.3)")))

# p
   


  })

}




    


  # output$plot2 <- renderPlot({

  
  #   dfplot2 <- subsetted()

  #   print(colnames(dfplot2))

  #     # gate_xvar <- resultdf_reactive()$Gate[resultdf_reactive()$Marker == input$xvar]
  #     # gate_yvar <- resultdf_reactive()$Gate[resultdf_reactive()$Marker == input$yvar]

  # #     resultdf_reactive()$Gate[
  # #   resultdf_reactive()$Marker == chosen_marker &
  # #   resultdf_reactive()$Patient == chosen_patient
  # # ]


  #   xvar <- input$xvar
  #   yvar <- input$yvar

  #   # gate_xvar <- input$gate_xvar_update
  #   # gate_yvar <- input$gate_xvar_update

  #   gate_vliner <- dfplot2[[xvar]]
  #   gate_hliner <- dfplot2[[yvar]]

  #       # dfplot2$density <- get_density(gate_vliner, gate_hliner, n = 100)


  #   gate_vliner <- remove_outliers2(gate_vliner)
  #   gate_hliner <- remove_outliers2(gate_hliner)

  #   gate_vliner <- gate_vliner[!is.infinite(gate_vliner)]
  #   gate_hliner <- gate_hliner[!is.infinite(gate_hliner)]

  #   gate_xvar <- skew_gate(gate_vliner, alpha=0.01)
  #   gate_yvar <- skew_gate(gate_hliner, alpha=0.01)


  #    # Calculate proportions
  #   num_points <- nrow(dfplot2)
  #   num_pp <- sum(dfplot2[[xvar]] > gate_xvar$cutoff & dfplot2[[yvar]] > gate_yvar$cutoff)
  #   num_pm <- sum(dfplot2[[xvar]] > gate_xvar$cutoff & dfplot2[[yvar]] <= gate_yvar$cutoff)
  #   num_mp <- sum(dfplot2[[xvar]] <= gate_xvar$cutoff & dfplot2[[yvar]] > gate_yvar$cutoff)
  #   num_mm <- sum(dfplot2[[xvar]] <= gate_xvar$cutoff & dfplot2[[yvar]] <= gate_yvar$cutoff)

  #   prop_pp <- num_pp / num_points
  #   prop_pm <- num_pm / num_points
  #   prop_mp <- num_mp / num_points
  #   prop_mm <- num_mm / num_points

  #         output$prop_summary <- renderPrint({
  #     cat("Proportion of ++:", prop_pp, "\n")
  #     cat("Proportion of -+:", prop_pm, "\n")
  #     cat("Proportion of +-:", prop_mp, "\n")
  #     cat("Proportion of --:", prop_mm, "\n")
  #   })

  #   # message("Proportion of ++:", prop_pp, "\n")
  #   # message("Proportion of -+:", prop_pm, "\n")
  #   # message("Proportion of +-:", prop_mp, "\n")
  #   # message("Proportion of --:", prop_mm, "\n")


  #   # dfplot2$density <- get_density(dfplot2$input$xvar, dfplot2$input$yvar, n = 100)


  #   p <- ggplot(dfplot2, aes(!!input$xvar,!!input$yvar)) + 
  #   geom_vline(xintercept = gate_xvar$cutoff, color = "red") +
  #   geom_hline(yintercept = gate_yvar$cutoff, color = "blue") +
  #   # scale_color_viridis(option="turbo") + 

    
  #   list(
  #     theme(legend.position = "bottom"),
  #     if (input$by_species)
  #       aes(color = imageid),
  #     geom_point(),

  #     if (input$smooth)
  #       geom_smooth()
  #   )
  #   if (input$show_margins) {
  #     margin_type <- if (input$by_species)
  #       "density"
  #     else
  #       "histogram"
  #     p <-
  #       p |> ggExtra::ggMarginal(
  #         type = margin_type,
  #         margins = "both",
  #         size = 8,
  #         groupColour = input$by_species,
  #         groupFill = input$by_species
  #       )
  #   }
    


  # print(p)


  # })



#   observeEvent(input$update_gates_bivariate, {

#    react_xgate(input$gate_xvar_update)
#    react_ygate(input$gate_yvar_update)


#     # # Check if react_xgate or react_ygate is defined
#     # if (!is.null(react_xgate())) {
#     #   # Code to execute when react_xgate is defined
#     #   gate_xvar$cutoff <- react_xgate()
#     # }


#     # if (!is.null(react_ygate())) {
#     #   # Code to execute when react_ygate is defined
#     #   gate_yvar$cutoff <- react_ygate()
#     # }


   

#    generatePlotBivar()

# })

# observeEvent(input$update_gates_bivariate, {

#   react_xgate_value <- input$gate_xvar_update
#   react_ygate_value <- input$gate_yvar_update

#   print(react_xgate_value)
#   print(react_ygate_value)

#   # Check if input is blank and set react_xgate accordingly
#   if (react_xgate_value == "") {
#     react_xgate(NULL)
#   } else {
#     react_xgate(react_xgate_value)
#   }

#   # Check if input is blank and set react_ygate accordingly
#   if (react_ygate_value == "") {
#     react_ygate(NULL)
#   } else {
#     react_ygate(react_ygate_value)
#   }

#   generatePlotBivar()

#   # input$gate_xvar_update <- ""
#   # input$gate_yvar_update <- ""
# })


# observeEvent(input$update_gates_bivariate, {
#   # Check if the input values are not blank
#   if (input$gate_xvar_update != "") {
#     react_xgate(as.numeric(input$gate_xvar_update))  # Assuming gate values are numeric
#   } else {
#     # If input is blank, reset to NULL or original value
#     react_xgate(NULL)  # Set to NULL or original value
#   }

#   if (input$gate_yvar_update != "") {
#     react_ygate(as.numeric(input$gate_yvar_update))  # Assuming gate values are numeric
#   } else {
#     # If input is blank, reset to NULL or original value
#     react_ygate(NULL)  # Set to NULL or original value
#   }

#   generatePlotBivar()
# })

# observeEvent(input$update_gates_bivariate, {
#   # Check if the input values are not blank
#   if (!is.null(input$gate_xvar_update) && input$gate_xvar_update != "") {
#     react_xgate(as.numeric(input$gate_xvar_update))  # Assuming gate values are numeric
#   } else {
#     # If input is blank, reset to NULL or original value
#     react_xgate(NULL)  # Set to NULL or original value
#   }

#   if (!is.null(input$gate_yvar_update) && input$gate_yvar_update != "") {
#     react_ygate(as.numeric(input$gate_yvar_update))  # Assuming gate values are numeric
#   } else {
#     # If input is blank, reset to NULL or original value
#     react_ygate(NULL)  # Set to NULL or original value
#   }

#   print("gate x and y var updates")
#   print(input$gate_xvar_update)

#   print(input$gate_yvar_update)

#   generatePlotBivar()
# })


# observeEvent(input$update_gates_bivariate, {
#   print("Before react_xgate update")
#   print(input$gate_xvar_update)
#   if (!is.null(input$gate_xvar_update) && input$gate_xvar_update != "" && !is.na(input$gate_xvar_update)) {
#     react_xgate(as.numeric(input$gate_xvar_update))
#   } else {
#     react_xgate(NULL)
#   }
#   print("After react_xgate update")

#   print("Before react_ygate update")
#   if (!is.null(input$gate_yvar_update) && input$gate_yvar_update != "") {
#     react_ygate(as.numeric(input$gate_yvar_update))
#   } else {
#     react_ygate(NULL)
#   }
#   print("After react_ygate update")

#   generatePlotBivar()
# })


observeEvent(input$update_gates_bivariate, {
  print("Before react_xgate update")
  if (!is.na(input$gate_xvar_update) && !is.null(input$gate_xvar_update) && input$gate_xvar_update != "") {
    react_xgate(as.numeric(input$gate_xvar_update))
  } else {
    react_xgate(NULL)
  }
  print("After react_xgate update")

  print("Before react_ygate update")
  if (!is.na(input$gate_yvar_update) && !is.null(input$gate_yvar_update) && input$gate_yvar_update != "") {
    react_ygate(as.numeric(input$gate_yvar_update))
  } else {
    react_ygate(NULL)
  }
  print("After react_ygate update")

  generatePlotBivar()
})






#   output$prop_summary <- renderPrint({
#   cat("Proportion of ++:", prop_pp, "\n")
#   cat("Proportion of -+:", prop_pm, "\n")
#   cat("Proportion of +-:", prop_mp, "\n")
#   cat("Proportion of --:", prop_mm, "\n")
# })


#   output$plot2 <- renderPlot({
#   dfplot2 <- subsetted()

#   xvar <- input$xvar
#   yvar <- input$yvar

#   gate_vliner <- dfplot2[[xvar]]
#   gate_hliner <- dfplot2[[yvar]]

#   gate_vliner <- remove_outliers2(gate_vliner)
#   gate_hliner <- remove_outliers2(gate_hliner)

#   gate_vliner <- gate_vliner[!is.infinite(gate_vliner)]
#   gate_hliner <- gate_hliner[!is.infinite(gate_hliner)]

#   gate_xvar <- skew_gate(gate_vliner, alpha=0.01)
#   gate_yvar <- skew_gate(gate_hliner, alpha=0.01)

#   print(gate_xvar$cutoff)
#   print(gate_yvar$cutoff)

#   # Calculate proportions
#   num_points <- nrow(dfplot2)
#   num_pp <- sum(dfplot2[[xvar]] > gate_xvar$cutoff & dfplot2[[yvar]] > gate_yvar$cutoff)
#   num_pm <- sum(dfplot2[[xvar]] > gate_xvar$cutoff & dfplot2[[yvar]] <= gate_yvar$cutoff)
#   num_mp <- sum(dfplot2[[xvar]] <= gate_xvar$cutoff & dfplot2[[yvar]] > gate_yvar$cutoff)
#   num_mm <- sum(dfplot2[[xvar]] <= gate_xvar$cutoff & dfplot2[[yvar]] <= gate_yvar$cutoff)

#   prop_pp <- num_pp / num_points
#   prop_pm <- num_pm / num_points
#   prop_mp <- num_mp / num_points
#   prop_mm <- num_mm / num_points

#   cat("Proportion of ++:", prop_pp, "\n")
#   cat("Proportion of -+:", prop_pm, "\n")
#   cat("Proportion of +-:", prop_mp, "\n")
#   cat("Proportion of --:", prop_mm, "\n")

#   # Rest of your ggplot code...
# })


#   output$plot2 <- renderPlot({
#   dfplot2 <- subsetted()

#   p <- ggplot(dfplot2, aes(!!input$xvar, !!input$yvar)) +
#     theme(legend.position = "bottom") +
#     if (input$by_species) aes(color = imageid) +
#     geom_point() +
#     if (input$smooth) geom_smooth()

#   if (input$show_margins) {
#     margin_type <- if (input$by_species) "density" else "histogram"
#     p <- p + ggExtra::ggMarginal(
#       type = margin_type,
#       margins = "both",
#       size = 8,
#       groupColour = input$by_species,
#       groupFill = input$by_species
#     )
#   }

#   # Add geom_vlines for xvar and yvar
#   gate_xvar <- resultdf_reactive()$Gate[resultdf_reactive()$Marker == input$xvar]
#   gate_yvar <- resultdf_reactive()$Gate[resultdf_reactive()$Marker == input$yvar]

#   p <- p +
#     geom_vline(xintercept = gate_xvar, color = "red") +
#     geom_hline(yintercept = gate_yvar, color = "blue")

#   p
# })

# output$plot2 <- renderPlot({

#   dfplot2 <- subsetted()

#   # ...

#   gate_xvar <- input$gate_xvar
#   gate_yvar <- input$gate_yvar

#   p <- ggplot(dfplot2, aes(!!input$xvar,!!input$yvar)) + 
#     geom_vline(xintercept = gate_xvar, color = "red") +
#     geom_hline(yintercept = gate_yvar, color = "blue") 

#   # ...

#   return(p)
# })

# output$proportions <- renderPrint({
#   dfplot2 <- subsetted()

#   gate_xvar <- input$gate_xvar
#   gate_yvar <- input$gate_yvar

#   # ...

#   cat("Proportion of ++:", prop_pp, "\n")
#   cat("Proportion of -+:", prop_pm, "\n")
#   cat("Proportion of +-:", prop_mp, "\n")
#   cat("Proportion of --:", prop_mm, "\n")
# })

# observeEvent(input$update_gates, {
#   # Triggered when the button is clicked

#   # Recalculate proportions and redraw the plot
#   invalidateLater(0)

# })



  


  output$plot3 <- renderPlotly({
  dfplot3 <- subsetted()
  max_vals <- sapply(dfplot3, max)


  if (any(max_vals > 20)) {
    numeric_cols <- sapply(dfplot3, is.numeric)
    dfplot3[numeric_cols] <- log(dfplot3[numeric_cols])
  }

  xvar <- input$xvar
  yvar <- input$yvar

  gate_vliner <- dfplot3[[xvar]]
  gate_hliner <- dfplot3[[yvar]]


  result_to_plot_x <- skew_gate(gate_vliner)
  result_to_plot_y <- skew_gate(gate_hliner)

  p <- ggplot(dfplot3, aes(!!input$xvar, !!input$yvar)) +
    theme(legend.position = "bottom") +
    geom_point(aes(color = imageid)) 
   
  p <- p + 
    geom_vline(aes(xintercept = result_to_plot_x$cutoff), color = "red") +
    geom_hline(aes(yintercept = result_to_plot_y$cutoff), color = "red")

 ggplotly(p)

})


  
  
 output$marker_checkboxes <- renderUI({
    marker_data <- input$selected_columns
    
    if (is.null(marker_data) || length(marker_data) == 0) {
      return(NULL) # If no markers are selected or available, return NULL
    }

    switch_list <- lapply(marker_data, function(marker) {

      fluidRow(
        column(2, checkboxInput(paste0("checkbox_", marker), marker)),
        # p("", style = "margin-bottom: -5px;"),
        br(),
        column(2, materialSwitch(paste0("switch_", marker), label = NULL, status = "success"))
      )
    })
    
    tagList(switch_list)

  })



  # Track selected markers and their statuses
  selected_markers_phenotype <- reactive({
    markers <- input$selected_columns
    markers_statuses <- sapply(markers, function(marker) {
      if (input[[paste0("checkbox_", marker)]] == TRUE) {
        if (input[[paste0("switch_", marker)]] == TRUE) {
          return(paste0(marker, "+"))
        } else {
          return(paste0(marker, "-"))
        }
      } else {
        return(NULL)  # Return NULL for unchecked markers
      }
    })
    return(markers_statuses)
  })
  
  # Track phenotype name
  phenotype_name <- reactive({
    input$phenotype_name
  })

    # Define the phenotype based on selected markers and phenotype name
  definePhenotype <- function(markers, phenotype) {
    markers <- markers[!is.null(markers)]  # Remove NULL entries
    print(markers)
    if (length(markers) > 0) {
      phenotype_str <- paste(phenotype, paste(markers, collapse = ", "), sep = "\n")
      return(phenotype_str)
    } else {
      return("")  # Return empty string if no markers are selected
    }
  }


  # Function to clean column names by removing extra characters
cleanColumnNames <- function(dataframe) {
  names(dataframe) <- gsub("\\.+|\\s+", "", names(dataframe))
  return(dataframe)
}
  

phenotype_wfl_reactive <- reactiveVal(NULL)
# Define the reactive dataframe for phenotypes
phenotype_df <- reactiveVal(data.frame(phenotype = character(), markers = character()))


# Update the reactive dataframe with the defined phenotype
updatePhenotypeDF <- function(phenotype, markers) {
  current_df <- phenotype_df()
  updated_df <- rbind(current_df, data.frame(phenotype = phenotype, markers = markers))

   # Clean column names before updating the reactive dataframe
  updated_df <- cleanColumnNames(updated_df)


  phenotype_df(updated_df)
}

# Function to update the phenotype dataframe
updatePhenotypeDF2 <- function(phenotype, markers) {
  current_df <- phenotype_df()
  
  # Split the concatenated string by newline and clean the values
  cleaned_markers <- unlist(strsplit(markers, "\n"))
  cleaned_markers <- gsub("NULL,?\\s*", "", cleaned_markers)  # Remove "NULL"
  cleaned_markers <- gsub("^\\s+|\\s+$", "", cleaned_markers)  # Trim leading/trailing spaces
  
  # Bind the phenotype and markers into the dataframe
  updated_df <- rbind(current_df, data.frame(phenotype = phenotype, markers = paste(cleaned_markers, collapse = ", ")))
  phenotype_df(updated_df)
}

updatePhenotypeDF3 <- function(phenotype_name, phenotype_markers) {
  current_df <- phenotype_df()
  phenotype <- unlist(strsplit(phenotype_markers, "\n"))[2] # Extract markers excluding the phenotype name
  updated_df <- rbind(current_df, data.frame(phenotype = phenotype_name, markers = phenotype))
  phenotype_df(updated_df)
}



  cleanPhenotype <- function(phenotype_output) {
  # Split the output string by newline character
  phenotype_list <- strsplit(phenotype_output, "\n")[[1]]
  
  # Remove the NULL values and extract the relevant information
  cleaned_phenotype <- lapply(phenotype_list, function(phenotype) {
    parts <- unlist(strsplit(phenotype, ", "))  # Split by ", "
    cleaned_parts <- parts[parts != "NULL"]     # Remove "NULL"
    cleaned_phenotype <- paste(cleaned_parts, collapse = ", ")  # Recreate the string
    return(cleaned_phenotype)
  })
  
  return(cleaned_phenotype)
}

cleanPhenotype2 <- function(phenotype_output) {
  # Split the output string by newline character
  phenotype_list <- strsplit(phenotype_output, "\n")[[1]]
  
  # Remove the NULL values and extract the relevant information
  cleaned_phenotype <- lapply(phenotype_list, function(phenotype) {
    parts <- unlist(strsplit(phenotype, ", "))  # Split by ", "
    cleaned_parts <- parts[parts != "NULL"]     # Remove "NULL"
    cleaned_parts <- parts[-1]  # Exclude the first part before the first comma
    cleaned_phenotype <- paste(cleaned_parts, collapse = ", ")  # Recreate the string
    return(cleaned_phenotype)
  })
  
  return(cleaned_phenotype)
}

cleanPhenotype3 <- function(phenotype_output) {
  # Split the output string by newline character
  phenotype_list <- strsplit(phenotype_output, "\n")[[1]]
  
  # Remove the NULL values and extract the relevant information
  cleaned_phenotype <- lapply(phenotype_list, function(phenotype) {
    parts <- unlist(strsplit(phenotype, ", "))  # Split by ", "
    cleaned_parts <- parts[parts != "NULL"]     # Remove "NULL"
    cleaned_phenotype <- paste(cleaned_parts, collapse = ", ")  # Recreate the string
    return(cleaned_phenotype)
  })
  
  return(cleaned_phenotype)
}


# # Function to subset the dataframe based on phenotype-markers
# subsetAndCount <- function(df, phenotype_df) {
#   if (is.null(phenotype_df())) {
#     return(NULL)
#   }
  
#   counts <- lapply(phenotype_df()$markers, function(markers) {
#     subset_df <- df[rowSums(df[, markers]) == length(markers), ]
#     nrow(subset_df)
#   })
  
#   return(counts)
# }

# # Function to subset df_positivity based on markers and return counts
# subsetAndCount <- function(df, phenotype_df) {
#   phenotype_counts <- lapply(1:nrow(phenotype_df), function(i) {
#     markers <- unlist(strsplit(phenotype_df[i, "markers"], ", "))
#     markers <- markers[markers != ""]  # Remove empty elements
#     subset_df <- df[rowSums(df[, paste0(markers, "_positivity")] == "+") == length(markers), ]
#     nrow(subset_df)
#   })
  
#   names(phenotype_counts) <- phenotype_df$phenotype
#   return(phenotype_counts)
# }

# subsetAndCount <- function(df, phenotype_df) {
#   phenotype_counts <- lapply(1:nrow(phenotype_df), function(i) {
#     markers <- unlist(strsplit(phenotype_df[i, "markers"], ", "))
#     markers <- markers[markers != ""]  # Remove empty elements
#     # Construct the column names for positivity markers
#     marker_columns <- paste0(markers, "_positivity")
    
#     # Subset the dataframe
#     subset_df <- df[rowSums(df[, marker_columns] == "+") == length(markers), ]
#     nrow(subset_df)
#   })
  
#   names(phenotype_counts) <- phenotype_df$phenotype
#   return(phenotype_counts)
# }
# Subset rows where specific columns have values equal to "+"
# subset_data <- subset(df, KRT14_positivity == "+" & CD207_positivity == "+" & SOX10_positivity == "+")



# subsetAndCount <- function(df, phenotype_df) {
#   phenotype_counts <- lapply(1:nrow(phenotype_df), function(i) {
#     markers <- unlist(strsplit(phenotype_df[i, "markers"], ", "))
#     markers <- markers[markers != ""]  # Remove empty elements
#     condition <- ""
#     for (marker in markers) {
#       if (grepl("\\+", marker)) {
#         condition <- paste0(condition, "& ", marker)
#       } else if (grepl("-", marker)) {
#         condition <- paste0(condition, "& ", marker, "== '-'")
#       }
#     }
#     condition <- substr(condition, 3, nchar(condition))  # Remove the leading "& "
#     subset_df <- subset(df, eval(parse(text = condition)))
#     nrow(subset_df)
#   })
  
#   names(phenotype_counts) <- phenotype_df$phenotype
#   return(phenotype_counts)
# }

# subsetAndCount <- function(df, phenotype_df) {
#   phenotype_counts <- lapply(1:nrow(phenotype_df), function(i) {
#     markers <- unlist(strsplit(phenotype_df[i, "markers"], ", "))
#     markers <- markers[markers != ""]  # Remove empty elements
#     positive_markers <- markers[grepl("\\+", markers)]
#     negative_markers <- markers[grepl("-", markers)]
    
#     conditions <- lapply(positive_markers, function(marker) {
#       marker <- gsub("[+,]", "", marker)
#       paste0(marker, "_positivity == '+'")
#     })
    
#     neg_conditions <- lapply(negative_markers, function(marker) {
#       marker <- gsub("[-,]", "", marker)
#       paste0(marker, "_positivity == '-'")
#     })
    
#     all_conditions <- c(conditions, neg_conditions)
#     combined_conditions <- paste(all_conditions, collapse = " & ")
    
#     subset_df <- subset(df, subset = eval(parse(text = combined_conditions)))
#     nrow(subset_df)
#   })
  
#   names(phenotype_counts) <- phenotype_df$phenotype
#   return(phenotype_counts)
# }

subsetAndCount <- function(df, phenotype_df) {
  original_df <- uploaded_df()
  phenotype_counts <- lapply(1:nrow(phenotype_df), function(i) {
    markers <- unlist(strsplit(phenotype_df[i, "markers"], ", "))
    markers <- markers[markers != ""]  # Remove empty elements
    positive_markers <- markers[grepl("\\+", markers)]
    negative_markers <- markers[grepl("-", markers)]
    
    conditions <- lapply(positive_markers, function(marker) {
      marker <- gsub("[+,]", "", marker)
      paste0(marker, "_positivity == '+'")
    })
    
    neg_conditions <- lapply(negative_markers, function(marker) {
      marker <- gsub("[-,]", "", marker)
      paste0(marker, "_positivity == '-'")
    })
    
    all_conditions <- c(conditions, neg_conditions)
    combined_conditions <- paste(all_conditions, collapse = " & ")

    print("combined_conditions below")
    print(combined_conditions)
    
    subset_df <- subset(df, subset = eval(parse(text = combined_conditions))) 
    
    #  Assign phenotype to the filtered cells in the original dataset
    if (nrow(subset_df) > 0) {
      original_df[rownames(subset_df), 'phenotype'] <- phenotype_df[i, 'phenotype']
      original_df[, 'phenotype'] <- ifelse(is.na(original_df[, 'phenotype']), 'other', original_df[, 'phenotype'])
    }

    uploaded_df(original_df)
    nrow(subset_df)
  })

  # Find rows that don't satisfy any defined conditions and count them as 'other'
  other_count <- nrow(df) - sum(unlist(phenotype_counts))
  phenotype_counts <- c(phenotype_counts, other = other_count)
  
  names(phenotype_counts) <- c(phenotype_df$phenotype, 'other')
  print(phenotype_counts)
  return(phenotype_counts)
}

# subsetAndCount <- function(df, phenotype_df, any_switch) {
#   original_df <- uploaded_df()
#   phenotype_counts <- lapply(1:nrow(phenotype_df), function(i) {
#     markers <- unlist(strsplit(phenotype_df[i, "markers"], ", "))
#     markers <- markers[markers != ""]  # Remove empty elements
#     positive_markers <- markers[grepl("\\+", markers)]
#     negative_markers <- markers[grepl("-", markers)]
    
#     # Modify the condition based on the 'ANY' switch
#     if (any_switch) {
#       positive_markers <- paste0(positive_markers, "_positivity == '+'")
#       negative_markers <- paste0(negative_markers, "_positivity == '-'")
#       combined_conditions <- paste(c(positive_markers, negative_markers), collapse = " | ")
#     } else {
#       conditions <- lapply(positive_markers, function(marker) {
#         marker <- gsub("[+,]", "", marker)
#         paste0(marker, "_positivity == '+'")
#       })
      
#       neg_conditions <- lapply(negative_markers, function(marker) {
#         marker <- gsub("[-,]", "", marker)
#         paste0(marker, "_positivity == '-'")
#       })
      
#       all_conditions <- c(conditions, neg_conditions)
#       combined_conditions <- paste(all_conditions, collapse = " & ")
#     }

#     print("combined_conditions below")
#     print(combined_conditions)
    
#     subset_df <- subset(df, subset = eval(parse(text = combined_conditions))) 
    
#     #  Assign phenotype to the filtered cells in the original dataset
#     if (nrow(subset_df) > 0) {
#       original_df[rownames(subset_df), 'phenotype'] <- phenotype_df[i, 'phenotype']
#       original_df[, 'phenotype'] <- ifelse(is.na(original_df[, 'phenotype']), 'other', original_df[, 'phenotype'])
#     }

#     uploaded_df(original_df)
#     nrow(subset_df)
#   })

#   # Find rows that don't satisfy any defined conditions and count them as 'other'
#   other_count <- nrow(df) - sum(unlist(phenotype_counts))
#   phenotype_counts <- c(phenotype_counts, other = other_count)
  
#   names(phenotype_counts) <- c(phenotype_df$phenotype, 'other')
#   print(phenotype_counts)
#   return(phenotype_counts)
# }







# observe({
#   counts <- subsetAndCount(df_positivity, phenotype_df)
#   print(counts)
#   # Here 'counts' will contain the counts of rows for each defined phenotype
# })

# Subset rows where specific columns have values equal to "+"
# subset_data <- subset(df, KRT14_positivity == "+" & CD207_positivity == "+" & SOX10_positivity == "+")


    # Render the phenotype string when the Define Phenotype button is clicked
  observeEvent(input$define_phenotype, {
    if (input$define_phenotype > 0) {
      print(selected_markers_phenotype())
      phenotype <- definePhenotype(selected_markers_phenotype(), phenotype_name())
      print(phenotype)
      print(class(phenotype))

      phenotype_list_marker_indicators <- phenotype

      #remove first occurrence of "e" from vector
      # phenotype <- str_remove_all(phenotype, "NULL")

      phenotype <- cleanPhenotype3(phenotype)

      print(phenotype)

      output_text_phenotype <- paste(phenotype, collapse = "\n")

      print(output_text_phenotype)


      # output$phenotype_output <- renderText({
      #   return(output_text_phenotype)
      # })
    }
      # Update the reactive dataframe with the defined phenotype
      # updatePhenotypeDF(phenotype)

    # Update the reactive dataframe with the defined phenotype and markers
    updatePhenotypeDF3(input$phenotype_name, output_text_phenotype)

    print(input$phenotype_name)

    print(phenotype)

    print("my phenotype df is below")

      print(phenotype_df())
      print(phenotype_df()$markers)
  })


   observeEvent(input$phen_wfl, {
     phenotype_wfl_reactive(read.csv(input$phen_wfl$datapath))
     print(phenotype_wfl_reactive())
  })

  output$phenotype_output <- renderText({
  return("")  # Placeholder or empty text
})


# Observe the Fetch Subsets button click
  observeEvent(input$define_phenotype_AUTO, {
    print(colnames(uploaded_df()))
      counts <- subsetAndCount(uploaded_df(), phenotype_df())
    print(counts)
  # Here 'counts' will contain the counts of rows for each defined phenotype
    bar_data <- data.frame(
          phenotype = names(counts),
          count = unlist(counts)
        )

  #  barchart <- ggplot(bar_data, aes(x = phenotype, y = count, fill=phenotype)) +
  #     geom_bar(stat = "identity", color = "black") +
  #     labs(title = "Phenotype Counts", x = "Phenotype", y = "Count") +
  #     theme_minimal() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
        # Reshape data to wide format for stacked bar chart
    # bar_data_wide <- spread(bar_data, key = phenotype, value = count)
    
    # barchart <- ggplot(bar_data_wide, aes(x = 1, y = value, fill = variable)) +
    #   geom_col(position = "stack") +
    #   labs(title = "Phenotype Composition", x = "", y = "Count") +
    #   theme_minimal() +
    #   scale_fill_brewer(palette = "Set3") +  # Use color palette
    #   guides(fill = guide_legend(title = "Phenotype")) 

          # Assuming 'bar_data' has 'phenotype' and 'count' columns
      total_counts <- sum(bar_data$count)

      # Calculate percentages
      bar_data$percentage <- (bar_data$count / total_counts) * 100

      print(bar_data)

      # Plot stacked bar chart
      barchart <- ggplot(bar_data, aes(fill = phenotype, x = "", y = percentage)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(x = NULL, y = "Percentage") +
        ggtitle("Phenotype Composition") +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.text  = element_text(size = 14),
          axis.title = element_text(size = 15),
          legend.text  = element_text(size = 13),
          legend.title = element_text(size = 14)
        )
        # barchart <- ggplot(bar_data, aes(x = phenotype, y = count, fill = phenotype)) +
        # geom_col() 

      # geom_bar(stat = "identity", position = "stack") +
      # labs(title = "Phenotype Composition", x = "Phenotype", y = "Count") +
      # theme_minimal() +
      # theme(axis.text.x = element_text(angle = 45, hjust = 1))

      # bar_data %>% pivot_longer(cols = !bar_data, names_to = "phenotype", values_to = "percentage") %>% 
      #   ggplot(aes(fill = phenotype, x = bar_data, y = percentage)) + 
      #   geom_bar(position = "stack", stat = "identity")

        population_cells <- uploaded_df()
        species <- population_cells$phenotype

        MLEP <- MLEp(abundance(species))

         output$post_statistics <- renderText({
            # Create a message for statistics
            message <- "Partition Diversity Estimate:"
            mle_statistics <- round(MLEP, 5)

            # Combine the message and MLE statistics
            paste(message, mle_statistics, sep = "\n")  # Use <br> to insert a line break

          }) 

        # bootstrap


        print(round(MLEP, 5))

        # print(species)
        # MLEestimate <- as.numeric(MLEp.bsci(as.numeric(species),rounds = 100, level = 0.99))

        # print("MLEestimate")
        # print(MLEestimate)





  
  

      output$pheno_bar <- renderPlot({
              barchart
        #       bar_data %>% pivot_longer(cols = !bar_data, names_to = "phenotype", values_to = "percentage") %>% 
        # ggplot(aes(fill = phenotype, x = bar_data, y = percentage)) + 
        # geom_bar(position = "stack", stat = "identity")
            })


  })

  output$phenotypeTable <- renderTable({
    phenotype_df()
  })


  #  observeEvent(input$define_phenotype_AUTO, {
  #   # Load the phenotype indicators from the CSV file
  #   phenotype_indicators <- phenotype_wfl_reactive()

  #   df_positivity <- uploaded_df()

  #   df_rules <- phenotype_wfl_reactive()

  #         # Loop through rows in df_positivity
  #     for (row in 1:nrow(df_positivity)) {
  #       # # Extract the marker positivity for this row
  #       # markers_positivity <- df_positivity[row]  # Exclude non-marker columns like imageid, CELLID
  #       markers_positivity <- df_positivity[row]  # Exclude non-marker columns like imageid, CELLID

        
  #       # Loop through rows in df_rules (phenotype definitions)
  #       for (i in 1:nrow(df_rules)) {
  #         phenotype <- df_rules[i, ]

  #         print(phenotype)
          
  #         # Filter markers based on the current phenotype's workflow
  #         # markers <- relevant_markers
          
  #         # Check and apply the rules
  #         if (phenotype$phenotype_workflow == "allpos" && all(markers_positivity[markers] == "+")) {
  #           detected_phenotype <- phenotype$phenotype_name
  #           print(detected_phenotype)
  #           break
  #         } else if (phenotype$phenotype_workflow == "allneg" && all(markers_positivity[markers] == "-")) {
  #           detected_phenotype <- phenotype$phenotype_name
  #           break
  #         } else if (phenotype$phenotype_workflow == "anypos" && any(markers_positivity[markers] == "+")) {
  #           detected_phenotype <- phenotype$phenotype_name
  #           break
  #         } else if (phenotype$phenotype_workflow == "anyneg" && any(markers_positivity[markers] == "-")) {
  #           detected_phenotype <- phenotype$phenotype_name
  #           break
  #         } else if (phenotype$phenotype_workflow == "pos" && any(markers_positivity[markers] == "+")) {
  #           detected_phenotype <- phenotype$phenotype_name
  #           break
  #         } else if (phenotype$phenotype_workflow == "neg" && any(markers_positivity[markers] == "-")) {
  #           detected_phenotype <- phenotype$phenotype_name
  #           break
  #         }
  #       }
        
  #       # Assign the detected phenotype to the corresponding row in df_positivity
  #       df_positivity[row, "Phenotype"] <- ifelse(is.null(detected_phenotype), "other", detected_phenotype)
  #     }


  #     uploaded_df(df_positivity)
    
  # })


     
   # Update Phenotype Table when the file is uploaded
  observeEvent(input$phen_wfl, {
    inFile <- input$phen_wfl

    df <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
    
    # Check for valid columns in the uploaded file
    if ("phenotype" %in% names(df) && "markers" %in% names(df)) {
      # phenotype_df_LOADED <<- bind_rows(phenotype_df_LOADED, df) %>%
      #   distinct()  # Remove duplicate rows if any
      
      phenotype_df(df)
    }
    
  output$phenotypeTable <- renderTable({
    phenotype_df()
  })

  }
  )
 




  
# In the UI, use conditionalPanel to conditionally render the plot
conditionalPanel(
  condition = "input.patient_number_plot2.length > 0",
  plotOutput("plot2")
)  

conditionalPanel(
  condition = "input.patient_number_plot2.length > 0",
  plotlyOutput("plot3")
)  



      })

shinyApp(ui = ui, server = server)
