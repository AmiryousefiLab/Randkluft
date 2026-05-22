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

df_example = read.csv(file.path("data", "exemplar-001--unmicst_cell.csv"))

ground_truth_gates_loaded = read.csv(file.path("data", "tuulia_data_GT.csv"))

df_example_PHENOTYPE = read.csv(file.path("data", "phenotype_table_help.csv"))

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


# options(shiny.maxRequestSize=100*1024^2)
options(shiny.maxRequestSize=100*1024^3)
server <- shinyServer(function(input, output, session) {

  remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
  low_threshold <- quantile(target, low_percentile / 100)
  high_threshold <- quantile(target, high_percentile / 100)
  target <- target[target >= low_threshold & target <= high_threshold]
  return(target) }

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
      # ggplot2 3.5+: geom_vline(xintercept = numeric(0)) is a hard error
      if (length(GroundTruthInput) > 0 && !is.na(GroundTruthInput[1])) {
        plot <- plot + geom_vline(xintercept = GroundTruthInput, color = "blue")
      }
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

  # outputinterceptreactive change automatically invalidates the renderPlot
  # above -- no direct generatePlot() call needed here (calling it from an
  # observer tries to open a system graphics device and fails on servers).

  selected_columns_man <- input$selected_columns
  selected_patients_man <- unique_patients_gating()

  chosen_patient <- selected_patients_man[current_patient()]
  chosen_marker  <- selected_columns_man[current_marker()]

  updated_df <- csv_save_file_path()
  if (!is.null(updated_df)) {
    updated_df$Gate[
      updated_df$Marker == chosen_marker &
      updated_df$Patient == chosen_patient
    ] <- new_gate_value
    csv_save_file_path(updated_df)
    resultdf_reactive(updated_df)
  }

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

    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "Cell", "Row", "X", "Y", "DNA2", "Hoechst1", "Hoechst2", "Hoechst3",
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
    runjs("document.getElementById('message_loading_data').style.display = 'none';")

  })

# Define a reactiveVal for plot height
#   plotHeight <- reactiveVal("2060px")


  reactive_markers_selected <- reactiveVal(NULL)

  uploaded_df <- reactiveVal(NULL)
# Define a reactive value for the GMM gate switch
  gmm_gate_switch <- reactiveVal(TRUE)  # Initialize with the desired default value
  # You can set the initial value based on the user's preference or the default value

  histplot_gate_switch <- reactiveVal(TRUE)  # Initialize with the desired default value
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
    shinyjs::hide("gated_histogram_on_page")
    shinyjs::hide("patient_number")
    shinyjs::hide("histogram_plot")

    shinyjs::hide("summary_output")
    shinyjs::hide("image_garage_output")
    shinyjs::hide("selected_columns_phenotyping")


    updateSelectInput(session, "yvar", choices = character(0), selected = character(0))
    updateSelectInput(session, "xvar", choices = character(0), selected = character(0))
    updateSelectInput(session, "marker", choices = character(0), selected = character(0))


  })

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


    for(i in seq_along(seq_length)){
      progress$set(value = i)
      if(i == no_batches) chunk_size = -1 else chunk_size = seq_length[i+1] - seq_length[i]

      df_temp = read.csv.sql(file_path, sql = query, dbname = tempfile())
      colnames(df_temp) = col_names
      df = rbind(df,df_temp)
    }

    progress$close()
    return(df)
  }


  observeEvent(input$cell_file, {


    n_rows = length(count.fields(input$cell_file$datapath))
    csv_columns <- names(read.csv(input$cell_file$datapath, nrows = 0, check.names = FALSE))


    # df_out = read_batch_with_progress(input$file1$datapath,n_rows,10)
    # Generate a sample of 10 numbers from 1 to 20 without replacement
    if (n_rows > 70000) {
      sample_numbers <- sample(1:n_rows, 70000, replace = FALSE)
    } else {
      sample_numbers <- seq(1:n_rows)
    }
    if ("CellID" %in% csv_columns) {
      # Create the SQL query string
      query <- sprintf("SELECT * FROM file WHERE CellID IN (%s)", paste(sample_numbers, collapse = ", ")) #the first row should be called CellID
      uploaded_df(read.csv.sql(input$cell_file$datapath, sql = query, dbname = tempfile()))
    } else {
      uploaded_intermediate <- read.csv(input$cell_file$datapath, header = TRUE, check.names = FALSE)
      if (nrow(uploaded_intermediate) > 70000) {
        uploaded_intermediate <- uploaded_intermediate[sample(seq_len(nrow(uploaded_intermediate)), 70000), , drop = FALSE]
      }
      uploaded_df(uploaded_intermediate)
    }


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
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "Cell", "Row", "X", "Y", "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
    "MinorAxisLength", "Extent", "Solidity", "Orientation", "", "DNA6a")  # List the columns to exclude
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

      uploaded_intermediate$imageid <- rep("sample", nrow(uploaded_intermediate))
      uploaded_df(uploaded_intermediate)
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

      if (!"X_centroid" %in% colnames(uploaded_intermediate)) {
        grid_width <- max(1, ceiling(sqrt(nrow(uploaded_intermediate))))
        uploaded_intermediate$X_centroid <- ((seq_len(nrow(uploaded_intermediate)) - 1) %% grid_width) + 1
        uploaded_df(uploaded_intermediate)
      }

      if (!"Y_centroid" %in% colnames(uploaded_intermediate)) {
        grid_width <- max(1, ceiling(sqrt(nrow(uploaded_intermediate))))
        uploaded_intermediate$Y_centroid <- ((seq_len(nrow(uploaded_intermediate)) - 1) %/% grid_width) + 1
        uploaded_df(uploaded_intermediate)
      }


          print("DONE LOADING CSV FINAL after POST")


          print(colnames(uploaded_df()))

    showNotification("File Ready for Use", duration = 10, id = "message")


    # Show the checkboxGroupInput
    shinyjs::show("selected_columns")
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

    unique_patients <- unique_patients_gating()

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


    if (length(selected_columns)==1){
      print("hi")
      print(colnames(sub_data))
      resultdf_to_save = get_gates_csv_single(sub_data, csv_save_file_name)
    }
    else {
      resultdf_to_save = get_gates_csv(sub_data, csv_save_file_name)
    }
    resultdf_reactive(resultdf_to_save)


    print(resultdf_to_save)

    csv_save_file_path(resultdf_to_save)

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

get_density_3d <- function(x, y, z, grid_size = 24) {
  marker_values <- cbind(
    suppressWarnings(as.numeric(x)),
    suppressWarnings(as.numeric(y)),
    suppressWarnings(as.numeric(z))
  )
  density_values <- rep(0, nrow(marker_values))
  finite_rows <- complete.cases(marker_values) & apply(marker_values, 1, function(row) all(is.finite(row)))
  marker_values <- marker_values[finite_rows, , drop = FALSE]

  if (!requireNamespace("ks", quietly = TRUE) ||
      nrow(marker_values) < 3 ||
      any(apply(marker_values, 2, function(values) diff(range(values)) == 0))) {
    return(density_values)
  }

  density_grid <- tryCatch(
    ks::kde(
      x = marker_values,
      H = ks::Hns(marker_values),
      binned = TRUE,
      bgridsize = rep(grid_size, 3),
      compute.cont = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(density_grid)) {
    return(density_values)
  }

  grid_indices <- lapply(seq_len(3), function(axis) {
    grid_axis <- density_grid$eval.points[[axis]]
    pmax(1, pmin(length(grid_axis), findInterval(marker_values[, axis], grid_axis)))
  })
  density_values[finite_rows] <- density_grid$estimate[
    cbind(grid_indices[[1]], grid_indices[[2]], grid_indices[[3]])
  ]
  density_values
}


  # Always register at server level so reactive to outputinterceptreactive()
  # at all times; visibility is controlled by the toggle below.
  output$gated_histogram_on_page <- renderPlot({
    req(uploaded_df(), resultdf_reactive(), histplot_gate_switch())
    generatePlot()
  })

  observeEvent(input$gen_hist_plots_on_off, {
    if (histplot_gate_switch()) {
      shinyjs::show("gated_histogram_on_page")
    } else {
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


      selected_columns_man <- input$selected_columns
      selected_patients_man <- unique_patients_gating()

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

histogram_plot <- make_marker_histogram_plot(df, chosen_patient, chosen_marker, gate_value, title_size = 18)

temp_df <- filtered_data_xy
temp_df[[chosen_marker]][temp_df[[chosen_marker]] < gate_value] <- 0

density_values <- get_density(filtered_data_xy$X_centroid, filtered_data_xy$Y_centroid, n = 100)

pos_subset_gp <- filtered_data_xy[filtered_data_xy[[chosen_marker]] > gate_value, ]
contour_plot <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
  geom_point(
    aes(color = ifelse(filtered_data_xy[[chosen_marker]] > gate_value,
                       density_values, 0)),
    size = 0.3,
    alpha = 0.35
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
if (nrow(pos_subset_gp) >= 2) {
  contour_plot <- contour_plot + geom_density_2d(data = pos_subset_gp, color = "black")
}
overlay_plot2 <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
  geom_point(aes(color = ifelse(filtered_data_xy[[chosen_marker]] > gate_value, "Positive", "Negative")), size = 0.3, alpha = 0.35) +
  scale_color_manual(
    guide = guide_legend(title = "", override.aes = list(size = 7, alpha = 1)),
    values = c("Positive" = "red", "Negative" = "grey")
  ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16, family = "Arial"),
        legend.title = element_text(size = 16, face = "bold", family = "Arial"),
        legend.key.size = grid::unit(1.2, "cm"),
        legend.key.width = grid::unit(1.4, "cm"),
        legend.key.height = grid::unit(0.9, "cm"),
        legend.spacing.x = grid::unit(0.35, "cm"),
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

  runjs("document.getElementById('message_gen_hist').style.display = 'none';")
  shinyjs::show("gated_histogram_on_page")

}


  # downloadall: one arranged PDF page per selected marker
  output$downloadall <- downloadHandler(
    filename = function() { "all_marker_plots.pdf" },
    contentType = "application/pdf",
    content = function(file) {
      tryCatch({
        data <- uploaded_df()
        if (is.null(data) || nrow(data) == 0) {
          return(write_message_pdf(file, "No plots available", "Upload data and run Randkluft before downloading plots."))
        }

        selected_columns_save <- input$selected_columns
        if (is.null(selected_columns_save) || length(selected_columns_save) == 0) {
          selected_columns_save <- setdiff(names(data), c("imageid", "X_centroid", "Y_centroid"))
        }

        if (!"imageid" %in% names(data)) {
          data$imageid <- "sample"
        }

        selected_patients_save <- unique_patients_gating()
        selected_patients_save <- selected_patients_save[selected_patients_save %in% unique(data$imageid)]
        if (length(selected_patients_save) == 0) {
          selected_patients_save <- unique(data$imageid)
        }

        sub_data <- data[data$imageid %in% selected_patients_save, , drop = FALSE]

        generate_histogram_pdf(
          subdata_to_plot = sub_data,
          pdf_file_name = file,
          markers = selected_columns_save,
          patient_ids = selected_patients_save
        )
      }, error = function(e) {
        write_message_pdf(file, "Could not generate all plots", conditionMessage(e))
      })
    }
  )

  # downloadcurrent: one arranged PDF page for the currently selected marker
  output$downloadcurrent <- downloadHandler(
    filename = function() { "current_plot.pdf" },
    contentType = "application/pdf",
    content = function(file) {
      tryCatch({
        data <- uploaded_df()
        if (is.null(data) || nrow(data) == 0) {
          return(write_message_pdf(file, "No current plot", "Upload data and run Randkluft before downloading the current plot."))
        }

        selected_columns_current <- input$selected_columns
        if (is.null(selected_columns_current) || length(selected_columns_current) == 0) {
          return(write_message_pdf(file, "No current plot", "Select at least one marker before downloading the current plot."))
        }

        current_marker_index <- min(current_marker(), length(selected_columns_current))
        selected_marker <- selected_columns_current[current_marker_index]

        if (!"imageid" %in% names(data)) {
          data$imageid <- "sample"
        }

        selected_patients_current <- unique_patients_gating()
        selected_patients_current <- selected_patients_current[selected_patients_current %in% unique(data$imageid)]
        if (length(selected_patients_current) == 0) {
          selected_patients_current <- unique(data$imageid)
        }
        current_patient_index <- min(current_patient(), length(selected_patients_current))
        selected_patient <- selected_patients_current[current_patient_index]

        sub_data <- data[data$imageid == selected_patient, , drop = FALSE]

        generate_histogram_pdf(
          subdata_to_plot = sub_data,
          pdf_file_name = file,
          markers = selected_marker,
          patient_ids = selected_patient
        )
      }, error = function(e) {
        write_message_pdf(file, "Could not generate current plot", conditionMessage(e))
      })
    }
  )


make_pdf_message_grob <- function(title, message) {
  gridExtra::arrangeGrob(
    grid::textGrob(
      message,
      gp = grid::gpar(fontsize = 14, col = "grey25"),
      x = 0.5,
      y = 0.5
    ),
    top = grid::textGrob(title, gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )
}

write_message_pdf <- function(file, title, message, width = 10, height = 10) {
  pdf(file, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  grid::grid.draw(make_pdf_message_grob(title, message))
  invisible(file)
}

resolve_gate_value <- function(patient_id, marker, data) {
  gates <- csv_save_file_path()
  if (is.null(gates)) {
    gates <- resultdf_reactive()
  }

  gate_value <- numeric(0)
  if (!is.null(gates) && all(c("Patient", "Marker", "Gate") %in% names(gates))) {
    gate_value <- gates$Gate[
      gates$Marker == marker &
        gates$Patient == patient_id
    ]
  }

  if (length(gate_value) == 0 || is.na(gate_value[1]) || !is.finite(gate_value[1])) {
    target <- suppressWarnings(as.numeric(data[data$imageid == patient_id, marker]))
    target <- target[is.finite(target)]
    gate_value <- tryCatch(skew_gate(target, 0.01)$cutoff, error = function(e) NA_real_)
  }

  gate_value <- as.numeric(gate_value[1])
  if (!is.finite(gate_value)) {
    gate_value <- NA_real_
  }
  gate_value
}

make_marker_histogram_plot <- function(data, patient_id, marker, gate_value, title_size = 15) {
  target <- suppressWarnings(as.numeric(data[data$imageid == patient_id, marker]))
  target <- target[is.finite(target)]

  if (length(target) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No finite marker values", size = 5) +
        theme_void() +
        labs(title = "Gate Histogram")
    )
  }

  target <- remove_outliers2(target)
  target <- target[is.finite(target)]
  plot_df <- data.frame(x = target)
  n_positive <- if (is.na(gate_value)) NA_integer_ else sum(target > gate_value)
  positive_rate <- if (is.na(gate_value)) NA_real_ else round(n_positive / length(target), 3)

  histogram_plot <- ggplot(plot_df, aes(x = x)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 80,
      fill = "#d8e9f7",
      color = "#34495e",
      size = 0.2
    ) +
    labs(
      title = "Gate Histogram",
      subtitle = if (is.na(gate_value)) {
        "Gate unavailable"
      } else {
        paste0("N+ = ", n_positive, "   +R = ", positive_rate, "   Gate = ", round(gate_value, 2))
      },
      x = marker,
      y = "Density"
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = title_size),
      plot.subtitle = element_text(hjust = 0.5, size = title_size),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    )

  if (length(unique(target)) > 1) {
    histogram_plot <- histogram_plot + geom_density(color = "#2c3e50", size = 0.6)
  }

  if (!is.na(gate_value)) {
    histogram_plot <- histogram_plot +
      geom_vline(xintercept = gate_value, color = "#c0392b", size = 0.8) +
      annotate(
        "label",
        x = gate_value,
        y = Inf,
        label = paste0("Gate: ", round(gate_value, 2)),
        hjust = -0.05,
        vjust = 1.2,
        size = 3.2,
        fill = "white",
        color = "#c0392b"
      )
  }

  gt_gate <- ground_truth_gates_loaded$Gate[
    ground_truth_gates_loaded$Patient == patient_id &
      ground_truth_gates_loaded$Marker == marker
  ]
  if (length(gt_gate) > 0 && !is.na(gt_gate[1])) {
    histogram_plot <- histogram_plot + geom_vline(xintercept = gt_gate[1], color = "#2c7fb8", size = 0.7)
  }

  histogram_plot
}

generate_four_panel_plot <- function(data, patient_id, marker) {
  filtered_data_xy <- data[data$imageid == patient_id, , drop = FALSE]
  if (nrow(filtered_data_xy) == 0 || !marker %in% names(filtered_data_xy)) {
    return(make_pdf_message_grob(marker, "No data available for this marker."))
  }

  if (!all(c("X_centroid", "Y_centroid") %in% names(filtered_data_xy))) {
    grid_width <- max(1, ceiling(sqrt(nrow(filtered_data_xy))))
    filtered_data_xy$X_centroid <- ((seq_len(nrow(filtered_data_xy)) - 1) %% grid_width) + 1
    filtered_data_xy$Y_centroid <- ((seq_len(nrow(filtered_data_xy)) - 1) %/% grid_width) + 1
  }

  filtered_data_xy[[marker]] <- suppressWarnings(as.numeric(filtered_data_xy[[marker]]))
  gate_value <- resolve_gate_value(patient_id, marker, data)
  histogram_plot <- make_marker_histogram_plot(data, patient_id, marker, gate_value)

  density_values <- tryCatch(
    get_density(filtered_data_xy$X_centroid, filtered_data_xy$Y_centroid, n = 100),
    error = function(e) rep(0, nrow(filtered_data_xy))
  )
  is_positive <- if (is.na(gate_value)) {
    rep(FALSE, nrow(filtered_data_xy))
  } else {
    filtered_data_xy[[marker]] > gate_value
  }
  is_positive[is.na(is_positive)] <- FALSE
  filtered_data_xy$pdf_density <- ifelse(is_positive, density_values, 0)
  filtered_data_xy$gate_status <- ifelse(is_positive, "Positive", "Negative")

  spatial_theme <- theme_minimal(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = grid::unit(0.9, "cm"),
      legend.key.width = grid::unit(1.1, "cm"),
      legend.key.height = grid::unit(0.75, "cm"),
      legend.spacing.x = grid::unit(0.35, "cm"),
      aspect.ratio = 1
    )

  digital_theme <- theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

  density_inside_legend_theme <- theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      color = "grey70"
    ),
    legend.key.height = grid::unit(12, "pt"),
    legend.key.width = grid::unit(6, "pt"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

  digrepresentation <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid, color = .data[[marker]])) +
    geom_point(shape = 20, size = 0.16, alpha = 0.5) +
    scale_color_gradient(low = "grey30", high = "white", guide = "none") +
    labs(title = "Digital Representation", x = "X Centroid", y = "Y Centroid") +
    digital_theme

  pos_subset <- filtered_data_xy[is_positive, , drop = FALSE]
  contour_plot <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid)) +
    geom_point(aes(color = pdf_density), size = 0.16, alpha = 0.35) +
    scale_color_viridis(option = "turbo", name = "Intensity") +
    labs(title = "Positive Density", x = "X Centroid", y = "Y Centroid") +
    spatial_theme +
    density_inside_legend_theme
  if (nrow(pos_subset) >= 2) {
    contour_plot <- contour_plot +
      geom_density_2d(data = pos_subset, aes(x = X_centroid, y = Y_centroid), inherit.aes = FALSE, color = "black", size = 0.35)
  }

  overlay_plot2 <- ggplot(filtered_data_xy, aes(x = X_centroid, y = Y_centroid, color = gate_status)) +
    geom_point(size = 0.16, alpha = 0.35) +
    scale_color_manual(
      guide = guide_legend(title = "", override.aes = list(size = 4.5, alpha = 1)),
      values = c("Positive" = "red", "Negative" = "grey")
    ) +
    labs(title = "Positive Cells", x = "X Centroid", y = "Y Centroid") +
    spatial_theme

  gridExtra::arrangeGrob(
    histogram_plot,
    digrepresentation,
    overlay_plot2,
    contour_plot,
    ncol = 2,
    top = grid::textGrob(
      paste0("Marker: ", marker),
      gp = grid::gpar(fontsize = 17, fontface = "bold")
    )
  )
}


generate_histogram_pdf <- function(subdata_to_plot, pdf_file_name, markers = NULL, patient_ids = NULL) {
  data <- subdata_to_plot
  if (!"imageid" %in% names(data)) {
    data$imageid <- "sample"
  }

  if (is.null(markers)) {
    markers <- setdiff(names(data), c("imageid", "X_centroid", "Y_centroid"))
  }
  markers <- markers[markers %in% names(data)]

  if (is.null(patient_ids)) {
    patient_ids <- unique(data$imageid)
  }
  patient_ids <- patient_ids[patient_ids %in% unique(data$imageid)]

  pdf(file = pdf_file_name, width = 10, height = 10, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  if (length(markers) == 0 || length(patient_ids) == 0) {
    grid::grid.draw(make_pdf_message_grob("No plots available", "Run Randkluft and select at least one marker before downloading all plots."))
    pdf_file_path(pdf_file_name)
    return(invisible(pdf_file_name))
  }

  first_page <- TRUE
  for (patient_id in patient_ids) {
    for (marker in markers) {
      page_grob <- tryCatch(
        generate_four_panel_plot(data, patient_id, marker),
        error = function(e) {
          make_pdf_message_grob(
            paste0("Marker: ", marker),
            paste("Could not generate this page:", conditionMessage(e))
          )
        }
      )
      if (first_page) {
        first_page <- FALSE
      } else {
        grid::grid.newpage()
      }
      grid::grid.draw(page_grob)
    }
  }

  pdf_file_path(pdf_file_name)
  invisible(pdf_file_name)
}


  output$icaAnalysis <- renderPlot({

    req(uploaded_df())
    colnames(uploaded_df())
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "Cell", "Row", "X", "Y", "DNA2", "Hoechst1", "Hoechst2", "Hoechst3",
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
        scale_color_discrete(guide = "none") +
        labs(title = "UMAP Plot of Selected Columns")
  })

  unique_patients <- reactive({
  req(uploaded_df())
  unique(uploaded_df()$imageid)
  })

  observe({
  shinyjs::show("patient_number")
  updateRadioButtons(session, "patient_number", choices = unique_patients())
})

 observe({
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
  shinyjs::show("precrev_marker")
  updateRadioButtons(session, "precrev_marker", choices = unique_markers_w_DNA())
})

  unique_markers_w_DNA <- reactive({
    req(uploaded_df())
    colnames(uploaded_df())
    # Define the column names you want to exclude
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "Cell", "Row", "X", "Y",
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
    columns_to_exclude <- c("imageid", "phenotype", "ROI_major_category", "CellID", "Cell", "Row", "X", "Y", "DNA2", "Hoechst1", "Hoechst2", "Hoechst3",
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

    exclude_pattern2 <- "_positivity"
    exclude_columns2 <- grep(exclude_pattern2, colnames(uploaded_df()), value = TRUE, ignore.case = TRUE)
    columns_to_exclude <- c(columns_to_exclude, exclude_columns2)

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
    unique_patients()
  })

  inputreactiveintercept_USERDEF <- reactive({
    input$intercept
  })


 selected_opacity <- reactive({
    input$opacity_slider
  })


  observe({
    updateSelectInput(session, "marker", choices = unique_markers())
})

 observe({
    updateSelectInput(session, "marker_im_garage", choices = unique_markers())
})


 observe({
    updateSelectInput(session, "xvar", choices = unique_markers())
})


 selected_var_name <- function(value) {
    if (is.null(value) || length(value) == 0) {
      return(NULL)
    }
    value <- as.character(value)
    if (length(value) == 0 || is.na(value[1]) || value[1] == "") {
      return(NULL)
    }
    value[1]
 }

 observe({
    markers <- unique_markers()
    current_yvar <- selected_var_name(input$yvar)
    selected_yvar <- if (!is.null(current_yvar) && current_yvar %in% markers) {
      current_yvar
    } else {
      markers[min(2, length(markers))]
    }
    updateSelectInput(session, "yvar", choices = markers, selected = selected_yvar)
})


 trivariate_marker_default <- function(markers, current_marker, position) {
    current_marker <- selected_var_name(current_marker)
    if (!is.null(current_marker) && current_marker %in% markers) {
      return(current_marker)
    }
    if (length(markers) == 0) {
      return(character(0))
    }
    if (length(markers) < 3) {
      return(markers[1])
    }
    markers[position]
 }

 observe({
    markers <- unique_markers()
    updateSelectInput(
      session,
      "xvarTri",
      choices = markers,
      selected = trivariate_marker_default(markers, input$xvarTri, 1)
    )
    updateSelectInput(
      session,
      "yvarTri",
      choices = markers,
      selected = trivariate_marker_default(markers, input$yvarTri, 2)
    )
    updateSelectInput(
      session,
      "zvar",
      choices = markers,
      selected = trivariate_marker_default(markers, input$zvar, 3)
    )
})


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


      find_mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }

        # Extract histogram data and find modes
        modes <- list()
        kldivs <- list()
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


          mean_KL_divergences <- mean(kl_df$KL_Divergence)


          normalized_value <- 1 - mean_KL_divergences

          first_cycle <- modes_df[1,2]

          # Get the minimum value from the second column of modes_df
          min_value_cycle <- min(modes_df[, 2], na.rm = TRUE)

          print(min_value_cycle)

          # Get the index of the minimum value from the second column of modes_df
          min_index <- which.min(modes_df[, 2])

          # Get the name of the row corresponding to the minimum value
          min_row_name <- rownames(modes_df)[min_index]

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

                })

 # Render the plot
        output$modes_plot_output <- renderPlot({
          regression_mode_react()
        })

         output$klPLOToutput <- renderPlot({
          klPlotReact()
        })


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
    uploaded_df(filtered_data_reactive())


  # showNotification("Data subset loaded for patient!", duration = NULL, id = "loadedsubset")
    # Display a completion message to the user
    shinyalert::shinyalert(
      title = "Subset Loaded",
      text = "You may use this subsetted data in the Randkluft Tab for analysis.",
      type = "success"
    )

})

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
    req(uploaded_df())
    uploaded_df() |> filter(imageid %in% unique_patients_gating())

  })


  subsetted_precrev_ts <- reactive({
    req(input$patients_ts)
    uploaded_df() |> filter(imageid %in% input$patients_ts)

  })

    subsetted_tri <- reactive({
    req(uploaded_df())
    uploaded_df() |> filter(imageid %in% unique_patients_gating())

  })

  react_xgate <- reactiveVal(NULL)
  react_ygate <- reactiveVal(NULL)

  resolve_bivariate_gate_value <- function(data, patient_id, marker, override = NULL) {
    override_value <- suppressWarnings(as.numeric(override[1]))
    if (!is.null(override) && length(override) > 0 && is.finite(override_value)) {
      return(override_value)
    }

    gate_value <- resolve_gate_value(patient_id, marker, data)
    if (!is.finite(gate_value)) {
      marker_values <- suppressWarnings(as.numeric(data[data$imageid == patient_id, marker]))
      marker_values <- marker_values[is.finite(marker_values)]
      gate_value <- if (length(marker_values) > 0) median(marker_values, na.rm = TRUE) else NA_real_
    }

    gate_value
  }

  prepare_bivariate_patient_data <- function(data, patient_id, markers) {
    if (!"imageid" %in% names(data)) {
      data$imageid <- "sample"
    }

    patient_data <- data[data$imageid == patient_id, , drop = FALSE]
    if (nrow(patient_data) == 0) {
      return(patient_data)
    }

    if (!all(c("X_centroid", "Y_centroid") %in% names(patient_data))) {
      grid_width <- max(1, ceiling(sqrt(nrow(patient_data))))
      patient_data$X_centroid <- ((seq_len(nrow(patient_data)) - 1) %% grid_width) + 1
      patient_data$Y_centroid <- ((seq_len(nrow(patient_data)) - 1) %/% grid_width) + 1
    }

    numeric_columns <- unique(c("X_centroid", "Y_centroid", markers))
    numeric_columns <- numeric_columns[numeric_columns %in% names(patient_data)]
    patient_data[numeric_columns] <- lapply(patient_data[numeric_columns], function(values) {
      suppressWarnings(as.numeric(values))
    })
    patient_data
  }

  precompute_bivariate_gates <- function(data, patient_id, markers, gate_overrides = list()) {
    markers <- unique(markers[markers %in% names(data)])
    stats::setNames(
      vapply(markers, function(marker) {
        resolve_bivariate_gate_value(data, patient_id, marker, gate_overrides[[marker]])
      }, numeric(1)),
      markers
    )
  }

  get_precomputed_bivariate_gate <- function(gate_values, marker, data, patient_id, override = NULL) {
    if (!is.null(gate_values) && marker %in% names(gate_values)) {
      return(as.numeric(gate_values[[marker]][1]))
    }
    resolve_bivariate_gate_value(data, patient_id, marker, override)
  }

  bivariate_pair_error <- function(xvar, yvar, message) {
    list(xvar = xvar, yvar = yvar, error = message)
  }

  prepare_bivariate_pair_data <- function(patient_data, patient_id, xvar, yvar,
                                          gate_values = NULL,
                                          gate_x_override = NULL,
                                          gate_y_override = NULL) {
    if (!all(c(xvar, yvar) %in% names(patient_data))) {
      return(bivariate_pair_error(
        xvar,
        yvar,
        "One or both selected markers are not available in the current data."
      ))
    }
    if (nrow(patient_data) == 0) {
      return(bivariate_pair_error(xvar, yvar, "No data available for this plot."))
    }

    complete_rows <- complete.cases(
      patient_data[[xvar]],
      patient_data[[yvar]],
      patient_data$X_centroid,
      patient_data$Y_centroid
    )
    finite_rows <- is.finite(patient_data[[xvar]]) &
      is.finite(patient_data[[yvar]]) &
      is.finite(patient_data$X_centroid) &
      is.finite(patient_data$Y_centroid)
    dfplot2 <- patient_data[complete_rows & finite_rows, , drop = FALSE]
    if (nrow(dfplot2) == 0) {
      return(bivariate_pair_error(xvar, yvar, "No finite values available for this plot."))
    }

    gate_xvar <- get_precomputed_bivariate_gate(
      gate_values, xvar, patient_data, patient_id, gate_x_override
    )
    gate_yvar <- get_precomputed_bivariate_gate(
      gate_values, yvar, patient_data, patient_id, gate_y_override
    )
    if (xvar == yvar) {
      if (!is.null(gate_x_override)) {
        gate_yvar <- gate_xvar
      } else if (!is.null(gate_y_override)) {
        gate_xvar <- gate_yvar
      }
    }

    is_x_pos <- if (is.finite(gate_xvar)) dfplot2[[xvar]] > gate_xvar else rep(FALSE, nrow(dfplot2))
    is_y_pos <- if (is.finite(gate_yvar)) dfplot2[[yvar]] > gate_yvar else rep(FALSE, nrow(dfplot2))
    is_x_pos[is.na(is_x_pos)] <- FALSE
    is_y_pos[is.na(is_y_pos)] <- FALSE

    dfplot2$bivariate_density <- tryCatch(
      get_density(dfplot2[[xvar]], dfplot2[[yvar]], n = 100),
      error = function(e) rep(0, nrow(dfplot2))
    )
    dfplot2$bivariate_gate_status <- dplyr::case_when(
      is_x_pos & is_y_pos ~ "+/+",
      is_x_pos & !is_y_pos ~ "+/-",
      !is_x_pos & is_y_pos ~ "-/+",
      TRUE ~ "-/-"
    )
    dfplot2$x_status <- ifelse(is_x_pos, "Positive", "Negative")
    dfplot2$y_status <- ifelse(is_y_pos, "Positive", "Negative")

    x_range <- range(dfplot2[[xvar]], na.rm = TRUE)
    y_range <- range(dfplot2[[yvar]], na.rm = TRUE)
    x_pad <- diff(x_range) * 0.04
    y_pad <- diff(y_range) * 0.04
    if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
    if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

    list(
      data = dfplot2,
      patient_id = patient_id,
      xvar = xvar,
      yvar = yvar,
      gate_xvar = gate_xvar,
      gate_yvar = gate_yvar,
      x_positive = is_x_pos,
      y_positive = is_y_pos,
      x_range = x_range,
      y_range = y_range,
      x_pad = x_pad,
      y_pad = y_pad,
      proportions = c(
        pp = sum(is_x_pos & is_y_pos) / nrow(dfplot2),
        pm = sum(is_x_pos & !is_y_pos) / nrow(dfplot2),
        mp = sum(!is_x_pos & is_y_pos) / nrow(dfplot2),
        mm = sum(!is_x_pos & !is_y_pos) / nrow(dfplot2)
      )
    )
  }

  bivariate_pdf_theme <- function() {
    theme_minimal(base_family = "sans") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.size = grid::unit(0.9, "cm"),
        legend.key.width = grid::unit(1.1, "cm"),
        legend.key.height = grid::unit(0.75, "cm"),
        legend.spacing.x = grid::unit(0.35, "cm"),
        aspect.ratio = 1
      )
  }

  bivariate_density_pdf_theme <- function() {
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.background = element_rect(fill = "#f8f8f8"),
      plot.background = element_rect(fill = "#f8f8f8"),
      aspect.ratio = 1
    )
  }

  make_bivariate_marker_pdf_plot <- function(data, marker, gate_value, positive_color,
                                             pdf_theme, pdf_point_size = 0.16) {
    if (!all(c(marker, "X_centroid", "Y_centroid") %in% names(data))) {
      return(make_pdf_message_grob(marker, "Marker data are not available for this plot."))
    }

    marker_rows <- complete.cases(data[[marker]], data$X_centroid, data$Y_centroid) &
      is.finite(data[[marker]]) &
      is.finite(data$X_centroid) &
      is.finite(data$Y_centroid)
    marker_data <- data[marker_rows, , drop = FALSE]
    if (nrow(marker_data) == 0) {
      return(make_pdf_message_grob(marker, "No finite marker values available for this plot."))
    }

    is_positive <- if (is.finite(gate_value)) marker_data[[marker]] > gate_value else rep(FALSE, nrow(marker_data))
    is_positive[is.na(is_positive)] <- FALSE
    marker_data$marker_status <- ifelse(is_positive, "Positive", "Negative")
    prop_positive <- sum(is_positive) / nrow(marker_data)

    marker_plot <- ggplot(marker_data, aes(x = X_centroid, y = Y_centroid, color = marker_status)) +
      geom_point(size = pdf_point_size, alpha = 0.35) +
      scale_color_manual(
        guide = guide_legend(title = "", override.aes = list(size = 4.5, alpha = 1)),
        values = c("Positive" = positive_color, "Negative" = "grey70")
      ) +
      labs(title = paste0(marker, "+ cells = ", round(prop_positive, 3)), x = "X Centroid", y = "Y Centroid") +
      pdf_theme

    pos_subset <- marker_data[is_positive, , drop = FALSE]
    if (nrow(pos_subset) >= 2 &&
        length(unique(pos_subset$X_centroid)) >= 2 &&
        length(unique(pos_subset$Y_centroid)) >= 2) {
      marker_plot <- marker_plot +
        geom_density_2d(data = pos_subset, aes(x = X_centroid, y = Y_centroid), inherit.aes = FALSE, color = "black", size = 0.35)
    }
    marker_plot
  }

  precompute_bivariate_marker_pdf_panels <- function(data, markers, gate_values,
                                                     pdf_theme, pdf_point_size = 0.16) {
    as_pdf_grob <- function(plot) {
      if (inherits(plot, "ggplot")) {
        return(ggplot2::ggplotGrob(plot))
      }
      plot
    }

    stats::setNames(
      lapply(markers, function(marker) {
        gate_value <- gate_values[[marker]]
        list(
          x_plot = as_pdf_grob(make_bivariate_marker_pdf_plot(data, marker, gate_value, "green4", pdf_theme, pdf_point_size)),
          y_plot = as_pdf_grob(make_bivariate_marker_pdf_plot(data, marker, gate_value, "blue", pdf_theme, pdf_point_size))
        )
      }),
      markers
    )
  }

  get_bivariate_marker_pdf_panel <- function(marker_panels, marker, axis, pair_data,
                                             pdf_theme, pdf_point_size = 0.16) {
    if (!is.null(marker_panels) && marker %in% names(marker_panels)) {
      return(marker_panels[[marker]][[paste0(axis, "_plot")]])
    }

    gate_value <- if (axis == "x") pair_data$gate_xvar else pair_data$gate_yvar
    positive_color <- if (axis == "x") "green4" else "blue"
    make_bivariate_marker_pdf_plot(pair_data$data, marker, gate_value, positive_color, pdf_theme, pdf_point_size)
  }

  generate_bivariate_plot_grob <- function(pair_data, marker_panels = NULL) {
    if (!is.null(pair_data$error)) {
      return(make_pdf_message_grob(paste0(pair_data$xvar, " vs ", pair_data$yvar), pair_data$error))
    }

    dfplot2 <- pair_data$data
    xvar <- pair_data$xvar
    yvar <- pair_data$yvar
    pdf_point_size <- 0.16
    pdf_theme <- bivariate_pdf_theme()

    density_plot <- ggplot(dfplot2, aes(x = .data[[xvar]], y = .data[[yvar]])) +
      geom_point(aes(color = bivariate_density), size = pdf_point_size, alpha = 0.35) +
      scale_color_viridis(option = "turbo") +
      labs(title = "Bivariate Density", x = xvar, y = yvar) +
      bivariate_density_pdf_theme()
    if (is.finite(pair_data$gate_xvar)) {
      density_plot <- density_plot + geom_vline(xintercept = pair_data$gate_xvar, color = "orange", size = 0.7)
    }
    if (is.finite(pair_data$gate_yvar)) {
      density_plot <- density_plot + geom_hline(yintercept = pair_data$gate_yvar, color = "orange", size = 0.7)
    }
    density_plot <- density_plot +
      geom_text(x = pair_data$x_range[2] - pair_data$x_pad, y = pair_data$y_range[2] - pair_data$y_pad, label = round(pair_data$proportions[["pp"]], 3), color = "red", size = 4) +
      geom_text(x = pair_data$x_range[2] - pair_data$x_pad, y = pair_data$y_range[1] + pair_data$y_pad, label = round(pair_data$proportions[["pm"]], 3), color = "green4", size = 4) +
      geom_text(x = pair_data$x_range[1] + pair_data$x_pad, y = pair_data$y_range[2] - pair_data$y_pad, label = round(pair_data$proportions[["mp"]], 3), color = "blue", size = 4) +
      geom_text(x = pair_data$x_range[1] + pair_data$x_pad, y = pair_data$y_range[1] + pair_data$y_pad, label = round(pair_data$proportions[["mm"]], 3), color = "black", size = 4)

    overlay_plot2 <- ggplot(dfplot2, aes(x = X_centroid, y = Y_centroid, color = bivariate_gate_status)) +
      geom_point(size = pdf_point_size, alpha = 0.55) +
      scale_color_manual(
        guide = guide_legend(title = "", override.aes = list(size = 4.5, alpha = 1)),
        values = c("+/+" = "red", "-/-" = "grey70", "+/-" = "green4", "-/+" = "blue")
      ) +
      labs(title = "Bivariate Gating", x = "X Centroid", y = "Y Centroid") +
      pdf_theme

    gridExtra::arrangeGrob(
      density_plot,
      overlay_plot2,
      get_bivariate_marker_pdf_panel(marker_panels, xvar, "x", pair_data, pdf_theme, pdf_point_size),
      get_bivariate_marker_pdf_panel(marker_panels, yvar, "y", pair_data, pdf_theme, pdf_point_size),
      ncol = 2,
      top = grid::textGrob(
        paste0("Bivariate: ", xvar, " vs ", yvar),
        gp = grid::gpar(fontsize = 17, fontface = "bold")
      )
    )
  }

  generate_bivariate_pdf <- function(data, pdf_file_name, marker_pairs, patient_id, gate_overrides = list()) {
    if (!"imageid" %in% names(data)) {
      data$imageid <- "sample"
    }

    pdf(pdf_file_name, width = 10, height = 10, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)

    if (length(marker_pairs) == 0) {
      grid::grid.draw(make_pdf_message_grob("No bivariate plots available", "Select at least two markers before downloading all bivariate plots."))
      return(invisible(pdf_file_name))
    }

    markers <- unique(unlist(marker_pairs, use.names = FALSE))
    patient_data <- prepare_bivariate_patient_data(data, patient_id, markers)
    gate_values <- precompute_bivariate_gates(patient_data, patient_id, markers, gate_overrides)
    marker_panels <- precompute_bivariate_marker_pdf_panels(
      patient_data,
      markers,
      gate_values,
      bivariate_pdf_theme()
    )

    first_page <- TRUE
    for (pair in marker_pairs) {
      page_grob <- tryCatch(
        generate_bivariate_plot_grob(
          prepare_bivariate_pair_data(
            patient_data = patient_data,
            patient_id = patient_id,
            xvar = pair[[1]],
            yvar = pair[[2]],
            gate_values = gate_values,
            gate_x_override = gate_overrides[[pair[[1]]]],
            gate_y_override = gate_overrides[[pair[[2]]]]
          ),
          marker_panels = marker_panels
        ),
        error = function(e) {
          make_pdf_message_grob(
            paste0(pair[[1]], " vs ", pair[[2]]),
            paste("Could not generate this bivariate plot:", conditionMessage(e))
          )
        }
      )
      if (first_page) {
        first_page <- FALSE
      } else {
        grid::grid.newpage()
      }
      grid::grid.draw(page_grob)
    }

    invisible(pdf_file_name)
  }

  generate_current_bivariate_pdf <- function(pair_data, pdf_file_name) {
    pdf(pdf_file_name, width = 10, height = 10, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)
    grid::grid.draw(generate_bivariate_plot_grob(pair_data))
    invisible(pdf_file_name)
  }

  output$download_bivariate_current <- downloadHandler(
    filename = function() { "current_bivariate_plot.pdf" },
    contentType = "application/pdf",
    content = function(file) {
      tryCatch({
        data <- uploaded_df()
        if (is.null(data) || nrow(data) == 0) {
          return(write_message_pdf(file, "No bivariate plot", "Upload data and run Randkluft before downloading the current bivariate plot."))
        }

        generate_current_bivariate_pdf(current_bivariate_pair_data(), file)
      }, error = function(e) {
        write_message_pdf(file, "Could not generate current bivariate plot", conditionMessage(e))
      })
    }
  )

  output$download_bivariate_all <- downloadHandler(
    filename = function() { "all_bivariate_plots.pdf" },
    contentType = "application/pdf",
    content = function(file) {
      tryCatch({
        data <- uploaded_df()
        if (is.null(data) || nrow(data) == 0) {
          return(write_message_pdf(file, "No bivariate plots", "Upload data and run Randkluft before downloading bivariate plots."))
        }

        markers <- input$selected_columns
        if (is.null(markers) || length(markers) < 2) {
          markers <- unique_markers()
        }
        markers <- markers[markers %in% names(data)]
        marker_pairs <- if (length(markers) >= 2) {
          utils::combn(markers, 2, simplify = FALSE)
        } else {
          list()
        }

        if (!"imageid" %in% names(data)) {
          data$imageid <- "sample"
        }

        patient_ids <- unique_patients_gating()
        patient_ids <- patient_ids[patient_ids %in% unique(data$imageid)]
        if (length(patient_ids) == 0) {
          patient_ids <- unique(data$imageid)
        }

        generate_bivariate_pdf(
          data = data[data$imageid == patient_ids[1], , drop = FALSE],
          pdf_file_name = file,
          marker_pairs = marker_pairs,
          patient_id = patient_ids[1]
        )
      }, error = function(e) {
        write_message_pdf(file, "Could not generate bivariate plots", conditionMessage(e))
      })
    }
  )

  current_bivariate_pair_data <- reactive({
    data <- subsetted()
    patient_id <- unique_patients_gating()[1]
    xvar <- selected_var_name(input$xvar)
    yvar <- selected_var_name(input$yvar)

    req(patient_id, !is.null(xvar), !is.null(yvar))
    markers <- unique(c(xvar, yvar))
    patient_data <- prepare_bivariate_patient_data(data, patient_id, markers)

    gate_overrides <- list()
    if (!is.null(react_xgate())) gate_overrides[[xvar]] <- react_xgate()
    if (!is.null(react_ygate())) gate_overrides[[yvar]] <- react_ygate()
    gate_values <- precompute_bivariate_gates(patient_data, patient_id, markers, gate_overrides)

    prepare_bivariate_pair_data(
      patient_data = patient_data,
      patient_id = patient_id,
      xvar = xvar,
      yvar = yvar,
      gate_values = gate_values,
      gate_x_override = react_xgate(),
      gate_y_override = react_ygate()
    )
  })

  bivariate_app_spatial_theme <- function() {
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
      axis.title = element_text(size = 16, family = "Arial"),
      axis.text = element_text(size = 14, family = "Arial"),
      legend.position = "bottom",
      legend.text = element_text(size = 16, family = "Arial"),
      legend.title = element_text(size = 16, face = "bold", family = "Arial"),
      legend.key.size = grid::unit(1.2, "cm"),
      legend.key.width = grid::unit(1.4, "cm"),
      legend.key.height = grid::unit(0.9, "cm"),
      legend.spacing.x = grid::unit(0.35, "cm"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  }

  make_bivariate_app_marker_plot <- function(pair_data, axis, positive_color) {
    dfplot2 <- pair_data$data
    marker <- if (axis == "x") pair_data$xvar else pair_data$yvar
    status_column <- paste0(axis, "_status")
    positive_rows <- pair_data[[paste0(axis, "_positive")]]
    prop_positive <- sum(positive_rows) / nrow(dfplot2)

    marker_plot <- ggplot(dfplot2, aes(x = X_centroid, y = Y_centroid, color = .data[[status_column]])) +
      geom_point(size = 0.3, alpha = 0.35) +
      scale_color_manual(
        guide = guide_legend(title = "", override.aes = list(size = 7, alpha = 1)),
        values = c("Positive" = positive_color, "Negative" = "grey")
      ) +
      labs(title = paste(marker, "+ cell=", round(prop_positive, 3), sep = ""), x = "X Centroid", y = "Y Centroid") +
      bivariate_app_spatial_theme()

    pos_subset <- dfplot2[positive_rows, , drop = FALSE]
    if (nrow(pos_subset) >= 2 &&
        length(unique(pos_subset$X_centroid)) >= 2 &&
        length(unique(pos_subset$Y_centroid)) >= 2) {
      marker_plot <- marker_plot +
        geom_density_2d(data = pos_subset, aes(x = X_centroid, y = Y_centroid), inherit.aes = FALSE, color = "black")
    }
    marker_plot
  }

  make_bivariate_app_grob <- function(pair_data) {
    dfplot2 <- pair_data$data
    xvar <- pair_data$xvar
    yvar <- pair_data$yvar

    density_plot <- ggplot(dfplot2, aes(x = .data[[xvar]], y = .data[[yvar]])) +
      geom_point(aes(color = bivariate_density), size = 0.3, alpha = 0.35) +
      scale_color_viridis(option = "turbo") +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial"),
        axis.title = element_text(size = 16, family = "Arial"),
        axis.text = element_text(size = 14, family = "Arial"),
        panel.background = element_rect(fill = "#f8f8f8"),
        plot.background = element_rect(fill = "#f8f8f8")
      ) +
      labs(title = "Bivariate Density")
    if (is.finite(pair_data$gate_xvar)) {
      density_plot <- density_plot + geom_vline(xintercept = pair_data$gate_xvar, color = "orange")
    }
    if (is.finite(pair_data$gate_yvar)) {
      density_plot <- density_plot + geom_hline(yintercept = pair_data$gate_yvar, color = "orange")
    }
    density_plot <- density_plot +
      geom_text(x = pair_data$x_range[2] - pair_data$x_pad, y = pair_data$y_range[2] - pair_data$y_pad, label = round(pair_data$proportions[["pp"]], 3), color = "red", size = 6) +
      geom_text(x = pair_data$x_range[2] - pair_data$x_pad, y = pair_data$y_range[1] + pair_data$y_pad, label = round(pair_data$proportions[["pm"]], 3), color = "green", size = 6) +
      geom_text(x = pair_data$x_range[1] + pair_data$x_pad, y = pair_data$y_range[2] - pair_data$y_pad, label = round(pair_data$proportions[["mp"]], 3), color = "blue", size = 6) +
      geom_text(x = pair_data$x_range[1] + pair_data$x_pad, y = pair_data$y_range[1] + pair_data$y_pad, label = round(pair_data$proportions[["mm"]], 3), color = "black", size = 6)

    overlay_plot2 <- ggplot(dfplot2, aes(x = X_centroid, y = Y_centroid, color = bivariate_gate_status)) +
      geom_point(size = 0.3, alpha = 0.7) +
      scale_color_manual(
        guide = guide_legend(title = "", override.aes = list(size = 7, alpha = 1)),
        values = c("+/+" = "red", "-/-" = "grey", "+/-" = "green", "-/+" = "blue", "Other" = "black")
      ) +
      labs(title = "Bivariate Gating", x = "X Centroid", y = "Y Centroid") +
      bivariate_app_spatial_theme()

    gridExtra::arrangeGrob(
      density_plot,
      overlay_plot2,
      make_bivariate_app_marker_plot(pair_data, "x", "green"),
      make_bivariate_app_marker_plot(pair_data, "y", "blue"),
      ncol = 2
    )
  }

  output$plot2 <- renderPlot({
    pair_data <- current_bivariate_pair_data()
    validate(need(is.null(pair_data$error), pair_data$error))
    grid::grid.draw(make_bivariate_app_grob(pair_data))
  })

  output$plot_trivariate <- renderPlotly({
    generatePlotTrivar()
  })


add_trivariate_gate_plane <- function(plot, axis, gate_value, plot_ranges, marker_label) {
  if (!is.finite(gate_value)) {
    return(plot)
  }

  plane <- switch(
    axis,
    x = list(
      x = rep(gate_value, 4),
      y = c(plot_ranges$y[1], plot_ranges$y[2], plot_ranges$y[2], plot_ranges$y[1]),
      z = c(plot_ranges$z[1], plot_ranges$z[1], plot_ranges$z[2], plot_ranges$z[2])
    ),
    y = list(
      x = c(plot_ranges$x[1], plot_ranges$x[2], plot_ranges$x[2], plot_ranges$x[1]),
      y = rep(gate_value, 4),
      z = c(plot_ranges$z[1], plot_ranges$z[1], plot_ranges$z[2], plot_ranges$z[2])
    ),
    z = list(
      x = c(plot_ranges$x[1], plot_ranges$x[2], plot_ranges$x[2], plot_ranges$x[1]),
      y = c(plot_ranges$y[1], plot_ranges$y[1], plot_ranges$y[2], plot_ranges$y[2]),
      z = rep(gate_value, 4)
    )
  )

  plotly::add_trace(
    plot,
    x = plane$x,
    y = plane$y,
    z = plane$z,
    i = c(0, 0),
    j = c(1, 2),
    k = c(2, 3),
    type = "mesh3d",
    # facecolor = rep("#FFA500", 2),
    color = "orange",
    opacity = 0.06,
    inherit = FALSE,
    showlegend = TRUE,
    showscale = FALSE,
    name = paste(marker_label, "gate"),
    hoverinfo = "text",
    text = rep(paste(marker_label, "gate =", round(gate_value, 3)), 4)
  )
}

generatePlotTrivar <- function() {
  dfplot2 <- subsetted_tri()
  patient_selected <- unique_patients_gating()[1]
  xvar <- selected_var_name(input$xvarTri)
  yvar <- selected_var_name(input$yvarTri)
  zvar <- selected_var_name(input$zvar)

  req(!is.null(xvar), !is.null(yvar), !is.null(zvar))
  validate(need(
    all(c(xvar, yvar, zvar) %in% names(dfplot2)),
    "Select three available markers for the trivariate plot."
  ))

  dfplot2[[xvar]] <- suppressWarnings(as.numeric(dfplot2[[xvar]]))
  dfplot2[[yvar]] <- suppressWarnings(as.numeric(dfplot2[[yvar]]))
  dfplot2[[zvar]] <- suppressWarnings(as.numeric(dfplot2[[zvar]]))
  finite_rows <- complete.cases(dfplot2[[xvar]], dfplot2[[yvar]], dfplot2[[zvar]]) &
    is.finite(dfplot2[[xvar]]) &
    is.finite(dfplot2[[yvar]]) &
    is.finite(dfplot2[[zvar]])
  dfplot2 <- dfplot2[finite_rows, , drop = FALSE]
  validate(need(nrow(dfplot2) > 0, "No finite marker values are available for the trivariate plot."))

  dfplot2$trivariate_density <- get_density_3d(
    dfplot2[[xvar]],
    dfplot2[[yvar]],
    dfplot2[[zvar]],
    grid_size = 50
  )

  gate_values <- c(
    x = resolve_gate_value(patient_selected, xvar, dfplot2),
    y = resolve_gate_value(patient_selected, yvar, dfplot2),
    z = resolve_gate_value(patient_selected, zvar, dfplot2)
  )
  plot_ranges <- list(
    x = range(dfplot2[[xvar]], na.rm = TRUE),
    y = range(dfplot2[[yvar]], na.rm = TRUE),
    z = range(dfplot2[[zvar]], na.rm = TRUE)
  )
  trivariate_density_range <- range(dfplot2$trivariate_density, finite = TRUE)
  if (!all(is.finite(trivariate_density_range)) || diff(trivariate_density_range) == 0) {
    trivariate_density_range <- c(0, 1)
  }
  trivariate_colors <- plotly::toRGB(viridisLite::viridis(256, option = "turbo"))
  trivariate_colorscale <- cbind(
    seq(0, 1, length.out = length(trivariate_colors)),
    trivariate_colors
  )

  p <- plot_ly(
    dfplot2,
    x = ~dfplot2[[xvar]],
    y = ~dfplot2[[yvar]],
    z = ~dfplot2[[zvar]],
    type = "scatter3d",
    mode = "markers",
    showlegend = FALSE,
    marker = list(
      size = 2,
      # alpha=0.9,
      color = dfplot2$trivariate_density,
      colorscale = trivariate_colorscale,
      cmin = trivariate_density_range[1],
      cmax = trivariate_density_range[2],
      showscale = FALSE
    )
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = paste(xvar)),
        yaxis = list(title = paste(yvar)),
        zaxis = list(title = paste(zvar))
      )
    )

  p <- add_trivariate_gate_plane(p, "x", gate_values[["x"]], plot_ranges, xvar)
  p <- add_trivariate_gate_plane(p, "y", gate_values[["y"]], plot_ranges, yvar)
  add_trivariate_gate_plane(p, "z", gate_values[["z"]], plot_ranges, zvar)
}


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
})


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

      output$pheno_bar <- renderPlot({
              barchart
            })


  })

  output$phenotypeTable <- renderTable({
    phenotype_df()
  })


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

      })

shinyApp(ui = ui, server = server)
