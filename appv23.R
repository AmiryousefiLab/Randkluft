# ==============================================================================
# Randkluft – appv23.R
# Improved version: reactivity, stability, performance, UX
# Changes vs appv22:
#   - Removed shiny.trace (was flooding logs)
#   - Font family switched from "Arial" to "sans" (Arial crashes on Linux servers)
#   - All print() calls removed / replaced with a DEBUG flag
#   - columns_to_exclude and columns_to_exclude_dna consolidated to one helper
#   - rbind-in-loop replaced with lapply + bind_rows (O(n) vs O(n²))
#   - determine_positivity vectorised with ifelse
#   - resultdf_reactive updated in both single and multi-column paths
#   - observe() for UI updates collapsed where possible / guarded with req()
#   - downloadHandler for downloadcurrent defined once, not inside generatePlot()
#   - ggsave() removed from inside renderPlot
#   - UMAP cbind column-name collision fixed
#   - plotly_selected event registered to suppress warning
#   - Disconnection guards: tryCatch around skew_gate calls
#   - Removed ~500 lines of commented-out dead code
# ==============================================================================

# ---- Libraries ---------------------------------------------------------------
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
library(fastICA)
library(tiff)
library(datasets)
data(iris)
library(stringr)
library(tidyr)
library(PEkit)
library(sqldf)
library(ggplot2)
library(gridExtra)

source("utils.R")

# ---- Static data -------------------------------------------------------------
df_example               <- read.csv("exemplar-001--unmicst_cell.csv")
ground_truth_gates_loaded <- read.csv("tuulia_data_GT.csv")
df_example_PHENOTYPE     <- read.csv("phenotype_table_help.csv")

# ---- Global options ----------------------------------------------------------
# shiny.trace removed – it was writing every WebSocket frame to the log and
# contributing to the 200 KB+ log bloat observed in Randkluft-logs.txt
options(shiny.maxRequestSize = 30 * 1024^3)

# Set to TRUE locally for verbose output; FALSE for production
DEBUG <- FALSE
dlog <- function(...) { if (DEBUG) message(...) }

# ---- Shared helpers ----------------------------------------------------------

# Columns that are never marker targets (exact names)
NON_MARKER_COLS <- c(
  "imageid", "phenotype", "ROI_major_category", "CellID", "X", "Y",
  "ROI_minor_category", "phenotype_v2", "X_centroid", "Y_centroid",
  "Eccentricity", "Area", "MajorAxisLength", "MinorAxisLength",
  "Extent", "Solidity", "Orientation", "",
  # Explicit DNA/Hoechst/DAPI numbered variants
  "DNA1","DNA2","DNA3","DNA4","DNA5","DNA6","DNA7","DNA8","DNA9","DNA10",
  "DNA11","DNA12","DNA13",
  "DNA_1","DNA_2","DNA_3","DNA_4","DNA_5","DNA_6","DNA_7","DNA_8","DNA_9",
  "DNA_10","DNA_11","DNA_12","DNA_13","DNA6a",
  "DAPI1","DAPI2","DAPI3","DAPI4","DAPI5","DAPI6","DAPI7","DAPI8","DAPI9","DAPI10",
  "DAPI_1","DAPI_2","DAPI_3","DAPI_4","DAPI_5","DAPI_6","DAPI_7","DAPI_8","DAPI_9","DAPI_10",
  "Hoechst1","Hoechst2","Hoechst3","Hoechst4","Hoechst5","Hoechst6",
  "Hoechst7","Hoechst8","Hoechst9","Hoechst10",
  "Hoechst_1","Hoechst_2","Hoechst_3","Hoechst_4","Hoechst_5",
  "Hoechst_6","Hoechst_7","Hoechst_8","Hoechst_9","Hoechst_10"
)
# Pattern catches any numbered variant not covered by the exact list above
NON_MARKER_PATTERN <- "DNA|DAPI|Hoechst"
POSITIVITY_PATTERN <- "_positivity"

get_marker_cols <- function(df, include_dna = FALSE) {
  excl <- NON_MARKER_COLS
  if (!include_dna) {
    dna_cols <- grep(NON_MARKER_PATTERN, colnames(df), value = TRUE, ignore.case = TRUE)
    excl     <- c(excl, dna_cols)
  }
  # Exclude positivity flag columns added after gating
  pos_cols <- grep(POSITIVITY_PATTERN, colnames(df), value = TRUE, ignore.case = TRUE)
  excl     <- c(excl, pos_cols)
  # Keep only numeric columns (true marker intensity values)
  candidate    <- setdiff(names(df), excl)
  numeric_cols <- candidate[sapply(candidate, function(cn) is.numeric(df[[cn]]))]
  numeric_cols
}

get_dna_cols <- function(df) {
  grep(NON_MARKER_PATTERN, colnames(df), value = TRUE, ignore.case = TRUE)
}

# ---- UI ----------------------------------------------------------------------
ui <- shinyUI(fluidPage(
  shinyjs::useShinyjs(),
  tags$header(HTML(html_code)),

  tags$style(HTML("
    .pagination-button      { float: right; margin-right: 10px; }
    .pagination-button-back { float: left;  margin-left:  10px; }
    .current-label          { font-size: 14px; color: #555; padding: 6px 0; }
  ")),

  tags$head(tags$link(
    rel  = "stylesheet",
    href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"
  )),

  tags$div(
    id = "title",
    tags$h1(id = "tool_name", "Randk\uft", style = "color:navy;"),
    tags$h4(
      id  = "tool_exp",
      tags$em("Unitary gating of the CyCIF markers"),
      style = "color:black;"
    )
  ),

  tags$div(
    id = "tool",
    tabsetPanel(

      # ---- Home --------------------------------------------------------------
      tabPanel(
        strong("Home"),
        br(),
        tags$div(includeMarkdown("./documents/home.md"), style = "max-width:800px;")
      ),

      # ---- Randkluft main tab ------------------------------------------------
      tabPanel(
        strong("Randkluft"),
        br(),
        sidebarLayout(
          sidebarPanel(
            width = 4,
            br(),
            tabsetPanel(
              id = "sidebartab",

              # ---- Upload panel --------------------------------------------
              tabPanel(
                strong("Upload file"),
                value = 1,
                br(),
                p("Please upload a .csv file with cell-marker data following the format described in",
                  strong("Help,"), " or download the example dataset found in",
                  downloadLink("download_example_data", strong(" here."))),

                fluidRow(
                  column(8,
                    div(id = "file_input_div",
                        fileInput("cell_file", label = "Upload a CSV file",
                                  accept = ".csv", buttonLabel = "Browse"))
                  ),
                  tags$div(id = "message_loading_data",
                           style = "font-size:14px; position:fixed; bottom:0; right:calc(5% + 10px);"),
                  column(4,
                    actionButton("remove_file", "Remove", icon = icon("trash"))
                  )
                )
              ),

              # ---- Randkluft analysis panel --------------------------------
              tabPanel(
                strong("Randkluft"),
                value = 3,
                tabsetPanel(
                  id = "sidebar2",

                  # Essential sub-panel
                  tabPanel(
                    strong("Essential"),
                    value = 1,
                    br(),
                    actionButton(inputId = "run_gate", "Randkluft"),
                    br(), br(),
                    checkboxGroupInput("selected_columns",       "Select Markers",  choices = NULL, selected = NULL),
                    checkboxGroupInput("unique_patients_gating", "Select Patients", choices = NULL, selected = NULL),
                    br(),
                    materialSwitch(inputId = "gen_hist_plots_on_off",
                                   label   = "Show gates",
                                   status  = "danger",
                                   right   = TRUE,
                                   value   = FALSE),
                    tags$div(id = "message",
                             style = "font-size:14px; position:fixed; bottom:0; right:calc(5% + 10px);"),
                    tags$div(id = "message_gen_hist",
                             style = "font-size:14px; position:fixed; bottom:0; right:calc(5% + 10px);"),
                    br(),
                    p("Download gate estimations as .csv:"),
                    downloadButton("downloadEstimations", "Download Gate Estimates"),
                    br(), br(),
                    p("Download the current 4-panel plot:"),
                    downloadButton("downloadcurrent",     "Download Current Plot"),
                    br(), br(),
                    p("Download all plots as PDF:"),
                    downloadButton("downloadall",         "Download All Plots")
                  ),

                  # Bivariate sub-panel
                  tabPanel(
                    strong("Bivariate"),
                    value = 2,
                    br(),
                    varSelectInput("xvar", "X variable", NULL, selected = NULL),
                    varSelectInput("yvar", "Y variable", NULL, selected = NULL),
                    radioButtons("patient_number_plot2", "Select Patients",
                                 choices = character(0), selected = character(0)),
                    hr()
                  ),

                  # Trivariate sub-panel
                  tabPanel(
                    strong("Trivariate"),
                    value = 3,
                    br(),
                    varSelectInput("xvarTri", "X variable", NULL, selected = NULL),
                    varSelectInput("yvarTri", "Y variable", NULL, selected = NULL),
                    varSelectInput("zvar",    "Z variable", NULL, selected = NULL),
                    radioButtons("trivar_patient", "Select Patients",
                                 choices = character(0), selected = character(0)),
                    hr()
                  )
                )
              ),

              # ---- Extra / phenotyping panel --------------------------------
              tabPanel(
                strong("Extra"),
                value = 4,
                tabsetPanel(
                  id = "sidebar_post",
                  tabPanel(
                    strong("Phenotyping"),
                    value = 1,
                    br(),
                    p("Please upload your phenotyping workflow following the format described in",
                      strong("Help,"), " or download the example found in",
                      downloadLink("download_example_data2", strong(" here."))),
                    fluidRow(
                      column(8,
                        div(id = "file_input_div2",
                            fileInput("phen_wfl", label = "Upload a CSV file",
                                      accept = ".csv", buttonLabel = "Browse"))
                      ),
                      tags$div(id = "message_loading_data2",
                               style = "font-size:14px; position:fixed; bottom:0; right:calc(5% + 10px);"),
                      br(),
                      column(4,
                        actionButton("remove_file2", "Remove", icon = icon("trash"))
                      )
                    ),
                    actionButton("define_phenotype_AUTO", "Phenotype my data"),
                    br(), br(),
                    uiOutput("marker_checkboxes"),
                    materialSwitch(inputId = "any_indicator",
                                   label   = "Any positive",
                                   status  = "info",
                                   right   = TRUE,
                                   value   = FALSE),
                    textInput("phenotype_name", "Phenotype Name"),
                    actionButton("define_phenotype", "Add phenotype definition"),
                    br(), br(),
                    p("Download your CSV with phenotypes, or the workflow you defined."),
                    downloadButton("downloadPhenotypes",    "Download Phenotyped Data"),
                    downloadButton("downloadPhenotypeTable","Download Workflow"),
                    br()
                  )
                )
              )

            ) # end sidebartab tabsetPanel
          ), # end sidebarPanel

          # ---- Main panel --------------------------------------------------
          mainPanel(
            width = 8,
            tabsetPanel(

              # Essential main panel
              conditionalPanel(
                condition = "input.sidebar2 == 1 && input.sidebartab == 3",
                br(),
                numericInput("intercept", "Type Gate Value", value = NA),
                actionButton("updateGates", "Update Gate"),
                br(),
                # Current marker / patient label
                uiOutput("current_position_label"),
                br(),
                plotOutput("gated_histogram_on_page", height = "1000px", width = "1000px"),
                br(),
                div(class = "pagination-button-back",
                    actionButton("prevMarker",  "Previous Marker ",  icon("arrow-left"))),
                div(class = "pagination-button",
                    actionButton("nextMarker",  "Next Marker ",      icon("arrow-right"))),
                div(class = "pagination-button-back",
                    actionButton("prevPatient", "Previous Patient ", icon("arrow-left"))),
                div(class = "pagination-button",
                    actionButton("nextPatient", "Next Patient ",     icon("arrow-right")))
              ),

              # Bivariate main panel
              conditionalPanel(
                condition = "input.sidebartab == 3 && input.sidebar2 == 2",
                numericInput("gate_xvar_update", "Enter Gate for X Variable:", value = NA),
                numericInput("gate_yvar_update", "Enter Gate for Y Variable:", value = NA),
                actionButton("update_gates_bivariate", "Update Gates"),
                plotOutput("plot2", height = "1000px", width = "1000px"),
                verbatimTextOutput("prop_summary")
              ),

              # Trivariate main panel
              conditionalPanel(
                condition = "input.sidebartab == 3 && input.sidebar2 == 3",
                numericInput("gate_xvar_updateTri", "Enter Gate for X Variable:", value = NA),
                numericInput("gate_yvar_updateTri", "Enter Gate for Y Variable:", value = NA),
                numericInput("gate_zvar_update",    "Enter Gate for Z Variable:", value = NA),
                actionButton("update_gates_trivariate", "Update Gates"),
                plotlyOutput("plot_trivariate", height = "1000px", width = "1000px")
              ),

              # Phenotyping main panel
              conditionalPanel(
                condition = "input.sidebartab == 4 && input.sidebar_post == 1",
                column(6, tableOutput("phenotypeTable")),
                column(6, plotOutput("pheno_bar", height = "500px", width = "500px")),
                verbatimTextOutput("post_statistics")
              ),

              # UMAP main panel
              conditionalPanel(
                condition = "input.sidebartab == 4 && input.sidebar_post == 2",
                plotOutput("icaAnalysis2")
              )

            ) # end inner tabsetPanel
          ) # end mainPanel
        ) # end sidebarLayout
      ),   # end Randkluft tabPanel

      # ---- Help, FAQ, Contact -----------------------------------------------
      tabPanel(
        strong("Help"),
        br(),
        tags$div(includeMarkdown("./documents/help.md"), style = "max-width:800px;")
      ),
      tabPanel(
        strong("FAQ"),
        br(),
        column(1, ""),
        br(),
        column(6,
          h4(strong("Q:"), tags$em(strong("Why are some proportions or total sample numbers zero?"))),
          p(strong("A:"),
            "Randkluft searches for a positively skewed signal emerging from background noise.
  In cases where the marker distribution is already negatively skewed or lacks a discernible
  positive tail, the algorithm terminates early without estimating a gate. In these situations,
  we recommend visual inspection of the distribution and manual gating."),
          br(),
          h4(strong("Q:"), tags$em(strong("Why do I get 'Disconnected from the server' after uploading my data?"))),
          p(strong("A:"),
            "This typically indicates that the uploaded file does not conform to the expected input
  format. Please consult the Help section and ensure that column names, data types, and required
  fields strictly follow the documented input structure before re-uploading."),
          br(),
          h4(strong("Q:"), tags$em(strong("Why are upload and analysis slow?"))),
          p(strong("A:"),
            "Randkluft treats all numeric columns as potential gating targets. If your input file
  contains many columns not required for analysis, removing them before upload can significantly
  improve performance."),
          br(),
          h4(strong("Q:"), tags$em(strong("Which data are used to detect the gates?"))),
          p(strong("A:"),
            "At each step of the workflow, Randkluft internally stores the active dataset associated
  with the selected panel. All subsequent analyses use this updated data. Both the modified datasets
  and the resulting gate estimates can be downloaded at each stage of the analysis.")
        ),
        br()
      ),
      tabPanel(
        strong("Contact"),
        br(),
        tags$div(includeMarkdown("./documents/contact.md"), style = "max-width:800px;")
      )

    ) # end top-level tabsetPanel
  )   # end #tool div
))    # end fluidPage / shinyUI


# ==============================================================================
# SERVER
# ==============================================================================
options(shiny.maxRequestSize = 100 * 1024^3)

server <- shinyServer(function(input, output, session) {

  # ---- Utility functions ---------------------------------------------------

  remove_outliers2 <- function(target, low_percentile = 1, high_percentile = 99) {
    lo <- quantile(target, low_percentile  / 100, na.rm = TRUE)
    hi <- quantile(target, high_percentile / 100, na.rm = TRUE)
    target[target >= lo & target <= hi]
  }

  find_mode2 <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix   <- findInterval(x, dens$x)
    iy   <- findInterval(y, dens$y)
    ii   <- cbind(ix, iy)
    dens$z[ii]
  }

  determine_positivity <- function(marker_intensity, gate_value) {
    # Vectorised – no sapply loop needed
    ifelse(marker_intensity > gate_value, "+", "-")
  }

  # ---- Core gating algorithm (UNCHANGED) -----------------------------------

  skew_gate <- function(x, alpha = 0.01) {
    sk <- moments::skewness(x)
    n  <- length(x)
    a  <- locmodes(x)$locations[1]
    b  <- max(x)

    if (sk < 0) {
      message("The skewness is negative – falling back to GMM gate.")
      outputGMM  <- Mclust(x, G = 2)
      gmm_gate   <- mean(outputGMM$parameters$mean)
      n_removed  <- sum(x > gmm_gate)
      perc_rem   <- round(n_removed / n, 3)
      return(list(
        skewness           = sk,
        cutoff             = gmm_gate,
        N                  = n,
        N_removed          = n_removed,
        percentage_removed = perc_rem,
        returnvalplot      = x[x > a + (b - a) / 2]
      ))
    }

    iteration <- 0
    while (abs(sk) > alpha && iteration <= 100) {
      if (sk >= 0) b <- a + (b - a) / 2
      else         a <- a + (b - a) / 2
      a  <- min(a, b)
      b  <- max(a, b)
      sk <- skewness(x[x < b])
      if (is.nan(sk)) { message("Warning: skewness is NaN"); break }
      iteration <- iteration + 1
    }

    n_removed <- sum(x > b)
    perc_rem  <- round(n_removed / n, 3)
    list(
      skewness           = sk,
      cutoff             = b,
      N                  = n,
      N_removed          = n_removed,
      percentage_removed = perc_rem,
      returnvalplot      = x[x > a + (b - a) / 2]
    )
  }

  # ---- Gate computation (performance-improved) -----------------------------

  get_gates_csv <- function(dataframe) {
    data          <- dataframe
    unique_ids    <- unique(data$imageid)
    alpha         <- 0.01
    total_iters   <- length(unique_ids) * (ncol(data) - 1)

    progress <- Progress$new(session, min = 1, max = total_iters)
    on.exit(progress$close(), add = TRUE)
    progress$set(message = "Randkluft in action…", value = 0)

    iter <- 0
    # Use lapply instead of rbind-in-loop (O(n) vs O(n²))
    rows <- lapply(unique_ids, function(imageid) {
      sub <- data[data$imageid == imageid, -1, drop = FALSE]
      lapply(seq_len(ncol(sub)), function(j) {
        iter <<- iter + 1
        progress$set(value = iter)
        col_name <- colnames(sub)[j]
        values   <- remove_outliers2(sub[, j])
        values   <- values[is.finite(values)]
        gate     <- tryCatch(skew_gate(values, alpha)$cutoff,
                             error = function(e) NA_real_)
        dlog("imageid:", imageid, "| marker:", col_name, "| gate:", gate)
        data.frame(Patient = imageid, Marker = col_name, Gate = gate,
                   stringsAsFactors = FALSE)
      })
    })
    bind_rows(unlist(rows, recursive = FALSE))
  }

  # ---- Plot: single 4-panel ------------------------------------------------

  outputinterceptreactive <- reactiveVal(NULL)

  plot_gating_individual <- function(target, gate_result, marker_label, gt_gate = NULL) {
    cutoff <- gate_result$cutoff

    # Compute a sensible y position for the annotation (density peak)
    dens_y <- tryCatch({
      d <- density(target, na.rm = TRUE)
      d$y[which.max(d$y)]
    }, error = function(e) 0.5)

    p <- ggplot(data.frame(x = target), aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 100,
                     fill = "lightblue", color = "black") +
      geom_density() +
      geom_vline(xintercept = cutoff, color = "red", linewidth = 0.8)

    if (!is.null(gt_gate) && length(gt_gate) == 1 && is.finite(gt_gate))
      p <- p + geom_vline(xintercept = gt_gate, color = "blue", linewidth = 0.8)

    p + geom_text(
      x     = cutoff,
      y     = dens_y,
      label = paste0("Gate: ", round(cutoff, 2)),
      hjust = -0.05,
      vjust = 1.2,
      size  = 4,
      color = "darkred"
    ) +
    labs(
      x     = marker_label,
      y     = "Density",
      title = paste0(
        "N+ = ", gate_result$N_removed,
        "   +R = ", round(gate_result$percentage_removed, 3),
        "   Gate = ", round(cutoff, 2)
      )
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", family = "sans"),
      axis.title = element_text(size = 13, family = "sans"),
      axis.text  = element_text(size = 11, family = "sans")
    )
  }

  plot_gating_grid <- function(df, unique_id, col_idx) {
    target <- df[df$imageid == unique_id, col_idx]
    target <- remove_outliers2(target)
    target <- target[is.finite(target)]
    marker <- paste0(c(unique_id, names(df)[col_idx]), collapse = "; ")

    gate_result <- tryCatch(skew_gate(target, 0.01),
                            error = function(e) list(cutoff = NA, N_removed = 0,
                                                     percentage_removed = 0))
    if (!is.null(outputinterceptreactive()))
      gate_result$cutoff <- outputinterceptreactive()

    gt_gate <- tryCatch(
      ground_truth_gates_loaded$Gate[
        ground_truth_gates_loaded$Patient == unique_id &
        ground_truth_gates_loaded$Marker  == names(df)[col_idx]
      ],
      error = function(e) NULL
    )
    plot_gating_individual(target, gate_result, marker, gt_gate)
  }

  # ---- Reactive state ------------------------------------------------------

  uploaded_df             <- reactiveVal(NULL)
  reactive_markers_sel    <- reactiveVal(NULL)
  pdf_file_path           <- reactiveVal(NULL)
  csv_save_file_path      <- reactiveVal(NULL)
  resultdf_reactive       <- reactiveVal(NULL)
  current_marker          <- reactiveVal(1)
  current_patient         <- reactiveVal(1)
  current_4panel          <- reactiveValues(p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL)
  histplot_gate_switch    <- reactiveVal(FALSE)
  gmm_gate_switch         <- reactiveVal(TRUE)
  brushed_data            <- reactiveVal(NULL)
  filtered_data_original  <- reactiveVal(NULL)
  filtered_data_reactive  <- reactiveVal(NULL)
  react_xgate             <- reactiveVal(NULL)
  react_ygate             <- reactiveVal(NULL)
  regression_mode_react   <- reactiveVal(NULL)
  klPlotReact             <- reactiveVal(NULL)
  filtered_data_updated   <- reactiveVal(NULL)
  phenotype_wfl_reactive  <- reactiveVal(NULL)
  phenotype_df            <- reactiveVal(data.frame(
                               phenotype = character(),
                               markers   = character()
                             ))

  # ---- Derived reactives ---------------------------------------------------

  unique_patients <- reactive({
    req(uploaded_df())
    unique(uploaded_df()$imageid)
  })

  unique_markers <- reactive({
    req(uploaded_df())
    get_marker_cols(uploaded_df(), include_dna = FALSE)
  })

  unique_markers_w_dna <- reactive({
    req(uploaded_df())
    get_marker_cols(uploaded_df(), include_dna = TRUE)
  })

  cycle_detected <- reactive({
    req(uploaded_df())
    get_dna_cols(uploaded_df())
  })

  selected_columns <- reactive({ input$selected_columns })
  unique_patients_gating <- reactive({ input$unique_patients_gating })

  # ---- Hide elements on startup -------------------------------------------

  shinyjs::hide("patient_number")

  # ---- File upload ---------------------------------------------------------

  observeEvent(input$cell_file, {
    req(input$cell_file)

    runjs("
      (function() {
        var el = document.getElementById('message_loading_data');
        el.innerText = 'Loading…';
        el.style.display = 'block';
      })();
    ")

    n_rows <- length(count.fields(input$cell_file$datapath))

    # Sample large files to 70k rows via SQL for speed
    if (n_rows > 70000) {
      samp  <- sample(seq_len(n_rows), 70000, replace = FALSE)
      query <- sprintf("SELECT * FROM file WHERE CellID IN (%s)",
                       paste(samp, collapse = ", "))
      df <- tryCatch(
        read.csv.sql(input$cell_file$datapath, sql = query, dbname = tempfile()),
        error = function(e) {
          showNotification(paste("SQL sampling failed:", e$message,
                                 "– loading full file."), type = "warning")
          read.csv(input$cell_file$datapath, header = TRUE, check.names = FALSE)
        }
      )
    } else {
      df <- read.csv(input$cell_file$datapath, header = TRUE, check.names = FALSE)
    }

    # Normalise column names
    if ("imageID" %in% colnames(df))
      colnames(df)[colnames(df) == "imageID"] <- "imageid"
    if ("X" %in% colnames(df))
      colnames(df)[colnames(df) == "X"] <- "X_centroid"
    if ("Y" %in% colnames(df))
      colnames(df)[colnames(df) == "Y"] <- "Y_centroid"

    # Drop unnamed columns
    df <- df[, colnames(df) != "", drop = FALSE]

    # Identify marker columns for log-transform check
    marker_cols <- get_marker_cols(df, include_dna = FALSE)
    num_vals    <- df[, marker_cols, drop = FALSE]
    num_vals    <- data.frame(lapply(num_vals, as.numeric))

    # Floor very small positives (< 10) to 10 before log
    num_vals <- data.frame(lapply(num_vals, function(x) ifelse(x >= 0 & x < 10, 10, x)))

    if (any(sapply(num_vals, function(x) max(x, na.rm = TRUE)) > 20)) {
      cat("Working with raw data – log-transforming marker columns.\n")
      num_vals <- data.frame(lapply(num_vals, log))
    } else {
      cat("Data appears already log-transformed.\n")
    }
    df[, marker_cols] <- num_vals

    # Ensure imageid column exists
    if (!"imageid" %in% colnames(df))
      df$imageid <- "image"

    uploaded_df(df)

    dlog("Loaded columns: ", paste(colnames(df), collapse = ", "))

    # Update UI
    col_choices <- get_marker_cols(df, include_dna = FALSE)
    reactive_markers_sel(col_choices)

    updateCheckboxGroupInput(session, "selected_columns",
                             choices  = col_choices,
                             selected = col_choices)
    updateCheckboxGroupInput(session, "unique_patients_gating",
                             choices  = unique(df$imageid),
                             selected = unique(df$imageid))

    shinyjs::show("selected_columns")
    shinyjs::show("unique_patients_gating")

    runjs("document.getElementById('message_loading_data').style.display = 'none';")
    showNotification("File ready for use.", duration = 8, id = "file_ready")
  })

  # ---- Remove file ---------------------------------------------------------

  observeEvent(input$remove_file, {
    shinyjs::reset("cell_file")
    shinyjs::hide("selected_columns")
    shinyjs::hide("unique_patients_gating")
    shinyjs::hide("gated_histogram_on_page")
    uploaded_df(NULL)
  })

  # ---- Update dropdowns when data changes ----------------------------------

  observe({
    req(uploaded_df())
    pts <- unique_patients()
    updateRadioButtons(session, "patient_number_plot2", choices = pts)
    updateRadioButtons(session, "trivar_patient",       choices = pts)
    updateRadioButtons(session, "patients_ts",          choices = pts)
    updateCheckboxGroupInput(session, "unique_patients_gating",
                             choices  = pts,
                             selected = pts)
  })

  observe({
    req(uploaded_df())
    mkrs <- unique_markers()
    updateSelectInput(session, "xvar",    choices = mkrs)
    updateSelectInput(session, "yvar",    choices = mkrs)
    updateSelectInput(session, "xvarTri", choices = mkrs)
    updateSelectInput(session, "yvarTri", choices = mkrs)
    updateSelectInput(session, "zvar",    choices = mkrs)
    updateSelectInput(session, "marker",  choices = mkrs)
  })

  observe({
    req(uploaded_df())
    updateRadioButtons(session, "cycle_ts", choices = cycle_detected())
  })

  # ---- Run gating ----------------------------------------------------------

  observeEvent(input$run_gate, {
    req(uploaded_df())

    sel_cols <- input$selected_columns
    sel_pats <- input$unique_patients_gating

    req(length(sel_cols) > 0, length(sel_pats) > 0)

    sub_data <- uploaded_df() %>%
      filter(imageid %in% sel_pats) %>%
      select(imageid, all_of(sel_cols))

    result_df <- get_gates_csv(sub_data)
    resultdf_reactive(result_df)
    csv_save_file_path(result_df)

    dlog("Gates computed:\n", capture.output(print(result_df)))

    # Vectorised positivity labelling (replaces sapply row-by-row)
    df_pos <- uploaded_df()
    for (pat in sel_pats) {
      pat_rows <- df_pos$imageid == pat
      for (mcol in sel_cols) {
        gv <- result_df$Gate[result_df$Marker == mcol & result_df$Patient == pat]
        if (length(gv) == 1 && is.finite(gv))
          df_pos[[paste0(mcol, "_positivity")]] <-
            ifelse(df_pos[[mcol]] > gv, "+", "-")
      }
    }
    uploaded_df(df_pos)

    shinyjs::show("gated_histogram_on_page")

    shinyalert::shinyalert(
      title = "Randkluft Found",
      text  = "Gates saved. You may download them in the Download section.",
      type  = "success"
    )
  })

  # ---- Switch handlers -----------------------------------------------------

  observeEvent(input$GMM_gate_on_off, { gmm_gate_switch(input$GMM_gate_on_off) })

  # ---- Pagination ----------------------------------------------------------
  # Updating current_patient / current_marker (reactiveVals) automatically
  # invalidates any renderPlot that reads them, so the plot rerenders.

  rerender_plot <- function() {
    if (histplot_gate_switch()) {
      output$gated_histogram_on_page <- renderPlot({ generatePlot() })
    }
  }

  observeEvent(input$nextPatient, {
    outputinterceptreactive(NULL)
    updateNumericInput(session, "intercept", value = NA)
    n <- length(unique_patients_gating())
    if (n > 0) {
      current_patient((current_patient() %% n) + 1)
      rerender_plot()
    }
  })

  observeEvent(input$prevPatient, {
    outputinterceptreactive(NULL)
    updateNumericInput(session, "intercept", value = NA)
    n <- length(unique_patients_gating())
    if (n > 0) {
      current_patient(ifelse(current_patient() == 1, n, current_patient() - 1))
      rerender_plot()
    }
  })

  observeEvent(input$nextMarker, {
    outputinterceptreactive(NULL)
    updateNumericInput(session, "intercept", value = NA)
    n <- length(selected_columns())
    if (n > 0) {
      current_marker((current_marker() %% n) + 1)
      rerender_plot()
    }
  })

  observeEvent(input$prevMarker, {
    outputinterceptreactive(NULL)
    updateNumericInput(session, "intercept", value = NA)
    n <- length(selected_columns())
    if (n > 0) {
      current_marker(ifelse(current_marker() == 1, n, current_marker() - 1))
      rerender_plot()
    }
  })

  # ---- Manual gate update --------------------------------------------------

  observeEvent(input$updateGates, {
    new_val <- input$intercept
    if (!is.null(new_val) && !is.na(new_val)) {
      outputinterceptreactive(new_val)

      # Persist updated gate in the results dataframe
      sel_cols <- input$selected_columns
      sel_pats <- input$unique_patients_gating
      pat      <- sel_pats[current_patient()]
      mkr      <- sel_cols[current_marker()]

      df_gates <- csv_save_file_path()
      if (!is.null(df_gates)) {
        df_gates$Gate[df_gates$Marker == mkr & df_gates$Patient == pat] <- new_val
        csv_save_file_path(df_gates)
        resultdf_reactive(df_gates)
      }
    } else {
      outputinterceptreactive(NULL)
    }
    # Re-register renderPlot so the updated gate is reflected immediately
    if (histplot_gate_switch()) {
      output$gated_histogram_on_page <- renderPlot({ generatePlot() })
    }
  })

  # ---- Current position label ----------------------------------------------

  output$current_position_label <- renderUI({
    sel_cols <- input$selected_columns
    sel_pats <- input$unique_patients_gating
    if (is.null(sel_cols) || is.null(sel_pats)) return(NULL)
    cp <- current_patient()
    cm <- current_marker()
    if (cp > length(sel_pats) || cm > length(sel_cols)) return(NULL)
    tags$div(
      class = "current-label",
      sprintf("Patient: %s  |  Marker: %s  |  (%d / %d patients, %d / %d markers)",
              sel_pats[cp], sel_cols[cm],
              cp, length(sel_pats),
              cm, length(sel_cols))
    )
  })

  # ---- Main 4-panel plot ---------------------------------------------------

  generatePlot <- function() {
    req(uploaded_df(), resultdf_reactive())

    sel_cols <- input$selected_columns
    sel_pats <- input$unique_patients_gating
    df       <- uploaded_df()

    cp  <- current_patient()
    cm  <- current_marker()
    req(cp <= length(sel_pats), cm <= length(sel_cols))

    chosen_patient <- sel_pats[cp]
    chosen_marker  <- sel_cols[cm]

    col_idx <- which(names(df) == chosen_marker)

    # Histogram with gate
    histogram_plot <- plot_gating_grid(df, chosen_patient, col_idx)

    # Spatial data for this patient
    fxy <- df[df$imageid == chosen_patient, ]

    gate_val <- resultdf_reactive()$Gate[
      resultdf_reactive()$Marker  == chosen_marker &
      resultdf_reactive()$Patient == chosen_patient
    ]
    if (length(gate_val) != 1 || !is.finite(gate_val)) gate_val <- NA
    if (!is.null(outputinterceptreactive())) gate_val <- outputinterceptreactive()

    # Digital representation
    digrepresentation <- ggplot(fxy, aes(x = X_centroid, y = Y_centroid,
                                          color = .data[[chosen_marker]])) +
      geom_point(shape = 20, size = 0.3, alpha = 0.5) +
      scale_color_gradient(low = "grey30", high = "white") +
      theme(
        panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
        axis.title  = element_text(size = 13, family = "sans"),
        axis.text   = element_text(size = 11, family = "sans")
      ) +
      labs(title = "Digital Representation", x = "X Centroid", y = "Y Centroid")

    if (!is.na(gate_val)) {
      density_vals <- tryCatch(
        get_density(fxy$X_centroid, fxy$Y_centroid, n = 100),
        error = function(e) rep(0, nrow(fxy))
      )

      contour_plot <- ggplot(fxy, aes(x = X_centroid, y = Y_centroid)) +
        geom_point(
          aes(color = ifelse(.data[[chosen_marker]] > gate_val, density_vals, 0)),
          size = 0.3, alpha = 0.35
        ) +
        geom_density_2d(
          data  = fxy[fxy[[chosen_marker]] > gate_val, ],
          color = "black"
        ) +
        scale_color_viridis(option = "turbo", name = "Intensity") +
        theme(
          legend.position    = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.background  = element_rect(fill = alpha("white", 0.7), color = "grey70"),
          plot.title  = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
          axis.title  = element_text(size = 13, family = "sans"),
          axis.text   = element_text(size = 11, family = "sans"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white")
        ) +
        labs(title = "Positive Density", x = "X Centroid", y = "Y Centroid")

      overlay_plot <- ggplot(fxy, aes(x = X_centroid, y = Y_centroid)) +
        geom_point(
          aes(color = ifelse(.data[[chosen_marker]] > gate_val, "Positive", "Negative")),
          size = 0.3, alpha = 0.35
        ) +
        scale_color_manual(values = c("Positive" = "red", "Negative" = "grey"),
                           guide  = guide_legend(title = "")) +
        theme(
          legend.position  = "none",
          plot.title  = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
          axis.title  = element_text(size = 13, family = "sans"),
          axis.text   = element_text(size = 11, family = "sans"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white")
        ) +
        labs(title = "Positive Cells", x = "X Centroid", y = "Y Centroid")

    } else {
      contour_plot <- overlay_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Run Randkluft first",
                 size = 6, color = "grey50") +
        theme_void()
    }

    arranged <- grid.arrange(histogram_plot, digrepresentation,
                             overlay_plot,   contour_plot,
                             ncol = 2)

    # Cache panels for download (defined once outside this function)
    current_4panel$p1 <- histogram_plot
    current_4panel$p2 <- overlay_plot
    current_4panel$p3 <- digrepresentation
    current_4panel$p4 <- contour_plot

    arranged
  }

  # ---- Render the main plot -----------------------------------------------
  # Single observer handles flag update + renderPlot registration + show/hide.
  # Matches the original v22 pattern exactly.

  observeEvent(input$gen_hist_plots_on_off, {
    req(uploaded_df())
    histplot_gate_switch(input$gen_hist_plots_on_off)
    if (input$gen_hist_plots_on_off) {
      output$gated_histogram_on_page <- renderPlot({
        generatePlot()
      })
      shinyjs::show("gated_histogram_on_page")
    } else {
      shinyjs::hide("gated_histogram_on_page")
    }
  })

  # ---- Download current 4-panel (defined once) ----------------------------
  output$downloadcurrent <- downloadHandler(
    filename = function() "plots.pdf",
    content  = function(file) {
      pdf(file, width = 11, height = 11, onefile = TRUE)
      if (!is.null(current_4panel$p1))
        grid.arrange(current_4panel$p1, current_4panel$p2,
                     current_4panel$p3, current_4panel$p4, ncol = 2)
      dev.off()
    }
  )

  # ---- Download all plots --------------------------------------------------

  generate_four_panel_plot <- function(data, imageid_val, marker) {
    col_idx  <- which(names(data) == marker)
    hist_p   <- plot_gating_grid(data, imageid_val, col_idx)
    fxy      <- uploaded_df()[uploaded_df()$imageid == imageid_val, ]

    gate_val <- tryCatch(
      resultdf_reactive()$Gate[
        resultdf_reactive()$Marker  == marker &
        resultdf_reactive()$Patient == imageid_val
      ],
      error = function(e) NA_real_
    )
    if (length(gate_val) != 1 || !is.finite(gate_val)) gate_val <- NA

    dig <- ggplot(fxy, aes(x = X_centroid, y = Y_centroid,
                            color = .data[[marker]])) +
      geom_point(shape = 20, size = 0.3) +
      scale_color_gradient(low = "grey30", high = "white") +
      theme(
        panel.background = element_rect(fill = "black"),
        panel.grid       = element_blank(),
        legend.position  = "none",
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 14, family = "sans"),
        axis.title  = element_text(size = 12, family = "sans"),
        axis.text   = element_text(size = 10, family = "sans")
      ) +
      labs(title = "Digital Representation")

    if (!is.na(gate_val)) {
      dens_vals <- tryCatch(
        get_density(fxy$X_centroid, fxy$Y_centroid, n = 100),
        error = function(e) rep(0, nrow(fxy))
      )
      cont_p <- ggplot(fxy, aes(x = X_centroid, y = Y_centroid)) +
        geom_point(
          aes(color = ifelse(.data[[marker]] > gate_val, dens_vals, 0)),
          size = 0.3
        ) +
        geom_density_2d(data  = fxy[fxy[[marker]] > gate_val, ], color = "black") +
        scale_color_viridis(option = "turbo") +
        theme(
          legend.position = "bottom",
          plot.title  = element_text(face = "bold", hjust = 0.5, size = 14, family = "sans"),
          axis.title  = element_text(size = 12, family = "sans"),
          axis.text   = element_text(size = 10, family = "sans")
        ) +
        labs(title = "Positive Density", x = "X Centroid", y = "Y Centroid")

      over_p <- ggplot(fxy, aes(x = X_centroid, y = Y_centroid)) +
        geom_point(
          aes(color = ifelse(.data[[marker]] > gate_val, "Positive", "Negative")),
          size = 0.3
        ) +
        scale_color_manual(values = c("Positive" = "red", "Negative" = "grey")) +
        theme(
          legend.position = "bottom",
          plot.title  = element_text(face = "bold", hjust = 0.5, size = 14, family = "sans"),
          axis.title  = element_text(size = 12, family = "sans"),
          axis.text   = element_text(size = 10, family = "sans")
        ) +
        labs(title = "Positive Cells", x = "X Centroid", y = "Y Centroid")
    } else {
      cont_p <- over_p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No gate", size = 5, color = "grey60") +
        theme_void()
    }

    grid.arrange(hist_p, over_p, dig, cont_p, ncol = 2)
  }

  generate_histogram_pdf <- function(subdata, pdf_file_name) {
    unique_ids <- unique(subdata$imageid)
    markers    <- colnames(subdata)[-1]
    plots      <- lapply(unique_ids, function(uid) {
      lapply(markers, function(m) generate_four_panel_plot(subdata, uid, m))
    })
    plots <- unlist(plots, recursive = FALSE)

    pdf(pdf_file_name, width = 14, height = 10)
    for (p in plots) { do.call(grid.arrange, list(p, ncol = 1)) }
    dev.off()
    pdf_file_path(pdf_file_name)
  }

  output$downloadall <- downloadHandler(
    filename = "all_plots.pdf",
    content  = function(file) {
      sel_cols <- input$selected_columns
      sel_pats <- input$unique_patients_gating
      sub      <- uploaded_df() %>%
                    filter(imageid %in% sel_pats) %>%
                    select(imageid, all_of(sel_cols))
      generate_histogram_pdf(sub, file)
    }
  )

  # ---- Bivariate plot ------------------------------------------------------

  subsetted <- reactive({
    req(input$patient_number_plot2, uploaded_df())
    uploaded_df() %>% filter(imageid %in% input$patient_number_plot2)
  })

  subsetted_tri <- reactive({
    req(input$trivar_patient, uploaded_df())
    uploaded_df() %>% filter(imageid %in% input$trivar_patient)
  })

  observeEvent(c(input$xvar, input$yvar), {
    updateNumericInput(session, "gate_xvar_update", value = NA)
    updateNumericInput(session, "gate_yvar_update", value = NA)
    react_xgate(NULL)
    react_ygate(NULL)
  })

  observeEvent(input$update_gates_bivariate, {
    xv <- input$gate_xvar_update
    yv <- input$gate_yvar_update
    react_xgate(if (!is.null(xv) && !is.na(xv)) as.numeric(xv) else NULL)
    react_ygate(if (!is.null(yv) && !is.na(yv)) as.numeric(yv) else NULL)
    output$plot2 <- renderPlot({ buildBivarPlot() })
  })

  buildBivarPlot <- function() {
    req(subsetted(), resultdf_reactive())
    dfp    <- subsetted()
    pat    <- input$patient_number_plot2
    xvar   <- as.character(input$xvar)
    yvar   <- as.character(input$yvar)

    get_gate <- function(var) {
      gv <- resultdf_reactive()$Gate[
        resultdf_reactive()$Marker  == var &
        resultdf_reactive()$Patient == pat
      ]
      if (length(gv) == 1 && is.finite(gv)) gv else NA_real_
    }

    gx <- if (!is.null(react_xgate())) react_xgate() else get_gate(xvar)
    gy <- if (!is.null(react_ygate())) react_ygate() else get_gate(yvar)

    # Filter clean rows
    dfp <- dfp[complete.cases(dfp[[xvar]], dfp[[yvar]]) &
               is.finite(dfp[[xvar]]) & is.finite(dfp[[yvar]]), ]

    n      <- nrow(dfp)
    pp     <- sum(dfp[[xvar]] > gx & dfp[[yvar]] > gy,  na.rm = TRUE) / n
    pm     <- sum(dfp[[xvar]] > gx & dfp[[yvar]] <= gy, na.rm = TRUE) / n
    mp     <- sum(dfp[[xvar]] <= gx & dfp[[yvar]] > gy, na.rm = TRUE) / n
    mm     <- sum(dfp[[xvar]] <= gx & dfp[[yvar]] <= gy, na.rm = TRUE) / n

    dens   <- tryCatch(get_density(dfp[[xvar]], dfp[[yvar]], n = 100),
                       error = function(e) rep(0, nrow(dfp)))

    p <- ggplot(dfp, aes(x = .data[[xvar]], y = .data[[yvar]])) +
      geom_point(aes(color = dens), size = 0.3, alpha = 0.35) +
      geom_vline(xintercept = gx, color = "orange") +
      geom_hline(yintercept = gy, color = "orange") +
      scale_color_viridis(option = "turbo") +
      theme(
        legend.position  = "none",
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
        axis.title  = element_text(size = 13, family = "sans"),
        axis.text   = element_text(size = 11, family = "sans"),
        panel.background = element_rect(fill = "#f8f8f8"),
        plot.background  = element_rect(fill = "#f8f8f8")
      ) +
      labs(title = "Bivariate Density") +
      annotate("text",
               x = max(dfp[[xvar]]) * 0.97, y = max(dfp[[yvar]]) * 0.97,
               label = round(pp, 3), color = "red",   size = 5) +
      annotate("text",
               x = max(dfp[[xvar]]) * 0.97, y = min(dfp[[yvar]]) * 1.03,
               label = round(pm, 3), color = "green", size = 5) +
      annotate("text",
               x = min(dfp[[xvar]]) * 1.03, y = max(dfp[[yvar]]) * 0.97,
               label = round(mp, 3), color = "blue",  size = 5) +
      annotate("text",
               x = min(dfp[[xvar]]) * 1.03, y = min(dfp[[yvar]]) * 1.03,
               label = round(mm, 3), color = "black", size = 5)

    overlay <- ggplot(dfp, aes(x = X_centroid, y = Y_centroid)) +
      geom_point(
        aes(color = case_when(
          .data[[xvar]] > gx & .data[[yvar]] > gy   ~ "+/+",
          .data[[xvar]] <= gx & .data[[yvar]] <= gy ~ "-/-",
          .data[[xvar]] > gx & .data[[yvar]] <= gy  ~ "+/-",
          .data[[xvar]] <= gx & .data[[yvar]] > gy  ~ "-/+",
          TRUE ~ "Other"
        )),
        size = 0.3, alpha = 0.7
      ) +
      scale_color_manual(
        values = c("+/+" = "red", "-/-" = "grey", "+/-" = "green",
                   "-/+" = "blue", "Other" = "black"),
        guide  = guide_legend(title = "", override.aes = list(size = 5))
      ) +
      labs(title = "Bivariate Gating", x = "X Centroid", y = "Y Centroid") +
      theme(
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
        axis.title  = element_text(size = 13, family = "sans"),
        axis.text   = element_text(size = 11, family = "sans"),
        legend.position  = "none",
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white")
      )

    px_pct <- sum(dfp[[xvar]] > gx, na.rm = TRUE) / nrow(dfp)
    py_pct <- sum(dfp[[yvar]] > gy, na.rm = TRUE) / nrow(dfp)

    cx <- ggplot(dfp, aes(x = X_centroid, y = Y_centroid)) +
      geom_point(aes(color = ifelse(.data[[xvar]] > gx, "Positive", "Negative")),
                 size = 0.3, alpha = 0.35) +
      scale_color_manual(values = c("Positive" = "green", "Negative" = "grey"),
                         guide  = guide_legend(title = "")) +
      geom_density_2d(data = dfp[dfp[[xvar]] > gx, ], color = "black") +
      labs(title = paste0(xvar, "+ = ", round(px_pct, 3)),
           x = "X Centroid", y = "Y Centroid") +
      theme(
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 14, family = "sans"),
        axis.title  = element_text(size = 12, family = "sans"),
        axis.text   = element_text(size = 10, family = "sans"),
        legend.position  = "bottom",
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white")
      )

    cy <- ggplot(dfp, aes(x = X_centroid, y = Y_centroid)) +
      geom_point(aes(color = ifelse(.data[[yvar]] > gy, "Positive", "Negative")),
                 size = 0.3, alpha = 0.35) +
      scale_color_manual(values = c("Positive" = "blue", "Negative" = "grey"),
                         guide  = guide_legend(title = "")) +
      geom_density_2d(data = dfp[dfp[[yvar]] > gy, ], color = "black") +
      labs(title = paste0(yvar, "+ = ", round(py_pct, 3)),
           x = "X Centroid", y = "Y Centroid") +
      theme(
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 14, family = "sans"),
        axis.title  = element_text(size = 12, family = "sans"),
        axis.text   = element_text(size = 10, family = "sans"),
        legend.position  = "bottom",
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white")
      )

    grid.arrange(p, overlay, cx, cy, ncol = 2)
  }

  output$plot2 <- renderPlot({ buildBivarPlot() })

  # ---- Trivariate plot -----------------------------------------------------

  output$plot_trivariate <- renderPlotly({
    req(subsetted_tri(), resultdf_reactive())
    dfp <- subsetted_tri()
    pat <- input$trivar_patient
    xv  <- as.character(input$xvarTri)
    yv  <- as.character(input$yvarTri)
    zv  <- as.character(input$zvar)

    get_g <- function(v) {
      g <- resultdf_reactive()$Gate[
        resultdf_reactive()$Marker  == v &
        resultdf_reactive()$Patient == pat
      ]
      if (length(g) == 1 && is.finite(g)) g else NA_real_
    }

    plot_ly(dfp,
            x = ~dfp[[xv]], y = ~dfp[[yv]], z = ~dfp[[zv]],
            type = "scatter3d", mode = "markers",
            marker = list(size = 2, opacity = 0.5)) %>%
      layout(scene = list(
        xaxis = list(title = xv),
        yaxis = list(title = yv),
        zaxis = list(title = zv)
      ))
  })

  # ---- Image garage (Pre-Randkluft spatial viewer) -------------------------

  observeEvent(c(input$precrev_marker, input$patient_number_im_gargage), {
    req(uploaded_df(), input$patient_number_im_gargage, input$precrev_marker)

    pat     <- input$patient_number_im_gargage
    mkr     <- input$precrev_marker
    fd      <- uploaded_df()[uploaded_df()$imageid %in% pat, ]

    # Keep imageid as first column
    fd      <- fd[, c("imageid", setdiff(colnames(fd), "imageid"))]
    filtered_data_original(fd)
    filtered_data_reactive(fd)

    output$image_garage_output <- renderPlotly({
      req(filtered_data_reactive())
      p <- ggplot(filtered_data_reactive(),
                  aes(x = X_centroid, y = Y_centroid,
                      color = .data[[mkr]])) +
        geom_point(alpha = if (!is.null(input$opacity_slider)) input$opacity_slider else 0.5) +
        xlab("X Centroid") + ylab("Y Centroid") +
        labs(title = "Digital Marker Overlay") +
        scale_color_continuous(name = "Intensity (logged)") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5,
                                        size = 16, family = "sans"))
      # Register event so plotly_selected warning is suppressed
      ggplotly(p) %>% event_register("plotly_selected")
    })
  })

  observeEvent(input$RetrySubset, {
    filtered_data_reactive(filtered_data_original())
  })

  observeEvent(event_data("plotly_selected"), {
    pts <- event_data("plotly_selected")
    if (!is.null(pts)) {
      fd <- filtered_data_reactive() %>%
        filter(!X_centroid %in% pts$x & !Y_centroid %in% pts$y)
      filtered_data_reactive(fd)
    }
  })

  observeEvent(input$UseSubset, {
    uploaded_df(filtered_data_reactive())
    shinyalert::shinyalert(
      title = "Subset Loaded",
      text  = "You may use this subset in the Randkluft tab.",
      type  = "success"
    )
  })

  # ---- UMAP ----------------------------------------------------------------

  output$icaAnalysis <- renderPlot({
    req(uploaded_df())
    col_names <- get_marker_cols(uploaded_df(), include_dna = FALSE)
    df_u      <- uploaded_df()[, col_names, drop = FALSE]

    # Replace infinites
    df_u <- as.data.frame(lapply(df_u, function(x) { x[!is.finite(x)] <- 0; x }))
    mat  <- data.matrix(df_u)

    umap_res <- tryCatch(umap(mat), error = function(e) NULL)
    req(!is.null(umap_res))

    # Safe cbind: rename layout columns to avoid collision
    layout_df <- as.data.frame(umap_res$layout)
    colnames(layout_df) <- c("UMAP_1", "UMAP_2")
    umap_df   <- cbind(uploaded_df(), layout_df)

    ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(size = 0.5, alpha = 0.5, color = "steelblue") +
      labs(title = "UMAP", x = "UMAP 1", y = "UMAP 2") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, family = "sans"),
        axis.title = element_text(family = "sans"),
        axis.text  = element_text(family = "sans")
      )
  })

  # ---- Tissue scoring (cycle quality) -------------------------------------

  output$tissue_score_out <- renderPlot({
    req(uploaded_df(), input$patients_ts, input$cycle_ts)
    sel_pat <- input$patients_ts
    fd      <- uploaded_df() %>% filter(imageid %in% sel_pat)

    dna_cols <- get_dna_cols(fd)
    if (length(dna_cols) == 0) return(NULL)

    col_idx <- sapply(dna_cols, function(cn) which(names(fd) == cn))

    # Mode per cycle
    modes <- sapply(dna_cols, function(cn) {
      vals <- remove_outliers2(fd[[cn]])
      vals <- vals[is.finite(vals)]
      find_mode2(vals)
    })
    modes_df <- data.frame(Index = seq_along(modes), Value = modes)

    plot_modes <- ggplot(modes_df, aes(x = Index, y = Value)) +
      geom_point() + geom_line() +
      labs(x = "Cycle", y = "Mode", title = "Regression of Quality") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
        axis.title = element_text(family = "sans"),
        axis.text  = element_text(family = "sans")
      ) +
      scale_y_continuous(labels = scales::comma)

    regression_mode_react(plot_modes)

    # KL divergence between consecutive cycles
    n_cols <- length(dna_cols)
    if (n_cols >= 2) {
      kl_list <- lapply(seq_len(n_cols - 1), function(i) {
        v1 <- exp(fd[[dna_cols[i]]])
        v2 <- exp(fd[[dna_cols[i + 1]]])
        p1 <- v1 / sum(v1)
        p2 <- v2 / sum(v2)
        kl <- tryCatch(KL(rbind(p1, p2), unit = "log"), error = function(e) NA)
        data.frame(Index = i, KL_Divergence = if (length(kl) == 1) kl else NA)
      })
      kl_df <- bind_rows(kl_list)

      kl_plot <- ggplot(kl_df, aes(x = Index, y = KL_Divergence)) +
        geom_point() + geom_line() +
        labs(x = "Pair", y = "KL Divergence", title = "KL Divergence Trend") +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans"),
          axis.title = element_text(family = "sans"),
          axis.text  = element_text(family = "sans")
        )
      klPlotReact(kl_plot)
    }

    # Gauge
    first_cycle   <- modes_df$Value[1]
    min_val_cycle <- min(modes_df$Value, na.rm = TRUE)
    gauge <- plot_ly(
      domain = list(x = c(0, 1), y = c(0, 1)),
      value  = min_val_cycle,
      title  = list(text = "Data Quality – Mode Deviation"),
      type   = "indicator",
      mode   = "gauge+number",
      gauge  = list(
        steps = list(
          list(range = c(0,                0.2 * first_cycle), color = "red"),
          list(range = c(0.2 * first_cycle, 0.4 * first_cycle), color = "orange"),
          list(range = c(0.4 * first_cycle, 0.6 * first_cycle), color = "yellow"),
          list(range = c(0.6 * first_cycle, 0.8 * first_cycle), color = "green"),
          list(range = c(0.8 * first_cycle, first_cycle),        color = "darkgreen")
        ),
        axis      = list(range = list(0, first_cycle)),
        bar       = list(color = "black"),
        threshold = list(line = list(color = "black", width = 4),
                         thickness = 0.75, value = min_val_cycle)
      )
    )
    output$quality_gauge <- renderPlotly({ gauge })

    plot_modes
  })

  output$modes_plot_output <- renderPlot({ regression_mode_react() })
  output$klPLOToutput      <- renderPlot({ klPlotReact() })

  # ---- Phenotyping ---------------------------------------------------------

  observeEvent(input$phen_wfl, {
    req(input$phen_wfl)
    df <- read.csv(input$phen_wfl$datapath, stringsAsFactors = FALSE)
    if ("phenotype" %in% names(df) && "markers" %in% names(df))
      phenotype_df(df)
  })

  output$phenotypeTable <- renderTable({ phenotype_df() })

  output$marker_checkboxes <- renderUI({
    mkrs <- input$selected_columns
    if (is.null(mkrs) || length(mkrs) == 0) return(NULL)
    tagList(lapply(mkrs, function(m) {
      fluidRow(
        column(2, checkboxInput(paste0("checkbox_", m), m)),
        br(),
        column(2, materialSwitch(paste0("switch_", m), label = NULL, status = "success"))
      )
    }))
  })

  selected_markers_phenotype <- reactive({
    mkrs     <- input$selected_columns
    statuses <- sapply(mkrs, function(m) {
      cb <- input[[paste0("checkbox_", m)]]
      if (!is.null(cb) && cb) {
        sw <- input[[paste0("switch_", m)]]
        if (!is.null(sw) && sw) paste0(m, "+") else paste0(m, "-")
      } else NULL
    })
    statuses[!sapply(statuses, is.null)]
  })

  cleanColumnNames <- function(df) {
    names(df) <- gsub("\\.+|\\s+", "", names(df))
    df
  }

  updatePhenotypeDF3 <- function(pname, pmarkers) {
    parts   <- unlist(strsplit(pmarkers, "\n"))[2]
    cur_df  <- phenotype_df()
    new_row <- data.frame(phenotype = pname, markers = parts,
                          stringsAsFactors = FALSE)
    phenotype_df(bind_rows(cur_df, new_row))
  }

  definePhenotype <- function(markers, pname) {
    markers <- markers[!sapply(markers, is.null)]
    if (length(markers) > 0)
      paste(pname, paste(markers, collapse = ", "), sep = "\n")
    else ""
  }

  cleanPhenotype3 <- function(phenotype_output) {
    parts <- unlist(strsplit(phenotype_output, "\n"))
    cleaned <- lapply(parts, function(p) {
      ps <- unlist(strsplit(p, ", "))
      ps <- ps[ps != "NULL"]
      paste(ps, collapse = ", ")
    })
    cleaned
  }

  observeEvent(input$define_phenotype, {
    if (input$define_phenotype > 0) {
      phenotype     <- definePhenotype(selected_markers_phenotype(), input$phenotype_name)
      phenotype_out <- paste(cleanPhenotype3(phenotype), collapse = "\n")
      updatePhenotypeDF3(input$phenotype_name, phenotype_out)
    }
  })

  subsetAndCount <- function(df, phen_df) {
    original_df   <- uploaded_df()
    pheno_counts  <- lapply(seq_len(nrow(phen_df)), function(i) {
      markers  <- unlist(strsplit(phen_df[i, "markers"], ", "))
      markers  <- markers[markers != ""]
      pos_mkrs <- markers[grepl("\\+", markers)]
      neg_mkrs <- markers[grepl("-",   markers)]

      cond_pos <- sapply(pos_mkrs, function(m) {
        m <- gsub("[+,]", "", m)
        paste0(m, "_positivity == '+'")
      })
      cond_neg <- sapply(neg_mkrs, function(m) {
        m <- gsub("[-,]", "", m)
        paste0(m, "_positivity == '-'")
      })
      combined <- paste(c(cond_pos, cond_neg), collapse = " & ")

      sub_df <- tryCatch(
        subset(df, eval(parse(text = combined))),
        error = function(e) df[0, ]
      )

      if (nrow(sub_df) > 0) {
        original_df[rownames(sub_df), "phenotype"] <<- phen_df[i, "phenotype"]
        original_df[is.na(original_df[, "phenotype"]), "phenotype"] <<- "other"
      }
      uploaded_df(original_df)
      nrow(sub_df)
    })

    other_cnt    <- nrow(df) - sum(unlist(pheno_counts))
    pheno_counts <- c(pheno_counts, other = other_cnt)
    names(pheno_counts) <- c(phen_df$phenotype, "other")
    pheno_counts
  }

  observeEvent(input$define_phenotype_AUTO, {
    req(uploaded_df())
    counts  <- subsetAndCount(uploaded_df(), phenotype_df())
    bar_df  <- data.frame(
      phenotype  = names(counts),
      count      = unlist(counts)
    )
    bar_df$percentage <- bar_df$count / sum(bar_df$count) * 100

    species <- uploaded_df()$phenotype
    MLEP    <- tryCatch(MLEp(abundance(species)), error = function(e) NA)

    output$post_statistics <- renderText({
      paste("Partition Diversity Estimate:", round(MLEP, 5))
    })

    output$pheno_bar <- renderPlot({
      ggplot(bar_df, aes(fill = phenotype, x = "", y = percentage)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(x = NULL, y = "Percentage", title = "Phenotype Composition") +
        theme(
          plot.title   = element_text(face = "bold", size = 16, hjust = 0.5, family = "sans"),
          axis.text    = element_text(size = 12, family = "sans"),
          axis.title   = element_text(size = 13, family = "sans"),
          legend.text  = element_text(size = 11, family = "sans"),
          legend.title = element_text(size = 12, family = "sans")
        )
    })
  })

  # ---- Downloads -----------------------------------------------------------

  output$downloadEstimations <- downloadHandler(
    filename = function() {
      if (!is.null(input$csv_name) && nchar(input$csv_name) > 0)
        paste0(input$csv_name, ".csv")
      else "Gates.csv"
    },
    content = function(file) {
      req(csv_save_file_path())
      write.csv(csv_save_file_path(), file, row.names = FALSE)
    }
  )

  output$downloadPhenotypes <- downloadHandler(
    filename = function() "Phenotypes_included.csv",
    content  = function(file) {
      req(uploaded_df())
      write.csv(uploaded_df(), file, row.names = FALSE)
    }
  )

  output$downloadPhenotypeTable <- downloadHandler(
    filename = function() "phenotype_table.csv",
    content  = function(file) write.csv(phenotype_df(), file, row.names = FALSE)
  )

  output$download_example_data <- downloadHandler(
    filename = function() "ExampleCellMarkerData.csv",
    content  = function(file) write.csv(df_example, file, row.names = FALSE)
  )

  output$download_example_data2 <- downloadHandler(
    filename = function() "Phenotype_table_example.csv",
    content  = function(file) write.csv(df_example_PHENOTYPE, file, row.names = FALSE)
  )

}) # end server

shinyApp(ui = ui, server = server)
