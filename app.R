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
library(stringr)
library(tidyr)
library(PEkit)
library(sqldf)
library(ggplot2)
library(gridExtra)

# source codes
source("utils.R")
source("R/gating.R")
source("R/phenotyping.R")
source("R/diagnostics.R")
source("R/data-input.R")
source("R/app-ui.R")

df_example = read.csv(file.path("data", "exemplar-001--unmicst_cell.csv"))

ground_truth_gates_loaded = read.csv(file.path("data", "tuulia_data_GT.csv"))

df_example_PHENOTYPE = read.csv(file.path("data", "phenotype_table_help.csv"))

options(shiny.trace = TRUE, shiny.maxRequestSize = 100 * 1024 ^ 3)


ui <- build_app_ui(html_code)


server <- shinyServer(function(input, output, session) {

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
      separation <- marker_separation_diagnostic(target)
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
          subtitle = separation$label,
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
          plot.subtitle = element_text(size = 12, family = "Arial"),
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
        target <- target[is.finite(target)]
        plot_target <- remove_outliers2(target)

        marker <- paste0(c(unique_id, names_columns[j]), collapse = "; ")
        output <- estimate_gate(target, 0.01)
        gtGate <- ground_truth_gates_loaded$Gate[ground_truth_gates_loaded$Patient == unique_id & ground_truth_gates_loaded$Marker == names_columns[j]]
        p_temp <- plot_gating_individual(plot_target, output, marker, gtGate)
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

  updated_df <- gate_results()
  if (!is.null(updated_df)) {
    updated_df$Gate[
      updated_df$Marker == chosen_marker &
      updated_df$Patient == chosen_patient
    ] <- new_gate_value
    gate_results(updated_df)
  }

})


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
    updateCheckboxGroupInput(session, "selected_columns", choices = column_names, selected = column_names)
    updateCheckboxGroupInput(session, "selected_columns_phenotyping", choices = column_names, selected = column_names)
    runjs("document.getElementById('message_loading_data').style.display = 'none';")

  })

  uploaded_df <- reactiveVal(NULL)

  histplot_gate_switch <- reactiveVal(TRUE)  # Initialize with the desired default value

  # Observe the histogram visibility switch.
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
    shinyjs::hide("histogram_plot")

    shinyjs::hide("summary_output")
    shinyjs::hide("image_garage_output")
    shinyjs::hide("selected_columns_phenotyping")


    updateSelectInput(session, "yvar", choices = character(0), selected = character(0))
    updateSelectInput(session, "xvar", choices = character(0), selected = character(0))


  })

  observeEvent(input$cell_file, {
    uploaded_df(prepare_cell_data(read_cell_data(input$cell_file$datapath)))
    showNotification("File Ready for Use", duration = 10, id = "message")
    shinyjs::show("selected_columns")
    shinyjs::show("selected_columns_phenotyping")
  })


  gate_results <- reactiveVal(NULL)


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

    # Show the progress bar


    if (length(selected_columns)==1){
      print("hi")
      print(colnames(sub_data))
      resultdf_to_save = get_gates_csv_single(sub_data)
    }
    else {
      resultdf_to_save = get_gates_csv(sub_data, session)
    }
    print(resultdf_to_save)
    gate_results(resultdf_to_save)

    dataframe_pos <- label_cells_with_gates(uploaded_df(), resultdf_to_save)

    # Print first few rows to check the changes
    print(head(dataframe_pos))

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

  # Define reactive values
current_marker <- reactiveVal(1)
current_patient <- reactiveVal(1)


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
    req(uploaded_df(), gate_results(), histplot_gate_switch())
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

gate_value <- gate_results()$Gate[
  gate_results()$Marker == chosen_marker &
    gate_results()$Patient == chosen_patient
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
  gates <- gate_results()

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
    gate_value <- tryCatch(estimate_gate(target, 0.01)$cutoff, error = function(e) NA_real_)
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

  n_positive <- if (is.na(gate_value)) NA_integer_ else sum(target > gate_value)
  positive_rate <- if (is.na(gate_value)) NA_real_ else round(n_positive / length(target), 3)
  plot_values <- remove_outliers2(target)
  separation <- marker_separation_diagnostic(plot_values)
  plot_df <- data.frame(x = plot_values)

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
        paste("Gate unavailable", separation$label, sep = "\n")
      } else {
        paste0("N+ = ", n_positive, "   +R = ", positive_rate,
               "   Gate = ", round(gate_value, 2), "\n", separation$label)
      },
      x = marker,
      y = "Density"
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = title_size),
      plot.subtitle = element_text(hjust = 0.5, size = title_size * 0.75,
                                   lineheight = 1.1),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    )

  if (length(unique(plot_values)) > 1) {
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

 selected_columns <- reactive({
    input$selected_columns
  })


  unique_patients_gating <- reactive({
    unique_patients()
  })


 selected_opacity <- reactive({
    input$opacity_slider
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


regression_mode_react <- reactiveVal(NULL)
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
      selected_markers <- input$selected_columns
      required_markers <- colnames(normalize_phenotype_definitions(phenotype_df())$states)
      workflow <- export_phenotype_workflow(
        phenotype_df(), unique(c(selected_markers, required_markers))
      )
      write.csv(workflow, file, row.names = FALSE, na = "NA")
    }
  )


      # Functions that return statistics from the models in list, save into csv and download
  output$downloadPhenotypes <- downloadHandler(
    filename = function() {
      if (!is.null(input$csv_name) && nzchar(input$csv_name)) {
        paste0(input$csv_name, ".csv")
      } else {
        "Phenotypes_included.csv"
      }
    },
    content = function(file) {
      req(phenotype_result())
      write.csv(phenotype_result()$cells, file, row.names = FALSE)
    }
  )


      # Functions that return statistics from the models in list, save into csv and download
  output$downloadEstimations <- downloadHandler(
    filename = function() {
      if (!is.null(gate_results())) {
        if (!is.null(input$csv_name)) {
          return(paste0(input$csv_name, ".csv"))
        } else {

          return('Gates.csv')
        }

      }
      showNotification("Run Randkluft first!", duration = NULL, id = "warngate")

    },
    content = function(file) {
      mydftosave = gate_results()
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
    dfplot2 <- dfplot2[keep_marker_inlier_rows(dfplot2, c(xvar, yvar)), , drop = FALSE]
    if (nrow(dfplot2) == 0) {
      return(bivariate_pair_error(xvar, yvar, "No values remain after outlier removal."))
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
    marker_data <- marker_data[keep_marker_inlier_rows(marker_data, marker), , drop = FALSE]
    if (nrow(marker_data) == 0) {
      return(make_pdf_message_grob(marker, "No values remain after outlier removal."))
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
    opacity = 0.2,
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
  dfplot2 <- dfplot2[keep_marker_inlier_rows(dfplot2, c(xvar, yvar, zvar)), , drop = FALSE]
  validate(need(nrow(dfplot2) > 0, "No marker values remain after outlier removal for the trivariate plot."))

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


  selected_markers_phenotype <- reactive({
    markers <- input$selected_columns
    stats::setNames(vapply(markers, function(marker) {
      if (!isTRUE(input[[paste0("checkbox_", marker)]])) return(NA_character_)
      if (isTRUE(input[[paste0("switch_", marker)]])) "+" else "-"
    }, character(1)), markers)
  })

  phenotype_df <- reactiveVal(data.frame(
    phenotype = character(), markers = character(), match_mode = character()
  ))
  phenotype_result <- reactiveVal(NULL)

  observeEvent(uploaded_df(), phenotype_result(NULL), ignoreInit = TRUE)
  observeEvent(phenotype_df(), phenotype_result(NULL), ignoreInit = TRUE)
  observeEvent(gate_results(), phenotype_result(NULL), ignoreInit = TRUE)

  observeEvent(input$define_phenotype, {
    tryCatch({
      selected <- selected_markers_phenotype()
      selected <- selected[!is.na(selected)]
      candidate <- rbind(phenotype_df(), data.frame(
        phenotype = trimws(input$phenotype_name),
        markers = paste0(names(selected), selected, collapse = ", "),
        match_mode = if (isTRUE(input$any_indicator)) "any_positive" else "all"
      ))
      phenotype_df(normalize_phenotype_definitions(candidate)$definitions)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  observeEvent(input$phen_wfl, {
    tryCatch({
      workflow <- read.csv(input$phen_wfl$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      phenotype_df(normalize_phenotype_definitions(workflow)$definitions)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  observeEvent(input$remove_file2, {
    phenotype_df(data.frame(phenotype = character(), markers = character(), match_mode = character()))
    reset("phen_wfl")
  })

  observeEvent(input$define_phenotype_AUTO, {
    result <- tryCatch(
      phenotype_partition(uploaded_df(), phenotype_df(), gate_results()),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
        NULL
      }
    )
    phenotype_result(result)
  })

  output$phenotypeTable <- renderTable(phenotype_df())
  output$phenotypeSummary <- renderTable({
    result <- phenotype_result()
    req(result)
    within(result$summary, percentage <- round(percentage, 1))
  })

  output$pheno_bar_ui <- renderUI({
    result <- phenotype_result()
    req(result)
    label_lines <- vapply(result$summary$phenotype, function(label) {
      max(1L, length(strsplit(wrap_phenotype_label(label, 25L), "\n", fixed = TRUE)[[1]]))
    }, integer(1))
    plot_height <- max(380L, 110L + sum(28L * label_lines + 24L))
    plotOutput("pheno_bar", height = paste0(plot_height, "px"), width = "100%")
  })

  output$pheno_bar <- renderPlot({
    result <- phenotype_result()
    req(result)
    plot_phenotype_composition(result$summary)
  })

  output$post_statistics <- renderText({
    result <- phenotype_result()
    req(result)
    estimate <- if (is.na(result$diversity)) "not estimable (no defined groups)" else
      if (is.infinite(result$diversity)) "unbounded (all included cells have distinct groups)" else
        format(round(result$diversity, 5))
    paste0("Partition Diversity Estimate (PEkit): ", estimate,
           "\nIncluded cells: ", result$n_index,
           "; observed intersection groups: ", result$k_index,
           "\nOther and Unresolved cells are excluded from this estimate.")
  })


      })

shinyApp(ui = ui, server = server)
