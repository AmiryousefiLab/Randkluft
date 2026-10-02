# Read and prepare an uploaded cell table. These functions keep the existing
# sampling and transformation rules; they do not impose new validation policy.

read_cell_data <- function(file_path) {
  n_rows <- length(count.fields(file_path))
  csv_columns <- names(read.csv(file_path, nrows = 0, check.names = FALSE))

  if (n_rows > 70000) {
    sample_numbers <- sample(1:n_rows, 70000, replace = FALSE)
  } else {
    sample_numbers <- seq(1:n_rows)
  }

  if ("CellID" %in% csv_columns) {
    query <- sprintf(
      "SELECT * FROM file WHERE CellID IN (%s)",
      paste(sample_numbers, collapse = ", ")
    )
    data <- read.csv.sql(file_path, sql = query, dbname = tempfile())
  } else {
    data <- read.csv(file_path, header = TRUE, check.names = FALSE)
    if (nrow(data) > 70000) {
      data <- data[sample(seq_len(nrow(data)), 70000), , drop = FALSE]
    }
  }

  print("DONE LOADING CSV INITIAL")
  data
}

prepare_cell_data <- function(data) {
  if ("imageID" %in% colnames(data)) {
    colnames(data)[colnames(data) == "imageID"] <- "imageid"
  } else {
    print("Column not found")
  }
  print(colnames(data))

  empty_string_cols <- which(colnames(data) == "")
  if (length(empty_string_cols) > 0) {
    data <- data[, -empty_string_cols, drop = FALSE]
  }

  print("colnames after drop blank")
  print(colnames(data))
  columns_to_exclude <- c(
    "imageid", "phenotype", "ROI_major_category", "CellID", "Cell", "Row",
    "X", "Y", "ROI_minor_category", "phenotype_v2", "X_centroid",
    "Y_centroid", "Eccentricity", "Area", "MajorAxisLength",
    "MinorAxisLength", "Extent", "Solidity", "Orientation", "", "DNA6a"
  )
  column_names <- setdiff(names(data), columns_to_exclude)
  print(column_names)

  print("pre")
  max_values <- sapply(data[column_names], max)
  print("post")

  data[column_names] <- sapply(data[column_names], as.numeric)
  data[column_names] <- data.frame(lapply(data[column_names], function(x) {
    ifelse(x >= 0 & x < 10, 10, x)
  }))

  if (any(max_values > 20)) {
    cat("Working with raw data\n")
    data[column_names] <- sapply(data[column_names], log)
  } else {
    cat("Working with logged values\n")
  }

  data$imageid <- rep("sample", nrow(data))
  if (any(c("DNA1", "DNA_1", "DAPI1", "DAPI_1", "Hoechst1", "Hoechst_1") %in% colnames(data))) {
    data$DNA1 <- rep(1, nrow(data))
  }

  if ("X" %in% colnames(data)) {
    colnames(data)[which(names(data) == "X")] <- "X_centroid"
  }
  if ("Y" %in% colnames(data)) {
    colnames(data)[which(names(data) == "Y")] <- "Y_centroid"
  }

  if (!"X_centroid" %in% colnames(data)) {
    grid_width <- max(1, ceiling(sqrt(nrow(data))))
    data$X_centroid <- ((seq_len(nrow(data)) - 1) %% grid_width) + 1
  }
  if (!"Y_centroid" %in% colnames(data)) {
    grid_width <- max(1, ceiling(sqrt(nrow(data))))
    data$Y_centroid <- ((seq_len(nrow(data)) - 1) %/% grid_width) + 1
  }

  print("DONE LOADING CSV FINAL after POST")
  print(colnames(data))
  data
}
