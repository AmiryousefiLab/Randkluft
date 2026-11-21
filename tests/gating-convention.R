# Run from the app directory: Rscript --vanilla tests/gating-convention.R
# Stub only the statistical search; verify its input and the resulting labels.
source("R/gating.R")

seen_inputs <- list()
skew_gate <- function(x, alpha = 0.01) {
  seen_inputs[[length(seen_inputs) + 1L]] <<- x
  list(cutoff = median(x))
}
Progress <- list(new = function(...) {
  list(set = function(...) invisible(NULL), close = function(...) invisible(NULL))
})

values <- c(seq_len(100), -1000, 1000, NA_real_, Inf)
expected <- remove_outliers2(values[is.finite(values)])
single_data <- data.frame(imageid = rep("sample", length(values)), marker = values)
single <- get_gates_csv_single(single_data)
multi_data <- single_data
multi_data$other <- values + 1
multiple <- get_gates_csv(multi_data, session = NULL)

stopifnot(
  identical(seen_inputs[[1]], expected),
  identical(seen_inputs[[2]], expected),
  identical(seen_inputs[[3]], remove_outliers2(multi_data$other[is.finite(multi_data$other)])),
  identical(single$Gate, multiple$Gate[1])
)

cells <- data.frame(
  imageid = c("sample", "sample", "sample", "other"),
  marker = c(0, 1000, NA_real_, 1000)
)
gates <- data.frame(Patient = c("sample", "other"), Marker = "marker", Gate = c(50, 2000))
labeled <- label_cells_with_gates(cells, gates)
stopifnot(
  nrow(labeled) == nrow(cells),
  identical(labeled$marker, cells$marker),
  identical(labeled$marker_positivity, c("-", "+", NA_character_, "-"))
)

cat("Gating convention checks passed\n")
