# Run from the app directory: Rscript --vanilla tests/refactor-smoke.R
# Uses only packages already required by the app; no testthat dependency.
library(moments)
library(multimode)
library(mclust)

source("R/gating.R")
source("R/data-input.R")

fixture <- data.frame(
  imageID = c("one", "two", "three"),
  X = c(1, 2, 3),
  MarkerA = c(1, 20, 30)
)
fixture_path <- tempfile(fileext = ".csv")
write.csv(fixture, fixture_path, row.names = FALSE)

prepared <- prepare_cell_data(read_cell_data(fixture_path))
stopifnot(
  nrow(prepared) == 3L,
  identical(prepared$imageid, rep("sample", 3)),
  isTRUE(all.equal(prepared$X_centroid, fixture$X)),
  isTRUE(all.equal(prepared$Y_centroid, c(1, 1, 2))),
  isTRUE(all.equal(as.numeric(prepared$MarkerA), log(c(10, 20, 30))))
)

logged_fixture <- data.frame(imageid = c("one", "two", "three"), MarkerA = c(1, 11, 19))
logged_prepared <- prepare_cell_data(logged_fixture)
stopifnot(isTRUE(all.equal(as.numeric(logged_prepared$MarkerA), c(10, 11, 19))))

set.seed(7)
positive_tail <- c(rnorm(400), rnorm(40, 4, 0.2))
negative_tail <- -positive_tail
positive_gate <- skew_gate(positive_tail)
negative_gate <- skew_gate(negative_tail)
stopifnot(
  moments::skewness(positive_tail) > 0,
  moments::skewness(negative_tail) < 0,
  is.finite(positive_gate$cutoff),
  isTRUE(all.equal(
    negative_gate$cutoff,
    mean(mclust::Mclust(negative_tail, G = 2)$parameters$mean)
  )),
  identical(determine_positivity(2, 1), "+"),
  identical(determine_positivity(1, 2), "-")
)

# Stub only the progress UI; the gate calculation itself is unchanged.
Progress <- list(new = function(...) {
  list(set = function(...) invisible(NULL), close = function() invisible(NULL))
})
gate_data <- data.frame(imageid = rep("sample", length(positive_tail)), marker = positive_tail)
single <- get_gates_csv_single(gate_data)
gate_data$other <- positive_tail + 1
multiple <- get_gates_csv(gate_data, session = NULL)
stopifnot(
  nrow(single) == 1L,
  nrow(multiple) == 2L,
  isTRUE(all.equal(single$Gate, skew_gate(remove_outliers2(positive_tail))$cutoff)),
  isTRUE(all.equal(multiple$Gate[1], skew_gate(remove_outliers2(positive_tail))$cutoff))
)

with_nonfinite <- data.frame(
  imageid = rep("sample", length(positive_tail) + 2L),
  marker = c(positive_tail, NA_real_, Inf)
)
stopifnot(isTRUE(all.equal(
  get_gates_csv_single(with_nonfinite)$Gate,
  single$Gate
)))

cells <- data.frame(
  imageid = c("one", "one", "one", "two", "two", "other"),
  marker = c(0, 100, NA_real_, 0, Inf, 100)
)
gates <- data.frame(Patient = c("one", "two"), Marker = "marker", Gate = c(1, -1))
labeled <- label_cells_with_gates(cells, gates)
stopifnot(
  nrow(labeled) == nrow(cells),
  identical(labeled$marker, cells$marker),
  identical(labeled$marker_positivity, c("-", "+", NA_character_, "+", NA_character_, NA_character_))
)

unlink(fixture_path)
cat("Refactor smoke tests passed\n")
