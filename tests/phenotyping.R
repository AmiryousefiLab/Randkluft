# Run from the app directory: Rscript --vanilla tests/phenotyping.R
source("R/gating.R")
source("R/phenotyping.R")

definitions <- data.frame(
  phenotype = c("T cells", "CD4+ helper T cells", "B cells"),
  markers = c("T+", "T+, H+", "B+"),
  stringsAsFactors = FALSE
)
cells <- data.frame(
  imageid = rep("sample", 11),
  T = c(1, 1, 1, 1, rep(0, 6), 1),
  H = c(1, 1, rep(0, 8), NA_real_),
  B = c(rep(0, 4), 1, rep(0, 6)),
  phenotype = rep("original label", 11)
)
gates <- data.frame(
  Patient = rep("sample", 3), Marker = c("T", "H", "B"), Gate = rep(0.5, 3)
)

seen_labels <- NULL
estimate <- function(labels) {
  seen_labels <<- labels
  42
}
result <- phenotype_partition(cells, definitions, gates, estimate)
stopifnot(
  identical(result$cells$phenotype, c(
    rep("T cells and CD4+ helper T cells", 2),
    rep("T cells only", 2), "B cells only", rep("Other", 5), "Unresolved"
  )),
  identical(result$cells$matched_phenotypes[1], "T cells; CD4+ helper T cells"),
  identical(result$cells$matched_phenotypes[11], "T cells"),
  identical(result$cells$input_phenotype, cells$phenotype),
  identical(result$summary$count, c(5L, 2L, 2L, 1L, 1L)),
  isTRUE(all.equal(sum(result$summary$percentage), 100)),
  identical(result$n_index, 5L),
  identical(result$k_index, 3L),
  identical(result$diversity, 42),
  identical(seen_labels, result$cells$phenotype[1:5]),
  is.na(result$marker_states["T cells", "H"]),
  is.na(result$marker_states["B cells", "T"])
)

# CSV workflows with two columns remain valid; unchecked markers are unconstrained.
round_trip <- normalize_phenotype_definitions(result$definitions, c("T", "H", "B"))
stopifnot(identical(round_trip$states, result$marker_states))
workflow <- export_phenotype_workflow(result$definitions, c("T", "H", "B"))
stopifnot(
  is.na(workflow$status__H[1]),
  identical(workflow$status__T[1], "+"),
  is.na(workflow$status__T[3]),
  identical(normalize_phenotype_definitions(workflow, c("T", "H", "B"))$states,
            result$marker_states)
)
workflow_file <- tempfile(fileext = ".csv")
write.csv(workflow, workflow_file, row.names = FALSE, na = "NA")
reloaded <- read.csv(workflow_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(normalize_phenotype_definitions(reloaded)$definitions,
                    result$definitions))
unlink(workflow_file)

cells_file <- tempfile(fileext = ".csv")
write.csv(result$cells, cells_file, row.names = FALSE)
exported_cells <- read.csv(cells_file, stringsAsFactors = FALSE)
stopifnot(identical(exported_cells$phenotype, result$cells$phenotype))
unlink(cells_file)

any_positive <- normalize_phenotype_definitions(data.frame(
  phenotype = "Either signal", markers = "T+, H+, B-", match_mode = "any_positive"
), c("T", "H", "B"))
rule_values <- data.frame(T_positivity = c("-", "-", NA),
                          H_positivity = c("+", "-", "-"),
                          B_positivity = c("-", "-", "-"))
stopifnot(identical(
  evaluate_phenotype_rule(rule_values, any_positive$states[1, ], "any_positive"),
  c(TRUE, FALSE, NA)
))

long_summary <- result$summary
long_summary$phenotype[2] <- paste(rep("CD4+ helper T cells with a very long intersection", 3), collapse = " and ")
plot_file <- tempfile(fileext = ".png")
png(plot_file, width = 1000, height = 900, res = 110)
plot_phenotype_composition(long_summary)
dev.off()
stopifnot(file.info(plot_file)$size > 1000)
unlink(plot_file)

cat("Phenotyping checks passed\n")
