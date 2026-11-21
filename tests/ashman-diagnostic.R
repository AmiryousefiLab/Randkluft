# Run from the app directory: Rscript --vanilla tests/ashman-diagnostic.R
library(mclust)
source("R/diagnostics.R")

set.seed(23)
separated <- c(rnorm(500, 0, 0.5), rnorm(500, 4, 0.5))
random_state <- .Random.seed
result <- marker_separation_diagnostic(separated)
stopifnot(identical(.Random.seed, random_state))

two_fit <- Mclust(separated, G = 2)
expected_d <- abs(diff(as.numeric(two_fit$parameters$mean))) /
  sqrt(mean(rep(as.numeric(two_fit$parameters$variance$sigmasq), length.out = 2)))
stopifnot(
  isTRUE(all.equal(result$d, unname(expected_d))),
  result$d > 2,
  isTRUE(result$two_components_favored),
  grepl("2G favored", result$label, fixed = TRUE),
  isTRUE(all.equal(marker_separation_diagnostic(c(separated, NA, Inf))$d, result$d))
)

set.seed(23)
overlapping <- marker_separation_diagnostic(rnorm(1000))
stopifnot(
  overlapping$d < 2,
  isFALSE(overlapping$two_components_favored),
  grepl("1G favored", overlapping$label, fixed = TRUE)
)

set.seed(77)
unequal_spreads <- c(rnorm(500, 0, 0.2), rnorm(500, 4, 1.2))
unequal_fit <- Mclust(unequal_spreads, G = 2)
stopifnot(identical(unequal_fit$modelName, "V"))
stopifnot(isTRUE(all.equal(
  marker_separation_diagnostic(unequal_spreads)$d,
  unname(abs(diff(unequal_fit$parameters$mean)) /
           sqrt(mean(unequal_fit$parameters$variance$sigmasq)))
)))

tiny <- marker_separation_diagnostic(separated, min_component_share = 0.6)
stopifnot(grepl("tiny component", tiny$label, fixed = TRUE))
stopifnot(is.na(marker_separation_diagnostic(rep(1, 50))$d))
stopifnot(is.na(marker_separation_diagnostic(1:10)$d))
cat("Ashman diagnostic checks passed\n")
