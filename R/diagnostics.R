# Descriptive diagnostics only; these fits never determine Randkluft's gate.
# D is the distance between two fitted means divided by their pooled SD.
# A positive BIC difference favors two components over one, but neither
# quantity establishes biological positivity or even two visible peaks.
marker_separation_diagnostic <- function(values, min_n = 20L,
                                         min_component_share = 0.01, seed = 1L) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]

  unavailable <- function(reason) {
    list(d = NA_real_, delta_bic = NA_real_, two_components_favored = NA,
         component_shares = numeric(0), label = paste0("D=n/a | ", reason))
  }
  spread <- stats::sd(values)
  if (length(values) < min_n || length(unique(values)) < 3L ||
      !is.finite(spread) || spread == 0) {
    return(unavailable("insufficient data"))
  }

  # mclust may sample observations when initializing large fits. Preserve the
  # caller's RNG state while making repeated diagnostics reproducible.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  fit <- function(components) {
    tryCatch(suppressWarnings(mclust::Mclust(values, G = components)),
             error = function(e) NULL)
  }
  one <- fit(1L)
  two <- fit(2L)
  if (is.null(one) || is.null(two)) return(unavailable("fit failed"))

  means <- as.numeric(two$parameters$mean)
  variances <- rep(as.numeric(two$parameters$variance$sigmasq), length.out = 2L)
  shares <- as.numeric(two$parameters$pro)
  bic_difference <- as.numeric(two$bic - one$bic)
  if (length(means) != 2L || length(shares) != 2L ||
      length(bic_difference) != 1L ||
      any(!is.finite(c(means, variances, shares, bic_difference))) ||
      any(variances <= 0) || any(shares <= 0)) {
    return(unavailable("fit failed"))
  }

  d <- abs(diff(means)) / sqrt(mean(variances))
  favored <- bic_difference > 0
  status <- if (min(shares) < min_component_share) {
    "tiny component"
  } else if (favored) {
    "2G favored"
  } else {
    "1G favored"
  }
  list(d = unname(d), delta_bic = bic_difference,
       two_components_favored = favored, component_shares = shares,
       label = sprintf("D=%.2f | %s", d, status))
}
