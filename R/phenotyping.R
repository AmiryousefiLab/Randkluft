# Evaluate phenotype rules once and use the resulting partition everywhere.

normalize_phenotype_definitions <- function(definitions, marker_names = NULL) {
  if (!is.data.frame(definitions) || !"phenotype" %in% names(definitions) ||
      nrow(definitions) == 0L) {
    stop("Add at least one phenotype definition.")
  }

  status_columns <- grep("^status__", names(definitions), value = TRUE)
  if ("markers" %in% names(definitions) && length(status_columns)) {
    stop("Use either a markers column or status__ marker columns in a workflow, not both.")
  }
  if (!"markers" %in% names(definitions)) {
    if (!length(status_columns)) stop("Workflow needs a markers column or status__ marker columns.")
    status_markers <- substring(status_columns, nchar("status__") + 1L)
    if (any(status_markers == "") || anyDuplicated(status_markers)) {
      stop("Workflow has invalid or repeated status__ marker columns.")
    }
    status_values <- as.data.frame(lapply(definitions[status_columns], function(x) {
      value <- trimws(as.character(x))
      value[is.na(value) | value == ""] <- NA_character_
      if (any(!is.na(value) & !value %in% c("+", "-"))) {
        stop("Workflow status columns must contain +, -, or NA.")
      }
      value
    }), check.names = FALSE)
    definitions$markers <- vapply(seq_len(nrow(definitions)), function(i) {
      states <- unlist(status_values[i, ], use.names = FALSE)
      paste0(status_markers[!is.na(states)], states[!is.na(states)], collapse = ", ")
    }, character(1))
  }

  rules <- definitions[, intersect(c("phenotype", "markers", "match_mode"), names(definitions)), drop = FALSE]
  rules$phenotype <- gsub("[[:space:]]+", " ", trimws(as.character(rules$phenotype)))
  rules$markers <- trimws(as.character(rules$markers))
  if (!"match_mode" %in% names(rules)) rules$match_mode <- "all"
  rules$match_mode <- trimws(as.character(rules$match_mode))
  rules$match_mode[is.na(rules$match_mode) | rules$match_mode == ""] <- "all"

  if (anyNA(rules$phenotype) || any(rules$phenotype == "") ||
      anyDuplicated(rules$phenotype) ||
      any(rules$phenotype %in% c("Other", "Unresolved"))) {
    stop("Phenotype names must be non-empty, unique, and different from Other and Unresolved.")
  }
  if (anyNA(rules$markers) || any(rules$markers == "")) {
    stop("Each phenotype must have at least one selected marker.")
  }
  if (any(!rules$match_mode %in% c("all", "any_positive"))) {
    stop("match_mode must be 'all' or 'any_positive'.")
  }

  criteria <- lapply(seq_len(nrow(rules)), function(i) {
    tokens <- trimws(strsplit(rules$markers[i], ",", fixed = TRUE)[[1]])
    if (any(tokens == "")) stop("Empty marker criterion in ", rules$phenotype[i], ".")
    states <- substring(tokens, nchar(tokens))
    markers <- trimws(substr(tokens, 1L, nchar(tokens) - 1L))
    if (any(!states %in% c("+", "-")) || any(markers == "") || anyDuplicated(markers)) {
      stop("Invalid or repeated marker criterion in ", rules$phenotype[i], ". Use Marker+ or Marker-.")
    }
    if (rules$match_mode[i] == "any_positive" && !any(states == "+")) {
      stop("Any positive requires a positive marker in ", rules$phenotype[i], ".")
    }
    stats::setNames(states, markers)
  })

  if (is.null(marker_names)) marker_names <- unique(unlist(lapply(criteria, names), use.names = FALSE))
  marker_names <- unique(as.character(marker_names))
  unknown <- setdiff(unique(unlist(lapply(criteria, names), use.names = FALSE)), marker_names)
  if (length(unknown)) stop("Gate these markers before phenotyping: ", paste(unknown, collapse = ", "))

  states <- matrix(NA_character_, nrow = nrow(rules), ncol = length(marker_names),
                   dimnames = list(rules$phenotype, marker_names))
  for (i in seq_along(criteria)) states[i, names(criteria[[i]])] <- criteria[[i]]
  rules <- rules[, c("phenotype", "markers", "match_mode"), drop = FALSE]
  list(definitions = rules, states = states)
}

export_phenotype_workflow <- function(definitions, marker_names = NULL) {
  normalized <- normalize_phenotype_definitions(definitions, marker_names)
  output <- normalized$definitions[c("phenotype", "match_mode")]
  states <- as.data.frame(normalized$states, check.names = FALSE)
  names(states) <- paste0("status__", names(states))
  cbind(output, states)
}

evaluate_phenotype_rule <- function(data, states, match_mode) {
  n <- nrow(data)
  positive <- names(states)[!is.na(states) & states == "+"]
  negative <- names(states)[!is.na(states) & states == "-"]

  negative_match <- rep(TRUE, n)
  for (marker in negative) {
    negative_match <- negative_match & (data[[paste0(marker, "_positivity")]] == "-")
  }

  positive_match <- rep(if (match_mode == "all") TRUE else FALSE, n)
  for (marker in positive) {
    observed <- data[[paste0(marker, "_positivity")]] == "+"
    positive_match <- if (match_mode == "all") positive_match & observed else positive_match | observed
  }
  negative_match & positive_match
}

phenotype_partition <- function(data, definitions, gates,
                                diversity_estimator = function(labels) PEkit::MLEp(PEkit::abundance(labels))) {
  if (!is.data.frame(data) || nrow(data) == 0L || !"imageid" %in% names(data)) {
    stop("Upload cell data before phenotyping.")
  }
  if (is.null(gates) || !all(c("Patient", "Marker", "Gate") %in% names(gates))) {
    stop("Run Randkluft to generate gates before phenotyping.")
  }

  normalized <- normalize_phenotype_definitions(definitions, unique(as.character(gates$Marker)))
  required <- colnames(normalized$states)[colSums(!is.na(normalized$states)) > 0L]
  for (marker in required) {
    for (patient in unique(data$imageid)) {
      gate <- gates$Gate[gates$Marker == marker & gates$Patient == patient]
      if (length(gate) != 1L || !is.finite(gate)) {
        stop("A finite gate is required for ", marker, " in ", patient, ".")
      }
    }
    if (!is.numeric(data[[marker]])) stop("Marker intensities must be numeric: ", marker)
  }

  cells <- label_cells_with_gates(data, gates[gates$Marker %in% required, , drop = FALSE])
  for (marker in required) {
    status <- cells[[paste0(marker, "_positivity")]]
    if (any(!is.na(status) & !status %in% c("+", "-"))) {
      stop("Invalid positivity labels for ", marker, ".")
    }
  }

  membership <- matrix(NA, nrow = nrow(cells), ncol = nrow(normalized$definitions),
                       dimnames = list(NULL, normalized$definitions$phenotype))
  for (i in seq_len(ncol(membership))) {
    membership[, i] <- evaluate_phenotype_rule(cells, normalized$states[i, ],
                                               normalized$definitions$match_mode[i])
  }

  known_matches <- lapply(seq_len(nrow(cells)), function(i) which(membership[i, ] %in% TRUE))
  unresolved <- rowSums(is.na(membership)) > 0L
  category_id <- vapply(seq_len(nrow(cells)), function(i) {
    if (unresolved[i]) return("Unresolved")
    if (!length(known_matches[[i]])) return("Other")
    paste(known_matches[[i]], collapse = "+")
  }, character(1))

  ids <- unique(category_id)
  labels <- vapply(ids, function(id) {
    if (id %in% c("Other", "Unresolved")) return(id)
    names <- normalized$definitions$phenotype[as.integer(strsplit(id, "+", fixed = TRUE)[[1]])]
    if (length(names) == 1L) paste0(names, " only") else paste(names, collapse = " and ")
  }, character(1))
  if (anyDuplicated(labels)) {
    duplicates <- labels %in% labels[duplicated(labels)]
    labels[duplicates] <- paste0(labels[duplicates], " [rules ", ids[duplicates], "]")
  }

  for (column in c("phenotype", "matched_phenotypes")) {
    if (column %in% names(cells)) {
      preserved_name <- paste0("input_", column)
      while (preserved_name %in% names(cells)) preserved_name <- paste0("input_", preserved_name)
      cells[[preserved_name]] <- cells[[column]]
    }
  }
  cells$phenotype <- unname(labels[match(category_id, ids)])
  cells$matched_phenotypes <- vapply(known_matches, function(indices) {
    paste(normalized$definitions$phenotype[indices], collapse = "; ")
  }, character(1))

  counts <- tabulate(match(category_id, ids), nbins = length(ids))
  summary <- data.frame(phenotype = unname(labels), count = counts,
                        percentage = 100 * counts / nrow(cells),
                        stringsAsFactors = FALSE)
  summary <- summary[order(-summary$count, summary$phenotype), , drop = FALSE]
  rownames(summary) <- NULL

  included <- !category_id %in% c("Other", "Unresolved")
  index_labels <- cells$phenotype[included]
  n_index <- length(index_labels)
  k_index <- length(unique(index_labels))
  diversity <- if (n_index == 0L) {
    NA_real_
  } else if (k_index == 1L) {
    0
  } else if (k_index == n_index) {
    Inf
  } else {
    diversity_estimator(index_labels)
  }

  list(cells = cells, summary = summary, definitions = normalized$definitions,
       marker_states = normalized$states, diversity = diversity,
       n_index = n_index, k_index = k_index)
}

wrap_phenotype_label <- function(label, width) {
  words <- strsplit(label, " ", fixed = TRUE)[[1]]
  pieces <- unlist(lapply(words, function(word) {
    starts <- seq.int(1L, nchar(word), by = width)
    substring(word, starts, pmin(starts + width - 1L, nchar(word)))
  }), use.names = FALSE)
  paste(strwrap(paste(pieces, collapse = " "), width = width), collapse = "\n")
}

plot_phenotype_composition <- function(summary) {
  if (nrow(summary) == 0L) {
    plot.new()
    text(0.5, 0.5, "No cells to display")
    return(invisible(NULL))
  }

  device_width <- grDevices::dev.size("in")[1]
  label_width <- max(18L, min(40L, floor(device_width * 0.40 / strwidth("M", units = "in", cex = 0.8))))
  wrapped <- vapply(summary$phenotype, wrap_phenotype_label, character(1), width = label_width)
  left_margin <- max(strwidth(wrapped, units = "in", cex = 0.8)) + 0.25

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mai = c(1.05, left_margin, 0.55, 0.25))
  ordered <- rev(seq_len(nrow(summary)))
  graphics::barplot(summary$percentage[ordered], horiz = TRUE,
                    names.arg = wrapped[ordered], las = 1, cex.names = 0.8,
                    col = grDevices::hcl.colors(nrow(summary), "Dark 3")[ordered],
                    border = NA, xlim = c(0, 100), axes = FALSE,
                    main = "Phenotype Composition", xlab = "Percentage of all cells")
  graphics::axis(1, at = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"))
  invisible(NULL)
}
