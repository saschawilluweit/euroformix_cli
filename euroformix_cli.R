#!/usr/bin/env Rscript

# nicely load required packages
for (package_name in c("docopt", "euroformix")) {
  if (!require(package_name, character.only = TRUE)) {
    print(paste0("Packege '", package_name, "' is needed. Please install using: R CMD INSTALL ", package_name))
    exit(1)
  }
  suppressPackageStartupMessages(library(package_name, quietly = TRUE, character.only = TRUE))
}

options(datatable.fread.datatable = FALSE)

###################### BEGIN HELPER FUNCTIONS ######################
print_nice <- function(s, file = stdout(), append = TRUE) {
  if (!arguments$quiet) {
    if (arguments$verbose) {
      write(s, stdout())
    }
    write(s, file, append = append)
  }
}
print_nice_sep_line <- function(prefix = "", ...) {
  print_nice(paste0(prefix, paste0(rep("*", 80), collapse = "")), ...)
}
print_nice_key_value <- function(key, value, prefix = "", sep = ": ", ...) {
  print_nice(paste0(prefix, key, sep, ifelse(is.null(value) || is.na(value), "", ifelse(typeof(value) == "double", signif(value, 3), value))), ...)
}
print_nice_list <- function(list_obj, list_names = NULL, ...) {
  if (length(list_obj) > 0) {
    if (is.null(list_names)) {
      list_names <- names(list_obj)
    }
    for (list_item_name in list_names) {
      print_nice_key_value(list_item_name, list_obj[[list_item_name]], ...)
    }
  }
}
print_nice_marker_sorted_list <- function(marker_sorted_list, markers = NULL, ...) {
  sample_names <- names(marker_sorted_list)
  for (sample_names_idx in seq_len(length(sample_names))) {
    marker_names <- names(marker_sorted_list[[sample_names[sample_names_idx]]])
    if(!is.null(markers)) {
      marker_names  <- markers
    }
    for (marker_names_idx in seq_len(length(marker_names))) {
      print_nice_key_value(paste0(sample_names[sample_names_idx], ": ", marker_names[marker_names_idx]), paste0(paste(marker_sorted_list[[sample_names[sample_names_idx]]][[ marker_names[marker_names_idx]]]$adata, collapse = ","), " [", paste(marker_sorted_list[[sample_names[sample_names_idx]]][[marker_names[marker_names_idx]]]$hdata, collapse = ","), "]"), ...)
    }
  }
  print_nice_sep_line(...)
}

# build a consistent (markers, allele and peak height) list out of a table matrix (got from tableReader())
table_to_marker_sorted_list <- function(src_table, allele_height_threshold = 0, markers = NULL) {
  column_names <- colnames(src_table)
  current_analysis_threshold <- ifelse(allele_height_threshold <= 0, 999999, allele_height_threshold)

  marker_column_index <- grep("marker", tolower(column_names), fixed = TRUE)[1]
  if (is.null(marker_column_index) || is.na(marker_column_index)) {
    stop("No column named like \"*marker*\" found at your table.")
  }

  marker_names <- unique(toupper(src_table[, marker_column_index]))
  if(!is.null(markers)) {
    marker_names <- markers
  }

  sample_name_column_index <- grep("sample", tolower(column_names), fixed = TRUE)[1]
  if (is.null(sample_name_column_index) || is.na(sample_name_column_index)) {
    sample_name_column_index <- grep("name", tolower(column_names), fixed = TRUE)[1]
    if (is.null(sample_name_column_index) || is.na(sample_name_column_index)) {
      stop("No column named like \"*sample*\" or \"*name*\" found at your table.")
    }
  }
  sample_names <- unique(as.character(src_table[, sample_name_column_index]))

  allele_column_indices <- grep("allele", tolower(column_names), fixed = TRUE)
  if (length(allele_column_indices) <= 0) {
    stop("No column named like \"*allele*\"e.g. \"Allele1\", \"Allele2\", ... or \"AlleleN\" found at your table.")
  }

  height_column_indices <- grep("height", tolower(column_names), fixed = TRUE)

  dest_list <- list()
  for (sample_names_idx in seq_len(length(sample_names))) {
    dest_list[[sample_names[sample_names_idx]]] <- list()
    for (marker_names_idx in seq_len(length(marker_names))) {
      logic_rows_to_use <- src_table[, sample_name_column_index] == sample_names[sample_names_idx] & toupper(src_table[, marker_column_index]) == marker_names[marker_names_idx]
      if (any(logic_rows_to_use)) {
        keep_row_indices <- which(!is.na(src_table[logic_rows_to_use, allele_column_indices]) & src_table[logic_rows_to_use, allele_column_indices] != "")
        if (length(height_column_indices) > 0) {
          peak_heights <- as.numeric(as.character(src_table[logic_rows_to_use, height_column_indices][keep_row_indices]))
          if (allele_height_threshold <= 0) {
            current_analysis_threshold <- min(c(current_analysis_threshold, peak_heights[which(peak_heights > 0)]))
          }
          keep_row_indices <- which(peak_heights >= allele_height_threshold)
          dest_list[[sample_names[sample_names_idx]]][[marker_names[marker_names_idx]]]$hdata <- peak_heights[keep_row_indices]
        }
        dest_list[[sample_names[sample_names_idx]]][[marker_names[marker_names_idx]]]$adata <- as.character(src_table[logic_rows_to_use, allele_column_indices][keep_row_indices])
      }
    }
  }
  names(dest_list) <- sample_names
  current_analysis_threshold <<- ifelse(allele_height_threshold <= 0, current_analysis_threshold, allele_height_threshold)
  return(dest_list)
}

# build a consistent list adn assure all alleles are known => use Qassignate function which converts non-observed alleles to "99" and adds missing allele frequencies, if any
euroformix_data_obj <- function(samples, references, allele_frequencies) {
  marker_names <- names(allele_frequencies)
  normalized_samples <- lapply(samples, function(x) {
    return(x[marker_names])
  })
  normalized_references <- list()
  for (marker_name in marker_names) {
    normalized_references[[marker_name]] <- lapply(references, function(x) {
      return(x[[marker_name]]$adata)
    })
  }
  # minF => the allele frequency of unobserved alleles
  # NB: NOTICE THE CHANGE HERE OF incS = FALSE even for stutter model (this has been updated in v2)
  qassignated <- euroformix::Qassignate(samples = normalized_samples, allele_frequencies, normalized_references, incS = FALSE, incR = FALSE)
  return(list(samples = normalized_samples, refData = qassignated$refData, popFreq = qassignated$popFreq))
}

estimate_dropin_lambda = function(src_table) {
  column_names <- colnames(src_table)
  height_column_indices <- grep("height", tolower(column_names), fixed = TRUE)

  height_values <- as.numeric(as.vector(unlist(src_table[, height_column_indices]))) # convert matrix to a numerical vector
  height_values <- height_values[!is.na(height_values)] # don't consider NAs
  height_values <- height_values[height_values >= current_analysis_threshold] # UPDATED in euroformix v2.1.1: consider dropinPH to be above threhsold

  lambda <- signif(length(height_values) / sum(height_values - current_analysis_threshold), 2) # estimated hyperparam with 2 significant numbers
 
  if (arguments$verbose) {
    print(paste("Height values:", paste0(x, collapse=",")))
    print(paste0("Estimated drop-in hyperparameter lambda:", lambda))
  }

  return(lambda)
}

cont_lik_mle_helper <- function(...) {
  if (n_cores > 1) {
    return(euroformix::contLikMLEpara(...))
  } else {
    return(euroformix::contLikMLE(...))
  }
}

cont_lik_mcmc_helper <- function(...) {
  if (n_cores > 1) {
    return(euroformix::contLikMCMCpara(...))
  } else {
    return(euroformix::contLikMCMC(...))
  }
}
###################### END HELPER FUNCTIONS ######################


###################### BEGIN ARGUMENTS PARSING ######################
"EuroForMix CLI.

Usage:
  euroformix_cli.R [options]
  euroformix_cli.R (-h | --help)
  euroformix_cli.R (-V | --version)

Options:
  -h --help                       Show this screen.
  -V --version                    Show version.
  -v --verbose                    Verbose output.
  --db_file <file>                Allele frequencies table.
  --samples_path <dir>            Path of sample files (one CSV file per sample with or without replicates).
  --sample_file <file>            Single sample file with or without replicates.
  --queried_file <file>           Genotypes of reference profiles in question (for use under Hp).
  --known_file <file>             Genotypes of known reference profiles (for use under Hp and Hd).
  --output_path <dir>             Path to place output files [default: results].
  --kit_database_file <kit_file>  Kit database file with all known kits (if not provided an internal one is used).
  --kit <name>                    Kit used for all given samples/replicates (in auto mode a separate column need to indicate the appropriate kit) [default: auto].
  --unknowns <n>                  Nubmber of unknows under prosecution hypotheses [default: 1].
  --fst <float>                   F_ST subpopulation/coancestry correction [default: 0.03].
  --dropin                        Consider Drop-In.
  --dropin_probability <float>    Probability of Drop-In [default: 0.05].
  --dropin_lambda <float>         Drop-In hyperparameter [default: auto].
  --stutter                       Consider stutter.
  --stutter_function <code>       Density function of stutter proportion parameter x [default: 1].
  --degradation                   Consider degradation.
  --analysis_threshold <n>        Analysis threshold in RFU (auto for autodetection = minimal observed peak height) [default: auto].
  --seed <n>                      Seed to reproduce randomness.
  --mle_steps <n>                 Maximum number of random evaluations nlm-optimizing routing [default: 4].
  --mle_iterations <n>            Maximum number of iterations for the MLE optimization [default: 100].
  --mcmc                          Calculate conservative LR estimate based on MCMC.
  --mcmc_iterations <n>           Number of MCMC iterations [default: 1000].
  --mcmc_delta <n>                A numerical parameter to scale with the covariance function Sigma [default: 10].
  --validation                    Do validation of assumption under Hp/Hd.
  --validation_alpha <float>      Rejection significance of peak-height-model under Hp/Hd [default: 0.01].
  --deconvolution                 Perform deconvolution under Hp/Hd.
  --parallel <n>                  Perfom calculatations in parallel at n cores (auto for autodetection) [default: auto].
  -Q --quiet                      Don't write any file or output anything except MLE LR (disables verbosity, MCMC, validation and deconvolution).
" -> doc

arguments <- docopt(doc, version = "EuroForMix CLI 1.0")

if ((is.null(arguments$samples_path) || is.na(arguments$samples_path) || !file.exists(arguments$samples_path)) && (is.null(arguments$sample_file) || is.na(arguments$sample_file) || !file.exists(arguments$sample_file))) {
  stop("No sample path (--sample_path <dir>) or sample file (--sample_file <file>) were given or the given one does not exists.")
} else {
  if (!is.null(arguments$samples_path) && !is.na(arguments$samples_path) && file.exists(arguments$samples_path)){
    sample_files <- list.files(arguments$samples_path, full.names = TRUE, pattern = ".(csv|txt)$")
    if (length(sample_files) <= 0) {
      stop(paste0("Sample path \"", arguments$samples_path, "\" does not contain any .csv or .txt file."))
    }
  } else {
    sample_files <- c(arguments$sample_file)
  }
}

if (is.null(arguments$queried_file) || is.na(arguments$queried_file) || !file.exists(arguments$queried_file)) {
  stop("No queried refererences file (--queried_file <file>) was given or given file does not exists.")
}

if (is.null(arguments$db_file) || is.na(arguments$db_file) || !file.exists(arguments$db_file)) {
  stop("No allele frequencies file (--db_file <file>) was given or given file does not exists.")
}
if (!arguments$quiet) {
  if (is.null(arguments$output_path) || is.na(arguments$output_path)) {
    stop("No output path (--output_path <dir>) was given.")
  } else {
    dir.create(arguments$output_path, showWarnings = FALSE, recursive = TRUE)
  }
}
if (is.null(arguments$kit_database_file) || is.na(arguments$kit_database_file)) {
  arguments$kit_database_file = system.file("extdata", "kit_database.txt", package = "euroformix", mustWork = TRUE)
} else if(!file.exists(arguments$kit_database_file)) {
  stop(paste0("kit_database '", , "' not found."))
}

if (arguments$quiet) {
  new_stdout <- vector("character")
  new_stdout_con <- textConnection("new_stdout", "wr")
  sink(new_stdout_con)

  options(warn = -1)

  arguments$verbose <- FALSE
  arguments$mcmc <- FALSE
  arguments$validation <- FALSE
  arguments$deconvolution <- FALSE
}

# seed for reproducting MLE and MCMC results
seed <- runif(1, min = 1, max = 1e6)
if (!is.null(arguments$seed) && !is.na(arguments$seed) && as.numeric(arguments$seed) > 0) {
  seed <- as.numeric(arguments$seed)
  print("Unfortunately, fixing the randomness via seed is not implemented yet.")
}
seed <- floor(seed)
arguments$seed <- seed

n_cores <- 1
if (!is.null(arguments$parallel) && !is.na(arguments$parallel)) {
  if (arguments$parallel == "auto") {
    n_cores <- parallel::detectCores()
  } else {
    n_cores <- max(c(1, as.numeric(arguments$parallel)))
  }
  if (n_cores > 1) {
    new_detect_cores_func <- function() {
      return(n_cores)
    }
    unlockBinding("detectCores", getNamespace("parallel"))
    assign("detectCores", new_detect_cores_func)
    assign("detectCores", new_detect_cores_func, getNamespace("parallel"))
  }
}

current_analysis_threshold <- 0
if (arguments$analysis_threshold == "auto" || (typeof(arguments$analysis_threshold) == "double" && as.numeric(arguments$analysis_threshold) < 1)) {
  arguments$analysis_threshold <- 0
} else {
  current_analysis_threshold <- as.numeric(arguments$analysis_threshold)
}

arguments_file <- stdout()
if (arguments$quiet) {
  arguments_file <- new_stdout_con
} else {
  arguments_file <- file(file.path(arguments$output_path, "arguments.txt"), "w")
}
print_nice_key_value("Start Time", format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), file = arguments_file, append = FALSE)
print_nice_sep_line(file = arguments_file)

print_nice_list(arguments, list_names = names(arguments[!startsWith(names(arguments), "--")]), file = arguments_file)
print_nice_key_value("Number of cores", n_cores, file = arguments_file)

print_nice_key_value("Euroformix version", paste(packageVersion("euroformix"), sep = "."), file = arguments_file)

print_nice_sep_line(file = arguments_file)
###################### END ARGUMENTS PARSING ######################


###################### BEGIN EUROFORMIX PREPARATIONS ######################
# kit database
kit_database = euroformix::tableReader(arguments$kit_database_file)
print_nice_key_value("Kit database", arguments$kit_database_file, file = arguments_file)
print_nice_key_value("Kit database avaiable kits", paste(unique(kit_database$Short.Name), collapse = ", "), file = arguments_file)

# allele frequencies
allele_frequencies <- euroformix::freqImport(arguments$db)[[1]]
print_nice_key_value("Loci to consider", paste(names(allele_frequencies), collapse = ", "), file = arguments_file)

# queried references (and set markers in use globally)
queried_references <- table_to_marker_sorted_list(euroformix::tableReader(arguments$queried))
golbal_marker_names <- names(queried_references[[names(queried_references)[1]]])

print_nice_marker_sorted_list(queried_references, prefix = "Queried ", file = arguments_file, markers = golbal_marker_names)

# known references
known_references <- list()
if (!is.null(arguments$known)) {
  known_references <- table_to_marker_sorted_list(euroformix::tableReader(arguments$known), markers = golbal_marker_names)
  print_nice_marker_sorted_list(known_references, prefix = "Known ", file = arguments_file, markers = golbal_marker_names)
}
# combine both list to a global references list
references <- append(queried_references, known_references)


# model setup
fst <- as.numeric(arguments$fst)
number_of_additional_unknowns <- as.numeric(arguments$unknowns)
number_of_contributors <- 1 + length(known_references) + number_of_additional_unknowns

xi <- 0 # if xi == 0, there is no stutter model applied
pXi <- function(x){1}
if (arguments$stutter) {
  xi <- NULL
  pXi <- eval(parse(text = paste0("function(x){", arguments$stutter_function, "}")))
}

#pXi <- function(x) {
#  dbeta(x, 1, 1)
#}

dropin_probability <- dropin_lambda <- 0
if (arguments$dropin) {
  dropin_probability <- as.numeric(arguments$dropin_probability)
}
###################### END EUROFORMIX PREPARATIONS ######################


###################### BEGIN MAIN CALCULATION(S) ######################
for (sample_file in sample_files) {
  sample_file_table_data <- euroformix::tableReader(sample_file)

  samples <- table_to_marker_sorted_list(sample_file_table_data, allele_height_threshold = as.numeric(arguments$analysis_threshold), markers = golbal_marker_names)
  sample_name <- sub(".csv", "", sub(".txt", "", basename(sample_file)))

  kit <- NULL # if kit == NULL, degradation is not estimated
  if (arguments$degradation) {
    kit_name <- NULL
    if (is.null(arguments$kit) || is.na(arguments$kit) || arguments$kit == "auto") {
      kit_column_index <- grep("kit", tolower(colnames(sample_file_table_data)), fixed = TRUE)[1]
      if (is.null(kit_column_index) || is.na(sample_file_table_data)) {
        kit_name <- unique(sample_file_table_data[sample_file_table_data, ])[1]
      } else {
        stop(paste0("Kit 'auto' was selected, but now kit column could be found at sample file '", sample_file, "'."))
      }
    } else {
      kit_name <- arguments$kit
    }
    kit <- euroformix::getKit(arguments$kit_database_file, kit_name)
    print_nice_key_value(paste0("Sample ", sample_name, " Kit"), unique(kit$Short.Name)[1])
  }

  print_nice_marker_sorted_list(samples, prefix = paste0("Sample ", sample_name, ": Replicates: "), file = arguments_file, markers = golbal_marker_names)

  data <- euroformix_data_obj(samples, references, allele_frequencies)

  if (!is.na(arguments$dropin_lambda)) {
    if (arguments$dropin_lambda == "auto") {
      dropin_lambda <- estimate_dropin_lambda(sample_file_table_data);
    } else {
      dropin_lambda <- as.numeric(arguments$dropin_lambda)
    }
  }

  for (queried_reference_name in names(queried_references)) {
    out_file <- stdout()
    if (arguments$quiet) {
      out_file <- new_stdout_con
    } else {
      out_file <- file(file.path(arguments$output_path, paste0("results_", sample_name, "_", queried_reference_name, ".txt")), "w")
    }

    print_nice_list(list("Sample" = sample_name, "Threshold" = current_analysis_threshold, "Queried" = queried_reference_name, "Known" = paste(names(known_references), collapse = ", "), "Unknown" = number_of_additional_unknowns, "Number of contributors" = number_of_contributors, "Drop-In lambda" = dropin_lambda), file = out_file)

    ########### Hp ############
    # generate condOrder
    # @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
    hp_cond_order <- c()
    hp_text <- c()
    hp_cond_order_index <- 1
    for (reference_name in names(references)) {
      if (reference_name == queried_reference_name || reference_name %in% names(known_references)) {
        hp_cond_order <- c(hp_cond_order, hp_cond_order_index)
        hp_cond_order_index <- hp_cond_order_index + 1
        hp_text <- c(hp_text, reference_name)
      } else {
        hp_cond_order <- c(hp_cond_order, 0)
      }
    }
    print_nice_key_value("Hp", paste(c(hp_text, rep("UK", number_of_additional_unknowns)), collapse = " + "), file = out_file)

    set.seed(seed)
    # refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=10,kit=NULL,verbose=TRUE,maxIter=100,knownRel=NULL,ibd=c(1,0,0)
    hp_fit <- cont_lik_mle_helper(number_of_contributors, data$samples, popFreq = data$popFreq, threshT = current_analysis_threshold, nDone = as.numeric(arguments$mle_steps), xi = xi, pXi = pXi, refData = data$refData, prC = dropin_probability, lambda = dropin_lambda, fst = fst, condOrder = hp_cond_order, kit = kit, verbose = arguments$verbose, maxIter = as.numeric(arguments$mle_iterations))
    print_nice_list(hp_fit$fit$thetahat2, prefix = "Hp ", file = out_file)

    if (arguments$deconvolution && number_of_additional_unknowns > 0) {
      hp_deconvolution <- euroformix::deconvolve(hp_fit)
      write.table(hp_deconvolution$table2, file = file.path(arguments$output_path, paste0("hp_deconvolution_", sample_name, "_", queried_reference_name, ".txt")), , sep = ";")
    }

    if (arguments$validation) {
      hp_validation_mle <- euroformix::validMLEmodel(hp_fit, kit = arguments$kit, plottitle = "Hp", alpha = as.numeric(arguments$validation_alpha))
      print_nice_list(as.list(hp_validation_mle), prefix = "Hp ", file = out_file)
    }

    if (arguments$mcmc) {
      set.seed(seed)
      hp_fit_mcmc <- cont_lik_mcmc_helper(hp_fit, niter = as.numeric(arguments$mcmc_iterations), delta = as.numeric(arguments$mcmc_delta))
    }

    ########### Hd ############
    # generate condOrder
    # @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
    hd_cond_order <- c(rep(0, length(queried_references)), seq_len(length(known_references)))
    # generate knownRef
    # @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
    hd_known_ref <- which(names(queried_references) == queried_reference_name)
    print_nice_key_value("Hd", paste(c(paste0("not(", queried_reference_name, ")"), names(known_references), rep("UK", number_of_additional_unknowns)), collapse = " + "), file = out_file)

    set.seed(seed)
    hd_fit <- cont_lik_mle_helper(number_of_contributors, data$samples, popFreq = data$popFreq, threshT = current_analysis_threshold, nDone = as.numeric(arguments$mle_steps), xi = xi, pXi = pXi, refData = data$refData, prC = dropin_probability, lambda = dropin_lambda, fst = fst, condOrder = hd_cond_order, knownRef = hd_known_ref, kit = kit, verbose = arguments$verbose, maxIter = as.numeric(arguments$mle_iterations))
    print_nice_list(hd_fit$fit$thetahat2, prefix = "Hd ", file = out_file)

    if (arguments$deconvolution) {
      hd_deconvolution <- euroformix::deconvolve(hd_fit)
      write.table(hd_deconvolution$table2, file = file.path(arguments$output_path, paste0("hd_deconvolution_", sample_name, "_", queried_reference_name, ".txt")), , sep = ";")
    }

    if (arguments$validation) {
      hd_validation_mle <- euroformix::validMLEmodel(hd_fit, kit = arguments$kit, plottitle = "Hd", alpha = as.numeric(arguments$validation_alpha))
      print_nice_list(as.list(hd_validation_mle), prefix = "Hd ", file = out_file)
    }

    if (arguments$mcmc) {
      set.seed(seed)
      hd_fit_mcmc <- cont_lik_mcmc_helper(hd_fit, niter = as.numeric(arguments$mcmc_iterations), delta = as.numeric(arguments$mcmc_delta))
    }


    ########### LR ############
    print_nice_key_value("LR", exp(hp_fit$fit$loglik - hd_fit$fit$loglik), file = out_file)
    print_nice_list(as.list(exp(euroformix::logLiki(hp_fit) - euroformix::logLiki(hd_fit))), prefix = "LR ", file = out_file)

    # Calculate conservatibe LR based on MCMC
    if (arguments$mcmc) {
      print_nice_key_value("conservative LR (MCMC)", (hp_fit_mcmc$margL / hd_fit_mcmc$margL), file = out_file)
    }

    if (arguments$quiet) {
      sink()
#      cat(paste0(sample_name, "\t", queried_reference_name, "\t", exp(hp_fit$fit$loglik - hd_fit$fit$loglik), "\n"), file = stdout(), append = TRUE)
      cat(paste0(exp(hp_fit$fit$loglik - hd_fit$fit$loglik), "\n"), file = stdout(), append = TRUE)
      sink(new_stdout_con)
    } else {
      close(out_file)
    }
  }
}

if (arguments$validation) {
  for (rplot_file in Sys.glob(file.path(arguments$output_path, "Rplots[0-9]*.pdf"))) {
    unlink(rplot_file)
  }
}

print_nice_key_value("End Time", format(Sys.time(), format = "%Y-%m-%d %H:%M"), file = arguments_file)
if (!arguments$quiet) {
  close(arguments_file)
}
