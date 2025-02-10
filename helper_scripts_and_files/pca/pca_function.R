spread_abundance_by <- function(abund, var, which_order) {
  # var <- lazyeval::lazy(var)
  abund <- data.table::as.data.table(abund)
  var_spread <- data.table::dcast(abund, target_id ~ sample, value.var = var)
  # there is a discrepancy between data table's sorting of character vectors
  # and how tidyr previously (or the order function) sorts character vectors
  # so next step is needed to make sure the order is correct
  var_spread <- var_spread[order(var_spread$target_id), ]
  var_spread <- as.data.frame(var_spread, stringsAsFactors = FALSE)
  rownames(var_spread) <- var_spread$target_id
  var_spread["target_id"] <- NULL
  result <- as.matrix(var_spread)
  
  result[, which_order, drop = FALSE]
}

check_norm_status <- function(obj) {
  if (!is(obj$norm_fun_counts, 'function') || !is(obj$norm_fun_tpm, 'function')) {
    stop("This sleuth object was prepared without normalization. If you wish to do this step,",
         " repeat 'sleuth_prep' with 'normalize' set to TRUE (the default).")
  } else {
    return(TRUE)
  }
}

check_quant_mode <- function(obj, units) {
  stopifnot( is(obj, 'sleuth') )
  if (obj$gene_mode & units == 'est_counts') {
    warning(paste("your sleuth object is in gene mode,",
                  "but you selected 'est_counts'. Selecting 'scaled_reads_per_base'..."))
    units <- 'scaled_reads_per_base'
  } else if (!obj$gene_mode & units == 'scaled_reads_per_base') {
    warning(paste("your sleuth object is not in gene mode,",
                  "but you selected 'scaled_reads_per_base'. Selecting 'est_counts'..."))
    units <- 'scaled_reads_per_base'
  }
  
  units
}

as_df <- function(x, ...) {
  as.data.frame(x, stringsAsFactors = FALSE, ...)
}

plot_pca <- function(obj,
                     pc_x = 1L,
                     pc_y = 2L,
                     use_filtered = TRUE,
                     units = 'est_counts',
                     text_labels = FALSE,
                     color_by = NULL,
                     point_size = 3,
                     point_alpha = 0.8,
                     ...) {
  stopifnot( is(obj, 'sleuth') )
  stopifnot( check_norm_status(obj) )
  units <- check_quant_mode(obj, units)
  
  mat <- NULL
  if (use_filtered) {
    mat <- spread_abundance_by(obj$obs_norm_filt, units,
                               obj$sample_to_covariates$sample)
  } else {
    mat <- spread_abundance_by(obj$obs_norm, units,
                               obj$sample_to_covariates$sample)
  }
  
  pca_res <- prcomp(t(mat))
  pcs <- as_df(pca_res$x[, c(pc_x, pc_y)])
  pcs$sample <- rownames(pcs)
  rownames(pcs) <- NULL
  
  pc_x <- paste0('PC', pc_x)
  pc_y <- paste0('PC', pc_y)
  
  pcs <- dplyr::left_join(pcs, obj$sample_to_covariates,
                          by = 'sample')
  
  return(pcs)
}
