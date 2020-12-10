#' Run rpart on a named gene, return clinical data with appended high/low based on cutoff
#'
#' @param expr_df data.frame containing expression in log2tpm with a unique column represented by param gene_id, rows are samples
#' @param gene_id string in colnames(expr_df) relating to log2tpm data
#'
#' @param clin_df data.frame containing clinical data with unique columns represented by surv_event and surv_time, rows are samples
#' @param surv_event colnames(clin_df) relating to survival event
#' @param surv_time colnames(clin_df) relating to survival event
#' @param join_el colname on which to join expression and survival data (default: rownames)
#'
#' @return tibble of clin_df with two columns appended, 'gene_id'_group, 'gene_id'_log2tpm
#'
#' @examples
#'
#' expr_df <- readRDS(system.file("extdata", "expr_df.rds", package="rpartSurvivalClassifier"))
#' clin_df <- readRDS(system.file("extdata", "clin_df.rds", package="rpartSurvivalClassifier"))
#' clin_new_tb <- rpartSurvivalClassifier::run_rpart(expr_df = expr_df, gene_id = "CRABP2", clin_df = clin_df, surv_event = "OS", surv_time = "OS.time", join_el = "sample")
#'
#' @export

run_rpart <- function(expr_df, gene_id, clin_df, surv_event, surv_time, join_el = NULL){

  ##parse relevant columns from inputs
  if(is.null(join_el)){
    expr_tb <- tibble::as_tibble(expr_df[,gene_id], rownames = "sample")
    surv_tb <- tibble::as_tibble(clin_df[,surv_event], rownames = "sample")
  } else {
    expr_tb <- tibble::as_tibble(expr_df[,c(join_el, gene_id)])
    colnames(expr_tb)[colnames(expr_tb) == join_el] <- "sample"
    surv_tb <- tibble::as_tibble(clin_df[,c(join_el, surv_event)])
    colnames(surv_tb)[colnames(surv_tb) == join_el] <- "sample"
  }

  #join and fit-tree
  surv_expr_tb <- dplyr::left_join(surv_tb, expr_tb, by = "sample")
  colnames(surv_expr_tb) <- c("sample", "OS", "gene")
  os_expr <- as.data.frame(surv_expr_tb[,c("OS", "gene")])
  fit_tree <- rpart::rpart(os_expr, method = "anova")

  # graph showing how patients are dichotomised
  pdf(paste0("rpart_", gene_id, "_", surv_event, ".pdf"), onefile = FALSE)
    rattle::fancyRpartPlot(fit_tree)
  dev.off()

  decision_values <- fit_tree$splits
  cut_off <- decision_values[1,4]

  high_low <- ifelse(surv_expr_tb$gene >= cut_off, "High", "Low")
  clin_tb <- dplyr::mutate(.data = clin_df, "{gene_id}_group" := high_low,
                                            "{gene_id}_log2tpm" := as.numeric(os_expr$gene))

  return(clin_tb)
}
