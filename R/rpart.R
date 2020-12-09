#' Run rpart
#' this can be multi-line if we like
#'
#' @param expr_df data.frame containing expression in log2tpm with a unique column represented by param gene_id, rows are samples
#' @param gene_id string in colnames(expr_df) relating to log2tpm data
#'
#' @param clin_df data.frame containing clinical data with unique columns represented by surv_event and surv_time, rows are samples
#' @param surv_event colnames(clin_df) relating to survival event
#' @param surv_time colnames(clin_df) relating to survival event
#' @param join_el colname on which to join expression and survival data (default: rownames)
#' @return output the output object data structure and what it contains
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

  ##run rpart
  rattle::fancyRpartPlot(fit_tree)  # graph showing how patients are dichotomised
  decision_values <- fit_tree$splits
  cut_off <- decision_values[1,4]

  high_low <- ifelse(surv_expr_tb$gene >= cut_off, "High", "Low")
  clin_tb <- dplyr::mutate(.data = clin_df, "{gene_id}_group" := high_low,
                                            "{gene_id}_log2tpm" := as.numeric(os_expr$gene))

  return(clin_tb)
}

result_df <- run_rpart(expr_df, "RARB", clin_df, "OS", "OS.time", "sample")




##combine expr and surv
rpart::fancyRpartPlot(fit_tree)  # graph showing how patients are dichotomised
decision_values <- fit_tree$splits
cut_off <- decision_values[1,4]

# grouping patients into their high and low cohorts
high_low <- ifelse(CRABP2_values$`8330` >= cut_off,"High","Low")
clinical_data$CRABP2_Cohort <- as.vector(high_low)

# Survival analysis using the Kaplan-Meier method for CRABP2
surv_object <- Surv(time = clinical_data$OS.time, event = clinical_data$OS)
fit1 <- survfit(surv_object ~ CRABP2_Cohort, data = clinical_data)
survdiff(surv_object ~ CRABP2_Cohort, data = clinical_data)  #log rank test
ggsurvplot(fit1, data = clinical_data, pval = TRUE,
           legend = "bottom",
           xlab = "Time(Days)",
           ylab = "Survival Probability",
           legend.title = "CRABP2 expression",
           legend.labs = c("High Expression (n=60)", "Low Expression (n=468)"),
           pval.size = 5,
           font.legend = c(10, "plain", "black"))

clin_df <- readRDS("inst/extdata/clin_df.rds")
expr_df <- readRDS("inst/extdata/expr_df.rds")
