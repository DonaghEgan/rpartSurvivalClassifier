#' Run survival analysis and plot using rpart output from run_rpart()
#'
#' @param clin_tb tibbe created by run_rpart()
#' @param gene_id string in colnames(expr_df) relating to log2tpm data
#'
#' @param surv_event colnames(clin_df) relating to survival event
#' @param surv_time colnames(clin_df) relating to survival event
#'
#' @return tibble of clin_df with two columns appended, 'gene_id'_group, 'gene_id'_log2tpm
#'
#' @examples
#'
#' expr_df <- readRDS(system.file("extdata", "expr_df.rds", package="rpartSurvivalClassifier"))
#' clin_df <- readRDS(system.file("extdata", "clin_df.rds", package="rpartSurvivalClassifier"))
#' clin_new_tb <- run_rpart(expr_df, "CRABP2", clin_df, "OS", "OS.time", "sample")
#' log_rank_res <- run_surv_plot(clin_new_tb, "CRABP2", "OS", "OS.time")
#'
#' @export

run_surv_plot <- function(clin_tb, gene_id, surv_event, surv_time){

  surv_object <- survival::Surv(time = unlist(clin_tb[,surv_time]),
                                event = unlist(clin_tb[,surv_event]))
  fit1 <- survival::survfit(as.formula(paste0("surv_object ~ ", gene_id, "_group")),
                            data = clin_tb)
  lrt <- survival::survdiff(as.formula(paste0("surv_object ~ ", gene_id, "_group")),
                            data = clin_tb)  #log rank test
  ntab <- table(clin_tb[paste0(gene_id, "_group")])
  ggs <- survminer::ggsurvplot(fit1, data = clin_tb,
                               pval = TRUE,
                               legend = "bottom",
                               xlab = surv_time,
                               ylab = paste0(surv_event, " Probability"),
                               legend.title = paste0(gene_id, " log2tpm: "),
                               legend.labs = c(paste0("High Expression (n = ", ntab["High"], ")"), paste0("Low Expression (n = ", ntab["Low"], ")")),
                               pval.size = 5,
                               font.legend = c(10, "plain", "black"))

  ##outputs
  pdf(paste0("ggsurvplot_", gene_id, "_", surv_event, ".pdf"), onefile = FALSE)
    print(ggs)
  dev.off()
  return(lrt)
}
