#' Run rpart
#' this can be multi-line if we like
#'
#' @param rpart_df a data.frame containing data of format:
#'
#' @param second the second  param, the type of data structure, format
#' @return output the output object data structure and what it contains
#' @export

run_rpart <- function(rpart_df, gene_col){
  #CRABP2_values -> ?
  #clinical_data -> ?
}

# regression tree analysis
fit_tree <- rpart(OS ~ CRABP2_values, data= clinical_data, method="anova")
fancyRpartPlot(fit_tree)  # graph showing how patients are dichotomised 
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
