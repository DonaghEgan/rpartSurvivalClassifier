#' Run rpart on a named gene, return clinical data with appended high/low based on cutoff
#'
#' @param expr_df data.frame containing expression in log2tpm with a unique column represented by param gene_id, rows are samples
#' @param gene_ids vector of strings in colnames(expr_df) relating to log2tpm data for one or more gene
#'
#' @param clin_df data.frame containing clinical data with unique columns represented by surv_event and surv_time, rows are samples
#' @param surv_event colnames(clin_df) relating to survival event
#' @param surv_time colnames(clin_df) relating to survival event
#' @param join_el colname on which to join expression and survival data (default: rownames)
#'
#' @return tibble of clin_df with two columns appended, 'gene_id'_group, 'gene_id'_log2tpm; rpart PDF printed
#'
#' @examples
#'
#' expr_df <- readRDS(system.file("extdata", "expr_df.rds", package="rpartSurvivalClassifier"))
#' clin_df <- readRDS(system.file("extdata", "clin_df.rds", package="rpartSurvivalClassifier"))
#' clin_new_tb <- rpartSurvivalClassifier::run_rpart(expr_df = expr_df, gene_id = "CRABP2", clin_df = clin_df, surv_event = "OS", surv_time = "OS.time", join_el = "sample")
#'
#' @export

run_rpart <- function(expr_df, gene_ids, clin_df, surv_event, surv_time, join_el = NULL){

  ##parse relevant columns from inputs
  if(is.null(join_el)){
    expr_tb <- tibble::as_tibble(expr_df[,gene_ids], rownames = "sample")
    surv_tb <- tibble::as_tibble(clin_df[,surv_event], rownames = "sample")
  } else {
    expr_tb <- tibble::as_tibble(expr_df[,c(join_el, gene_ids)])
    colnames(expr_tb)[colnames(expr_tb) == join_el] <- "sample"
    surv_tb <- tibble::as_tibble(clin_df[,c(join_el, surv_event)])
    colnames(surv_tb)[colnames(surv_tb) == join_el] <- "sample"
  }

  #join and fit-tree
  surv_expr_tb <- dplyr::left_join(surv_tb, expr_tb, by = "sample")

  colnames(surv_expr_tb)[1:2] <- c("sample", surv_event)

  rpart_tb_list <- lapply(seq_along(gene_ids), function(x){
    gene_id <- gene_ids[x]
    print(paste0("Working on: ", gene_id))
    s_expr <- as.data.frame(surv_expr_tb[,c(surv_event, gene_id)])
    fit_tree <- rpart::rpart(s_expr, method = "anova")

    if(is.null(fit_tree$splits)){
      return(NULL)
    } else {
    # graph showing how patients are dichotomised
      pdf(paste0("rpart_", x, "_", surv_event, ".pdf"), onefile = FALSE)
        rattle::fancyRpartPlot(fit_tree)
      dev.off()

      decision_values <- fit_tree$splits
      cut_off <- decision_values[1,4]

      high_low <- ifelse(unlist(surv_expr_tb[gene_id]) >= cut_off, "High", "Low")
      gene_tb <- dplyr::mutate(.data = clin_df[],
                               "{gene_id}_group" := high_low,
                               "{gene_id}_log2tpm" := as.numeric(unlist(s_expr[x]))) %>%
                 dplyr::select(patient, barcode,
                               tidyselect::starts_with(!!as.vector(gene_id)))
     return(gene_tb)
    }
  })

  names(rpart_tb_list) <- gene_ids

  ##remove NULL
  rpart_tb_list <- rpart_tb_list[!sapply(rpart_tb_list, is.null)]

  if(!is.null(dim(rpart_tb_list))){
    rpart_tb <- rpart_tb_list %>% purrr::reduce(dplyr::left_join) %>%
                                  dplyr::left_join(., clin_df)

    return(rpart_tb)
  }
}
