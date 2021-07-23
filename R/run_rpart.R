#' Run rpart on a named gene, return clinical data with appended high/low based on cutoff
#'
#' @param expr_df data.frame containing expression in log2tpm with a unique column represented by param gene_id, rows are samples
#' @param gene_ids vector of strings in colnames(expr_df) relating to log2tpm data for one or more gene
#'
#' @param clin_df data.frame containing clinical data with unique columns represented by surv_event and surv_time, rows are samples
#' @param surv_event colnames(clin_df) relating to survival event
#' @param surv_time colnames(clin_df) relating to survival event
#' @param join_el colname on which to join expression and survival data (default: rownames)
#' @param print_pdf print PDF to file (else return in output list)
#' @param title_text title text for plot
#'
#' @return list containing: [[1]] tibble of clin_df with columns appended: 'gene_id'_group, 'gene_id'_log2tpm; [[2]] named list of rpart plots by gene_ids
#'
#' @examples
#'
#' expr_df <- readRDS(system.file("extdata", "expr_df.rds", package="rpartSurvivalClassifier"))
#' clin_df <- readRDS(system.file("extdata", "clin_df.rds", package="rpartSurvivalClassifier"))
#' clin_new_tb <- rpartSurvivalClassifier::run_rpart(expr_df = expr_df, gene_id = "CRABP2", clin_df = clin_df, surv_event = "OS", surv_time = "OS.time", join_el = "sample")
#'
#' @export

run_rpart <- function(expr_df, gene_ids, clin_df, surv_event, surv_time, join_el = NULL, print_pdf = NULL, title_text = NULL){

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

  if(is.null(title_text)){
    title_text = ""
  } else {
    title_text <- paste0(title_text, " - ")
  }

  plot_catch <- function(fit_tree, gene_id, surv_event, title_text){
    rattle::fancyRpartPlot(fit_tree, main = paste0(title_text, gene_id, " - ", surv_event))
  }

  rpart_tb_list <- lapply(seq_along(gene_ids), function(x){
    gene_id <- gene_ids[x]
    print(paste0("Working on: ", gene_id))
    s_expr <- as.data.frame(surv_expr_tb[,c(surv_event, gene_id)])
    fit_tree <- rpart::rpart(s_expr, method = "anova")

    if(is.null(fit_tree$splits)){
      return(NULL)
    } else {
    # graph showing how patients are dichotomised
      if(!is.null(print_pdf)){
        pdf(paste0("rpart_", gene_id, "_", surv_event, ".pdf"), onefile = FALSE)
          rattle::fancyRpartPlot(fit_tree, main = paste0(title_text,  gene_id, " - ", surv_event))
        dev.off()
      }

      fit_tree_plot <- function(){plot_catch(fit_tree, gene_id, surv_event, title_text)}

      decision_values <- fit_tree$splits
      cut_off <- decision_values[1,4]

      high_low <- ifelse(unlist(surv_expr_tb[gene_id]) >= cut_off, "High", "Low")
      gene_tb <- dplyr::mutate(.data = clin_df[],
                               "{gene_id}_rpart_group" := high_low,
                               "{gene_id}_log2tpm" := as.numeric(unlist(s_expr[gene_id]))) %>%
                 dplyr::select(patient, barcode,
                               tidyselect::starts_with(!!as.vector(gene_id)))
     return(list(fit_tree_plot, gene_tb))
    }
  })

  names(rpart_tb_list) <- gene_ids

  ##split plots and gene_tbs
  rpart_plot_list <- lapply(rpart_tb_list, function(f){
      return(f[[1]])
    })
  rpart_gene_list <- lapply(rpart_tb_list, function(f){
      return(f[[2]])
    })

  ##remove NULL
  rpart_plot_list <- rpart_plot_list[!sapply(rpart_plot_list, is.null)]
  rpart_gene_list <- rpart_gene_list[!sapply(rpart_gene_list, is.null)]

  if(length(rpart_gene_list)>0){
    rpart_tb <- rpart_gene_list %>% purrr::reduce(dplyr::left_join) %>%
                                  dplyr::left_join(., clin_df)

    return(list(rpart_tb = rpart_tb, rpart_split_plot_list = rpart_plot_list))
  }
}
