#' Run survival analysis and plot using rpart output from run_rpart()
#'
#' @param clin_tb tibble created by run_rpart(), or a tibble containing a columns named '{gene_ids[1]}_{group_name}' and {surv_event, surv_time} as below
#' @param gene_ids vector of strings used in run_rpart to define gene used in classification; if named vector, uses names and ids for plotting
#' @param surv_event colnames(clin_tb) relating to survival event
#' @param surv_time colnames(clin_tb) relating to survival event
#' @param group_name character string of the 'group_name' column, therefore the gene_id concatenated with _{group_name} suffix; default: group
#' @param expr_unit unit of expression in clin_tb; default - log2tpm
#' @param col_palette colours to use in plotting (vector, high -> low expression; think palette in ggsurvplot is alphanum sorted...)
#' @param print_pdf print PDF to file (else return in output list)
#' @param print_png print PDF to file (else return in output list)
#' @param title_text title text for plot
#' @param sub_text sub text for plot
#' @param plot_prefix prefix for plot filename

#'
#' @return table from survival::survdiff (log rank test), ggsurvplot PDF printed
#'
#' @examples
#'
#' expr_df <- readRDS(system.file("extdata", "expr_df.rds", package="rpartSurvivalClassifier"))
#' clin_df <- readRDS(system.file("extdata", "clin_df.rds", package="rpartSurvivalClassifier"))
#' clin_new_tb <- rpartSurvivalClassifier::run_rpart(expr_df = expr_df, gene_id = "CRABP2", clin_df = clin_df, surv_event = "OS", surv_time = "OS.time", join_el = "sample")
#' log_rank_res <- rpartSurvivalClassifier::run_surv_plot(clin_tb = clin_new_tb, gene_id = "CRABP2", surv_event = "OS", surv_time = "OS.time")
#'
#' @export

run_surv_plot <- function(clin_tb, gene_ids, surv_event, surv_time, expr_unit = "log2tpm", group_name = "_group", col_palette = NULL, print_pdf = NULL, print_png = NULL, title_text = "", sub_text = "", plot_prefix = "rpart"){

  surv_object <- survival::Surv(time = unlist(clin_tb[,surv_time]),
                                event = unlist(clin_tb[,surv_event]))

  gene_lrts <- lapply(seq_along(gene_ids), function(x){

    gene_id <- gene_ids[x]
    ##include names if named vector
    if(!is.null(names(gene_ids))){
      gene_id <- paste(gene_id, names(gene_ids)[x], sep = "_")
    }

    if(paste0(gene_id, group_name) %in% colnames(clin_tb)){
      print(paste0("Data available for: ", gene_id))
      form1 <- paste0("surv_object ~ ", gene_id, group_name)
      form1<- as.formula(form1)
      fit1 <- survminer::surv_fit(form1,
                                data = clin_tb)
      lrt <- survival::survdiff(form1,
                                data = clin_tb)  #log rank test
      ntab <- table(clin_tb[paste0(gene_id, group_name)])

      ##colouring
      if(is.null(col_palette)){
        col_palette <- c("red", "dodgerblue")
      }

      ggs <- survminer::ggsurvplot(fit1, data = clin_tb,
                                   pval = TRUE,
                                   legend = "bottom",
                                   xlab = surv_time,
                                   ylab = paste0(surv_event, " Probability"),
                                   legend.title = paste0(gene_id, " ", expr_unit, ": "),
                                   legend.labs = c(paste0("High Expression (n = ", ntab["High"], ")"), paste0("Low Expression (n = ", ntab["Low"], ")")),
                                   pval.size = 5,
                                   font.legend = c(10, "plain", "black"),
                                   palette = col_palette,
                                   title = paste0(title_text, " - ", gene_id, "\n", sub_text))

      ##outputs
      if(!is.null(print_pdf)){
        pdf(paste0(print_pdf, "/", plot_prefix, ".ggsurvplot_", gene_id, "_", surv_event, ".pdf"), onefile = FALSE)
          print(ggs)
        dev.off()
      }
      if(!is.null(print_png)){
        png(paste0(print_png, "/", plot_prefix, ".ggsurvplot_", gene_id, "_", surv_event, ".png"), width = 800, height = 800)
          print(ggs)
        dev.off()
      }
      return(log_rank_test = lrt)
    } else {
      print(paste0("Data not available for: ", gene_id))
    }
  })

  if(!is.null(names(gene_ids))){
    gene_ids <- paste(gene_ids, names(gene_ids), sep = "_")
  }

  names(gene_lrts) <- gene_ids
  return(gene_lrts)
}
