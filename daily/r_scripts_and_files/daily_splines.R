library(sleuth, splines)
library(biomaRt)
ensembl <- biomaRt::useEnsembl(biomart = "genes",
                               dataset = "hsapiens_gene_ensembl",
                               mirror = "useast")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "external_gene_name"),
                      mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id,
                     ext_gene = external_gene_name)

setwd("~/longsaliva")
s2c <- read.delim("s2c_daily.txt", sep=" ", header=TRUE)
time <- seq(from=1, to=length(s2c$sample), by=1)
s2c <- dplyr::mutate(s2c, time=time)
colnames(s2c) <- c("path", "sample", "time")

new_filter <- function(row, min_reads = 5, min_prop = 0.67) {
  mean(row >= min_reads) >= min_prop
}

so <- sleuth_prep(s2c, target_mapping = t2g, 
                  aggregation_column = "ens_gene",
                  extra_bootstrap_summary = TRUE,
                  filter_fun=new_filter)
X <- splines::ns(s2c$time, df=5)
full_design <- model.matrix(formula(~ X))
so <- sleuth_fit(so, formula = full_design, fit_name = "full")
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_de <- dplyr::filter(sleuth_table, qval <= 0.05)
#write.csv(sleuth_de, 'daily_new_filter.csv')
