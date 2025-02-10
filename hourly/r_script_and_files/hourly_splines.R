library(sleuth, splines)
library(biomaRt)
ensembl <- biomaRt::useEnsembl(biomart = "genes",
                               dataset = "hsapiens_gene_ensembl",
                               mirror='www')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "external_gene_name"),
                      mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id,
                     ext_gene = external_gene_name)

setwd("~/longsaliva")
s2c <- read.delim("s2c_hourly.txt", sep=" ", header=TRUE)
time <- rep(seq(from=1, to=length(s2c$sample)/2, by=1), times=2)
s2c <- dplyr::mutate(s2c, time=time)
sample <- paste0(rep(c('pre', 'post'), each=20),'_',s2c$time)
s2c$sample <- sample

colnames(s2c) <- c("path", "sample", "condition", "time")

group <- relevel(factor(s2c$condition), ref='pre_vaccination')

new_filter <- function(row, min_reads = 5, min_prop = 0.67) {
  mean(row >= min_reads) >= min_prop
}

so <- sleuth_prep(s2c, target_mapping = t2g,
                  aggregation_column = "ens_gene",
                  extra_bootstrap_summary = TRUE,
                  filter_fun = new_filter)

X <- splines::ns(s2c$time, df=5)
full_design <- model.matrix(formula(~0 + group + group:X))
colnames(full_design)

so <- sleuth_fit(so, full_design, "full")
so <- sleuth_fit(so, ~0+group, "reduced")
so <- sleuth_lrt(so, "reduced", "full")

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_de <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_de)
#write.csv(sleuth_de, 'hourly_de_newfilter.csv')
