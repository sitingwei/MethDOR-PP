library(DSS)
library(bsseq)
library(tidyverse)
library(data.table)
library(grid)
library(VennDiagram)
library(magrittr)
library(ggrepel)
library(ggplot2)
library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggsci)
library(ggpubr)
library(patchwork)
library(BiocParallel)

# DSS input --------------------------
setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")

file <- list.files(path = "data/", pattern = "(.*).txt", full.names = T)
input <- lapply(file, function(x) fread(x))
names(input) <- word(file, 2, sep = "//") %>% word(., 1, sep = "\\.")
bs_data <- lapply(1:length(input), function(num) {
  tmp <- input[[num]] %>%
    transmute(
      chr = V1,
      pos = ifelse(V3 == "-", V2 - 1, V2),
      mc = V4,
      umc = V5
    ) %>%
    group_by(chr, pos) %>%
    summarise(mC = sum(mc), umC = sum(umc)) %>%
    as.data.frame()
  return(tmp)
})
names(bs_data) <- names(input)
dss_input <- lapply(1:length(bs_data), function(num) {
  tmp <- bs_data[[num]] %>%
    transmute(
      chr = paste0("chr", chr),
      pos = pos,
      N = mC + umC,
      X = mC
    )
  write.table(tmp, file = paste0("process_data/DSS_input_data/", names(bs_data)[num], ".txt"), sep = "\t", quote = F, row.names = F)
})

# DSS dmltest ######
file <- list.files(path = "process_data/DSS_input_data", pattern = "^F2_BS(.*).txt", full.names = T)
input <- lapply(file, function(x) fread(x))
names(input) <- word(file, 3, sep = "/") %>% word(., 1, sep = "\\.")

BSobj <- makeBSseqData(input, names(input))

mParam <- MulticoreParam(workers = 20, progressbar = TRUE)
dmlTest <- DMLtest(BSobj,
                   group1 = str_subset(names(input), "Con"), group2 = str_subset(names(input), "PrP"),
                   BPPARAM = mParam, smoothing = TRUE
)
save(dmlTest, file = "process_data/DSS_dmltest_data/F2_dmltest.RData")

file <- list.files(path = "process_data/DSS_input_data", pattern = "^F2_Oocy(.*).txt", full.names = T)
input <- lapply(file, function(x) fread(x))
names(input) <- word(file, 3, sep = "/") %>% word(., 1, sep = "\\.")

BSobj <- makeBSseqData(input, names(input))

mParam <- MulticoreParam(workers = 20, progressbar = TRUE)
dmlTest <- DMLtest(BSobj,
                   group1 = str_subset(names(input), "Con"), group2 = str_subset(names(input), "PrP"),
                   BPPARAM = mParam, smoothing = TRUE
)
save(dmlTest, file = "process_data/DSS_dmltest_data/F2Oocy_dmltest.RData")

file <- list.files(path = "process_data/DSS_input_data", pattern = "^F3_BS(.*).txt", full.names = T)
input <- lapply(file, function(x) fread(x))
names(input) <- word(file, 3, sep = "/") %>% word(., 1, sep = "\\.")

BSobj <- makeBSseqData(input, names(input))

mParam <- MulticoreParam(workers = 20, progressbar = TRUE)
dmlTest <- DMLtest(BSobj,
                   group1 = str_subset(names(input), "Con"), group2 = str_subset(names(input), "PrP"),
                   BPPARAM = mParam, smoothing = TRUE
)
save(dmlTest, file = "process_data/DSS_dmltest_data/F3_dmltest.RData")

# call DMR ------------------------
file <- list.files(
  path = "/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji/process_data/DSS_dmltest_data",
  pattern = ".RData$", full.names = T
)

dss_dmltest <- lapply(file, function(dat) {
  load(dat)
  tmp <- dmlTest
  return(tmp)
})
names(dss_dmltest) <- word(file, 10, sep = "/") %>% word(., 1, sep = "\\.")

save(dss_dmltest, file = "process_data/F1_to_F3_DSS_dmltest.RData")