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

##### call DMR #####
setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")
load("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji/DSS_dmltest_data/F3_dmltest.RData")

dmr <- callDMR(dmlTest, delta = 0.1, p.threshold = 0.05, minlen = 40, minCG = 2, dis.merge = 50, pct.sig = 0.5)
save(dmr, file = "process_data/Methy_callDMR_delta0.1_minCG2_minlen40.RData")

##### chipseeker annotation #####
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)

dmr_anno <- dmr %>%
    dplyr::select(chr, start, end) %>%
    as(., "GRanges") %>%
    annotatePeak(., tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db") %>%
    as.data.frame() %>%
    mutate(location = word(annotation, 1, sep = " \\(")) %>%
    rename(., chr = seqnames, length = width) %>%
    inner_join(dmr, ., by = c("chr", "start", "end", "length"))
dmr_anno_hyper <- dmr_anno %>% filter(diff.Methy < (-0.1))
dmr_anno_hypo <- dmr_anno %>% filter(diff.Methy > 0.1)
save(dmr_anno, file = "process_data/Methy_DMR_annotation_with_dalta0.1_minCG2_minlen40.RData")

##### DMR num summarise #####
DMR_num <- data.frame(
    type = c("hyper", "hypo", "all"),
    DMR = c(nrow(dmr_anno_hyper), nrow(dmr_anno_hypo), nrow(dmr_anno)),
    gene = c(length(dmr_anno_hyper$SYMBOL %>% unique()), length(dmr_anno_hypo$SYMBOL %>% unique()), length(dmr_anno$SYMBOL %>% unique()))
)
save(DMR_num, file = "process_data/DMR_num_summarise.RData")

##### DMR region summarise #####
region_table <- bind_rows(
    data.frame(location = names(table(dmr_anno$location)), count = as.numeric(table(dmr_anno$location)), type = "all"),
    data.frame(location = names(table(dmr_anno_hyper$location)), count = as.numeric(table(dmr_anno_hyper$location)), type = "hyper"),
    data.frame(location = names(table(dmr_anno_hypo$location)), count = as.numeric(table(dmr_anno_hypo$location)), type = "hypo")
)
save(region_table, file = "process_data/DMR_region_summarise.RData")

##### DMR circos #####
data <- dmr_anno %>% transmute(chr, start, end, diff.methy = -diff.Methy)
hyper_dat <- data %>% filter(diff.methy > 0.1)
hypo_dat <- data %>% filter(diff.methy < (-0.1))
save(data, hyper_dat, hypo_dat, file = "process_data/DMR_circos.RData")

##### dotplot #####
load("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji/process_data/RNA_deseq2_analysis_with_symbol_with_delC1P2_with_combat.RData")

methy <- dmr_anno %>%
    dplyr::select(chr:end, diff.Methy, SYMBOL) %>%
    rename(symbol = SYMBOL)
DEG <- sym_F3 %>%
    dplyr::select(symbol, log2FoldChange, pvalue)
dat <- inner_join(methy, DEG, by = "symbol") %>%
    na.omit() %>%
    filter(abs(diff.Methy) > 0.1) %>%
    mutate(type = case_when(
        diff.Methy > 0.1 & log2FoldChange > 1 ~ "hypo_up",
        diff.Methy > 0.1 & log2FoldChange < (-1) ~ "hypo_down",
        diff.Methy < (-0.1) & log2FoldChange > 1 ~ "hyper_up",
        diff.Methy < (-0.1) & log2FoldChange < (-1) ~ "hyper_down",
        T ~ "no"
    ))
save(dat, file = "process_data/methy_and_RNA_scatter_Plot.RData")

##### venn plot #####
DEG_gene <- filter(sym_F3, abs(log2FoldChange) > 1 & pvalue < 0.05) %>%
    .$symbol %>%
    na.omit() %>%
    unique()
hyper_gene <- filter(dmr_anno, diff.Methy < (-0.1)) %>%
    .$SYMBOL %>%
    na.omit() %>%
    unique()
hypo_gene <- filter(dmr_anno, diff.Methy > (0.1)) %>%
    .$SYMBOL %>%
    na.omit() %>%
    unique()
save(DEG_gene, hyper_gene, hypo_gene, file = "process_data/venn_plot.RData")
