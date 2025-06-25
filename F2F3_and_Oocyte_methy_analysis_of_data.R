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

setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")

# CpG analysis ----------------------
file <- list.files(
    path = "/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji/data/bulk_methy_mapping/DSS_input_data",
    pattern = "^F2|F3", full.names = T
)
file2 <- list.files(
    path = "/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji/data/scMethy_mapping/DSS_input_data",
    full.names = T
)
all_file <- c(file, file2)
methy_dat <- lapply(all_file, function(x) fread(x) %>% as.data.frame())
names(methy_dat) <- word(all_file, 10, sep = "/") %>%
    word(., 1, sep = "\\.") %>%
    str_replace(., "_BS", "")
save(methy_dat, file = "process_data/all_samples_cpg_location_methylation.RData")

cpg_dat <- lapply(methy_dat, function(x) {
    tmp <- x %>%
        unite("loc", c("chr", "pos"), sep = "-", remove = T) %>%
        rowwise() %>%
        mutate(methy = round(X / N * 100, 2)) %>%
        dplyr::select(loc, methy) %>%
        column_to_rownames(var = "loc")
}) %>% do.call(cbind, .)
colnames(cpg_dat) <- names(methy_dat)
save(cpg_dat, file = "process_data/all_samples_cpg_methy_ratio.RData")

# violin chart ---------------------
F2 <- cpg_dat %>%
    dplyr::select(starts_with("F2"))
F2 <- F2[is.finite(rowSums(F2)), ]
f2_samp <- sample(rownames(F2), 198177)
F2_sum <- F2 %>%
    rownames_to_column() %>%
    filter(rowname %in% f2_samp) %>%
    rowwise() %>%
    transmute(label = rowname, F2_Con = mean(c_across(contains("Con"))), F2_PrP = mean(c_across(contains("PrP"))))
F2_sum2 <- F2_sum %>%
    pivot_longer(cols = -label, names_to = "sample", values_to = "methy") %>%
    mutate(group = word(sample, 1, sep = "_"))

Oocyte <- cpg_dat %>%
    dplyr::select(starts_with("Oocyte"))
Oocyte <- Oocyte[is.finite(rowSums(Oocyte)), ]
Oocyte_samp <- sample(rownames(Oocyte), 7585)
Oocy_sum <- Oocyte %>%
    rownames_to_column() %>%
    filter(rowname %in% Oocyte_samp) %>%
    rowwise() %>%
    transmute(label = rowname, Oocy_Con = mean(c_across(contains("Con"))), Oocy_PrP = mean(c_across(contains("PrP"))))
Oocy_sum2 <- Oocy_sum %>%
    pivot_longer(cols = -label, names_to = "sample", values_to = "methy") %>%
    mutate(group = word(sample, 1, sep = "_"))

F3 <- cpg_dat %>%
    dplyr::select(starts_with("F3"))
F3 <- F3[is.finite(rowSums(F3)), ]
f3_samp <- sample(rownames(F3), 197860)
F3_sum <- F3 %>%
    rownames_to_column() %>%
    filter(rowname %in% f3_samp) %>%
    rowwise() %>%
    transmute(label = rowname, F3_Con = mean(c_across(contains("Con"))), F3_PrP = mean(c_across(contains("PrP"))))
F3_sum2 <- F3_sum %>%
    pivot_longer(cols = -label, names_to = "sample", values_to = "methy") %>%
    mutate(group = word(sample, 1, sep = "_"))
cpg_sap_dat <- bind_rows(F2_sum2, Oocy_sum2, F3_sum2)
save(cpg_sap_dat, file = "process_data/violin_chart_cpg_methy.RData")

# DMR annotation ------------------
load("process_data/dmr_methylation.RData")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
dmr_anno <- lapply(dmr, function(x) {
    tmp <- x %>%
        dplyr::select(chr, start, end) %>%
        as(., "GRanges") %>%
        annotatePeak(., tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db") %>%
        as.data.frame() %>%
        mutate(location = word(annotation, 1, sep = " \\(")) %>%
        rename(., chr = seqnames, length = width) %>%
        inner_join(x, ., by = c("chr", "start", "end", "length"))
    return(tmp)
})
save(dmr_anno, file = "process_data/dmr_anno_methylation.RData")

# DMR number summarize --------------
sep_dmr <- lapply(dmr_anno, function(x) {
    tmp <- x %>% mutate(change = case_when(
        diff.Methy < (-0.1) ~ "hyper",
        diff.Methy > 0.1 ~ "hypo",
        T ~ "no"
    ))
    return(tmp)
})
dmr_num <- lapply(1:length(sep_dmr), function(x) {
    tmp <- sep_dmr[[x]] %>%
        group_by(change) %>%
        summarise(sum = n()) %>%
        mutate(type = names(sep_dmr)[x])
}) %>% do.call(rbind, .)
dmr_total <- dmr_num %>%
    group_by(type) %>%
    summarise(sum = sum(sum)) %>%
    mutate(change = "total") %>%
    bind_rows(dmr_num) %>%
    filter(., change %in% c("total", "hyper", "hypo")) %>%
    as.data.frame()
save(dmr_total, file = "process_data/DMR_num_summarise.RData")

# DMR region summarize ---------------------------
reg_dat <- lapply(1:length(sep_dmr), function(x) {
    tmp1 <- sep_dmr[[x]] %>%
        filter(change == "hyper") %>%
        group_by(location) %>%
        summarise(sum = n()) %>%
        mutate(type = names(sep_dmr)[x], change = "hyper", per = round(.$sum / sum(.$sum) * 100, 2)) %>%
        rowwise() %>%
        mutate(label = ifelse(per %in% sort(.$per, decreasing = T)[1:3], paste0(per, "%"), "")) %>%
        dplyr::select(-per) %>%
        as.data.frame()
    tmp2 <- sep_dmr[[x]] %>%
        filter(change == "hypo") %>%
        group_by(location) %>%
        summarise(sum = n()) %>%
        mutate(type = names(sep_dmr)[x], change = "hypo", per = round(.$sum / sum(.$sum) * 100, 2)) %>%
        rowwise() %>%
        mutate(label = ifelse(per %in% sort(.$per, decreasing = T)[1:3], paste0(per, "%"), "")) %>%
        dplyr::select(-per) %>%
        as.data.frame()
    tmp <- bind_rows(tmp1, tmp2)
    return(tmp)
}) %>% do.call(rbind, .)
save(reg_dat, file = "process_data/DMR_region_summarise.RData")

# DMR change Sankey diagram ------------------------
dmr_exp <- lapply(dmr, function(x) {
    tmp <- x %>%
        transmute(chr = chr, start = start - 500, end = end + 500, diff = -diff.Methy)
    return(tmp)
})
CHR <- paste0("chr", c(1:19, "X", "Y"))
f2_cpg <- lapply(CHR, function(x) {
    cpg <- data.frame(chr = dmltest_dat[[1]]$chr, pos = dmltest_dat[[1]]$pos) %>% filter(chr == x)
    dmr <- dmr_exp[[1]]
    tmp <- dmr %>%
        arrange(chr, start) %>%
        filter(chr == x) %>%
        mutate(diff = case_when(diff > 0 ~ "hyper", diff < 0 ~ "hypo", T ~ "no"))
    ind <- tmp %>%
        mutate(start = start - 1) %>%
        pivot_longer(-c(chr, diff)) %>%
        .$value %>%
        c(., Inf)
    label <- tmp %>%
        pivot_longer(-c(chr, diff)) %>%
        .$diff
    label[seq(2, length(label), by = 2)] <- "NA"
    cpg$rank <- cut(cpg$pos, breaks = ind, labels = label)
    cpg1 <- cpg %>% filter(rank != "NA")
    return(cpg1)
}) %>%
    do.call(rbind, .) %>%
    mutate(type = "F2", rank = as.character(rank)) %>%
    unite(., "label", c("chr", "pos"), sep = "-", remove = F)
Oocy_cpg <- lapply(CHR, function(x) {
    cpg <- data.frame(chr = dmltest_dat[[2]]$chr, pos = dmltest_dat[[2]]$pos) %>% filter(chr == x)
    dmr <- dmr_exp[[2]]
    tmp <- dmr %>%
        arrange(chr, start) %>%
        filter(chr == x) %>%
        mutate(diff = case_when(diff > 0 ~ "hyper", diff < 0 ~ "hypo", T ~ "no"))
    ind <- tmp %>%
        mutate(start = start - 1) %>%
        pivot_longer(-c(chr, diff)) %>%
        .$value %>%
        c(., Inf)
    label <- tmp %>%
        pivot_longer(-c(chr, diff)) %>%
        .$diff
    label[seq(2, length(label), by = 2)] <- "NA"
    cpg$rank <- cut(cpg$pos, breaks = ind, labels = label)
    cpg1 <- cpg %>% filter(rank != "NA")
    return(cpg1)
}) %>%
    do.call(rbind, .) %>%
    mutate(type = "Oocyte", rank = as.character(rank)) %>%
    unite(., "label", c("chr", "pos"), sep = "-", remove = F)
f3_cpg <- lapply(CHR, function(x) {
    cpg <- data.frame(chr = dmltest_dat[[3]]$chr, pos = dmltest_dat[[3]]$pos) %>% filter(chr == x)
    dmr <- dmr_exp[[3]]
    tmp <- dmr %>%
        arrange(chr, start) %>%
        filter(chr == x) %>%
        mutate(diff = case_when(diff > 0 ~ "hyper", diff < 0 ~ "hypo", T ~ "no"))
    ind <- tmp %>%
        mutate(start = start - 1) %>%
        pivot_longer(-c(chr, diff)) %>%
        .$value %>%
        c(., Inf)
    label <- tmp %>%
        pivot_longer(-c(chr, diff)) %>%
        .$diff
    label[seq(2, length(label), by = 2)] <- "NA"
    cpg$rank <- cut(cpg$pos, breaks = ind, labels = label)
    cpg1 <- cpg %>% filter(rank != "NA")
    return(cpg1)
}) %>%
    do.call(rbind, .) %>%
    mutate(type = "F3", rank = as.character(rank)) %>%
    unite(., "label", c("chr", "pos"), sep = "-", remove = F)
all_cpg <- reduce(list(f2_cpg$label, f3_cpg$label), union)
tmp1 <- f2_cpg %>%
    dplyr::select(label, rank) %>%
    bind_rows(., data.frame(label = setdiff(all_cpg, f2_cpg$label), rank = "no")) %>%
    rename(F2 = rank)
tmp2 <- Oocy_cpg %>%
    dplyr::select(label, rank) %>%
    filter(label %in% intersect(all_cpg, Oocy_cpg$label)) %>%
    bind_rows(., data.frame(label = setdiff(all_cpg, Oocy_cpg$label), rank = "no")) %>%
    rename(Oocyte = rank)
tmp3 <- f3_cpg %>%
    dplyr::select(label, rank) %>%
    bind_rows(., data.frame(label = setdiff(all_cpg, f3_cpg$label), rank = "no")) %>%
    rename(F3 = rank)
all_dat <- Reduce(function(x, y) merge(x, y, by = "label"), list(tmp1, tmp2, tmp3), accumulate = FALSE) %>%
    group_by(F2, Oocyte, F3) %>%
    summarise(Freq = n()) %>%
    mutate(change = F2)
save(all_dat, file = "process_data/DMR_change_Sankey_diagram_tmp.RData")

# overlap gene ----------------
cpg_dat1 <- cpg_dat %>% rownames_to_column()
over_dat <- cpg_dat1 %>%
    filter(str_detect(rowname, "chr10-692")) %>%
    mutate(pos = word(rowname, 2, sep = "-") %>% as.numeric()) %>%
    filter(pos >= 69272415 & pos <= 69283309) %>%
    dplyr::select(-rowname) %>%
    pivot_longer(-pos) %>%
    mutate(group = gsub("(.*)\\d", "\\1", name))
over_dat[which(over_dat$value == "NaN"), "value"] <- NA

# ucsc gemone browers data
ucsc_dat <- over_dat %>%
    group_by(pos, group) %>%
    summarise(dist = mean(value, na.rm = T)) %>%
    as.data.frame()
ucsc_dat[which(ucsc_dat$dist == "NaN"), "dist"] <- NA
F2_Con <- ucsc_dat %>%
    filter(group == "F2_Con") %>%
    transmute(chr = "chr10", start = pos - 1, end = pos, methy = dist) %>%
    na.omit()
write.table(F2_Con, file = "process_data/F2_Con_over_gene.bedGraph", sep = "\t", col.names = F, row.names = F, quote = F)
F2_PrP <- ucsc_dat %>%
    filter(group == "F2_PrP") %>%
    transmute(chr = "chr10", start = pos - 1, end = pos, methy = dist) %>%
    na.omit()
write.table(F2_PrP, file = "process_data/F2_PrP_over_gene.bedGraph", sep = "\t", col.names = F, row.names = F, quote = F)

Oocy_con <- ucsc_dat %>%
    filter(group == "Oocyte_con") %>%
    transmute(chr = "chr10", start = pos - 1, end = pos, methy = dist) %>%
    na.omit()
write.table(Oocy_con, file = "process_data/Oocy_Con_over_gene.bedGraph", sep = "\t", col.names = F, row.names = F, quote = F)
Oocy_PrP <- ucsc_dat %>%
    filter(group == "Oocyte_PrP") %>%
    transmute(chr = "chr10", start = pos - 1, end = pos, methy = dist) %>%
    na.omit()
write.table(Oocy_PrP, file = "process_data/Oocy_PrP_over_gene.bedGraph", sep = "\t", col.names = F, row.names = F, quote = F)

F3_Con <- ucsc_dat %>%
    filter(group == "F3_Con") %>%
    transmute(chr = "chr10", start = pos - 1, end = pos, methy = dist) %>%
    na.omit()
write.table(F3_Con, file = "process_data/F3_Con_over_gene.bedGraph", sep = "\t", col.names = F, row.names = F, quote = F)
F3_PrP <- ucsc_dat %>%
    filter(group == "F2_PrP") %>%
    transmute(chr = "chr10", start = pos - 1, end = pos, methy = dist) %>%
    na.omit()
write.table(F3_PrP, file = "process_data/F3_PrP_over_gene.bedGraph", sep = "\t", col.names = F, row.names = F, quote = F)

over_dat <- cpg_dat1 %>%
    filter(str_detect(rowname, "chr10-692")) %>%
    mutate(pos = word(rowname, 2, sep = "-") %>% as.numeric()) %>%
    filter(pos >= 69272415 & pos <= 69283309) %>%
    dplyr::select(-rowname) %>%
    pivot_longer(-pos) %>%
    mutate(group = gsub("(.*)\\d", "\\1", name))
over_dat[which(over_dat$value == "NaN"), "value"] <- NA
save(over_dat, file = "process_data/overlap_gene_Rhobtb1.RData")
