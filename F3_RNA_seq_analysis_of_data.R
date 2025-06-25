library(tidyverse)
library(data.table)
library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)
library(sva)
library(ComplexHeatmap)

setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")

##### get counts #####
file <- list.files(path = "/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji/all_analysis_together/data/bulk_rna_with_own_mapping", pattern = "bulk_RNA_F3", full.names = T)
counts <- lapply(file, function(x) read.table(x, row.names = 1, skip = 4)[1]) %>% do.call(cbind, .)
colnames(counts) <- gsub("(.*)/(.*)/bulk_RNA_(.*)ReadsPerGene.out.tab", "\\3", file)
write.csv(counts, file = "process_data/RNA_row_counts_with_all_sample.csv")

##### deseq2 analysis #####
difGene <- function(mycounts, control_pattern, case_pattern, geneSymbol_col) {
    library(tidyverse)
    library(magrittr)
    library(DESeq2)
    mycounts <- mycounts %>%
        dplyr::select(
            str_subset(colnames(.), control_pattern),
            str_subset(colnames(.), case_pattern)
        )

    colData <- data.frame(
        row.names = colnames(mycounts),
        condition = c(
            stringr::str_extract(colnames(mycounts), control_pattern),
            stringr::str_extract(colnames(mycounts), case_pattern)
        ) %>% na.omit()
    )

    print(colData)

    colData$condition <- factor(colData$condition, levels = c(control_pattern, case_pattern))
    dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)
    dds <- DESeq(dds)

    print(resultsNames(dds)[2])

    tmp_res <- results(dds)

    print(summary(tmp_res))

    dds_counts <- counts(dds, normalize = TRUE) %>% as.data.frame() # library size normalize matrix
    colnames(dds_counts) <- paste0("dds_counts_", colnames(dds_counts))
    colnames(mycounts) <- paste0("raw_", colnames(mycounts))

    if (ncol(mycounts) > 50) {
        vsd <- assay(vst(dds, blind = TRUE)) %>% as.data.frame()
        colnames(vsd) <- paste0("vsd_", colnames(vsd))
        res <- bind_cols(as.data.frame(tmp_res), dds_counts, vsd, mycounts)
    } else {
        rld <- assay(rlogTransformation(dds)) %>% as.data.frame() ## rlogs 矩阵
        colnames(rld) <- paste0("rld_", colnames(rld))
        res <- bind_cols(as.data.frame(tmp_res), dds_counts, rld, mycounts)
    }

    res <- res %>%
        mutate(ensembl = rownames(.) %>% stringr::word(., 1, sep = "\\.")) %>%
        dplyr::select(ensembl, everything())

    return(res)
}
F3 <- difGene(mycounts = count_cob, control_pattern = "Con", case_pattern = "PrP")

gene <- rownames(count_cob) %>% stringr::word(., 1, sep = "\\.") 
gene2 <- bitr(gene, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db) %>%
    dplyr::select(symbol = SYMBOL, ensembl = ENSEMBL) %>%
    filter(symbol != "a")

sym_F3 <- merge(F3, gene2, by = "ensembl") %>% dplyr::select(symbol, everything())
DEG_up <- filter(sym_F3, log2FoldChange > 1 & pvalue < 0.05) %>%
    .$symbol %>%
    na.omit() %>%
    unique() # 120
DEG_down <- filter(sym_F3, log2FoldChange < (-1) & pvalue < 0.05) %>%
    .$symbol %>%
    na.omit() %>%
    unique() # 331
save(DEG_up, DEG_down, sym_F3, file = "process_data/RNA_deseq2_symgene_DEGgene_logfc1_with_up_and_down.RData")

##### volcano #####
top10gene <- sym_F3 %>%
    arrange(desc(-log10(.$pvalue))) %>%
    .[1:15, ] %>%
    .$symbol
top50gene <- sym_F3 %>%
    arrange(desc(-log10(.$pvalue))) %>%
    .[1:50, ] %>%
    .$symbol
gene <- c(
    "Cox7a2l", "Atp4b", "Cyct", "Ndufb4c", "Serpine1", "Cyct", "Rprm",
    "Adh7", "Mgst2", "Gm4846", "Ugt1a7c", "Hsd11b1", "Adh7", "Mgst2", "Ugt1a7c"
)
vol_dat <- sym_F3 %>%
    transmute(
        symbol = symbol,
        log2FoldChange = log2FoldChange,
        pvalue = pvalue,
        change = case_when(
            log2FoldChange > 1 & pvalue < 0.05 ~ "up",
            log2FoldChange < (-1) & pvalue < 0.05 ~ "down",
            T ~ "no"
        ),
        topgene = ifelse(.$symbol %in% top10gene, .$symbol, ""),
        size = ifelse(.$symbol %in% top50gene, "a", "b"),
        alpha = ifelse(.$symbol %in% gene, "a", "b")
    ) %>%
    na.omit()
save(vol_dat, file = "process_data/volcano_with_4V4_combat.RData")

##### KEGG enrichment #####
counts <- read.csv("process_data/RNA_row_counts_with_all_sample.csv", header = T) %>% column_to_rownames(var = "X")
F3 <- difGene(mycounts = counts, control_pattern = "Con", case_pattern = "PrP")

gene <- rownames(counts) %>% stringr::word(., 1, sep = "\\.") # word从字符串向量中提取字符
gene2 <- bitr(gene, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db) %>%
    dplyr::select(symbol = SYMBOL, ensembl = ENSEMBL) %>%
    filter(symbol != "a") # 32980

sym_F3 <- merge(F3, gene2, by = "ensembl") %>% dplyr::select(symbol, everything())
DEG_up <- filter(sym_F3, log2FoldChange > 1 & pvalue < 0.05) %>%
    .$symbol %>%
    na.omit() %>%
    unique() # 120
DEG_down <- filter(sym_F3, log2FoldChange < (-1) & pvalue < 0.05) %>%
    .$symbol %>%
    na.omit() %>%
    unique() # 331
save(DEG_up, DEG_down, sym_F3, file = "process_data/RNA_deseq2_symgene_DEGgene_logfc1_with_up_and_down_with_4V4.RData")

get_KEGG <- function(data) {
    tmp <- bitr(data, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db") %>%
        .$ENTREZID %>%
        enrichKEGG(., organism = "mmu", keyType = "kegg", pvalueCutoff = 1) %>%
        setReadable(., OrgDb = org.Mm.eg.db, keyType = "ENTREZID") %>%
        .@result %>%
        arrange(pvalue) %>%
        mutate(ratio = as.numeric(word(.$GeneRatio, 1, sep = "/")) / as.numeric(word(.$GeneRatio, 2, sep = "/")))
    return(tmp)
}
up <- get_KEGG(DEG_up)
down <- get_KEGG(DEG_down)
save(up, down, file = "process_data/KEGG_enrichment.RData")

