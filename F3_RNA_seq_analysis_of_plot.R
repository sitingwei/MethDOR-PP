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
library(ggsci)
library(ggpubr)
library(patchwork)

setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")

##### volcano #####
load("process_data/Fig2B_volcano_with_4V4.RData")
p1 <- vol_dat %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
    geom_point(aes(size = size, alpha = alpha)) +
    scale_size_manual(limits = c("a", "b"), values = c(6, 3)) +
    scale_alpha_manual(limits = c("a", "b"), values = c(1, 0.5)) +
    scale_color_manual(
        limits = c("up", "no", "down"),
        values = c("#EE5A5A", "#67696B", "#638bd1")
    ) +
    theme_bw() +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 23),
        plot.subtitle = element_text(hjust = 0),
        legend.title = element_blank()
    ) +
    geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2, size = 0.5) +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.5) +
    labs(subtitle = paste("up_gene:", nrow(vol_dat %>% filter(change == "up")), ";down_gene:", nrow(vol_dat %>% filter(change == "down")))) +
    labs(x = "log2 Fold Change", y = "-log10 p-value") +
    geom_text_repel(aes(x = log2FoldChange, y = -log10(pvalue), label = topgene),
        color = "black",
        size = 4, box.padding = unit(0.25, "lines"),
        point.padding = unit(0.4, "lines"),
        segment.color = "black",
        max.overlaps = 20000
    )
p1
pdf(file = "plot_out/volcano_chart.pdf", height = 6, width = 6.5)
p1
dev.off()

##### KEGG enrichment #####
load("process_data/KEGG_enrichment.RData")
up <- filter(up, Description != "Coronavirus disease - COVID-19")
up_label <- c("p53 signaling pathway", "Oxidative phosphorylation")
down_label <- c("Steroid biosynthesis", "Metabolism of xenobiotics by cytochrome P450")
tmp <- up %>%
    filter(pvalue < 0.05) %>%
    arrange(-log10(pvalue)) %>%
    rowwise() %>%
    mutate(., labels = ifelse(Description %in% up_label, "#D20A13", "black"))
label_color <- tmp$labels
names(label_color) <- tmp$Description
p1 <- tmp %>%
    ggplot(aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue), fill = ratio)) +
    geom_bar(stat = "identity", color = "black") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    scale_fill_continuous(low = "#fc8f83", high = "#D20A13", space = "rgb") +
    labs(title = "upregulated genes", fill = "ratio") +
    theme(
        axis.text.y = element_text(colour = label_color),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(color = "black", size = 16),
        plot.title = element_text(color = "black", size = 23, face = "bold"),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(color = "black", size = 16)
    ) +
    coord_flip() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))
p1
pdf(file = "plot_out/upregulated_gene_KEGG_barplot.pdf", height = 5, width = 8)
p1
dev.off()

tmp <- down %>%
    filter(pvalue < 0.05) %>%
    arrange(-log10(pvalue)) %>%
    rowwise() %>%
    mutate(., labels = ifelse(Description %in% down_label, "#045EAA", "black"))
label_color <- tmp$labels
names(label_color) <- tmp$Description
p2 <- tmp %>%
    ggplot(aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue), fill = ratio)) +
    geom_bar(stat = "identity", color = "black") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    scale_fill_continuous(low = "#8dc8f5", high = "#223D6C", space = "rgb") +
    labs(title = "hypomethylated genes", fill = "ratio") +
    theme(
        axis.text.y = element_text(colour = label_color),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(color = "black", size = 16),
        plot.title = element_text(color = "black", size = 23, face = "bold"),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(color = "black", size = 16)
    ) +
    coord_flip() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 62))
p2
pdf(file = "plot_out/downregulated_gene_KEGG_barplot.pdf", height = 4.2, width = 8.4)
p1
dev.off()
