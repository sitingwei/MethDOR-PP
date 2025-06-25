library(tidyverse)
library(data.table)
library(grid)
library(VennDiagram)
library(magrittr)
library(ggrepel)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(gg.gap)

setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")

# violin chart ----------------------
load("process_data/violin_chart_cpg_methy.RData")
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = .95, .drop = TRUE) {
    library(plyr)

    length2 <- function(x, na.rm = FALSE) {
        if (na.rm) {
            sum(!is.na(x))
        } else {
            length(x)
        }
    }
    
    datac <- ddply(data, groupvars,
        .drop = .drop,
        .fun = function(xx, col) {
            c(
                N = length2(xx[[col]], na.rm = na.rm),
                mean = mean(xx[[col]], na.rm = na.rm),
                sd = sd(xx[[col]], na.rm = na.rm)
            )
        },
        measurevar
    )

    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N) # Calculate standard error of the mean

    ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

dat_summary <- summarySE(cpg_sap_dat, measurevar = "methy", groupvars = c("group", "sample"))
p1 <- ggplot(cpg_sap_dat, aes(x = factor(group, levels = c("F3", "Oocy", "F2")), y = methy, fill = sample)) +
    geom_violin() +
    geom_point(data = dat_summary, aes(x = factor(group, levels = c("F3", "Oocy", "F2")), y = methy), pch = 19, position = position_dodge(0.9), size = 2.5, color = "#E25C49") +
    scale_fill_manual(
        values = c("#c3f1c4", "#59A95A", "#fbcda6", "#F7903D", "#b8d3ee", "#4D85BD"),
        limits = c("F2_Con", "F2_PrP", "Oocy_Con", "Oocy_PrP", "F3_Con", "F3_PrP")
    ) +
    theme_bw() +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 23),
        # axis.title.y.left = element_text(margin = margin(0, 0.7, 0, 0, "cm")),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 15)
    ) +
    labs(y = "CpG methylation(%)") +
    coord_flip() +
    scale_x_discrete(labels = c("F3", "Oocyte", "F2"))
p1
pdf(file = "plot_out/violin_chart_of_CpG_methy.pdf", height = 7.5, width = 6)
p1
dev.off()

# DMR number and gene summarize ---------------------
load("process_data/DMR_num_summarise.RData")

p <- dmr_total %>%
    ggplot(aes(x = factor(type, levels = c("F2", "Oocyte", "F3")), y = sum, color = change, group = change, shape = change)) +
    geom_point(size = 2.5) +
    geom_line(size = 0.7) +
    scale_color_manual(values = c("#F5724E", "#0E9CC5", "#C15A9D"), limits = c("hyper", "hypo", "total")) +
    theme_classic() +
    theme(
        # panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        # panel.border = element_rect(color = "black", size = 2),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.ticks.length = unit(5, "pt"),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 15),
        legend.title = element_blank(),
        # legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.8, "cm")
    ) +
    geom_text_repel(aes(label = sum), nudge_y = -1.5, nudge_x = -0.1, size = 2.5) +
    labs(y = "numbers of DMRs")
p
b <- dmr_gene_total %>%
    ggplot(aes(x = factor(type, levels = c("F2", "Oocyte", "F3")), y = gene_sum, fill = change, group = change)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#F5724E", "#0E9CC5", "#C15A9D"), limits = c("hyper", "hypo", "total")) +
    theme_classic() +
    theme(
        # panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        # panel.border = element_rect(color = "black", size = 2),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.ticks.length = unit(5, "pt"),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 15),
        legend.title = element_blank(),
        # legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.8, "cm")
    ) +
    geom_text(aes(label = gene_sum), position = position_dodge(0.9), vjust = 0, size = 2.5) +
    labs(y = "numbers of genes")
b
p1 <- p / b
pdf(file = "plot_out/line_chart_of_DMR_sum.pdf", height = 5, width = 4.7)
p1
dev.off()

##### 3.Fig7C DMR region summarise #####
load("process_data/DMR_region_summarise.RData")
hyper <- function(dat, a) {
    tmp <- dat %>%
        filter(change == "hyper" & type == a) %>%
        ggplot(aes(x = type, y = sum, fill = factor(location, levels = c("Downstream", "3' UTR", "Exon", "Intron", "Promoter", "5' UTR", "Distal Intergenic")))) +
        geom_bar(stat = "identity", position = "stack", color = "black", width = 1) +
        scale_fill_manual(
            values = c("#5e5f8d", "#CC3399", "#009966", "#F85F67", "#FF9933", "#deb573", "#87cde9"),
            limits = c("Downstream", "3' UTR", "Exon", "Intron", "Promoter", "5' UTR", "Distal Intergenic")
        ) +
        coord_polar(theta = "y") +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            line = element_line(color = "black", size = 0.1)
        )
    return(tmp)
}
p1 <- hyper(reg_dat, "F2")
p3 <- hyper(reg_dat, "Oocyte")
p5 <- hyper(reg_dat, "F3")

hypo <- function(dat, a) {
    tmp <- dat %>%
        filter(change == "hypo" & type == a) %>%
        ggplot(aes(x = type, y = sum, fill = factor(location, levels = c("Downstream", "3' UTR", "Exon", "Intron", "Promoter", "5' UTR", "Distal Intergenic")))) +
        geom_bar(stat = "identity", position = "stack", color = "black", width = 1) +
        scale_fill_manual(
            values = c("#5e5f8d", "#CC3399", "#009966", "#F85F67", "#FF9933", "#deb573", "#87cde9"),
            limits = c("Downstream", "3' UTR", "Exon", "Intron", "Promoter", "5' UTR", "Distal Intergenic")
        ) +
        coord_polar(theta = "y") +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 11),
            legend.key.size = unit(0.5, "cm"),
            line = element_line(color = "black", size = 0.1)
        )
}
p2 <- hypo(reg_dat, "F2")
p4 <- hypo(reg_dat, "Oocyte")
p6 <- hypo(reg_dat, "F3")

p <- (p1 / p3 / p5) | (p2 / p4 / p6) + plot_layout(heights = 1)
pdf(file = "plot_out/DMR_region_pie_chart.pdf", height = 9, width = 7.8)
p
dev.off()

# DMR change Sankey diagram --------------------
library(ggalluvial)
load("process_data/DMR_change_Sankey_diagram_tmp.RData")

data <- all_dat %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(cols = c("F2", "Oocyte", "F3"), names_to = "type", values_to = "stratum") %>%
    mutate(stratum = as.factor(stratum), type = as.factor(type))
p1 <- data %>%
    ggplot(aes(x = type, stratum = stratum, alluvium = rowname, y = Freq, label = stratum)) +
    scale_x_discrete(limits = c("F2", "Oocyte", "F3"), expand = c(.1, .1)) +
    geom_alluvium(aes(fill = change), knot.pos = 1 / 5) +
    geom_stratum(aes(fill = stratum)) +
    scale_fill_manual(values = c("#fa7d5a", "#0E9CC5", "#7671B5"), limits = c("hyper", "hypo", "no")) +
    geom_text(stat = "stratum", size = 6, color = "white", fontface = "bold") +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
    )
p1
pdf(file = "plot_out/DMR_change_Sankey_diagram.pdf", height = 7, width = 7)
p1
dev.off()

# overlap gene ----------------------
load("process_data/overlap_gene_Rhobtb1.RData")
dat <- data.frame(x = c(69276919, 69277810, 69276919, 69277810), y = c(0, 0, 100, 100))
p1 <- over_dat %>%
    filter(pos >= 69275415 & pos <= 69279309) %>%
    ggplot(aes(x = pos, y = value, color = group)) +
    geom_smooth(se = F) +
    geom_rug(sides = "b", color = "black", size = 0.1) +
    scale_color_manual(
        values = c("#c3f1c4", "#59A95A", "#fbcda6", "#f5862b", "#b8d3ee", "#4D85BD"),
        limits = c("F2_Con", "F2_PrP", "Oocyte_con", "Oocyte_PrP", "F3_Con", "F3_PrP"), na.value = NA
    ) +
    theme_bw() +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 13),
        # axis.title.y.left = element_text(margin = margin(0, 0.7, 0, 0, "cm")),
        legend.position = "none"
    ) +
    labs(y = "CpG methylation(%)") +
    coord_cartesian(ylim = c(50, 100)) +
    geom_vline(xintercept = 69276919, colour = "#bebebe", linetype = "dashed", size = 0.1) +
    geom_vline(xintercept = 69277810, colour = "#bebebe", linetype = "dashed", size = 0.1)
p1
pdf(file = "plot_out/overlap_gene_V3.pdf", height = 1.5, width = 5)
p1
dev.off()
