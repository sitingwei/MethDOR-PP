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

setwd("/mnt/memory2/weisiting/project_transgenerational_epigenetic_of_tongji")

##### DMR num summarise #####
load("process_data/DMR_num_summarise.RData")
p1 <- DMR_num %>%
    ggplot(aes(x = factor(type, levels = c("hyper", "hypo", "all")))) +
    scale_y_continuous(
        name = "number of DMRs", breaks = c(seq(500, 3000, 500)),
        sec.axis = sec_axis(~ . * 1, breaks = c(seq(500, 3000, 500)), name = "number of genes")
    ) +
    theme_bw() +
    geom_segment(aes(
        x = factor(type, levels = c("hyper", "hypo", "all")),
        y = 500,
        xend = factor(type, levels = c("hyper", "hypo", "all")),
        yend = DMR
    ), color = "gray", size = 1.5) +
    geom_point(aes(y = DMR), size = 15, color = "#334f65") +
    geom_point(aes(y = gene), size = 15, color = "#f05326") +
    geom_text(aes(y = DMR, label = DMR), color = "white", size = 4.5, fontface = "bold") +
    geom_text(aes(y = gene, label = gene), color = "white", size = 4.5, fontface = "bold") +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y.left = element_text(color = "#334f65"),
        axis.text.y.right = element_text(color = "#f05326"),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 23),
        axis.title.y.right = element_text(margin = margin(0, 0, 0, 0.7, "cm")),
        axis.title.y.left = element_text(margin = margin(0, 0.7, 0, 0, "cm"))
    )
p1
pdf(file = "plot_out/lollipop_chart_V2.pdf", height = 6, width = 6)
p1
dev.off()

##### DMR region summarise #####
load("process_data/DMR_region_summarise.RData")
p2 <- region_table %>%
    ggplot(aes(
        x = factor(type, levels = c("hyper", "hypo", "all")), y = count,
        fill = factor(location, levels = c("Downstream", "3' UTR", "Exon", "Intron", "Promoter", "5' UTR", "Distal Intergenic"))
    )) +
    geom_bar(stat = "identity", position = "fill", color = "black", width = 0.8) +
    scale_fill_manual(values = c("#5e5f8d", "#CC3399", "#009966", "#F85F67", "#FF9933", "#deb573", "#87cde9")) +
    theme_bw() +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 23),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 15)
    ) +
    labs(y = "relative frequency")
p2
pdf(file = "plot_out/barplot_V2.pdf", height = 6, width = 6)
p2
dev.off()

##### DMR circos #####
library(circlize)
load("process_data/DMR_circos.RData")
bed_list <- rbind(arrange(hyper_dat, desc(diff.methy))[1:50, ], arrange(hypo_dat, diff.methy)[1:50, ])
pdf(file = "plot_out/Fig_6C_DMR_circos_V3.pdf", height = 7, width = 7)
cytoband <- read.cytoband(species = "mm10")$df
circos.clear()
circos.par("gap.degree" = c(rep(c(1), 20), 10), "start.degree" = 82)
circos.initializeWithIdeogram(species = "mm10", plotType = NULL)
circos.track(
    ylim = c(0, 1), track.height = mm_h(1),
    panel.fun = function(x, y) {
        circos.text(
            CELL_META$xcenter,
            CELL_META$ylim[2] + mm_y(5.5),
            CELL_META$sector.index,
            cex = 0.8 * par("cex"),
            niceFacing = TRUE
        )
    },
    cell.padding = c(0, 0, 0, 0), bg.border = NA
)
circos.genomicIdeogram(cytoband)
circos.track(
    track.index = get.current.track.index(),
    panel.fun = function(x, y) {
        circos.genomicAxis(h = "top", labels.cex = 0.4 * par("cex"), labels.facing = "clockwise")
    }
)
circos.genomicTrack(
    data,
    ylim = c(-0.8, 0.8),
    panel.fun = function(region, value, ...) {
        # for (h in seq(-0.8, 0.8, by = 0.4)) {
        #     circos.lines(CELL_META$cell.xlim, c(h, h),
        #         lty = 3, col = "#AAAAAA"
        #     )
        # }
        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = "#888888")
        circos.genomicPoints(region, value,
            pch = 16, cex = 0.5,
            col = ifelse(value[[1]] > 0.1, "#E41A1C", "#377EB8")
        )
    }
)
circos.yaxis(
    side = "left",
    at = seq(-0.8, 0.8, by = 0.4),
    sector.index = get.all.sector.index()[1],
    labels.cex = 0.4
)
circos.genomicDensity(
    hyper_dat,
    col = c("#E41A1C"), track.height = 0.1
)
circos.yaxis(
    side = "left",
    at = seq(0, 0.003, by = 0.001),
    sector.index = get.all.sector.index()[1],
    labels.cex = 0.3
)
circos.genomicDensity(
    hypo_dat,
    col = c("#377EB8"), track.height = 0.1
)
circos.yaxis(
    side = "left",
    at = seq(0, 0.003, by = 0.001),
    sector.index = get.all.sector.index()[1],
    labels.cex = 0.3
)
circos.genomicTrack(bed_list,
    panel.fun = function(region, value, ...) {
        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "#00000040")
        x <- (region[[1]] + region[[2]]) / 2
        circos.segments(x, 0, x, value[[1]], col = ifelse(value[[1]] > 0.1, "#F4573A", "#00A6B1"))
    }
)
circos.yaxis(
    side = "left",
    at = seq(-0.8, 0.8, by = 0.4),
    sector.index = get.all.sector.index()[1],
    labels.cex = 0.3
)
dev.off()

##### dotplot #####
library(ggpointdensity)
library(ggExtra)

load("process_data/methy_and_RNA_scatter_Plot.RData")

p1 <- dat %>%
    ggplot(aes(x = -diff.Methy, y = log2FoldChange, color = type)) +
    geom_point(size = 2) +
    scale_color_manual(
        limits = c("hyper_up", "hyper_down", "hypo_up", "hypo_down", "no"),
        values = c("#F4573A", "#00A6B1", "#EF9136", "#437EB0", "#cecdcd")
    ) +
    labs(x = "diff.methy", y = "log2foldchange") +
    theme_bw() +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 22),
        legend.title = element_blank(),
        legend.key.size = unit(1.3, "cm"),
        legend.text = element_text(size = 16)
    )
p1
pdf(file = "plot_out/express_methy_relation_V2.pdf", height = 6, width = 8.2)
p1
dev.off()


##### venn plot #####
load("process_data/venn_plot.RData")
library(ggvenn)
dat <- list(
    "hypermethy" = hyper_gene,
    "hypomethy" = hypo_gene,
    "DEGs" = DEG_gene
)
p1 <- ggvenn(dat,
    fill_color = c("#FC8A61", "#67C2A3", "#8EA0C9"),
    fill_alpha = 0.8,
    stroke_color = "white",
    stroke_size = 0.5,
    show_percentage = F,
    text_size = 8,
    set_name_size = 8
)
p1
pdf(file = "plot_out/venn_plot_V2.pdf", height = 6, width = 6)
p1
dev.off()

