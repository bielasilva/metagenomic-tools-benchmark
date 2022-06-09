# Figura 1 (Parcial) - Facilidade de setup e uso, tempo de execução, uso de RAM e tamanho da database
library(ggplot2)
library(hrbrthemes)
library(viridis)

setwd("~/Documents/TCC/resultados/")
size_dbDf <- read.csv("size_db.csv")  # Data collected using command "du"

size_dbDf$values <- round(size_dbDf$values/(1024*1024), digits=2)



## Run time and RAM usage
time_ramDf <- read.csv("stats_time_ram.csv") # Data collected using GNU Time
time_ramDf$ram <- round(time_ramDf$ram/(1024*1024), digits = 2)

time_ramDf$tempo <- sapply(strsplit(time_ramDf$tempo,":"),
                           function(x) sum(as.numeric(x)*c(1,1/60,1/3600)))

time_ramDf$tempo <- round(time_ramDf$tempo, digits = 2)

# Heatmap Max RAM
ggplot(time_ramDf, aes(x = programa, y = etapa, fill = ram)) +
    geom_tile() +
    geom_text(aes(label=ram), vjust=-1, color="black", size=3.5) +
    scale_fill_viridis(option="viridis", direction = -1) +
    theme_light() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
    labs(title = "RAM")

# Heatmap Execution Time
ggplot(time_ramDf, aes(x = programa, y = etapa, fill = tempo)) +
    geom_tile() +
    geom_text(aes(label=tempo), vjust=-1, color="black", size=3.5) +
    scale_fill_viridis(option="viridis", direction = -1) +
    theme_light() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
    labs(title = "Tempo")

# Database size
ggplot(size_dbDf, aes(x=programa, y=values, color=condition)) +
    geom_point(size=3) + 
    geom_segment(aes(x=programa, xend=programa, y=0, yend=values), color="black") +
    geom_text(aes(label=values), vjust=-1, color="black", size=3.5) +
    coord_flip() +
    theme_light() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank())

# Preparation for Figures 2 and 3
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggpubr)

df_creator <- function(file) {
    # Import results
    results <- read.csv(file)
    # Create empty data frame
    dfPlot <- data.frame(
        program = character(),
        sample = character(),
        specificity = character(),
        precision = character()
    )
    
    # Calculate the specificity and precision
    for (number in row.names(results)) {
        program <- results[number, "program"]
        sample <- results[number, "sample"]
        correct <- results[number, "correct"]
        incorrect <- results[number, "incorrect"]
        unclassified <- results[number, "unclassified"]
        specificity <-
            (correct / (correct + incorrect + unclassified)) * 100
        precision <- (correct / (correct + incorrect)) * 100
        new_row <- data.frame(program, sample, specificity, precision)
        dfPlot <- rbind(dfPlot, new_row)
    }
    
    #Calculate the specificity and precision means
    for (program in unique(dfPlot$program)) {
        specificity <- mean(dfPlot[dfPlot$program == program, ]$specificity)
        precision <- mean(dfPlot[dfPlot$program == program, ]$precision)
        sample <- "mean"
        new_row <- data.frame(program, sample, specificity, precision)
        dfPlot <- rbind(dfPlot, new_row)
    }
    return(dfPlot)
}

setwd("~/Documents/TCC/resultados/stats/") 
speciesDf <- df_creator("statistics_result_species.csv") # Output from correct_sort.py
genusDf <- df_creator("statistics_result_genus.csv") # Output from correct_sort.py
familyDf <- df_creator("statistics_result_family.csv") # Output from correct_sort.py
orderDf <- df_creator("statistics_result_order.csv") # Output from correct_sort.py

diversityDf <- read.csv("count_diversity_cami.csv") # Output from diversity_count.py

samples <- c("s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9")

for (sample in samples) {
    speciesDf$diversity[speciesDf$sample == sample] <- diversityDf$species[diversityDf$sample == sample]
    genusDf$diversity[genusDf$sample == sample] <- diversityDf$genus[diversityDf$sample == sample]
    familyDf$diversity[familyDf$sample == sample] <- diversityDf$family[diversityDf$sample == sample]
    orderDf$diversity[orderDf$sample == sample] <- diversityDf$order[diversityDf$sample == sample]
}

speciesDf$taxonomy <- "Espécie"
genusDf$taxonomy <- "Gênero"
familyDf$taxonomy <- "Família"
orderDf$taxonomy <- "Ordem"

dfPlot <- rbind(speciesDf, genusDf, familyDf, orderDf)
dfPlot$taxonomy <- factor(dfPlot$taxonomy, levels = c("Espécie", "Gênero", "Família", "Ordem"))

dfPlot <- dfPlot[!dfPlot$sample == "mean", ]
dfPlot_no_means <-  dfPlot[!dfPlot$sample == "mean", ]

# Figura 2 - Boxplots Precisão e Especificidade
precision_plot <- ggplot(dfPlot_no_means, aes(x = taxonomy, y = precision, fill = program)) +
    geom_boxplot(position=position_dodge(.85)) +
    labs(title = "Precisão", y = "(%)", x = "Taxonomia", fill = "Program") +
    stat_summary(fun=mean, geom="point", aes(group=program), shape=23, fill="red4", size=4, color="white", position=position_dodge(.85)) +
    ylim(7,100) +
    theme_ipsum(
        base_family = "Arial",
        axis_title_just = "cc",
        axis_title_size = 10,
        base_size = 11,
        plot_title_size = 14,
        plot_title_margin = 3,
        plot_margin = margin(5, 5, 5, 5)
    ) +
    theme(plot.title = element_text(hjust = 0.5))

specificity_plot <- ggplot(dfPlot_no_means, aes(x = taxonomy, y = specificity, fill = program)) +
    geom_boxplot(position=position_dodge(.85)) +
    ylim(7,100) +
    labs(title = "Especificidade", y = "(%)", x = "Taxonomia", fill = "Programa") +
    stat_summary(fun=mean, geom="point", aes(group=program), shape=23, fill="red4", size=4, color="white", position=position_dodge(.85)) +
    theme_ipsum(
        base_family = "Arial",
        axis_title_just = "cc",
        axis_title_size = 10,
        base_size = 11,
        plot_title_size = 14,
        plot_title_margin = 3,
        plot_margin = margin(5, 5, 5, 5)
    ) +
    theme(plot.title = element_text(hjust = 0.5))

ggarrange(specificity_plot, precision_plot,
          nrow = 2, common.legend = TRUE)

# Figura 3 - Quantidade de espécies por amostra com Precisão e Especificidade
scater_diversity <- function(dfPlot, title) {
    # Plot the graph
    plot <- ggplot(dfPlot, aes(x = specificity, y = precision, color = diversity, shape = program)) +
        scale_shape_manual(values = c(15, 16, 17, 18)) +
        scale_color_viridis(option="viridis", direction = -1) +
        geom_point(size = 4) +
        xlim(7, 100) +
        ylim(15, 100) +
        theme_ipsum(
            base_family = "Arial",
            axis_title_just = "cc",
            axis_title_size = 10,
            base_size = 12,
            plot_title_size = 14,
            plot_title_margin = 3,
            plot_margin = margin(5, 5, 5, 5)
        ) +
        labs(
            title = title,
            y = "Precisão (%)",
            x = "Especificidade (%)",
            color = "Diversidade da amostra",
            shape = "Programa"
        ) +
        theme(plot.title = element_text(hjust = 0.5))
    return(plot)
}

ggarrange(
    scater_diversity(dfPlot[dfPlot$taxonomy == "Espécie", ], "Espécie"),
    scater_diversity(dfPlot[dfPlot$taxonomy == "Gênero", ], "Gênero"),
    scater_diversity(dfPlot[dfPlot$taxonomy == "Família", ], "Família"),
    scater_diversity(dfPlot[dfPlot$taxonomy == "Ordem", ], "Ordem"),
    nrow = 1,
    common.legend = F
)

# Figura 4
## UpsetPlot all CAMI
library(ggplot2)
library(svglite)
library(ComplexUpset)
setwd("~/Documents/TCC/resultados/upset_plots/")

for (tax_level in c("species", "genus", "family", "order")) {
  df_list = list()
  for (i in (0:9)) {
    file <- paste("cami/s", i, "/taxa_matrix_", tax_level, ".csv", sep = "")
    df <- as.data.frame(t(read.csv(file, row.names = 1)))
    df_list <- c(list(df), df_list)
  }
  programs <- c("cami", "centrifuge", "kaiju", "kraken2", "metacache")
  
  df_t <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(df_t) <- programs
  
  for (df in df_list) {
    for (taxon in row.names(df)) {
      if (any(taxon == row.names(df_t))) {
        df_t[taxon,] <- (((df_t[taxon,] + df[taxon,]) > 0) + 0)
      } else {
        df_t <- rbind(df_t, df[taxon,])
      }
    }
  }
  
  
  upset(df_t, programs, name = "Programs",
        #min_size = 7,
        #width_ratio = 0.25,
        #intersections='all',
        set_sizes = (
          upset_set_size()
          + geom_text(aes(label=..count..), hjust=-0.05, stat = 'count')
          + theme(axis.ticks.x = element_line())
        ) + ggtitle(paste("CAMI", tax_level)))
  
  ggsave(paste("CAMI", tax_level, ".svg"), path = "plots/", width = 40, height = 16, units = "cm", dpi = 500)
}

## UpsetPlot all Real
library(ggplot2)
library(svglite)
library(ComplexUpset)
setwd("~/Documents/TCC/resultados/upset_plots/")

for (tax_level in c("species", "genus", "family", "order")) {
  df_list = list()
  for (i in (37:48)) {
    file <- paste("real/S", i, "/taxa_matrix_", tax_level, ".csv", sep = "")
    df <- as.data.frame(t(read.csv(file, row.names = 1)))
    df_list <- c(list(df), df_list)
  }
  programs <- c("centrifuge", "kaiju", "kraken2", "metacache")
  
  df_t <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_t) <- programs
  
  for (df in df_list) {
    for (taxon in row.names(df)) {
      if (any(taxon == row.names(df_t))) {
        df_t[taxon,] <- (((df_t[taxon,] + df[taxon,]) > 0) + 0)
      } else {
        df_t <- rbind(df_t, df[taxon,])
      }
    }
  }
  
  
  upset(df_t, programs, name = "Programs",
        #min_size = 7,
        #width_ratio = 0.25,
        #intersections='all',
        set_sizes = (
          upset_set_size()
          + geom_text(aes(label=..count..), hjust=-0.05, stat = 'count')
          + theme(axis.ticks.x = element_line())
        ) + ggtitle(paste("Real", tax_level)))
  
  ggsave(paste("Real", tax_level, ".svg"), path = "plots/", width = 40, height = 16, units = "cm", dpi = 500)
}

## Heatmap CAMI
library(ggplot2)
library(viridis)
library(ggpubr)

df_heatmap <- as.data.frame(read.csv("upset_matrix_heatmap_cami.csv"))

level_order_x <- c("cen_kai_kra_met",	"cen_kai_kra",	"cam",	"cen",	"cen_kai_kra_met_cam",	"cen_kai",	"cen_kra",	"kai",	"cen_kai_met",	"cen_met",	"cen_kra_met")
level_order_y <- c("species", "genus", "family", "order")

species_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Species", fill = species)) +
  geom_tile() +
  geom_text(aes(label=species), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

genus_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Genus", fill = genus)) +
  geom_tile() +
  geom_text(aes(label=genus), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

family_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Family", fill = family)) +
  geom_tile() +
  geom_text(aes(label=family), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

order_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Order", fill = order)) +
  geom_tile() +
  geom_text(aes(label=order), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

ggarrange(species_plt, genus_plt, family_plt, order_plt,
          ncol = 1,
          common.legend = F)

## Heatmap real
library(ggplot2)
library(viridis)
library(ggpubr)

df_heatmap <- as.data.frame(read.csv("upset_matrix_heatmap_real.csv"))

level_order_x <- c("cen_kai_kra_met", "cen_kai_kra", "cen", "cen_kai", "kai", "cen_kai_met",
                   "cen_kra", "cen_met", "cen_kra_met", "met", "kai_met", "kai_kra")
level_order_y <- c("species", "genus", "family", "order")

species_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Species", fill = species)) +
  geom_tile() +
  geom_text(aes(label=species), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

genus_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Genus", fill = genus)) +
  geom_tile() +
  geom_text(aes(label=genus), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

family_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Family", fill = family)) +
  geom_tile() +
  geom_text(aes(label=family), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

order_plt <- ggplot(df_heatmap, aes(x = factor(programs, levels = level_order_x), y = "Order", fill = order)) +
  geom_tile() +
  geom_text(aes(label=order), vjust=-1, color="black", size=3.5) +
  scale_fill_viridis(option="viridis", direction = -1)

ggarrange(species_plt, genus_plt, family_plt, order_plt,
          ncol = 1,
          common.legend = F)

## Bar plot total
library(ggplot2)
library(ggbreak)

df_bar <- as.data.frame(read.csv("upset_matrix_bar.csv"))

level_order_y <- c("species", "genus", "family", "order")

ggplot(df_bar[df_bar$dataset == "cami", ], aes(x = value, y = programs, fill = tax_level)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_x_reverse(
    limits = c(6500, 0), 
    breaks = seq(6500, 0, by = -100)) +
  scale_x_break(c(3800, 6100), scales = "fixed") +
  scale_x_break(c(1700, 3500), scales = "fixed") +
  scale_x_break(c(450, 1250), scales = "fixed")


ggplot(df_bar[df_bar$dataset == "real", ], aes(x = value, y = programs, fill = tax_level)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_x_reverse(
    limits = c(6500, 0), 
    breaks = seq(6500, 0, by = -200)) +
  scale_x_break(c(4100, 6100), scales = "fixed") + 
  scale_x_break(c(1700, 3700), scales = "fixed")