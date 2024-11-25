
library(ggplot2)
library(ggalluvial)
library(dplyr)

l3_colors <- list(
  'ATL' = '#ffff00',
  'B' = '#1ce6ff',
  'C-TAL' = '#ff34ff',
  'PC' = '#aec7e8',
  'CCD-IC-A' = '#ff4a46',
  'CCD-PC' = '#008941',
  'CNT' = '#006fa6',
  'CNT-IC-A' = '#a30059',
  'CNT-PC' = '#ffdbe5',
  'DCT1' = '#7a4900',
  'DCT2' = '#0000a6',
  'DTL' = '#63ffac',
  'EC-AEA' = '#b79762',
  'EC-AVR' = '#004d43',
  'EC-DVR' = '#8fb0ff',
  'EC-GC' = '#997d87',
  'EC-LYM' = '#5a0007',
  'EC-PTC' = '#809693',
  'FIB' = '#6a3a4c',
  'IC-B' = '#1b4400',
  'IMCD' = '#4fc601',
  'M-FIB' = '#3b5dff',
  'M-TAL' = '#4a3b53',
  'MAC-M2' = '#ff2f80',
  'MAST' = '#61615a',
  'MC' = '#ba0900',
  'MD' = '#6b7900',
  'MDC' = '#00c2a0',
  'MYOF' = '#ffaa92',
  'N' = '#ff90c9',
  'NEU' = '#b903aa',
  'NKC/T' = '#d16100',
  'OMCD-IC-A' = '#ddefff',
  'OMCD-PC' = '#000035',
  'PEC' = '#7b4f4b',
  'PL' = '#a1c299',
  'POD' = '#300018',
  'PT-S1/2' = '#0aa6d8',
  'PT-S3' = '#013349',
  'PapE' = '#00846f',
  'REN' = '#372101',
  'T' = '#ffb500',
  'VSMC' = '#c2ffed',
  'VSMC/P' = '#a079bf',
  'aFIB' = '#cc0744',
  'aPT' = '#c0b9b2',
  'aTAL1' = '#c2ff99',
  'aTAL2' = '#001e09',
  'cDC' = '#00489c',
  'cycCNT' = '#6f0062',
  'cycDCT' = '#0cbd66',
  'cycEC' = '#eec3ff',
  'cycMNP' = '#456d75',
  'cycMYOF' = '#b77b68',
  'cycNKC/T' = '#7a87a1',
  'cycPT' = '#788d66',
  'dC-IC-A' = '#885578',
  'dC-TAL' = '#fad09f',
  'dCNT' = '#ff8a9a',
  'dDCT' = '#d157a0',
  'dEC' = '#bec459',
  'dEC-PTC' = '#456648',
  'dFIB' = '#0086ed',
  'dIMCD' = '#886f4c',
  'dM-FIB' = '#34362d',
  'dM-TAL' = '#b4a8bd',
  'dOMCD-PC' = '#00a6aa',
  'dPT' = '#452c2c',
  'dVSMC' = '#636375',
  'ncMON' = '#a3c8c9',
  'pDC' = '#ff913f',
  'tPC-IC' = '#938a81',
  'Unclassified' = '#d3d3d3'
)
l3_colors_un <- unlist(l3_colors)

l1_colors <- list(
  'ATL' = '#1f77b4',
  'PT_VCAM1' = '#c5b0d5',
  'CNT' = '#ff7f0e',
  'DCT' = '#279e68',
  'DCT1' = '#279e68',
  'DCT2' = '#ffb500',
  'DTL' = '#d62728',
  'EC' = '#aa40fc',
  'ENDO' = '#aa40fc',
  'FIB' = '#8c564b',
  'IC' = '#e377c2',
  'ICA' = '#e377c2',
  'ICB' = '#7b4f4b',
  'IMM' = '#b5bd61',
  'LEUK' = '#b5bd61',
  'NEU' = '#17becf',
  'PC' = '#aec7e8',
  'PEC' = '#17becf',
  'POD' = '#98df8a',
  'PODO' = '#98df8a',
  'PT' = '#ff9896',
  'PapE' = '#c5b0d5',
  'TAL' = '#c49c94',
  'MES' = '#f7b6d2',
  'VSM/P' = '#f7b6d2',
  'Unclassified' = '#d3d3d3',
  'scRNA' = '#1f77b4',
  'scRNA5p' = '#ff7f0e',
  'snRNA' = '#2ca02c'
)
l1_colors_un <- unlist(l1_colors)

techs_color <- list(
  'scRNA' = '#1f77b4',
  'scRNA5p' = '#ff7f0e',
  'snRNA' = '#2ca02c'
)

l3_colors_un <- unlist(l3_colors)
# Set working directory
setwd('/home/macera/Documentos/CZI/MANUSCRIPT_PREP/FIGURES')

# Read the data from CSV
df <- read.csv("objects/alluvial_data.csv")
df <- df %>% select('X', 'Deepscore_HCA_l1_Clean', 'Deepscore_HCA_l3_Clean')
colnames(df) <- c("X", "L1", "L2")
df$X <- as.factor(df$X)
df$L1 <- as.factor(df$L1)
df$L2 <- as.factor(df$L2)

# Calculate frequency
freq <- df %>% count(X, L1, L2)
colnames(freq) <- c("Methods", "L1", "L2", "freq")

# Check alluvia form
is_alluvia_form(as.data.frame(freq), axes = 1:3, silent = TRUE)

# Prepare color mapping
l2_levels <- unique(freq$L2)
filtered_colors <- l3_colors[names(l3_colors) %in% l2_levels]
filtered_colors <- unlist(filtered_colors)


# Remove rows where L1 is "NEU"
freq <- freq %>% filter(L1 != "NEU")



ggplot(as.data.frame(freq), aes(y = freq, axis1 = Methods, axis2 = L1, axis3 = L2)) +
  geom_alluvium(aes(fill = L1), width = 5/12) +
  geom_stratum(aes(fill = L1), width = c(rep(5/12, 3), rep(4/12, 12), rep(6/12, 39)), color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6, label.size = 0.5, 
             nudge_x = c(rep(0, 3), rep(c(-0.08, 0.08), 13/2), rep(c(-0.25, -0.125,0, 0.125, 0.25), 40/2))) +
  scale_x_discrete(limits = c("Methods", "L1", "L2"), expand = c(.15, .15)) +
  scale_fill_manual(values = unlist(l1_colors), guide = guide_legend(title = "Layer 1 Annotation")) +
  ggtitle("Frequency Alluvial") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background without border
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.border = element_blank(),  # No panel border
    axis.line = element_blank(),  # No axis lines
    panel.grid.major = element_blank(),  # No major grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  scale_y_continuous(breaks = seq(0, 2000, by = 2500)) +
  labs(y = "# Cells", title = "")

theme(plot.margin = unit(c(1,1,1,1), "lines"))


ggsave("../FIGURES/alluvial_plot_v2.png",dpi = 600)
ggsave("../FIGURES/alluvial_plot_v2.pdf", dpi = 600, device = "pdf")


ggplot(as.data.frame(freq), aes(y = freq, axis1 = Methods, axis2 = L1, axis3 = L2)) +
  geom_alluvium(aes(fill = L2), width = 5/12) +
  geom_stratum(aes(fill = L2), width = c(rep(5/12, 3), rep(4/12, 12), rep(6/12, 39)), color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6, label.size = 0.5, 
             nudge_x = c(rep(0, 3), rep(c(-0.08, 0.08), 13/2), rep(c(-0.25, -0.125,0, 0.125, 0.25), 40/2))) +
  scale_x_discrete(limits = c("Methods", "L1", "L2"), expand = c(.15, .15)) +
  scale_fill_manual(values = unlist(l3_colors), guide = guide_legend(title = "Layer 1 Annotation")) +
  ggtitle("Frequency Alluvial") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background without border
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.border = element_blank(),  # No panel border
    axis.line = element_blank(),  # No axis lines
    panel.grid.major = element_blank(),  # No major grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  scale_y_continuous(breaks = seq(0, 2000, by = 2500)) +
  labs(y = "# Cells", title = "")

theme(plot.margin = unit(c(1,1,1,1), "lines"))

ggsave("../FIGURES/figures/alluvial/alluvial_L3.png",dpi = 600)
ggsave("../FIGURES/figures/alluvial/alluvial_L3.pdf", dpi = 600, device = "pdf")




ggplot(as.data.frame(freq), aes(y = freq, axis1 = Methods, axis2 = L1, axis3 = L2)) +
  geom_alluvium(aes(fill = L1), width = 5/12) +
  geom_stratum(aes(fill = L1), width = c(rep(5/12, 3), rep(4/12, 12), rep(6/12, 39)), color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 4, label.size = 0.5, 
             nudge_x = c(rep(0, 3), rep(c(-0.08, 0.08), 13/2), rep(c(-0.25, -0.125,0, 0.125, 0.25), 40/2))) +
  scale_x_discrete(limits = c("Methods", "L1", "L2"), expand = c(.15, .15)) +
  scale_fill_manual(values = unlist(l1_colors), guide = guide_legend(title = "Layer 1 Annotation")) +
  ggtitle("Frequency Alluvial") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background without border
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.border = element_blank(),  # No panel border
    axis.line = element_blank(),  # No axis lines
    panel.grid.major = element_blank(),  # No major grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  scale_y_continuous(breaks = seq(0, 2000, by = 2500)) +
  labs(y = "# Cells", title = "")

theme(plot.margin = unit(c(1,1,1,1), "lines"))

ggsave("../FIGURES/figures/alluvial/alluvial_L1_caps.png",dpi = 600)
ggsave("../FIGURES/figures/alluvial/alluvial_L1_caps.pdf", dpi = 600, device = "pdf")




ggplot(as.data.frame(freq), aes(y = freq, axis1 = Methods, axis2 = L1, axis3 = L2)) +
  geom_alluvium(aes(fill = Methods), width = 7/12) +
  geom_stratum(aes(fill = Methods), width = c(rep(7/12, 3), rep(4/12, 12), rep(6/12, 39)), color = "black") +
  geom_label(
    stat = "stratum", 
    aes(label = sprintf("%s\n%d cells", after_stat(stratum), after_stat(count))),  # Combine name and count
    size = 5, 
    label.size = 0.5, 
    nudge_x = c(rep(0, 3), rep(c(-0.08, 0.08), 13/2), rep(c(-0.25, -0.125,0, 0.125, 0.25), 40/2))
  ) +
  scale_x_discrete(limits = c("Methods", "L1", "L2"), expand = c(.15, .15)) +
  scale_fill_manual(values = unlist(techs_color), guide = guide_legend(title = "Layer 1 Annotation")) +
  ggtitle("Frequency Alluvial") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background without border
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.border = element_blank(),  # No panel border
    axis.line = element_blank(),  # No axis lines
    panel.grid.major = element_blank(),  # No major grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  scale_y_continuous(breaks = seq(0, 2000, by = 2500)) +
  labs(y = "# Cells", title = "")

theme(plot.margin = unit(c(1,1,1,1), "lines"))

ggsave("../FIGURES/figures/alluvial/alluvial_Methods_cells.png",dpi = 600)
ggsave("../FIGURES/figures/alluvial/alluvial_Methods_cells.pdf", dpi = 600, device = "pdf")












