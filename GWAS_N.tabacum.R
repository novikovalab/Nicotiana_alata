rm(list=ls())

library(circlize)
library(RColorBrewer)

# Read data
RS_BS62 <- read.table('./BSA6.result', header = TRUE, sep = "\t")
RS_BS52 <- read.table('./BSA5.result', header = TRUE, sep = "\t")
RS_BS42 <- read.table('./BSA4.result', header = TRUE, sep = "\t")
RS_BS12 <- read.table('./BSA1.result', header = TRUE, sep = "\t")

# Read synteny data
file_path <- "D:/Google download/S_T.link.txt"
link_data <- read.table(file_path, header=FALSE, col.names=c("Chr1", "Start1", "End1", "Chr2", "Start2", "End2"))

# Create chromosome information
chr_info <- data.frame(
  Chr = c("Nt01", "Nt02", "Nt03", "Nt04", "Nt05", "Nt06", "Nt07", "Nt08", "Nt09", "Nt10", "Nt11", "Nt12", 
          "Nt13", "Nt14", "Nt15", "Nt16", "Nt17", "Nt18", "Nt19", "Nt20", "Nt21", "Nt22", "Nt23", "Nt24"),
  Source = c("S", "T", "S", "T", "S", "S", "S", "S", "T", "S", "S", "T", 
             "T", "T", "T", "S", "T", "S", "T", "S", "S", "T", "T", "T")
)

# Correct chromosome labels by appending source
RS_BS62$Chr <- ifelse(!is.na(match(RS_BS62$Chr, chr_info$Chr)), 
                      paste0(RS_BS62$Chr, "_", chr_info$Source[match(RS_BS62$Chr, chr_info$Chr)]), 
                      RS_BS62$Chr)
RS_BS52$Chr <- ifelse(!is.na(match(RS_BS52$Chr, chr_info$Chr)), 
                      paste0(RS_BS52$Chr, "_", chr_info$Source[match(RS_BS52$Chr, chr_info$Chr)]), 
                      RS_BS52$Chr)
RS_BS42$Chr <- ifelse(!is.na(match(RS_BS42$Chr, chr_info$Chr)), 
                      paste0(RS_BS42$Chr, "_", chr_info$Source[match(RS_BS42$Chr, chr_info$Chr)]), 
                      RS_BS42$Chr)
RS_BS12$Chr <- ifelse(!is.na(match(RS_BS12$Chr, chr_info$Chr)), 
                      paste0(RS_BS12$Chr, "_", chr_info$Source[match(RS_BS12$Chr, chr_info$Chr)]), 
                      RS_BS12$Chr)

# Process chromosome labels in synteny data
link_data$Chr1 <- ifelse(!is.na(match(link_data$Chr1, chr_info$Chr)), 
                         paste0(link_data$Chr1, "_", chr_info$Source[match(link_data$Chr1, chr_info$Chr)]), 
                         link_data$Chr1)
link_data$Chr2 <- ifelse(!is.na(match(link_data$Chr2, chr_info$Chr)), 
                         paste0(link_data$Chr2, "_", chr_info$Source[match(link_data$Chr2, chr_info$Chr)]), 
                         link_data$Chr2)

# Calculate chromosome lengths: determine based on the maximum end position in GWAS and synteny data
max_lengths <- function(df, chr_col, pos_col) {
  aggregate(df[[pos_col]] ~ df[[chr_col]], FUN = max)
}

chr_length_S <- max_lengths(RS_BS62, "Chr", "Pos")
chr_length_T <- max_lengths(link_data, "Chr2", "End2")

# Ensure column names match
names(chr_length_S) <- c("Chr", "Length")
names(chr_length_T) <- c("Chr", "Length")

# Remove duplicate chromosome labels and keep the maximum length
chr_length_S <- aggregate(Length ~ Chr, data = chr_length_S, FUN = max)
chr_length_T <- aggregate(Length ~ Chr, data = chr_length_T, FUN = max)

# Merge chromosome lengths
chr_length <- unique(rbind(chr_length_S, chr_length_T))

# Check again and remove duplicates
chr_length <- aggregate(Length ~ Chr, data = chr_length, FUN = max)

# Rearrange chromosome order with S on the left and T on the right
chr_length <- chr_length[order(factor(chr_length$Chr, levels = c(
  paste0("Nt", sprintf("%02d", 24:13), "_S"),
  paste0("Nt", sprintf("%02d", 12:1), "_S"),
  paste0("Nt", sprintf("%02d", 1:12), "_T"),
  paste0("Nt", sprintf("%02d", 13:24), "_T")
))), ]

# Separate chromosomes by S and T order
colors_S <- rep("yellow", sum(grepl("_S$", chr_length$Chr)))
colors_T <- rep("lightcoral", sum(grepl("_T$", chr_length$Chr)))

# Combine colors
colors <- setNames(c(colors_S, colors_T), chr_length$Chr)

# Generate enough colors
link_colors <- colorRampPalette(brewer.pal(8, "Set1"))(nrow(chr_length))

# Create color map
color_map <- setNames(link_colors, chr_length$Chr)

# Initialize circos plot
png("circos_plot_symmetric_updated2.png", width = 400, height = 400, res = 600)
#cairo_svg("circos_plot_symmetric_updated1.svg", width = 800, height = 800)
circos.clear()
circos.par("start.degree" = 90, "gap.degree" = 2)
circos.initialize(factors = factor(chr_length$Chr, levels = chr_length$Chr), xlim = cbind(rep(0, nrow(chr_length)), chr_length$Length))

# Plot chromosomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xcenter = get.cell.meta.data("xcenter")
  circos.text(xcenter, 0.5, sector.index, facing = "bending.inside", niceFacing = TRUE, cex = 0.6)
  circos.axis(h = "top", labels = FALSE)
}, bg.col = colors[chr_length$Chr], bg.border = NA, track.height = 0.1)

# Plot four GWAS results and add vertical labels
gwas_data <- list(RS_BS62, RS_BS52, RS_BS42, RS_BS12)
track_colors <- c("blue", "red", "brown", "purple")
track_height <- 0.1
track_spacing <- 0.05
track_labels <- c("R:S", "DEF:def", "DEF1:def", "DEF2:def")

for (i in 1:4) {
  bed <- data.frame(chr = gwas_data[[i]]$Chr,
                    start = gwas_data[[i]]$Pos, 
                    end = gwas_data[[i]]$Pos + 1, 
                    value = gwas_data[[i]]$ED4)
  
  # Sort data
  bed <- bed[order(bed$chr, bed$start), ]
  
  # Plot track
  circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, col = track_colors[i], pch = 20, cex = 0.2)
  }, ylim = c(min(bed$value, na.rm = TRUE), max(bed$value, na.rm = TRUE)), track.height = track_height, bg.border = NA, track.margin = c(track_spacing, 0))
  
  # Add track name with matching color, labels placed vertically
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + 0.5, labels = track_labels[i], cex = 0.8, facing = "inside", niceFacing = TRUE, col = track_colors[i])
}

# Add synteny links (enhanced color contrast)
for (i in 1:nrow(link_data)) {
  col = adjustcolor(color_map[as.character(link_data$Chr2[i])], alpha.f = 0.7)
  circos.link(sector.index1 = as.character(link_data$Chr1[i]), point1 = c(link_data$Start1[i], link_data$End1[i]),
              sector.index2 = as.character(link_data$Chr2[i]), point2 = c(link_data$Start2[i], link_data$End2[i]),
              col = col, border = col)
}

# Close the device to display the plot
dev.off()
