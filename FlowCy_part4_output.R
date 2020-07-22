rm(list = ls())

# Load packages
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
# library(cytofkit)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(Rphenograph)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory
# Define workingDirectory
wdName <- "Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")

setwd(workingDirectory)

# Load workspace and SCEobject
load("workspaceSCEDA.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "mean")

# Save .fcs per cluster
outputDirectory <- getwd()
outputDirectory <- paste(outputDirectory, "output", sep = "/")
dir.create(outputDirectory)
setwd(outputDirectory)
sce$Clusters <- cluster_ids(sce, "meta12")
# keep_dr = TRUE not all cells have DR
flowSet <- sce2fcs(sce, split_by = "Clusters", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(flowSet, outdir = outputDirectory, filename = "bycluster")
merged <- sce2fcs(sce, split_by = NULL, keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.FCS(merged, filename = "merged.fcs")

# Save final workspace
save(list = ls(), file = "workspaceFinal.rds")
# Make plots
getwd()
plotFolder <- paste(outputDirectory, "plots", sep = "/")
dir.create(plotFolder)
setwd(plotFolder)

display.brewer.all(colorblindFriendly = TRUE)
png(filename = "colorblindFriendly.png", bg = "white")
display.brewer.all(colorblindFriendly = TRUE)
dev.off()
png(filename = "brewerAll.png", bg = "white")
display.brewer.all(colorblindFriendly = FALSE)
dev.off()
# MDS plot
tiff(filename = "MDSplot.tiff", compression = "lzw", bg = "white")
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition")
dev.off()

# boxplot abundances per cluster
tiff(filename = "boxplot_clusters", compression = "lzw", bg = "white")
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
dev.off()
# Expression Heatmap - clusters

tiff(filename = "ExpressionHeatmap.tiff", compression = "lzw", bg = "white")
plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id")
dev.off()

# UMAP color_by = Clusters
tiff(filename = "UMAP_clusters_condition.tiff", compression = "lzw", bg = "white")
plotDR(sce, dr = "UMAP", color_by = "Clusters", facet_by = "condition") + scale_color_brewer(palette = "Dark2")
dev.off()
png(filename = "UMAP_clusters_conditionSet2.png", bg = "white")
plotDR(sce, dr = "UMAP", color_by = "Clusters", facet_by = "condition") + scale_color_brewer(palette = "Set2")
dev.off()
