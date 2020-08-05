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
load("workspaceSCE.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "mean")
# Run FlowSOM and ConsensusClusterPlus
seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 2500
n_cells < n_events
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
sce <- runDR(sce, dr = "TSNE", cells = n_cells, features = "type", theta = 0.5, max_iter = 1000, 
             distMethod = "euclidean",
             PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")
sce <- runDR(sce, dr = "DiffusionMap", cells = n_cells, features = "type", assay = "exprs")

save(list = ls(), file = "workspaceSCEclusterDR.rds")

# Plots
plotAbundances(sce, k = "meta20", by = "cluster_id", group_by = "condition")
plotExprHeatmap(sce, features = type_markers(sce), k = "meta20", by = "cluster_id", fun = "mean", bars = TRUE, perc = TRUE, 
                scale = "last")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta8",  fun = "mean", scale = "last")

plotDR(sce, dr = "UMAP", color_by = "condition") + scale_color_manual(values = c("blue", "red"))
plotDR(sce, dr = "UMAP", color_by = "condition")
# UMAP plot color_by = "meta8", facet_by = "condition"
plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "condition")

plot <- plotDR(sce, dr = "UMAP", color_by = "meta12", facet_by = "condition")
plot$facet$params$ncol <- 2
plot
plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "condition")
plotDR(sce, dr = "TSNE", color_by = "condition", facet_by = "condition")
plotDR(sce, dr = "DiffusionMap", color_by = "condition") + scale_color_manual(values = c("blue", "red"))
markers <- c("CD28", "TCF1", "PD1", "Ki67", "Tbet", "CD69", "CD127")
plotDR(sce, dr = "DiffusionMap", color_by = markers)
plotDR(sce, dr = "UMAP", color_by = markers)
plotDR(sce, dr = "TSNE", color_by = markers)
# UMAP plot color_by = "meta8", facet_by = "sample_id"
plot <- plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "condition")
plot$facet$params$ncol <- 2
plot
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
plot <- plotDR(sce, dr = "DiffusionMap", color_by = "meta8", facet_by = "sample_id")
plot$facet$params$ncol <- 3
plot
# UMAP color_by CD27 DNAM1

plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "sample_id")
