# adapted from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html

#remotes::install_version("Seurat", version = "3.2.0")
library(Seurat)
library(ggplot2)
library(patchwork)
packageVersion("Seurat")

pbmc.rna <- readRDS(file = "~/DiseaseLabelsMethod/data/bmcite_rna_counts.RDS")
pbmc.meta <- readRDS(file = "~/DiseaseLabelsMethod/data/bmcite_metadata.RDS")

# Clustering
pbmc <- CreateSeuratObject(counts = pbmc.rna, meta.data = pbmc.meta)
pbmc <- NormalizeData(pbmc)
# choose ~1k variable features
pbmc <- FindVariableFeatures(pbmc)
# standard scaling (no regression)
pbmc <- ScaleData(pbmc)
# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
ElbowPlot(pbmc, ndims = 50)

#remove.packages("Matrix")
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
## restart r studio

pbmc <- FindNeighbors(pbmc, dims = 1:25)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunTSNE(pbmc, dims = 1:25, method = "FIt-SNE")

# Find the markers that define each cluster, and use these to annotate the clusters, we use
# max.cells.per.ident to speed up the process
pbmc.rna.markers <- FindAllMarkers(pbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)

# Hierarchical clustering
Idents(pbmc) = "celltype.l2"
#install.packages('ape')
pbmc <- BuildClusterTree(object = pbmc)
PlotClusterTree(pbmc, type = 'c')

head(pbmc@meta.data)
DimPlot(pbmc, label = TRUE, group.by = "celltype.l1") + NoLegend()
DimPlot(pbmc, label = TRUE, group.by = "celltype.l2") + NoLegend()
table(pbmc@meta.data$celltype.l2)
DimPlot(pbmc, label = TRUE, group.by = "celltype.l2")
pl <- DimPlot(pbmc, label = FALSE,
              cols = c(rep("white", 5), "lightgray", rep("white", 10), "#980C13", rep("white", 10)), 
              label.size = 0) #+ NoLegend()
LabelClusters(plot = pl, id='ident', color = "red")
pl
# redo this plot by subsetting all Naive B and all Memory B and doing their own dim red
pbmc_naiveB_memoryB <- subset(pbmc, idents = c("Naive B", "Memory B"))
DimPlot(pbmc_naiveB_memoryB, label = TRUE, group.by = "celltype.l2") + NoLegend()
# recompute the embedding
pbmc_naiveB_memoryB <- NormalizeData(pbmc_naiveB_memoryB)
pbmc_naiveB_memoryB <- FindVariableFeatures(pbmc_naiveB_memoryB)
pbmc_naiveB_memoryB <- ScaleData(pbmc_naiveB_memoryB)
pbmc_naiveB_memoryB <- RunPCA(pbmc_naiveB_memoryB)
#pbmc_naiveB_memoryB <- FindNeighbors(pbmc_naiveB_memoryB, dims = 1:25)
#pbmc_naiveB_memoryB <- FindClusters(pbmc_naiveB_memoryB, resolution = 0.8)
pbmc_naiveB_memoryB <- RunTSNE(pbmc_naiveB_memoryB, method = "FIt-SNE")
DimPlot(pbmc_naiveB_memoryB, label = TRUE, group.by = "celltype.l2") + NoLegend()
pl <- DimPlot(pbmc_naiveB_memoryB, label = TRUE,
              cols = c("lightgray", "#980C13"), 
              label.size = 0) #+ NoLegend()
#LabelClusters(plot = pl, id='ident', color = "black")
ggsave('~/DiseaseLabelsMethod/data/Figure2/tsne_pbmc_naiveB_memoryB.png', dpi='print')
#install.packages('svglite')
library(svglite)
ggsave('~/DiseaseLabelsMethod/data/Figure2/tsne_pbmc_naiveB_memoryB.svg')
pl

# install.packages("RColorBrewer")
library(RColorBrewer)
DimPlot(pbmc, cols = c(brewer.pal(n = 5, name = "Set2"), "gray", 
                       brewer.pal(n = 10, name = "Spectral"), "#980C13", 
                       brewer.pal(n = 10, name = "PiYG"))) #+ NoLegend()

#brewer.pal(n = 11, name = "Spectral")
#palette(brewer.pal(n = 8, name = "Set2"))

# Remove all mouse genes
mouse_vargenes_bool <- startsWith(VariableFeatures(pbmc), "MOUSE-")
VariableFeatures(pbmc) = VariableFeatures(pbmc)[!mouse_vargenes_bool]

# Find markers
pbmc.naive.vs.memory.markers <- FindMarkers(pbmc, ident.1 = "Naive B", ident.2 = "Memory B", test.use = "wilcox", group.by = "celltype.l2")
# TCL1A, IGHD, B2M, CLECL1
#> dim(pbmc.naive.vs.memory.markers)
#[1] 140   5
pbmc.memory.vs.naive.markers <- FindMarkers(pbmc, ident.1 = "Memory B", ident.2 = "Naive B", test.use = "wilcox", group.by = "celltype.l2")

# Mean var plot
#install.packages("matrixStats")
library(matrixStats)
geneVar = rowVars(as.matrix(pbmc@assays$RNA@counts))
geneMean = rowMeans(as.matrix(pbmc@assays$RNA@counts))
ggplot(data.frame(logmean = log(geneMean), logvariance = log(geneVar)), aes(logmean, logvariance)) + geom_point()
ggplot(data.frame(mean = geneMean, variance = geneVar), aes(mean, variance)) + geom_point()

# Subset the celltypes of interest
Idents(pbmc) = "celltype.l2"
pbmc_naiveB <- subset(pbmc, idents = c("Naive B"))
pbmc_memoryB <- subset(pbmc, idents = c("Memory B"))

num_minority_cells_vec = rep(NA, length(seq(from=2.5, to=50, by=2.5)))
i = 1
for (a in seq(from=2.5, to=50, by=2.5)) {
  num_minority_cells_vec[i] = round(a*1900/(100-a), digits = 0)
  i = i+1
}
num_minority_cells_vec
# 49  100  154  211  271  335  403  475  552  633  721  814  915 1023 1140 1267 1404 1555 1719 1900

set.seed(9)
memoryB_49 <- subset(pbmc_memoryB, downsample = 49) #2.5%
memoryB_100 <- subset(pbmc_memoryB, downsample = 100) #5%
memoryB_154 <- subset(pbmc_memoryB, downsample = 154) #7.5%
memoryB_211 <- subset(pbmc_memoryB, downsample = 211) #10%
memoryB_271 <- subset(pbmc_memoryB, downsample = 271) #12.5%
memoryB_335 <- subset(pbmc_memoryB, downsample = 335) #15%
memoryB_403 <- subset(pbmc_memoryB, downsample = 403) #17.5%
memoryB_475 <- subset(pbmc_memoryB, downsample = 475) #20%
memoryB_552 <- subset(pbmc_memoryB, downsample = 552) #22.5%
memoryB_633 <- subset(pbmc_memoryB, downsample = 633) #25%
memoryB_721 <- subset(pbmc_memoryB, downsample = 721) #27.5%
memoryB_814 <- subset(pbmc_memoryB, downsample = 814) #30%
memoryB_915 <- subset(pbmc_memoryB, downsample = 915) #32.5%
memoryB_1023 <- subset(pbmc_memoryB, downsample = 1023) #35%
memoryB_1140 <- subset(pbmc_memoryB, downsample = 1140) #37.5%
memoryB_1267 <- subset(pbmc_memoryB, downsample = 1267) #40%
memoryB_1404 <- subset(pbmc_memoryB, downsample = 1404) #42.5%
memoryB_1555 <- subset(pbmc_memoryB, downsample = 1555) #45%

#set.seed(9)
#cd4_memory_49 <- subset(pbmc_cd4_memory, downsample = 49) #2.5%

#==== Create datasets of interest ====
merge_majority_minority <- function(x, y, cell_ids, project = "mix") {
  set.seed(8)
  half_cellnames_by_ident = sample(x=Cells(x), size=ceiling((length(Cells(x))+length(Cells(y)))/2))
  x$batch = as.numeric(!(Cells(pbmc_naiveB) %in% unlist(half_cellnames_by_ident)))
  merged_s = merge(x, y, add.cell.ids = cell_ids, project = project)
  merged_s@meta.data$batch[is.na(merged_s@meta.data$batch)] = '1'
  return(merged_s)
}

# 1 # Majority naiveB, minority memoryB
naiveB_1900_memoryB_49 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_49, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_100 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_100, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_154 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_154, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_211 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_211, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_271 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_271, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_335 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_335, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_403 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_403, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_475 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_475, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_552 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_552, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_633 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_633, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_721 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_721, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_814 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_814, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_915 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_915, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_1023 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_1023, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_1140 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_1140, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_1267 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_1267, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_1404 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_1404, cell_ids = c("naiveB", "memoryB"))
naiveB_1900_memoryB_1555 <- merge_majority_minority(x=pbmc_naiveB, y=memoryB_1555, cell_ids = c("naiveB", "memoryB"))

table(naiveB_1900_memoryB_49@meta.data$batch)

##################################
#==== Biological perturbation ====
##################################
col_perturb = c("lightgray", "#980C13")
col_batch = c("#1B3EF3", "#EC8233") #c("#2E4ED7", "#D28239") #old c("#1B9E77", "#D95F02")
minority_celltype = 'Memory B'

cluster_s_list <- list(naiveB_1900_memoryB_49, naiveB_1900_memoryB_100, naiveB_1900_memoryB_154, naiveB_1900_memoryB_211, 
                       naiveB_1900_memoryB_271, naiveB_1900_memoryB_335, naiveB_1900_memoryB_403, naiveB_1900_memoryB_475, 
                       naiveB_1900_memoryB_552, naiveB_1900_memoryB_633, naiveB_1900_memoryB_721, naiveB_1900_memoryB_814, 
                       naiveB_1900_memoryB_915, naiveB_1900_memoryB_1023, naiveB_1900_memoryB_1140, naiveB_1900_memoryB_1267, 
                       naiveB_1900_memoryB_1404, naiveB_1900_memoryB_1555)
names(cluster_s_list) <- c("naiveB_1900_memoryB_49", "naiveB_1900_memoryB_100", "naiveB_1900_memoryB_154", "naiveB_1900_memoryB_211", 
                           "naiveB_1900_memoryB_271", "naiveB_1900_memoryB_335", "naiveB_1900_memoryB_403", "naiveB_1900_memoryB_475", 
                           "naiveB_1900_memoryB_552", "naiveB_1900_memoryB_633", "naiveB_1900_memoryB_721", "naiveB_1900_memoryB_814", 
                           "naiveB_1900_memoryB_915", "naiveB_1900_memoryB_1023", "naiveB_1900_memoryB_1140", "naiveB_1900_memoryB_1267", 
                           "naiveB_1900_memoryB_1404", "naiveB_1900_memoryB_1555")

npcs = 60
n_dims = 5 
fixed_npcs = 5 

violin_plot_data_pcs = data.frame()
variable_genes_list <- list()
for (i in seq_along(cluster_s_list)){
  cluster_s = cluster_s_list[[i]]
  cluster_ident = names(cluster_s_list)[i]
  print(cluster_ident)
  dir.create(paste0("data/", cluster_ident, "/"))
  cluster_dir_path = paste0("data/", cluster_ident, "/")
  cluster_s <- NormalizeData(cluster_s)
  all.genes <- rownames(cluster_s)
  cluster_s <- ScaleData(cluster_s, features = all.genes)
  cluster_s <- FindVariableFeatures(cluster_s)
  mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
  VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
  #VariableFeatures(cluster_s) = VariableFeatures(pbmc)
  cluster_s <- RunPCA(cluster_s, npcs=npcs)
  #==== Can we find the structure with PCA? ====
  DimPlot(cluster_s, reduction='pca', group.by = 'batch', cols = col_batch) + ggtitle(cluster_ident)
  ggsave(paste0(cluster_dir_path, "pca2_batch.png"))
  DimPlot(cluster_s, reduction='pca', group.by = 'celltype.l2', cols = c("red", "gray")) + ggtitle(cluster_ident)
  ggsave(paste0(cluster_dir_path, "pca2_label.png"))
  cluster_s <- RunTSNE(cluster_s, dims = 1:n_dims)
  DimPlot(cluster_s, reduction = "tsne", group.by = 'batch', cols = col_batch) + ggtitle(cluster_ident)
  ggsave(paste0(cluster_dir_path, "tsne5_batch.png"))
  DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("red", "gray")) + ggtitle(n_dims)
  ggsave(paste0(cluster_dir_path, "tsne5_label.png"))
  #==== Neighbors in PC space ==== 
  # track the mixing in the latent space, compute for the perturbed.minority cells only
  num_minority_cells = sum(cluster_s@meta.data$celltype.l2==minority_celltype)
  n_dims_vec = seq(from=2, to=npcs, by=1)
  violin_plot_data = data.frame()
  for(n_dims in n_dims_vec) {
    print(n_dims)
    cluster_s <- FindNeighbors(cluster_s, compute.SNN = TRUE, dims = 1:n_dims)
    tmp = apply(cluster_s@graphs[[1]][cluster_s@meta.data$celltype.l2==minority_celltype,], 1, 
                function(x) {cluster_s@meta.data$celltype.l2[x == 1]})
    # vector of length number of minoruty cells containing the proportion of their neighborst might minority label
    prop_minority_neighbors = apply(tmp, 2, function(x) {mean(x==minority_celltype)})
    numPCs = rep(n_dims, num_minority_cells)
    violin_plot_data = rbind(violin_plot_data, cbind(prop_minority_neighbors, numPCs))
    print(paste0('done with ', n_dims, 'PCs!'))
  }
  # Plots
  # 1. across numPCs for the same num_minority cells
  violin_plot_data$numPCs = as.factor(violin_plot_data$numPCs)
  p <- ggplot(violin_plot_data, aes(x=numPCs, y=prop_minority_neighbors)) + geom_boxplot(width=0.1, outlier.size = 0.1)#+ geom_violin()
  p = p + stat_summary(fun.y=median, geom="point", size=2, color='red')
  p = p #+ theme_classic() #+ geom_jitter(size=1, position=position_jitter(0.1), alpha=0.1)  
  p = p + ggtitle(paste0(num_minority_cells, " ", minority_celltype, " cells"))
  p = p + xlab('Num PCs') + ylab(paste0("Fraction of ", minority_celltype, " neighbors"))
  p = p + theme(axis.text.x = element_text(angle = 90, hjust=0.95, size=7))
  p
  ggsave(paste0(cluster_dir_path, "prop_neighbors_acrosspcs.png"))
  # 2. across num_minority cells, fixed numPCs
  cluster_s <- FindNeighbors(cluster_s, compute.SNN = TRUE, dims = 1:fixed_npcs)
  tmp = apply(cluster_s@graphs[[1]][cluster_s@meta.data$celltype.l2==minority_celltype,], 1, 
              function(x) {cluster_s@meta.data$celltype.l2[x == 1]})
  # vector of length number of minoruty cells containing the proportion of their neighborst with minority label
  prop_minority_neighbors = apply(tmp, 2, function(x) {mean(x==minority_celltype)})
  violin_plot_data_pcs = rbind(violin_plot_data_pcs, cbind(prop_minority_neighbors, num_minority_cells))
  # Save
  variable_genes_list[[i]] <- VariableFeatures(cluster_s)
  counts_allgenes = t(as.matrix(cluster_s@assays$RNA@counts))
  counts_hvgs = t(as.matrix(cluster_s@assays$RNA@counts[VariableFeatures(cluster_s), ]))
  write.table(x=counts_allgenes, file = paste0(cluster_dir_path, "counts_allgenes.csv"), sep = ',', row.names=T, col.names=T)
  write.table(x=counts_hvgs, file = paste0(cluster_dir_path, "counts_hvgs.csv"), sep = ',', row.names=T, col.names=T)
  write.table(x=cluster_s$batch, file = paste0(cluster_dir_path, "batch.csv"), sep = ',', row.names=T, col.names=T)
  write.table(x=cluster_s$celltype.l2, file = paste0(cluster_dir_path, "celltype.csv"), sep = ',', row.names=T, col.names=T)
  write.table(x=VariableFeatures(cluster_s), file = paste0(cluster_dir_path, "hvgenes_names.csv"), sep = ',', row.names=F, col.names=F)
  print(paste0('done with ', cluster_ident, '! ', i, '/ ', length(cluster_s_list)))
}

violin_plot_data_pcs$num_minority_cells = as.factor(violin_plot_data_pcs$num_minority_cells)
tmpLabels = seq(5, 90, by=5)
names(tmpLabels) = levels(violin_plot_data_pcs$num_minority_cells)
violin_plot_data_pcs$pct_minority_cells = as.factor(tmpLabels[as.character(violin_plot_data_pcs$num_minority_cells)])
p <- ggplot(violin_plot_data_pcs, 
            aes(x=pct_minority_cells, 
                y=prop_minority_neighbors)) + geom_boxplot(width=0.1, outlier.size = 0.1)#+ geom_violin()
p = p + stat_summary(fun.y=median, geom="point", size=2, color='red')
p = p #+ theme_classic() #+ geom_jitter(size=1, position=position_jitter(0.1), alpha=0.1)  
p = p + ggtitle(paste0(fixed_npcs, " PCs, varying #", minority_celltype, " cells"))
p = p + xlab('Percent Memory B cells in case') + ylab(paste0("Fraction of ", minority_celltype, " neighbors"))
p = p + theme(axis.text.x = element_text(angle = 90, hjust=0.95, size=7))
p
ggsave(paste0("prop_neighbors_acrossnumminority.png"))

saveRDS(variable_genes_list, "variable_genes_list.RDS")

#####################
#==== Clustering ====
#####################

for (i in seq_along(cluster_s_list)){
  cluster_s = cluster_s_list[[i]]
  cluster_ident = names(cluster_s_list)[i]
  print(cluster_ident)
  #dir.create(paste0("data/", cluster_ident, "/"))
  cluster_dir_path = paste0("data/", cluster_ident, "/")
  cluster_s <- NormalizeData(cluster_s)
  all.genes <- rownames(cluster_s)
  cluster_s <- ScaleData(cluster_s, features = all.genes)
  cluster_s <- FindVariableFeatures(cluster_s)
  mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
  VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
  #VariableFeatures(cluster_s) = VariableFeatures(pbmc)
  cluster_s <- RunPCA(cluster_s, npcs=npcs)
  ElbowPlot(cluster_s)
  cluster_s <- FindNeighbors(cluster_s)
  cluster_s <- FindClusters(cluster_s, resolution = 0.05) # adjust resultion to get two clusters
  #cluster_s <- FindClusters(cluster_s)
  #head(Idents(cluster_s), 5)
  #levels(Idents(cluster_s))
  #cluster_s <- RunTSNE(cluster_s, dims = 1:fixed_npcs)
  #DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters')# + ggtitle(cluster_ident)
  print(table(cluster_s@meta.data[c("batch", "seurat_clusters")])) # original labels
  table(cluster_s@meta.data[c("celltype.l2", "seurat_clusters", "seurat_clusters")]) # true labels
  
  cluster_s <- RunTSNE(cluster_s, dims = 1:fixed_npcs)
  DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters') + ggtitle("Seurat clusters")
  ggsave(paste0(cluster_dir_path, "tsne5_seurat_", 
                length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters.png"))
  
  plot_seurat_clusetrs <- function(cluster_s, cluster_dir_path){
    
    
    prop_test <- prop.test(
      x = table(cluster_s@meta.data[c("batch", "seurat_clusters")])[2, ], # case counts
      n = colSums(table(cluster_s@meta.data[c("batch", "seurat_clusters")])), # total_counts
      p = rep(0.5, length(unique(cluster_s@meta.data$seurat_clusters)))
    )
    
    p_df = cluster_s@meta.data[, c("batch", "celltype.l2", "seurat_clusters")]
    # Stacked barplot with multiple groups
    p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=batch)) +
      geom_bar(stat="identity")
    p <- p + scale_fill_manual(values=col_batch)
    p <- p + ggtitle(paste0("Chi-squared prop test p-val=", round(prop_test$p.value, 4)))
    p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic()
    p
    ggsave(paste0(cluster_dir_path, "seurat_", 
                  length(unique(cluster_s@meta.data$seurat_clusters)), 
                  "clusters_props_origlabels.png"))
    
    # Stacked barplot with multiple groups
    p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=celltype.l2)) +
      geom_bar(stat="identity")#, position=position_dodge())
    p <- p + scale_fill_manual(values=c("red", "gray"))
    p <- p + ggtitle("Ground truth labels")
    p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic()
    p
    ggsave(paste0(cluster_dir_path, "seurat_", 
                  length(unique(cluster_s@meta.data$seurat_clusters)), 
                  "clusters_props_truthlabels.png"))
  }
  
  plot_seurat_clusetrs(cluster_s, cluster_dir_path)
  write.table(x=cluster_s$seurat_clusters, file = paste0(cluster_dir_path, "seurat_cluster_labs.csv"), sep = ',', row.names=T, col.names=T)
  #cluster_s <- FindClusters(cluster_s)
  #plot_seurat_clusetrs(cluster_s, cluster_dir_path)
}

## the above code for one dataset
# 49  100  154  211  271  335  403  475  552  633  721  814  915 1023 1140 1267 1404 1555 1719 1900
i = 1

# run with default values
cluster_s = cluster_s_list[[i]]
cluster_ident = names(cluster_s_list)[i]
print(cluster_ident)
dir.create(paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/"))
cluster_dir_path = paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/")
print(cluster_dir_path)
cluster_s <- NormalizeData(cluster_s)
all.genes <- rownames(cluster_s)
cluster_s <- FindVariableFeatures(cluster_s)
mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
# find the DE genes and use those instead of variable genes later to check for separability
Idents(cluster_s) = 'celltype.l2'

markers_Memory_B_res = FindMarkers(cluster_s, ident.1 = 'Memory B', only.pos = FALSE, test.use = 't')
head(markers_Memory_B_res)
marker_genes = rownames(markers_Memory_B_res)[markers_Memory_B_res$p_val_adj<0.05]
marker_genes
print(length(marker_genes))
#VariableFeatures(cluster_s) = marker_genes
VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
#VariableFeatures(cluster_s) = VariableFeatures(pbmc)
cluster_s <- ScaleData(cluster_s, features = all.genes)
cluster_s <- RunPCA(cluster_s)#, npcs=npcs) # hyperparameter npc # default npcs is 50, by default it is ran on the hvgs as features
ElbowPlot(cluster_s)
cluster_s <- FindNeighbors(cluster_s)#, dims = 1:10) # hyperparameter dims
cluster_s <- FindClusters(cluster_s)#, resolution = 0.1)# hyperparameter resolution #, resolution = 0.1) # adjust resultion to get two clusters
#cluster_s <- FindClusters(cluster_s)
#head(Idents(cluster_s), 5)
levels(Idents(cluster_s))
#cluster_s <- RunTSNE(cluster_s) # # hyperparameter dims #, dims = 1:fixed_npcs)
#DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters')# + ggtitle(cluster_ident)
print(table(cluster_s@meta.data[c("batch", "seurat_clusters")])) # original labels
table(cluster_s@meta.data[c("celltype.l2", "seurat_clusters")]) # true labels
write.table(x=cluster_s$seurat_clusters, 
            file = paste0(cluster_dir_path, "seurat_cluster_default.csv"), 
            sep = ',', row.names=T, col.names=T)
print(cluster_dir_path)

cluster_s <- RunTSNE(cluster_s)#, dims = 1:fixed_npcs)
DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters') + ggtitle("Seurat clusters")
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_hvgs.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction="tsne", group.by='batch', cols=alpha(col_batch, 0.8), shuffle=T) + ggtitle(cluster_ident) + NoLegend()
ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_hvgs.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("#980C13", "lightgray")) + ggtitle(cluster_ident) + NoLegend()
ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_hvgs.svg"), width = 5, height = 5)
plot_seurat_clusetrs <- function(cluster_s, cluster_dir_path){
  prop_test <- prop.test(
    x = table(cluster_s@meta.data[c("batch", "seurat_clusters")])[2, ], # case counts
    n = colSums(table(cluster_s@meta.data[c("batch", "seurat_clusters")])), # total_counts
    p = rep(0.5, length(unique(cluster_s@meta.data$seurat_clusters)))
  )
  p_df = cluster_s@meta.data[, c("batch", "celltype.l2", "seurat_clusters")]
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=batch)) + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values=col_batch)
  p <- p + ggtitle(paste0("Chi-squared prop test p-val=", round(prop_test$p.value, 4)))
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_hvgs.png"), dpi='print', width = 3, height = 6)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_hvgs.svg"), width = 3, height = 6)
  print(p)
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=celltype.l2)) + geom_bar(stat="identity")#, position=position_dodge())
  p <- p + scale_fill_manual(values=c("#980C13", "lightgray"))
  p <- p + ggtitle("Ground truth labels")
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_hvgs.png"), dpi='print', width = 3, height = 6)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_hvgs.svg"), width = 3, height = 6)
  print(p)
}
plot_seurat_clusetrs(cluster_s, cluster_dir_path)
#write.table(x=cluster_s$seurat_clusters, file = paste0(cluster_dir_path, "seurat_cluster_labs.csv"), sep = ',', row.names=T, col.names=T)
#cluster_s <- FindClusters(cluster_s)
#plot_seurat_clusetrs(cluster_s, cluster_dir_path)

# with everything else kept at default, and using the hvgs, force the resolution to yield two clusters
cluster_s = cluster_s_list[[i]]
cluster_ident = names(cluster_s_list)[i]
print(cluster_ident)
dir.create(paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/"))
cluster_dir_path = paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/")
cluster_s <- NormalizeData(cluster_s)
all.genes <- rownames(cluster_s)
cluster_s <- FindVariableFeatures(cluster_s)
mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
# find the DE genes and use those instead of variable genes later to check for separability
Idents(cluster_s) = 'celltype.l2'

markers_Memory_B_res = FindMarkers(cluster_s, ident.1 = 'Memory B', only.pos = FALSE, test.use = 't')
head(markers_Memory_B_res)
marker_genes = rownames(markers_Memory_B_res)[markers_Memory_B_res$p_val_adj<0.05]
marker_genes
print(length(marker_genes))
#VariableFeatures(cluster_s) = marker_genes

VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
#VariableFeatures(cluster_s) = VariableFeatures(pbmc)
cluster_s <- ScaleData(cluster_s, features = all.genes)
cluster_s <- RunPCA(cluster_s)#, npcs=npcs) # hyperparameter npc # default npcs is 50, by default it is ran on the hvgs as features
ElbowPlot(cluster_s)
cluster_s <- FindNeighbors(cluster_s)#, dims = 1:10) # hyperparameter dims
cluster_s <- FindClusters(cluster_s, resolution = 0.1)# hyperparameter resolution #, resolution = 0.1) # adjust resultion to get two clusters
#cluster_s <- FindClusters(cluster_s)
#head(Idents(cluster_s), 5)
levels(Idents(cluster_s))
#cluster_s <- RunTSNE(cluster_s) # # hyperparameter dims #, dims = 1:fixed_npcs)
#DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters')# + ggtitle(cluster_ident)
print(table(cluster_s@meta.data[c("batch", "seurat_clusters")])) # original labels
table(cluster_s@meta.data[c("celltype.l2", "seurat_clusters")]) # true labels

#fixed_npcs = 15

cluster_s <- RunTSNE(cluster_s)#, dims = 1:fixed_npcs)
DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters') + ggtitle("Seurat clusters")
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_resolution01_hvgs.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_resolution01_hvgs.svg"), width = 5, height = 5)
#DimPlot(cluster_s, reduction="tsne", group.by='batch', cols=alpha(col_batch, 0.8), shuffle=T) + ggtitle(cluster_ident) + NoLegend()
#ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
#ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_hvgs.svg"), width = 5, height = 5)
#DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("#980C13", "lightgray")) + ggtitle(cluster_ident) + NoLegend()
#ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
#ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_hvgs.svg"), width = 5, height = 5)
plot_seurat_clusetrs <- function(cluster_s, cluster_dir_path){
  prop_test <- prop.test(
    x = table(cluster_s@meta.data[c("batch", "seurat_clusters")])[2, ], # case counts
    n = colSums(table(cluster_s@meta.data[c("batch", "seurat_clusters")])), # total_counts
    p = rep(0.5, length(unique(cluster_s@meta.data$seurat_clusters)))
  )
  p_df = cluster_s@meta.data[, c("batch", "celltype.l2", "seurat_clusters")]
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=batch)) + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values=col_batch)
  p <- p + ggtitle(paste0("Chi-squared prop test p-val=", round(prop_test$p.value, 4)))
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_resolution01_hvgs.png"), dpi='print', width = 5, height = 5)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_resolution01_hvgs.svg"), width = 5, height = 5)
  print(p)
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=celltype.l2)) + geom_bar(stat="identity")#, position=position_dodge())
  p <- p + scale_fill_manual(values=c("#980C13", "lightgray"))
  p <- p + ggtitle("Ground truth labels")
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_resolution01_hvgs.png"), dpi='print', width = 5, height = 5)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_resolution01_hvgs.svg"), width = 5, height = 5)
  print(p)
}
plot_seurat_clusetrs(cluster_s, cluster_dir_path)

# with everything else kept at default, but using the marker genes instead of hvgs, force the resolution to yield two clusters
cluster_s = cluster_s_list[[i]]
cluster_ident = names(cluster_s_list)[i]
print(cluster_ident)
dir.create(paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/"))
cluster_dir_path = paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/")
cluster_s <- NormalizeData(cluster_s)
all.genes <- rownames(cluster_s)
cluster_s <- FindVariableFeatures(cluster_s)
mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
# find the DE genes and use those instead of variable genes later to check for separability
Idents(cluster_s) = 'celltype.l2'
markers_Memory_B_res = FindMarkers(cluster_s, ident.1 = 'Memory B', only.pos = FALSE, test.use = 't')
head(markers_Memory_B_res)
marker_genes = rownames(markers_Memory_B_res)[markers_Memory_B_res$p_val_adj<0.05]
marker_genes
VariableFeatures(cluster_s) = marker_genes
#VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
#VariableFeatures(cluster_s) = VariableFeatures(pbmc)
cluster_s <- ScaleData(cluster_s, features = all.genes)
cluster_s <- RunPCA(cluster_s)#, npcs=npcs) # hyperparameter npc # default npcs is 50, by default it is ran on the hvgs as features
ElbowPlot(cluster_s)
cluster_s <- FindNeighbors(cluster_s)#, dims = 1:10) # hyperparameter dims
cluster_s <- FindClusters(cluster_s, resolution = 0.1)# hyperparameter resolution #, resolution = 0.1) # adjust resultion to get two clusters
#cluster_s <- FindClusters(cluster_s)
#head(Idents(cluster_s), 5)
levels(Idents(cluster_s))
#cluster_s <- RunTSNE(cluster_s) # # hyperparameter dims #, dims = 1:fixed_npcs)
#DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters')# + ggtitle(cluster_ident)
print(table(cluster_s@meta.data[c("batch", "seurat_clusters")])) # original labels
table(cluster_s@meta.data[c("celltype.l2", "seurat_clusters")]) # true labels
#fixed_npcs = 15
cluster_s <- RunTSNE(cluster_s)#, dims = 1:fixed_npcs)
DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters') + ggtitle("Seurat clusters")
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_resolution02_markers.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_resolution02_markers.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction="tsne", group.by='batch', cols=alpha(col_batch, 0.8), shuffle=T) + ggtitle(cluster_ident) + NoLegend()
ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_resolution02_markers.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_resolution02_markers.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("#980C13", "lightgray")) + ggtitle(cluster_ident) + NoLegend()
ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_resolution02_markers.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_resolution02_markers.svg"), width = 5, height = 5)
plot_seurat_clusetrs <- function(cluster_s, cluster_dir_path){
  prop_test <- prop.test(
    x = table(cluster_s@meta.data[c("batch", "seurat_clusters")])[2, ], # case counts
    n = colSums(table(cluster_s@meta.data[c("batch", "seurat_clusters")])), # total_counts
    p = rep(0.5, length(unique(cluster_s@meta.data$seurat_clusters)))
  )
  p_df = cluster_s@meta.data[, c("batch", "celltype.l2", "seurat_clusters")]
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=batch)) + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values=col_batch)
  p <- p + ggtitle(paste0("Chi-squared prop test p-val=", round(prop_test$p.value, 4)))
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_resolution02_markers.png"), dpi='print', width = 5, height = 5)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_resolution02_markers.svg"), width = 5, height = 5)
  print(p)
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=celltype.l2)) + geom_bar(stat="identity")#, position=position_dodge())
  p <- p + scale_fill_manual(values=c("#980C13", "lightgray"))
  p <- p + ggtitle("Ground truth labels")
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_resolution02_markers.png"), dpi='print', width = 5, height = 5)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_resolution02_markers.svg"), width = 5, height = 5)
  print(p)
}
plot_seurat_clusetrs(cluster_s, cluster_dir_path)

# with everything else kept at default, but using the marker genes instead of hvgs, default resolution
cluster_s = cluster_s_list[[i]]
cluster_ident = names(cluster_s_list)[i]
print(cluster_ident)
dir.create(paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/"))
cluster_dir_path = paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/")
cluster_s <- NormalizeData(cluster_s)
all.genes <- rownames(cluster_s)
cluster_s <- FindVariableFeatures(cluster_s)
mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
# find the DE genes and use those instead of variable genes later to check for separability
Idents(cluster_s) = 'celltype.l2'
markers_Memory_B_res = FindMarkers(cluster_s, ident.1 = 'Memory B', only.pos = FALSE, test.use = 't')
head(markers_Memory_B_res)
marker_genes = rownames(markers_Memory_B_res)[markers_Memory_B_res$p_val_adj<0.05]
marker_genes
VariableFeatures(cluster_s) = marker_genes
#VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
#VariableFeatures(cluster_s) = VariableFeatures(pbmc)
cluster_s <- ScaleData(cluster_s, features = all.genes)
cluster_s <- RunPCA(cluster_s)#, npcs=npcs) # hyperparameter npc # default npcs is 50, by default it is ran on the hvgs as features
ElbowPlot(cluster_s)
cluster_s <- FindNeighbors(cluster_s)#, dims = 1:10) # hyperparameter dims
cluster_s <- FindClusters(cluster_s)# hyperparameter resolution #, resolution = 0.1) 
#cluster_s <- FindClusters(cluster_s)
#head(Idents(cluster_s), 5)
levels(Idents(cluster_s))
#cluster_s <- RunTSNE(cluster_s) # # hyperparameter dims #, dims = 1:fixed_npcs)
#DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters')# + ggtitle(cluster_ident)
print(table(cluster_s@meta.data[c("batch", "seurat_clusters")])) # original labels
table(cluster_s@meta.data[c("celltype.l2", "seurat_clusters")]) # true labels
#fixed_npcs = 15
cluster_s <- RunTSNE(cluster_s)#, dims = 1:fixed_npcs)
DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters') + ggtitle("Seurat clusters")
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_default_resolution_markers.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
              "clusters_default_parameters_default_resolution_markers.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction="tsne", group.by='batch', cols=alpha(col_batch, 0.8), shuffle=T) + ggtitle(cluster_ident) + NoLegend()
ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_default_resolution_markers.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_default_resolution_markers.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("#980C13", "lightgray")) + ggtitle(cluster_ident) + NoLegend()
ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_default_resolution_markers.png"), dpi='print', width = 5, height = 5)
ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_default_resolution_markers.svg"), width = 5, height = 5)
plot_seurat_clusetrs <- function(cluster_s, cluster_dir_path){
  prop_test <- prop.test(
    x = table(cluster_s@meta.data[c("batch", "seurat_clusters")])[2, ], # case counts
    n = colSums(table(cluster_s@meta.data[c("batch", "seurat_clusters")])), # total_counts
    p = rep(0.5, length(unique(cluster_s@meta.data$seurat_clusters)))
  )
  p_df = cluster_s@meta.data[, c("batch", "celltype.l2", "seurat_clusters")]
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=batch)) + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values=col_batch)
  p <- p + ggtitle(paste0("Chi-squared prop test p-val=", round(prop_test$p.value, 4)))
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_default_resolution_markers.png"), dpi='print', width = 5, height = 5)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_batch_default_parameters_default_resolution_markers.svg"), width = 5, height = 5)
  print(p)
  p <- ggplot(data=p_df, aes(x=seurat_clusters, y=seurat_clusters, fill=celltype.l2)) + geom_bar(stat="identity")#, position=position_dodge())
  p <- p + scale_fill_manual(values=c("#980C13", "lightgray"))
  p <- p + ggtitle("Ground truth labels")
  p <- p + labs(x="Seurat clusters", y = "Num cells") + theme_classic() + NoLegend()
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_default_resolution_markers.png"), dpi='print', width = 5, height = 5)
  ggsave(paste0(cluster_dir_path, "seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
                "clusters_props_truelabels_default_parameters_default_resolution_markers.svg"), width = 5, height = 5)
  print(p)
}
plot_seurat_clusetrs(cluster_s, cluster_dir_path)


# Separability of the latent space plots
#npcs = 60
n_dims = 5 #5, 10, 50
#fixed_npcs = 5 # fixed

#violin_plot_data_pcs = data.frame()
#variable_genes_list <- list()
#for (i in seq_along(cluster_s_list)){
i=1
cluster_s = cluster_s_list[[i]]
cluster_ident = names(cluster_s_list)[i]
print(cluster_ident)
cluster_dir_path = paste0("~/DiseaseLabelsMethod/data/Figure2/", cluster_ident, "/")
cluster_s <- NormalizeData(cluster_s)
cluster_s <- NormalizeData(cluster_s)
all.genes <- rownames(cluster_s)
cluster_s <- FindVariableFeatures(cluster_s)
mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]
cluster_s <- ScaleData(cluster_s, features = all.genes)
cluster_s <- RunPCA(cluster_s)#, npcs=npcs) # hyperparameter npc # default npcs is 50, by default it is ran on the hvgs as features
ElbowPlot(cluster_s)
cluster_s <- FindNeighbors(cluster_s)#, dims = 1:10) # hyperparameter dims
cluster_s <- FindClusters(cluster_s)#, resolution = 0.1)# hyperparameter resolution #, resolution = 0.1) # adjust resultion to get two clusters
levels(Idents(cluster_s))
cluster_s <- RunTSNE(cluster_s)#, dims = 1:fixed_npcs)
DimPlot(cluster_s, reduction = "tsne", group.by = 'seurat_clusters') + ggtitle("Seurat clusters")
#ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
#              "clusters_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
#ggsave(paste0(cluster_dir_path, "supp_tsne_seurat_", length(unique(cluster_s@meta.data$seurat_clusters)), 
#              "clusters_default_parameters_hvgs.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction="tsne", group.by='batch', cols=alpha(col_batch, 0.8), shuffle=T) + ggtitle(cluster_ident) + NoLegend()
#ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
#ggsave(paste0(cluster_dir_path, "tsne_batch_default_parameters_hvgs.svg"), width = 5, height = 5)
DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("#980C13", "lightgray")) + ggtitle(cluster_ident) + NoLegend()
#ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_hvgs.png"), dpi='print', width = 5, height = 5)
#ggsave(paste0(cluster_dir_path, "tsne_truelabels_default_parameters_hvgs.svg"), width = 5, height = 5)

#==== Neighbors in PC space ==== 
# track the mixing in the latent space, compute for the perturbed.minority cells only
num_minority_cells = sum(cluster_s@meta.data$celltype.l2==minority_celltype)
n_dims_vec = seq(from=2, to=50, by=1)
violin_plot_data = data.frame()
for(n_dims in n_dims_vec) {
  print(n_dims)
  cluster_s <- FindNeighbors(cluster_s, compute.SNN = TRUE, dims = 1:n_dims)
  tmp = apply(cluster_s@graphs[[1]][cluster_s@meta.data$celltype.l2==minority_celltype,], 1, 
              function(x) {cluster_s@meta.data$celltype.l2[x == 1]})
  # vector of length number of minoruty cells containing the proportion of their neighborst might minority label
  prop_minority_neighbors = apply(tmp, 2, function(x) {mean(x==minority_celltype)})
  numPCs = rep(n_dims, num_minority_cells)
  violin_plot_data = rbind(violin_plot_data, cbind(prop_minority_neighbors, numPCs))
  print(paste0('done with ', n_dims, 'PCs!'))
}
# Plots
# 1. across numPCs for the same num_minority cells
violin_plot_data$numPCs = as.factor(violin_plot_data$numPCs)
p <- ggplot(violin_plot_data, aes(x=numPCs, y=prop_minority_neighbors)) + geom_boxplot(width=0.1, outlier.size = 0.1)#+ geom_violin()
p = p + stat_summary(fun.y=median, geom="point", size=2, color='red')
p = p #+ theme_classic() #+ geom_jitter(size=1, position=position_jitter(0.1), alpha=0.1)  
p = p + ggtitle(paste0(num_minority_cells, " ", minority_celltype, " cells"))
p = p + xlab('Num PCs') + ylab(paste0("Fraction of ", minority_celltype, " neighbors"))
p = p + theme(axis.text.x = element_text(angle = 90, hjust=0.95, size=7))
ggsave(paste0(cluster_dir_path, "prop_neighbors_acrosspcs.png"), dpi='print', width = 10, height = 5)
ggsave(paste0(cluster_dir_path, "prop_neighbors_acrosspcs.svg"), width = 10, height = 5)
print(p)
#  # 2. across num_minority cells, fixed numPCs
#  cluster_s <- FindNeighbors(cluster_s, compute.SNN = TRUE, dims = 1:fixed_npcs)
#  tmp = apply(cluster_s@graphs[[1]][cluster_s@meta.data$celltype.l2==minority_celltype,], 1, 
#              function(x) {cluster_s@meta.data$celltype.l2[x == 1]})
#  # vector of length number of minoruty cells containing the proportion of their neighborst with minority label
#  prop_minority_neighbors = apply(tmp, 2, function(x) {mean(x==minority_celltype)})
#  violin_plot_data_pcs = rbind(violin_plot_data_pcs, cbind(prop_minority_neighbors, num_minority_cells))
#  # Save
#  variable_genes_list[[i]] <- VariableFeatures(cluster_s)
#  counts_allgenes = t(as.matrix(cluster_s@assays$RNA@counts))
#  counts_hvgs = t(as.matrix(cluster_s@assays$RNA@counts[VariableFeatures(cluster_s), ]))
#  write.table(x=counts_allgenes, file = paste0(cluster_dir_path, "counts_allgenes.csv"), sep = ',', row.names=T, col.names=T)
#  write.table(x=counts_hvgs, file = paste0(cluster_dir_path, "counts_hvgs.csv"), sep = ',', row.names=T, col.names=T)
#  write.table(x=cluster_s$batch, file = paste0(cluster_dir_path, "batch.csv"), sep = ',', row.names=T, col.names=T)
#  write.table(x=cluster_s$celltype.l2, file = paste0(cluster_dir_path, "celltype.csv"), sep = ',', row.names=T, col.names=T)
#  write.table(x=VariableFeatures(cluster_s), file = paste0(cluster_dir_path, "hvgenes_names.csv"), sep = ',', row.names=F, col.names=F)
#  print(paste0('done with ', cluster_ident, '! ', i, '/ ', length(cluster_s_list)))
#}
#violin_plot_data_pcs$num_minority_cells = as.factor(violin_plot_data_pcs$num_minority_cells)
#tmpLabels = seq(5, 90, by=5)
#names(tmpLabels) = levels(violin_plot_data_pcs$num_minority_cells)
#violin_plot_data_pcs$pct_minority_cells = as.factor(tmpLabels[as.character(violin_plot_data_pcs$num_minority_cells)])
#p <- ggplot(violin_plot_data_pcs, 
#            aes(x=pct_minority_cells, 
#                y=prop_minority_neighbors)) + geom_boxplot(width=0.1, outlier.size = 0.1)#+ geom_violin()
#p = p + stat_summary(fun.y=median, geom="point", size=2, color='red')
#p = p #+ theme_classic() #+ geom_jitter(size=1, position=position_jitter(0.1), alpha=0.1)  
#p = p + ggtitle(paste0(fixed_npcs, " PCs, varying #", minority_celltype, " cells"))
#p = p + xlab('Percent Memory B cells in case') + ylab(paste0("Fraction of ", minority_celltype, " neighbors"))
#p = p + theme(axis.text.x = element_text(angle = 90, hjust=0.95, size=7))
#p
#ggsave(paste0("prop_neighbors_acrossnumminority.png"))



#########################
#==== Variable Genes ====
#########################
col_perturb = c("gray", "red")
col_batch = c("#1B9E77", "#D95F02")
minority_celltype = 'Memory B'
variable_genes_list = readRDS("variable_genes_list.RDS")
i = 18
j = 1
# when looking at the genes unique to a harder dataset compared to an easy one, and plotting them in the easier dataset,
# the unique genes are less variable than the shared genes, if any such trend at all could be called significant
genes_int_ij = intersect(variable_genes_list_tmp[[i]], variable_genes_list_tmp[[j]])
genes_unique_to_i = setdiff(variable_genes_list_tmp[[i]], genes_int_ij)
genes_unique_to_j = setdiff(variable_genes_list_tmp[[j]], genes_int_ij)

cluster_s_i = cluster_s_list[[i]]

cells_majority_i = cluster_s_i$celltype.l2!=minority_celltype
cells_minority_i = cluster_s_i$celltype.l2==minority_celltype

#genes_int_ij = intersect(variable_genes_list_tmp[[i]], variable_genes_list_tmp[[j]])
#genes_unique_to_i = setdiff(variable_genes_list_tmp[[i]], genes_int_ij)
#genes_unique_to_j = setdiff(variable_genes_list_tmp[[j]], genes_int_ij)

type = rep(c("shared"), length(variable_genes_list[[j]]))
unique_gene_indx = variable_genes_list[[j]] %in% genes_unique_to_j
type[unique_gene_indx] = "unique"
PlotDF = data.frame(
  MeanNaive = rowMeans(as.matrix(cluster_s_i@assays$RNA@data[variable_genes_list[[j]], cells_majority_i])),
  MeanMemory = rowMeans(as.matrix(cluster_s_i@assays$RNA@data[variable_genes_list[[j]], cells_minority_i])),
  Type = type)
#PlotDF
ggplot(PlotDF, aes(x=MeanNaive, y=MeanMemory,shape=Type, color=Type)
) + geom_point() + geom_abline(slope=1, intercept = 0
) + ggtitle(paste0('genes unique to ', j, ' compared to ', i, ' dataset, in dataset ', i))

library(pheatmap)
library(RColorBrewer)
MyMatrix = matrix(NA, nrow=length(variable_genes_list),ncol=length(variable_genes_list))
for (i in 1:length(variable_genes_list)) {
  for (j in 1:length(variable_genes_list)){
    MyMatrix[i, j] = length(intersect(variable_genes_list_tmp[[i]], variable_genes_list_tmp[[j]]))
  }
}
dimnames(MyMatrix) = list(paste0("Cluster",1:nrow(MyMatrix)),paste0("Cluster",1:nrow(MyMatrix)))
pheatmap(
  mat = MyMatrix, cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100)
)

# use the hvgs from 18 (easy case), to try and separate the hard case
i=1
j=18
#genes_int_ij = intersect(variable_genes_list_tmp[[i]], variable_genes_list_tmp[[j]])
#genes_unique_to_i = setdiff(variable_genes_list_tmp[[i]], genes_int_ij)
#genes_unique_to_j = setdiff(variable_genes_list_tmp[[j]], genes_int_ij)


violin_plot_data_pcs = data.frame()
cluster_s = cluster_s_list[[i]]
cluster_ident = names(cluster_s_list)[i]

#cluster_s = naiveB_1900_cd4_memory_49
#cluster_ident = 'naiveB_1900_cd4_memory_49'
#minority_celltype = 'CD4 Memory'
#table(naiveB_1900_cd4_memory_49$celltype.l2)

print(cluster_ident)


#dir.create(paste0("data/", cluster_ident, "/"))
#cluster_dir_path = paste0("data/", cluster_ident, "/")
cluster_s <- NormalizeData(cluster_s)
all.genes <- rownames(cluster_s)
cluster_s <- ScaleData(cluster_s, features = all.genes)

#cluster_s <- FindVariableFeatures(cluster_s)
#mouse_vargenes_bool <- startsWith(VariableFeatures(cluster_s), "MOUSE-")
#VariableFeatures(cluster_s) = VariableFeatures(cluster_s)[!mouse_vargenes_bool]

VariableFeatures(cluster_s) = critical_genes2000#genes_unique_to_j#variable_genes_list[[j]]
cluster_s <- RunPCA(cluster_s, npcs=npcs)
#==== Can we find the structure with PCA? ====
DimPlot(cluster_s, reduction='pca', group.by = 'batch', cols = col_batch) + ggtitle(cluster_ident)
#ggsave(paste0(cluster_dir_path, "pca2_batch.png"))
DimPlot(cluster_s, reduction='pca', group.by = 'celltype.l2', cols = c("red", "gray")) + ggtitle(cluster_ident)
#ggsave(paste0(cluster_dir_path, "pca2_label.png"))
n_dims = 5
cluster_s <- RunTSNE(cluster_s, dims = 1:n_dims)
DimPlot(cluster_s, reduction = "tsne", group.by = 'batch', cols = col_batch) + ggtitle(cluster_ident)
#ggsave(paste0(cluster_dir_path, "tsne5_batch.png"))
DimPlot(cluster_s, reduction = "tsne", group.by = 'celltype.l2', cols = c("red", "gray")) + ggtitle(n_dims)
#ggsave(paste0(cluster_dir_path, "tsne5_label.png"))
#==== Neighbors in PC space ==== 
# track the mixing in the latent space, compute for the perturbed.minority cells only
num_minority_cells = sum(cluster_s@meta.data$celltype.l2==minority_celltype)
n_pcs_vec = seq(from=2, to=npcs, by=1)
violin_plot_data = data.frame()
for(n in n_pcs_vec) {
  print(n)
  cluster_s <- FindNeighbors(cluster_s, compute.SNN = TRUE, dims = 1:n)
  tmp = apply(cluster_s@graphs[[1]][cluster_s@meta.data$celltype.l2==minority_celltype,], 1, 
              function(x) {cluster_s@meta.data$celltype.l2[x == 1]})
  # vector of length number of minoruty cells containing the proportion of their neighborst might minority label
  prop_minority_neighbors = apply(tmp, 2, function(x) {mean(x==minority_celltype)})
  numPCs = rep(n, num_minority_cells)
  violin_plot_data = rbind(violin_plot_data, cbind(prop_minority_neighbors, numPCs))
  print(paste0('done with ', n, 'PCs!'))
}
# Plots
# 1. across numPCs for the same num_minority cells
violin_plot_data$numPCs = as.factor(violin_plot_data$numPCs)
p <- ggplot(violin_plot_data, aes(x=numPCs, y=prop_minority_neighbors)) + geom_boxplot(width=0.1)#+ geom_violin()
p = p + stat_summary(fun=median, geom="point", size=2, color='red')
p = p #+ theme_classic() #+ geom_jitter(size=1, position=position_jitter(0.1), alpha=0.1)  
p = p + ggtitle(paste0(num_minority_cells, " ", minority_celltype, " cells"))
p = p + xlab('Num PCs') + ylab(paste0("Fraction of ", minority_celltype, " neighbors"))
p
#ggsave(paste0(cluster_dir_path, "prop_neighbors_acrosspcs.png"))
# using the hvgs from the easy case maybe very slightly improves the separationn, but still definnitely doesn't solve the problem
# if we used specifically the unique genes to the easy case, the result is even worse than using just default hvgs for hard case

############################
##### Find ctirical pcs ####
############################
get_cor_labelsvec_allpcs <- function(perturbed_labels) {
  cor_vec = rep(NA, npcs)
  for (k in 1:npcs) {
    cor_vec[k] = cor(cluster_s@reductions$pca@cell.embeddings[,k], 
                     perturbed_labels)}
  return(cor_vec)
}

cor_vec = get_cor_labelsvec_allpcs(as.numeric(cluster_s$celltype.l2==minority_celltype))
plot(cluster_s@reductions$pca@stdev, cor_vec, pch=20)
abline(a=0, b=0, h=1, v=1)
order(cor_vec, decreasing = T)
DimPlot(cluster_s, group.by = 'celltype.l2', reduction = 'pca', 
        dims = c(1, 7),
        cols = c("red", "gray")) + ggtitle('After perturbation')

# get a null distribution by making other 'perturbed' vectors - sampling 5% disease cells 
# and then do an approx normal test
num_null_labelsvec = 100
cor_null_labels_pcs_matrix = matrix(data=NA, 
                                    nrow = num_null_labelsvec, 
                                    ncol = npcs)
for (i in 1:num_null_labelsvec) {
  set.seed(i+50)
  perturbed_cells_null = sample(x=Cells(cluster_s[, cluster_s$batch=='1']), 
                                size=sum(cluster_s$celltype.l2==minority_celltype))
  perturbed_null = as.numeric(Cells(cluster_s) %in% perturbed_cells_null)
  cor_null_labels_pcs_matrix[i, ] = get_cor_labelsvec_allpcs(perturbed_null)
}

null_means = colMeans(cor_null_labels_pcs_matrix)
null_std = matrixStats::colSds(cor_null_labels_pcs_matrix)
print(dim(cor_null_labels_pcs_matrix))
plot(cor_null_labels_pcs_matrix[2,], cor_vec, pch=20)
abline(a=0, b=0, h=1, v=1)

p_vals_cor = sapply(1:npcs,function(k){2*(1-pnorm(abs((cor_vec[k]-null_means[k])/null_std[k])))})
#mean(abs(cor_vec[1]) > abs(cor_null_labels_pcs_matrix[, 1]))
p_vals_cor_nonpar = sapply(1:npcs,function(k){mean(abs(cor_vec[k]) < abs(cor_null_labels_pcs_matrix[, k]))})

plot(p_vals_cor, p_vals_cor_nonpar, pch=20)
sum(p_vals_cor<(0.05/npcs))
which(p_vals_cor<(0.05/npcs))
head(order(p_vals_cor, decreasing = F), 10)

##############################
##### Find ctirical genes ####
##############################
n_genes = dim(cluster_s@assays$RNA@data)[1]
get_cor_labelsvec_allgenes <- function(perturbed_labels) {
  cor_vec = rep(NA, npcs)
  for (k in 1:n_genes) {
    cor_vec[k] = cor(cluster_s@assays$RNA@data[k,], 
                     perturbed_labels)}
  return(cor_vec)
}

cor_vec = get_cor_labelsvec_allgenes(as.numeric(cluster_s$celltype.l2==minority_celltype))
plot(seq(1:n_genes), cor_vec, pch=19, cex=0.1)
abline(a=0, b=0, col="red")
order(abs(cor_vec), decreasing = T)
critical_genes = rownames(cluster_s@assays$RNA@data)[head(order(abs(cor_vec), decreasing = T), 10)]
FeaturePlot(cluster_s, features=critical_genes, reduction = 'tsne')
Idents(cluster_s)
VlnPlot(cluster_s, features = critical_genes, cols = c("red", "gray"))

DimPlot(cluster_s, group.by = 'celltype.l2', reduction = 'pca', 
        dims = c(1, 7),
        cols = c("red", "gray")) + ggtitle('After perturbation')

critical_genes2000 = rownames(cluster_s@assays$RNA@data)[head(order(abs(cor_vec), decreasing = T), 2000)]
critical_genes100 = rownames(cluster_s@assays$RNA@data)[head(order(abs(cor_vec), decreasing = T), 100)]
table(abs(cor_vec[order(abs(cor_vec), decreasing = T)])>0.05)
critical_genes800 = rownames(cluster_s@assays$RNA@data)[head(order(abs(cor_vec), decreasing = T), 800)]
# what percent of the critical genes are in the hvgs?
cluster_s <- FindVariableFeatures(cluster_s)
hvgs = VariableFeatures(cluster_s)
length(intersect(hvgs, critical_genes2000))#/2000
# there is a total of 427 genes that are both critical and are hvgs
# are they the top critical genes?
plot(x=seq(1, 2000, 1), 
     y=critical_genes2000 %in% intersect(hvgs, critical_genes2000), 
     cex = 0.3)
# another look at the same question - as we dialate the set of critical genes,
# what is the rate of increase of the size oof the overlap with the 2000 hvgs?
size_g_int = rep(NA, 2000)
for (i in 1:2000) {
  critical_genes_v = rownames(cluster_s@assays$RNA@data)[
    head(order(abs(cor_vec), decreasing = T), i)]
  size_g_int[i] = length(intersect(hvgs, critical_genes_v))
}
plot(x=seq(1, 2000, 1), y=size_g_int, cex = 0.3,
     xlab = 'size of set of critical genes', ylab = 'size of intersection',
     main=paste0(cluster_ident, '\n', '2000 HVGs and x critical genes'), 
     ylim=c(0, 2000), xlim=c(0, 2000))
abline(a=0, b=1, col = 'red')