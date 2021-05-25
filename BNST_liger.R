
library(Seurat)

library(liger)

# https://singlecell.broadinstitute.org/single_cell/study/SCP466/analysis-of-human-substantia-nigra-sn-and-mouse-bed-nucleus-of-the-stria-terminalis-bnst-with-liger#study-download
load("~/DiseaseLabelsMethod/data/BNST_only_liger.robj")
str(a.bed.josh.clean)
a.bed.josh.clean
class(a.bed.josh.clean)
a.bed.josh.clean@raw.data
dim(a.bed.josh.clean@raw.data[[1]])
table(a.bed.josh.clean@clusters)
#"ac_Sln", "ac_Col5a3", "al_Chat", "al_Cyp26b1", "al_Sema3e", "al_Dpp4", "al_Impg1", "am_Slc17a8", "am_Trp73"       
#"am_Wif1", "a_Npffr2", "a_Esr1", "a_Gli3", "a_Zeb2", "a_Bmpr1b", "a_Synpo2", "a_Pax6", "a_Th", "a_Arhgef38"
#"a_4932435O22Rik", "a_Lrrc9", "ov_Sh3d21", "ov_Vipr2", "p_Ror1", "p_Nxph2", "p_Myo1h", "pr_St18", "p_Epsti1"
#"p_Haus4", "p_Bnc2", "p_Tac2", "p_Ebf1", "p_Angpt1", "p_Sst", "pr_Esr2", "p_Cplx3", "_Vip", "_Avp", "_C1ql3"         
#"_Piezo2", "_Cdh23" "pr_Esr2"
for (cluster_name in levels(a.bed.josh.clean@clusters)){
  print(cluster_name)
  bnstpr_esr2_cellnames = names(a.bed.josh.clean@clusters)[as.character(a.bed.josh.clean@clusters)==cluster_name]
  bnstpr_esr2_male = a.bed.josh.clean@raw.data$MALE[,
                                                    colnames(a.bed.josh.clean@raw.data$MALE) %in% bnstpr_esr2_cellnames]
  
  bnstpr_esr2_female = a.bed.josh.clean@raw.data$FEMALE[,
                                                        colnames(a.bed.josh.clean@raw.data$FEMALE) %in% bnstpr_esr2_cellnames]
  ###################
  #### Save Data ####
  ###################
  #variable_genes_list[[i]] <- VariableFeatures(cluster_s)
  all_genes = intersect(rownames(bnstpr_esr2_female), rownames(bnstpr_esr2_male))
  #length(all_genes) # intersection of Male and Female genes
  counts_allgenes = t(as.matrix(cbind(bnstpr_esr2_female[all_genes,], bnstpr_esr2_male[all_genes,])))
  dir.create(paste0("~/DiseaseLabelsMethod/data/", cluster_name, "/"))
  write.table(x=counts_allgenes, file = paste0("~/DiseaseLabelsMethod/data/", cluster_name, "/", 
                                               "counts_allgenes.csv"), sep = ',', row.names=T, col.names=T)
  batch = as.numeric(grepl("FEMALE", rownames(counts_allgenes)))
  # 1 is FEMALE, 0 is MALE
  sex = ifelse(batch, "FEMALE", "MALE")
  table(sex)
  write.table(x=sex, file = paste0("~/DiseaseLabelsMethod/data/", cluster_name, "/", 
                                   "sex.csv"), sep = ',', row.names=T, col.names=T)
  #write.table(x=cluster_s$celltype.l2, file = paste0(cluster_dir_path, "celltype.csv"), sep = ',', row.names=T, col.names=T)
  #write.table(x=VariableFeatures(cluster_s), file = paste0(cluster_dir_path, "hvgenes_names.csv"), sep = ',', row.names=F, col.names=F)

}
head(a.bed.josh.clean@clusters)
#levels(a.bed.josh.clean@clusters)
#length(a.bed.josh.clean@clusters)
#sum(a.bed.josh.clean@clusters == "BNSTpr_Esr2")
#dim(bnstpr_esr2_male)
#dim(bnstpr_esr2_female)
#'Xist' %in% (rownames(bnstpr_esr2_male))
#summary(bnstpr_esr2_male["Xist",])
#summary(bnstpr_esr2_female["Xist",])
#summary(bnstpr_esr2_male["Esr2", ])
#summary(bnstpr_esr2_female["Esr2",])





