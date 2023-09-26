library(cmapR)
library(rhdf5)
library(tidyverse)
library(AMARETTO)


col_meta1 <- read.delim(path_sig_info1, sep="\t", stringsAsFactors=F)
col_meta2 <- read.delim(path_sig_info2, sep="\t", stringsAsFactors=F)

table(col_meta1$cell_id)
table(col_meta1$pert_type)

col_meta_filtered_1<-dplyr::filter(col_meta1,cell_id %in% c("PHH","HEPG2"))%>%dplyr::filter(pert_type %in% c("trt_sh.cgs","trt_oe"))
dim(col_meta_filtered_1)
head(col_meta_filtered_1)

col_meta_filtered_2<-dplyr::filter(col_meta2,cell_id %in% c("PHH","HEPG2"))%>%dplyr::filter(pert_type %in% c("trt_sh.cgs","trt_oe"))
dim(col_meta_filtered_2)
head(col_meta_filtered_2)

table(col_meta2$pert_type)
table(col_meta2$cell_id)

sig_ids_1<-col_meta_filtered_1$sig_id
subset_gct_file1 <- parse.gctx(path_gct_file1, cid=sig_ids_1)

saveRDS(list(gct_file=subset_gct_file1,info_file=col_meta_filtered_1),"./subset_data/subset_gct_info_list.rds")
write_tsv(col_meta_filtered_1,"./subset_data/col_meta_filtered_1")
write.gctx(subset_gct_file1, "./subset_data/subset_gct_file1")

###################################### #start from here

library(annotate)
library(org.Hs.eg.db)

subset_all_data<-readRDS("./subset_data/subset_gct_info_list.rds")
subset_gct_file1<-subset_all_data$gct_file
col_meta_filtered_1<-subset_all_data$info_file
signature_matrix<-subset_gct_file1@mat
gene_names<-as.vector(getSYMBOL(rownames(signature_matrix),data='org.Hs.eg'))
rownames(signature_matrix)<-gene_names
pert_gene_name<-sapply(colnames(signature_matrix), function(x) col_meta_filtered_1$pert_iname[which(col_meta_filtered_1$sig_id==x)])
pert_type<-sapply(colnames(signature_matrix), function(x) col_meta_filtered_1$pert_type[which(col_meta_filtered_1$sig_id==x)])

########################################
# module based fgsea computations :


AMARETTOinit<-readRDS("AMARETTOinit_TCGA.rds")
AMARETTOresults<-readRDS("AMARETTOresults_TCGA.rds")

regulator_list<-list()
target_list<-list()
regulator_list_coefficients<-list()
for (ModuleNr in 1:AMARETTOresults$NrModules){
  ModuleRegulators <- names(which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0))
  ModuleRegulators_coefficients<-AMARETTOresults$RegulatoryPrograms[ModuleNr,ModuleRegulators]
  ModuleTargets <- names(AMARETTOresults$ModuleMembership[which(AMARETTOresults$ModuleMembership==ModuleNr),1])
  regulator_list[[ModuleNr]]<-ModuleRegulators
  regulator_list_coefficients[[ModuleNr]]<-ModuleRegulators_coefficients
  target_list[[ModuleNr]]<-ModuleTargets
}

overlapping_perts_indexes<-which(pert_gene_name%in%regulator_list[[ModuleNr]])
module_specific_signature_matrix<-signature_matrix[,overlapping_perts_indexes]
module_specific_pert_gene_names<-pert_gene_name[overlapping_perts_indexes]
module_specific_pert_types<-pert_type[overlapping_perts_indexes]

if(length(overlapping_perts_indexes)==0){
  print("No overlap between drivers and existing perturbations")
}
if(length(overlapping_perts_indexes)==1){
  module_specific_signature_matrix<-matrix(module_specific_signature_matrix,nrow = length(module_specific_signature_matrix),ncol = 1)
  colnames(module_specific_signature_matrix)<-colnames(signature_matrix)[overlapping_perts_indexes]
  rownames(module_specific_signature_matrix)<-rownames(signature_matrix)
}


dim(module_specific_signature_matrix)
head(module_specific_signature_matrix)

resultloop<-NULL
for(j in 1:length(colnames(module_specific_signature_matrix))){
  print(j)
  fgseaRes <- fgsea(target_list[ModuleNr],
                    stats = module_specific_signature_matrix[,j],
                    minSize=1,
                    maxSize=50000,
                    nperm=1000)
  fgseaRes$sig_id <- rep(colnames(module_specific_signature_matrix)[j],nrow(fgseaRes))
  fgseaRes$pert_iname <- rep(module_specific_pert_gene_names[j],nrow(fgseaRes))
  fgseaRes$pert_type <- rep(module_specific_pert_types[j],nrow(fgseaRes))
  fgseaRes$regulator_amaretto_value <- rep(regulator_list_coefficients[[ModuleNr]][module_specific_pert_gene_names[j]],nrow(fgseaRes))
  resultloop<-rbind(resultloop,fgseaRes)
} 

resultloop$padj <-p.adjust(resultloop$pval, method = "BH")
library(DT)
DT::datatable(resultloop)



##########

library(fgsea)
data(examplePathways)
data(exampleRanks)

AMARETTOinit<-readRDS("AMARETTOinit_TCGA.rds")
AMARETTOresults<-readRDS("AMARETTOresults_TCGA.rds")

regulator_list<-list()
target_list<-list()

for (ModuleNr in 1:AMARETTOresults$NrModules){
  ModuleRegulators <- AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
  ModuleTargets <- names(AMARETTOresults$ModuleMembership[which(AMARETTOresults$ModuleMembership==ModuleNr),1])
  regulator_list[[ModuleNr]]<-ModuleRegulators
  target_list[[ModuleNr]]<-ModuleTargets
  }

NrCores=1
cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
registerDoParallel(cluster,cores=NrCores)


resultloop2<-foreach(j=1:5500, .combine='rbind') %do%{
  print(j)
  fgseaRes <- fgsea(target_list,
                    stats = signature_matrix[,j],
                    minSize=5,
                    maxSize=5000,
                    nperm=10)
  fgseaRes$pert <- rep(colnames(signature_matrix)[j],nrow(fgseaRes))
  fgseaRes$pert_iname<-rep(colnames(signature_matrix)[j],nrow(fgseaRes))
}

resultloop3<-NULL

for(j in 1:5500){
  print(j)
  fgseaRes <- fgsea(target_list,
                    stats = signature_matrix[,j],
                    minSize=5,
                    maxSize=5000,
                    nperm=10)
  fgseaRes$pert <- rep(colnames(signature_matrix)[j],nrow(fgseaRes))
  resultloop3<-rbind(resultloop3,fgseaRes)
} 


# find genes associated with the Perturbation name
pert_gene_name<-sapply(resultloop$pert, function(x) col_meta_filtered_1$pert_iname[which(col_meta_filtered_1$sig_id==x)])
resultloop$pert_gene_name<-pert_gene_name



#write_csv(resultloop,'resultloop.csv')
#saveRDS(resultloop,"resultloop.RDS")
#resultloop$leadingEdge<-NULL
fwrite(resultloop,"resultloop.csv")
#saveRDS(resultloop,'resultloop.rds')
resultloop<-readRDS('resultloop.rds')

##########################################################################################
##########################################################################################














dim(resultloop)
head(resultloop)

#Nathalie's Experiment     change pathway number and cell.
sum("CES2"%in%resultloop$pert_name)
kk<-filter(resultloop,pert=="CGS001_HEPG2_96H:CES2:1.5")%>%filter(pathway==32)
ppp<-signature_matrix[,colnames(signature_matrix)=="CGS001_HEPG2_96H:CES2:1.5"]
ppp[names(ppp) %in% kk$leadingEdge[[1]] ]

AMARETTOresults$RegulatoryPrograms[1:5,1:5]

dim(AMARETTOresults$RegulatoryPrograms)

view(AMARETTOresults$RegulatoryPrograms)

AMARETTOresults$RegulatoryPrograms
dim(AMARETTOresults$RegulatoryPrograms)

for (i in 1:175){
print(i)
mm<-AMARETTOresults$RegulatoryPrograms[i,][AMARETTOresults$RegulatoryPrograms[i,]<0]
print(mm)
print(names(mm) %in%resultloop$pert_name )
}


ppp[names(ppp) %in% c("ERG") ]













sum(resultloop$pval<0.05)

head(resultloop)
resultloop



library(tidyverse)

module_number<-1
regulator_list[[module_number]]

drivers_interested<-regulator_list[[module_number]][regulator_list[[module_number]] %in% gene_names]

drivers_interested
mm<-dplyr::filter(resultloop,pathway==module_number)%>%dplyr::filter(pert_name %in% drivers_interested )
mm

dim(mm)
mm$padj
mm$pval

ModuleNr<-1
ModuleRegulators <- AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]

ModuleTargets <- names(AMARETTOresults$ModuleMembership[which(AMARETTOresults$ModuleMembership==ModuleNr),1])

cat(ModuleTargets)

resultloop$pert_name
