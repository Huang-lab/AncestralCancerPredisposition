##### analyze_path_burden_byAncestry.R #####
# Kuan-lin Huang 2018
# find non-cancer pathogenic variant in the ExAC cohort
# targeting specific ethnicity

# set work dir for testing in dev environ
bdir = "~/Box\ Sync/Huang_lab/manuscripts/germlineEthnicPower/analysis/burden_assoc"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")
source("../stat_functions.R")
#source("./label_onco_var_ExAC.R")
source("./TFT.R")

##### analyze #####

#all_cancer_stat_byCancer_all_ancestry = data.frame()

for (ethni in unique(clin$consensus_call))# some ethni are too lacking
#for (ethni in c("eur","asian","afr"))
{

#ethni = "afr" # use for testing
  
clin_ethni = clin[clin$consensus_call==ethni,] 
#clin_ethni = clin[clin$self_identified_race==ethni,] #for using self-reported ethnicity

pathVarP_ethni = pathVarP[pathVarP$bcr_patient_barcode %in% clin_ethni$bcr_patient_barcode,]
# compile by cancer type table
sample_size = data.frame(table(clin_ethni$type))
colnames(sample_size) = c("Cancer","Sample_size")
PCA_count_byCancer = data.frame(table(pathVarP_ethni$HUGO_Symbol,pathVarP_ethni$cancer))
colnames(PCA_count_byCancer) = c("Gene","Cancer","Count")
PCA_count_byCancer = merge(PCA_count_byCancer,sample_size,by="Cancer")
PCA_count_byCancer$Freq = PCA_count_byCancer$Count/PCA_count_byCancer$Sample_size #each cancer type
  
##### implement burden testing #####

# limit to cancers with at least 25 cases in that ancestry
test_cancer = unique(PCA_count_byCancer$Cancer[PCA_count_byCancer$Sample_size > 24])
# 2. Cancer vs. cancer
AF_thres = 0.0005
all_cancer_stat_byCancer = vector("list")
for (cancer in test_cancer){
  data_c = PCA_count_byCancer[PCA_count_byCancer$Cancer %in% cancer & PCA_count_byCancer$Count > 1,]
  
  if (nrow(data_c) > 0){
    cancer_stat = run_TFT_against_others(data_c,PCA_count_byCancer,NULL,AF_thres)
    cancer_stat$Cancer = cancer
    all_cancer_stat_byCancer[[cancer]]=cancer_stat
  }
}
all_cancer_stat_byCancer_m = do.call(rbind,all_cancer_stat_byCancer)

all_cancer_stat_byCancer_m$FDR = p.adjust(all_cancer_stat_byCancer_m$P, method="BH")
all_cancer_stat_byCancer_m = all_cancer_stat_byCancer_m[order(all_cancer_stat_byCancer_m$P),]

tn = paste("out/cancer_vs_otherCancer_pathogenic_variants_burden",ethni,"tsv",sep=".")
write.table(all_cancer_stat_byCancer_m, quote=F, sep="\t", file = tn, row.names = F)

all_cancer_stat_byCancer_m$cohort_AF = 100*all_cancer_stat_byCancer_m$cohort_AC/all_cancer_stat_byCancer_m$cohort_AN
all_cancer_stat_byCancer_m_sig = all_cancer_stat_byCancer_m[all_cancer_stat_byCancer_m$FDR < 0.05,]
all_cancer_stat_byCancer_m_suggest = all_cancer_stat_byCancer_m[all_cancer_stat_byCancer_m$FDR < 0.15,]

##### let's plot the results #####
### volcano plots ###
getPalette = colorRampPalette(c("#9e9ac8","#3f007d"))
p = ggplot(data=all_cancer_stat_byCancer_m_suggest)
#p = p + geom_point(aes(y=-log10(FDR),x= cohort_AF,color = Cancer),alpha=0.5)
#p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
p = p + geom_text_repel(aes(y=cohort_AF,x=Cancer,color = Cancer,label=ifelse(FDR<0.03, Gene,NA)))
p = p + getPCACancerColor() #+ scale_y_log10()
p = p + labs(y= "Cohort Carrier Frequency (%)") + expand_limits(y=0)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
#p = p + theme(legend.position = "bottom")
p
fn = paste('out/Gene_by_carrierFrequency_suggestive',ethni,'pdf',sep=".")
ggsave(fn,h=5,w=6,useDingbat=F)

# ### frequency ###
# # merge with ExAC first #
# colnames(variants_cancerP_gsum) = colnames(PCA_count) 
# all_counts = rbind(variants_cancerP_gsum,PCA_count,PCA_count_byCancer)
# 
# all_counts$Gene_Classification = "TSG"
# all_counts$Gene_Classification[all_counts$Gene %in% all_oncogenes] = "Oncogene"
# 
# sigColors=c("NA", "#000000")
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# 
# gene_order = PCA_count$Gene[order(PCA_count$Count,decreasing = T)]
# cancer_order = c("ExAC_nonTCGA","TCGA",as.character(unique(PCA_count_byCancer$Cancer)))
# 
# all_counts$Gene = factor(all_counts$Gene, levels = rev(gene_order))
# all_counts$Cancer = factor(all_counts$Cancer, levels = cancer_order)
# all_counts = all_counts[!is.na(all_counts$Gene),]
# #high_genes = unique(all_counts$Gene[all_counts$Count > 19])
# 
# all_counts$Freq_plot = 100*signif(all_counts$Freq, digits = 2)
# all_counts$Freq_plot_fill = all_counts$Freq_plot
# all_counts$Freq_plot_fill[all_counts$Freq_plot_fill>2] = 2
# all_counts$Freq_plot[all_counts$Freq_plot==0] = NA
# 
# p = ggplot(data=all_counts)
# #p = p + facet_grid(.~Cancer,drop=T,scales = "free_x", space = "free_x")
# p = p + geom_tile(data=all_counts,aes(x=Cancer, y=Gene, fill= Freq_plot_fill), linetype="blank") + scale_fill_gradientn(name= "% Carriers", colours=getPalette(100), na.value=NA, limit=c(0,NA))
# p = p + geom_text(data=all_counts,aes(x=Cancer, y=Gene, label = Freq_plot), color="black", size=3)
# p = p + geom_tile(data=all_cancer_stat_byCancer_m_suggest,aes(x=Cancer, y=Gene), color="grey",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
# p = p + geom_tile(data=all_cancer_stat_byCancer_m_sig,aes(x=Cancer, y=Gene), color="black",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
# p = p  + theme_bw() + theme_nogrid() +
#   theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# fn = 'out/ExAC_nonTCGA_vs_ExAC_geneFreq_heatmap_byCancers_sig.pdf'
# ggsave(fn,height = 20,width=15, useDingbat=F)
# 
# ### showing less genes
# genes = all_cancer_stat_byCancer_m_sig$Gene
# all_counts_g = all_counts[all_counts$Gene %in% genes,]
# all_cancer_stat_byCancer_m_suggest_p = all_cancer_stat_byCancer_m_suggest[all_cancer_stat_byCancer_m_suggest$Gene %in% genes,]
# 
# p = ggplot()
# #p = p + facet_grid(Gene_Classification~.,drop=T,scales = "free_y", space = "free_y")
# p = p + geom_tile(data=all_counts_g,aes(x=Cancer, y=Gene, fill= Freq_plot_fill), linetype="blank") + scale_fill_gradientn(name= "% Carriers", colours=getPalette(100), na.value=NA, limit=c(0,NA))
# p = p + geom_text(data=all_counts_g,aes(x=Cancer, y=Gene, label = Freq_plot), color="black", size=3)
# p = p + geom_tile(data=all_cancer_stat_byCancer_m_suggest_p,aes(x=Cancer, y=Gene), color="grey",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
# p = p + geom_tile(data=all_cancer_stat_byCancer_m_sig,aes(x=Cancer, y=Gene), color="black",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
# p = p  + theme_bw() + theme_nogrid() +
#   theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# 
# fn = 'out/ExAC_nonTCGA_vs_ExAC_geneFreq_heatmap_byCancers_sig_sele.pdf'
# ggsave(fn,height = 6,width=13, useDingbat=F)
# 
# write(genes, "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/PCA_feature_gene_list.txt", sep="\n")
# 
# ### using ExAC assoc as significance
# genes = all_cancer_stat_m_sig$Gene
# all_counts_g = all_counts[all_counts$Gene %in% genes,]
# 
# p = ggplot()
# p = p + geom_tile(data=all_counts_g,aes(x=Cancer, y=Gene, fill= Freq_plot_fill), linetype="blank") + scale_fill_gradientn(name= "% Carriers", colours=getPalette(100), na.value=NA, limit=c(0,NA))
# p = p + geom_text(data=all_counts_g,aes(x=Cancer, y=Gene, label = Freq_plot), color="black", size=3)
# p = p + geom_tile(data=all_cancer_stat_m_suggest,aes(x=Cancer, y=Gene), color="grey",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
# p = p + geom_tile(data=all_cancer_stat_m_sig,aes(x=Cancer, y=Gene), color="black",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
# p = p  + theme_bw() + theme_nogrid() +
#   theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# 
# fn = 'out/ExAC_nonTCGA_vs_ExAC_geneFreq_heatmap_byCancers_ExAC_sig_sele.pdf'
# ggsave(fn,height = 8,width=15, useDingbat=F)
}

