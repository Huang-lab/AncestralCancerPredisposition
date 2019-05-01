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


# prepare data files
tn = "~/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/burden_assoc/out/ExAC_vs_cancer_pathogenic_variants_burden.tsv"
all_cancer_stat_against_exac = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = tn)

# correct for sequencing center
fileName = "~/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/variant_QC/PCA.r1.TCGAbarcode.merge.exon.vcf.istats.tsv"
istat = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
istat$center = gsub(".*(..)","\\1",istat$ID)
# table(istat$center) #https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes
istat$center[istat$center=="01"] = "BI" # both BI
istat$center[istat$center=="08"] = "BI" # both BI
istat$center[istat$center=="09"] = "WUSM"
istat$center[istat$center=="10"] = "BCM"
istat$center[istat$center=="32"] = "SANGER"
istat$bcr_patient_barcode = substr(istat$ID,1,12)
istat = istat[!duplicated(istat$bcr_patient_barcode),]
colnames(PCs)[1] = "bcr_patient_barcode"
PCs = PCs[!duplicated(PCs$bcr_patient_barcode),]
clin_pc = merge(clin,PCs, by= "bcr_patient_barcode",all.x=T,all.y=F)
clin_pc_center = merge(clin_pc,istat[,c("bcr_patient_barcode","center")],by="bcr_patient_barcode",all.x=T,all.y=F)

##### analyze #####
results_list = as.list(NULL)
i=1
for (ethni in unique(clin$consensus_call)[!is.na(unique(clin$consensus_call))]){# some ethni are too lacking
  #for (ethni in c("eur","asian","afr"))
  #ethni = "afr" # use for testing
  
  clin_ethni = clin_pc_center[clin_pc_center$consensus_call==ethni,] 
  clin_ethni = clin_ethni[!is.na(clin_ethni$type),]
  
  for (gene in unique(pathVarP$HUGO_Symbol)){
    
    #merge into that gene's dataframe
    pathVarG = pathVarP[pathVarP$HUGO_Symbol == gene & !duplicated(pathVarP$bcr_patient_barcode),]
    clin_pathvarG = merge(clin_ethni, pathVarG[,-which(colnames(pathVarG) %in% colnames(clin_ethni)[-1])], by= "bcr_patient_barcode",all.x =T)
    clin_pathvarG$HUGO_Symbol[is.na(clin_pathvarG$HUGO_Symbol)] = "WT"
    
    # # exclude suggestive cancers per testing against ExAC
    # exclude_cancers = all_cancer_stat_against_exac$Cancer[all_cancer_stat_against_exac$Gene==gene & all_cancer_stat_against_exac$FDR < 0.15]
    # 
    for (cancer in unique(clin_pathvarG$type)){
      # # exclude the suggestive cancer types above; not implemented in the end as that's based on the pan-ancestry test
      # exclude_cancers = exclude_cancers[exclude_cancers!= cancer]# ensure include this original cancer type
      # clin_pathvarG = clin_pathvarG[!(clin_pathvarG$type %in% exclude_cancers),]
      
      clin_pathvarG$inCancer = clin_pathvarG$type== cancer
      subset_pathvar_count = sum(clin_pathvarG$HUGO_Symbol!="WT" & clin_pathvarG$inCancer)
      if (sum(clin_pathvarG$inCancer,na.rm=T ) < 25){next} # only test for cancer types with at least 25 cases 
      if (subset_pathvar_count < 2){ next} # only test for cancer types with at least 2 pathvar
      cat(paste("Processing: gene =", gene, " cancer =", cancer, " ancestry =", ethni, "\n") )
      model_results = run_glm(data = clin_pathvarG, yi = "inCancer", xi = "HUGO_Symbol", ytype = "Binary", covi = c("age_at_initial_pathologic_diagnosis","gender","PC1","PC2"))
      cancer_gene_stat = data.frame(cbind(ethni, cancer, gene,subset_pathvar_count,model_results))
      results_list[[i]] = cancer_gene_stat; i = i + 1
    }
  }
}

tt = do.call(rbind,results_list)


colnames(tt) = c("ancestry","cancer","gene","gene_path_count","y","y_type","gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "pvalue_chisq","xi_lvl1","beta_coefficient","StdErr","Zvalue","pvalue_Walds","covariants");
tt$OR = exp(-tt$beta_coefficient)
tt$OR_lower_bound95 = exp(-tt$beta_coefficient + qnorm(0.025) * tt$StdErr)
tt$OR_higher_bound95 = exp(-tt$beta_coefficient + qnorm(0.975) * tt$StdErr)
tt$conf_int95 = paste(tt$OR_lower_bound95,tt$OR_higher_bound95,sep="-")

tt$FDR_chisq = p.adjust(tt[,"pvalue_chisq"], method="fdr") 
tt=tt[order(tt$pvalue_chisq, decreasing=FALSE),]

# should use Wald's test; results seems more robust
tt$FDR_walds = p.adjust(tt[,"pvalue_Walds"], method="fdr")
tt=tt[order(tt$pvalue_Walds, decreasing=FALSE),]

tn = "out/within_ancestry_regression_stats.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

tt_sub=tt[tt$gene %in% featGenes,]
tt_sub$FDR_walds = p.adjust(tt_sub[,"pvalue_Walds"], method="fdr") 
tt_sub$FDR_chisq = p.adjust(tt_sub[,"pvalue_chisq"], method="fdr") 
tn = "out/within_ancestry_regression_stats_featGenes.txt"
write.table(tt_sub, quote=F, sep="\t", file = tn, row.names = F)
