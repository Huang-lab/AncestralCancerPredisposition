##### dependency_files.R #####
# Kuan-lin Huang @ WashU 2018 April
# dependent files for analysis in the PCA Germline project

# the new pathVar is already filtered, tn_swap adjusted
fn = "../../data/PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv"
pathVar = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)

##### subsets #####
pathVarP = pathVar[pathVar$Overall_Classification %in% c("Pathogenic","Likely Pathogenic"),]
pathVarOT = pathVar[!is.na(pathVar$Gene_Classification) & pathVar$Gene_Classification != "None",]
pathVarPOT = pathVarP[!is.na(pathVarP$Gene_Classification) & pathVarP$Gene_Classification != "None",]

##### clinical files #####
clin_f = "../../data/PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)

clin_complete_f = "../../data/clinical_PANCAN_patient_with_followup.tsv.gz"
clin_complete = read.table(header=T, quote = "", sep="\t", fill =T, file = gzfile(clin_complete_f), stringsAsFactors=FALSE)

# featured genes in Pan10389 manuscript
gene_fn = "../../data/PCA_feature_gene_list.txt"
glist_f = read.table(header=FALSE, stringsAsFactors = F, file = gene_fn)
featGenes = as.vector(t(glist_f))