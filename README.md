## fRNC: An R package to uncover RBP-ncRNA circuits from multi-omics data


## Description

The RNA binding protein (RBP) and non-coding RNA (ncRNA) interacting networks are increasingly recognized as main
mechanism in post-transcriptional regulation, and tightly associated with cellular malfunction and disease. Here,
we present fRNC, a systems biology tool to uncover dynamic spectrum of RBP-ncRNA circuits (RNC) by integrating 
transcriptomics, interactomics and clinical data. fRNC constructs the RBP-ncRNA network from experiment derived
CLIP-seq or PARE data. Given scoring on each node in the network, it finds a RNC containing global maximum significant 
genes. Alternatively, it can also search locally maximum RNCs according to user defined nodes. It enables users flexibly
to analyse and visualize the collective behaviors between a RBP and its interacting ncRNAs in a malfunctioned biological process.
The attachment and reference manual of this package is available in the docs directory 

## Getting Started
### Step 1. Install package dependencies
Enter the R function (install_dependpackages) in the R below and run it. some messages will appear to inform you whether or not any R packages dependencies have been installed.
```R
     install_dependpackages <- function(){
         metr_pkgs <- c("limma", "ggpubr", "XML", "igraph", "multtest","RBGL","edgeR")  
         list_installed <- installed.packages()
         new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"])) 
         if(length(new_pkgs)!=0){   
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
                BiocManager::install(new_pkgs)
                print(c(new_pkgs, " packages will be installed..."))
            }  
          if((length(new_pkgs)<1)){
                print("No new packages will be installed...")
            }
        }
```
```R
     install_dependpackages()
```
### Step 2. Install the package
fRNC is freely available from GitHub.
```R
# Step 1: Install devtools
install.packages("devtools")
library(devtools)
# Step 2: Install fRNC
devtools::install_github("leiming8886/fRNC",ref = "main")
library(fRNC)
```
## Example
-------
Below is a quick example of how the package can be used for case 1 from the ESCA dataset.
We downloaded miRNA, mRNA and circRNA expression data of esophageal carcinoma (ESCA) and normal samples from TCGA. 
The miRNA, mRNA-seq data consist of 161 tumors with 11 normal samples in the TCGA. The miRNA-lncRNA and miRNA-mRNA(RBP) interactions 
belonged to the category " transcriptome " were extracted from ENCORI with a high stringency, where the number of the supported CLIP
 experimental evidence is 3 or greater. The R package edgeR was used to analyze the differentially expressed miRNAs, 
 mRNAs and lncRNAs. 
```R
library(fRNC)

rm(list=ls())
data("case.exp_miRNA")
data("control.exp_miRNA")
result_miR <- DEGs(case.exp_miRNA,control.exp_miRNA, geneid= rownames(control.exp_miRNA), data_type = "RNAseq_counts")
data("case.exp_rna")
data("control.exp_rna")
result_rna <- DEGs(case.exp_rna,control.exp_rna, geneid= rownames(control.exp_rna), data_type = "RNAseq_counts")
print("1: Data Loaded and differentially expressed analyse")
dataNo <- rbind(result_miR$DEGs, result_rna$DEGs)
gene2weight <- combinp(dataNo[,c("type","logFC","PValue")], islog = T)
interac_rbp <- interStringency(type = "transcription", spec ="hg",stringency = "high")
print("2: interaction Loaded")
interac_rbp <- interac_rbp[,c("node_gene_ID","type","target_gene_ID")]
res.list_global <- runmodule(network = interac_rbp, gene2weight, method = "global", FDR = 1e-5)
print("3: global module detected")
res.list_local <- runmodule(network = interac_rbp, gene2weight, method = "local", maxsize=15, seletN = c("MIMAT0000089") )
print("4: local module detected with the seed")
```
Save global and local module results respectively. And, the result was saved as XGMML file and then observed it in the Cytoscape environment.
```R
saveNetwork(res.list_global$module,file="ceRNA_module_transcription",type = "XGMML")
savelocalM(res.list_local)
```
the RNC relevant to BRCA at protein level can be detected with the fRNC for case 2.
The datasets contain 852 cancer samples and 18 normal samples, simultaneously measured all mRNA, miRNA, and protein expression. To save running time, we processed in advance to perform 
differential expression analysis and store internal data.
```R
#load the result from differential expression analysis
data("brca_miRNA_re_se")
data("brca_RNA_re_se")
data("brca_pro_re_se")

#Expression data of RBP (mRNA) was removed from mRNA expression.
RNA_re_se_type <- RNA_re_se$DEGs[,c("genes","logFC","PValue","type")]
RNA_re_se_type_no_RBP <- RNA_re_se_type[which(RNA_re_se_type$type != "rbp"),]
dataN_age_all <- rbind(pro_re_se$DEGs[,c("genes","logFC","PValue","type")], RNA_re_se_type_no_RBP, miRNA_re_se$DEGs[,c("genes","logFC","PValue","type")])

#Load interaction data
gene2weight_age <- combinp(dataN_age_all[,c("type","logFC","PValue")], islog = T)
interac_RBP_age <- interStringency(type = "Protein", spec ="hg",stringency = "high")
interac_RBP_age <- interac_RBP_age[,c("node_gene_ID","type","target_gene_ID")]
#head(interac_RBP_age)
#global module detected with FDR = 1e-20
res.list_global_brca_e_20 <- runmodule(network = interac_RBP_age, gene2weight = gene2weight_age, method = "global", FDR =  1e-20)
```
## Contact
If you have an comments, suggestions, corrections or ideas to install ,use, improve or extend this package, feel free to contact me. 
my email is leiming8886\@163.com. You also submit a report on the  [Github issues page](https://github.com/leiming8886/fRNC/issues)