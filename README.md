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
devtools::install_github("leiming8886/fRNC",ref = "master")
library(fRNC)
```
## Example
-------
Below is a quick example of how the package can be used.
We downloaded miRNA, mRNA and circRNA expression data of esophageal carcinoma (ESCA) and normal samples from TCGA and MiOncoCirc. 
The miRNA, mRNA-seq data consist of 161 tumors with 11 normal samples in the TCGA, and circRNA data with 19 tumors and 25 normal samples
 in MiOncoCirc. The miRNA-lncRNA and miRNA-circRNA interactions were extracted from ENCORI with a high stringency, where the number of 
 the supported CLIP experimental evidence is 3 or greater. The R package edgeR was used to analyze the differentially expressed miRNAs, 
 mRNAs and lncRNAs. The result was saved with the data "dataN". 
```R
library(fRNC)
data("dataN")
gene2weight <- combinp(dataN[,c("type","logFC","PValue")])
interac <- interStringency(type = "ncRNA",stringency = "strict")
interac <- interac[,c("node_gene_ID","type","target_gene_ID")]
res.list_global <- runmodule(network = interac, gene2weight, method = "global",FDR = 1e-14)
res.list_local <- runmodule(network = interac, gene2weight, method = "local",
maxsize=15, seletN = "MIMAT0000461")
```
Save global and local module results respectively. And, the result was saved as XGMML file and then observed it in the Cytoscape environment.
```R
saveNetwork(res.list_global$module,file="ceRNA_module",type = "XGMML")
savelocalM(res.list_local)
```

## Contact
If you have an comments, suggestions, corrections or ideas to install ,use, improve or extend this package, feel free to contact me. 
my email is leiming8886\@163.com. You also submit a report on the  [Github issues page](https://github.com/leiming8886/fRNC/issues)