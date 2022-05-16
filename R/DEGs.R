#' @title Performe differential expression analysis
#'
#' @description Differential expression analysis using edgeR or limma for two group comparison
#'
#' @param case.exp the case expression matrix, in which the row name is gene id and the column name is sample id
#' @param control.exp the control expression matrix, in which the row name is gene id (the same with the case) and the column name is sample id
#' @param geneid  gene id in the case or control expression matrix
#' @param data_type a character string indicating which date type to deal with is to be choosed, One of "RNAseq_counts" , "fpkm" and "microarray": can be abbreviated
#'
#' @return DEGlist list contain four elements,
#'  DEGs the differential expression matrix,
#'  Nor_expr the normalized expression matrix,
#'  data_type a character string of the data type
#'  algorithm  a character string indicating which algorithm  was used, One of "edgeR" , "limma"
#'
#' @examples
#' \dontrun{
#' data("case.exp_miRNA")
#' data("control.exp_miRNA")
#' result_miR <- DEGs(case.exp_miRNA,control.exp_miRNA,
#' geneid= rownames(control.exp_miRNA), data_type = "RNAseq_counts")
#'
#'}
#' @export
DEGs <- function(case.exp,control.exp,geneid,data_type){
	## some ----------------------------------------------------
	# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
	fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
	normalize_rna <-function(counts,number = 3){
          keep.counts <- ( rowMeans(counts) >= 0 ) & ( rowSums(counts > 0) >= number )
          counts <- counts[keep.counts, ]
          return(counts)
	}
	if(data_type=="RNAseq_counts"){
		#library(edgeR)
		exp1=cbind(case.exp,control.exp)
		rownames(exp1)=geneid
		group=factor(c(rep("Tumor",ncol(case.exp)),rep("Normal",ncol(control.exp))))
		design <- model.matrix(~ 0+group)
		colnames(design) <- levels(group)
		contrast.matrix <- makeContrasts(contrasts = "Tumor-Normal",
                                 levels = design)
		#创建edgeR数据格式
		data1=DGEList(counts=exp1,genes=geneid,group=group)
		#过滤
		index1=rowSums(cpm(data1)>1)>=(ncol(exp1)/2)
		data1=data1[index1,]
		#标准化，默认为TMN
		data1=calcNormFactors(data1)
		nor_data<-cpm(data1)
    #data1=estimateDisp(data1,design)
		data1=estimateGLMCommonDisp(data1,design)
		data1=estimateGLMTrendedDisp(data1,design)
		data1=estimateGLMTagwiseDisp(data1,design)

		fit <- glmFit(data1, design)
		lrt <- glmLRT(fit, contrast = contrast.matrix)
		lrt <- topTags(lrt,nrow(lrt))
		DEG <- lrt$table
		one_ID <- rownames(nor_data)[1]
		if(grepl("ENS",one_ID, ignore.case = T)){
		  DEG_type <- merge(DEG,fRNC::IDsymbol[,c("genes","type")], by="genes", all.x=TRUE)
		  DEG_type <- DEG_type[!duplicated(DEG_type),]
		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		if(grepl("MIMAT",one_ID, ignore.case = T)){
		  type <- c(rep("miRNA",dim(DEG)[1]))
		  DEG_type <- cbind(DEG,type)

		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		if(grepl("circ",one_ID, ignore.case = T)){
		  type <- c(rep("circRNA",dim(DEG)[1]))
		  DEG_type <- cbind(DEG,type)

		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		DEGlist=list(DEGs = DEG_type,Nor_expr = nor_data, data_type = data_type, algorithm = "edgeR")
	}
	if(data_type=="microarray"){
		#library(limma)
		labe1 = factor(c(rep("Tumor",ncol(case.exp)),rep("Normal",ncol(control.exp))))
		exp1 = cbind(case.exp,control.exp)
		rownames(exp1) = geneid
		design1 = model.matrix(~0+labe1)
		colnames(design1) = levels(labe1)
		fit1 = lmFit(exp1,design1)
		contrast.matrix = makeContrasts(contrasts = "Tumor-Normal",levels=design1)
		fit11 = contrasts.fit(fit1,contrast.matrix)
		fit12 = eBayes(fit11)
		nor_data = exp1#exprs(exp1)
		DEG1 = topTable(fit12,adjust.method ="BH",number = nrow(exp1))
		DEG1 = cbind(as.character(rownames(DEG1)),DEG1)
		colnames(DEG1)[c(1,5,6)] = c("genes","PValue","fdr")
		one_ID <- rownames(nor_data)[1]
		if(grepl("ENS",one_ID, ignore.case = T)){
		  DEG_type <- merge(DEG1,fRNC::IDsymbol[,c("genes","type")], by="genes", all.x=TRUE)
		  DEG_type <- DEG_type[!duplicated(DEG_type),]
		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		if(grepl("MIMAT",one_ID, ignore.case = T)){
		  type <- c(rep("miRNA",dim(DEG1)[1]))
		  DEG_type <- cbind(DEG1,type)

		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		if(grepl("circ",one_ID, ignore.case = T)){
		  type <- c(rep("circRNA",dim(DEG1)[1]))
		  DEG_type <- cbind(DEG1,type)

		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		DEGlist = list(DEGs = DEG_type, Nor_expr = nor_data, data_type=data_type, algorithm="limma")
	}
	if(data_type=="TPM"|| data_type=="FPKM"){
		#比较靠谱的方法是limma-trend，以log2(TPM+1)替代logCPM以limma使用手册中15.4章。
		library(limma)
		labe1 = factor(c(rep("Tumor",ncol(case.exp)),rep("Normal",ncol(control.exp))))
		exp1 = cbind(case.exp,control.exp)
		#FPKM convert to TPM
		if(data_type == "FPKM"){
		exp1 <- apply(exp1,2,fpkmToTpm)
		}
		#log2(TPM+1)
		exp1 <- log2(exp1 + 1)
		rownames(exp1) = geneid
		### remove lowly expressed genes，mean of gene >0 and rowSum of (gene >0) > half of sample
		exp1 <- normalize_rna(exp1,number = (ncol(exp1)/2))
		#limma归一化
		exp2 <- normalizeBetweenArrays(exp1)
		design1 = model.matrix(~0+labe1)
		colnames(design1) = levels(labe1)
		fit1 = lmFit(exp2,design1)
		contrast.matrix = makeContrasts(contrasts = "Tumor-Normal",levels=design1)
		fit2 = contrasts.fit(fit1,contrast.matrix)
		#trend = TRUE
		fit3 = eBayes(fit2,trend = TRUE)
		nor_data = exp1#exprs(exp1)
		DEG1 = topTable(fit3,adjust.method ="BH",number = nrow(exp2))
		DEG1 = cbind(as.character(rownames(DEG1)),DEG1)
		colnames(DEG1)[c(1,5,6)] = c("genes","PValue","fdr")
		one_ID <- rownames(nor_data)[1]
		if(grepl("ENS",one_ID, ignore.case = T)){
		  DEG_type <- merge(DEG1,fRNC::IDsymbol[,c("genes","type")], by="genes", all.x=TRUE)
		  DEG_type <- DEG_type[!duplicated(DEG_type),]
		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		if(grepl("MIMAT",one_ID, ignore.case = T)){
		  type <- c(rep("miRNA",dim(DEG1)[1]))
		  DEG_type <- cbind(DEG1,type)

		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		if(grepl("circ",one_ID, ignore.case = T)){
		  type <- c(rep("circRNA",dim(DEG1)[1]))
		  DEG_type <- cbind(DEG1,type)

		  rownames(DEG_type) <- DEG_type$genes
		  DEG_type <- DEG_type[order(DEG_type$PValue),]
		}
		DEGlist = list(DEGs = DEG_type, Nor_expr = nor_data, data_type = data_type, algorithm="limma")
	}
	return(DEGlist)
}
