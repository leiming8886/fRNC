#' @title survival analysis
#'
#' @description The function calculate the P-values from a univariable Cox proportional hazards regression model
#' between survival data and corresponding expression of node. First, the overlap of the sample IDs among the survival data and
#' the two expression data is obtained, then based on the overlap sample IDs, the p-value of every node from a univariable Cox proportional
#' hazards regression model is calculated on the R package survival.
#'
#' @param miRNA_profile the miRNA expressin data matrix,in which the row name is gene id and the column name is sample id.
#' @param gene_profile the mRNA expressin data matrix,in which the row name is gene id and the column name is sample id.
#' @param clinData the data matrix of survival data,in which the column name is "ID"(sample IDs),  "Survival"(months) and "Status"(0,1)
#' @return surivallis list contain three elements,
#'  miR_p the p-value matrix of IDs,
#'  rna_p the p-value matrix of IDs,
#'  algorithm  a character string indicating which algorithm  was used
#'
#' @examples
#' \dontrun{
#' result_miR <- DEGs(case.exp_miRNA,control.exp_miRNA,
#' geneid= rownames(control.exp_miRNA),data_type = "RNAseq_counts")
#' result_rna <- DEGs(case.exp_rna,control.exp_rna,
#' geneid= rownames(control.exp_rna), data_type = "RNAseq_counts")
#' interac <- interac[,c("node_gene_ID","type","target_gene_ID")]
#' survival.miR_rna(miRNA_profile=result_miR$Nor_expr,
#' gene_profile = result_rna$Nor_expr, clinData = clinData)
#' }
#' @export

survival.miR_rna <- function(miRNA_profile= NULL, gene_profile = NULL, clinData = NULL){
	#library("survival")
	#library("survminer")
	intersection1<-intersect(colnames(gene_profile),clinData[,1])
	intersection<-intersect(intersection1,colnames(miRNA_profile))
	length(intersection)
	a=0;b=0;d=0
	for(i in 1:length(intersection)){
	a[i]=which(clinData[,1]==intersection[i])  ##clinData
	b[i]=which(colnames(miRNA_profile)==intersection[i])  #profileData
	d[i]=which(colnames(gene_profile)==intersection[i])  ##profileData
	}
	clinData<-clinData[a,]
	rownames(clinData)<-c(1:nrow(clinData))
	miRNA_profile<-miRNA_profile[,b]
	gene_profile<-gene_profile[,d]
	ab<-c()
	#file_miR<-paste(name_surP_mir,"_P.txt",sep = "")
	#file_rna<-paste(name_surP_rna,"_P.txt",sep = "")
	cox_gene<-c()
	cox_mir<-c()
	p_gene<-c()
	p_mir<-c()
	cox<-c()
	p<-c()
	cox_three<-c()
	p_three<-c()
	mirname <- rownames(miRNA_profile)
	for ( i in 1:length(mirname)){
		miRNA<-miRNA_profile[which(rownames(miRNA_profile) == mirname[i]),]
		miRNA=as.matrix(miRNA)
		n2<-grep(paste("^",0,"$",sep=""),perl=T,miRNA)
		miRNA[n2]<-0.00001
		miRNA_log<-log2(as.numeric(miRNA))
		miRNA_log_norm<-(miRNA_log-mean(miRNA_log))/var(miRNA_log)
		test2 <- list(time=clinData$Survival,
				status=clinData$Status,
				x=as.numeric(miRNA_log_norm)
		)
	t2 <- coxph(survival::Surv(time, status) ~ x , test2)
	cox_mir<-c(cox_mir,summary(t2)[[7]][1])
	p_mir<-rbind(p_mir,c(as.character(mirname[i]),summary(t2)[[7]][5]))
	}
	genename <- rownames(gene_profile)
	for ( i in 1:length(genename)){
		gene<-gene_profile[which(rownames(gene_profile)==genename[i]),]
		gene=as.matrix(gene)
		n1<-grep(paste("^",0,"$",sep=""),perl=T,gene)
		gene[n1]<-0.00001
		gene_log<-log2(as.numeric(gene))
		gene_log_norm<-(gene_log-mean(gene_log))/var(gene_log)
		test1 <- list(time=clinData$Survival,
				status=clinData$Status,
				x=as.numeric(gene_log_norm)
		)
	t1<-coxph(survival::Surv(time, status) ~ x , test1)
	cox_gene<-c(cox_gene,summary(t1)[[7]][1])
	p_gene <- rbind(p_gene,c(as.character(genename[i]),summary(t1)[[7]][5]))
	}
	#result = list(surP_mir = p_mir, surP_gene= p_gene)
	colnames(p_mir) <- c("genes","PValue")
	type_mir <- c(rep("miRNA",dim(p_mir)[1]))
	p_mir <- cbind(p_mir,type_mir)
	p_mir <- as.data.frame(p_mir)
	colnames(p_mir) <- c("genes","PValue","type")
	rownames(p_mir) <- p_mir$genes

	#p_gene
	colnames(p_gene) <- c("genes","PValue")
	p_gene_type <- merge(p_gene,fRNC::IDsymbol[,c("genes","type")], by="genes", all.x=TRUE)
	p_gene_type <- p_gene_type[!duplicated(p_gene_type),]
	rownames(p_gene_type) <- p_gene_type$genes
	p_gene_type <- p_gene_type[order(p_gene_type$PValue),]

	surivallist=list(miR_p = p_mir,rna_p = p_gene_type, algorithm = "univariable Cox")

	#write.table(p_mir,file = file_miR, sep ="\t", row.names = F,col.names = F, quote = F)
	#write.table(p_gene, file = file_rna, sep ="\t", row.names = F,col.names = F, quote = F)
	return(surivallist)
}




