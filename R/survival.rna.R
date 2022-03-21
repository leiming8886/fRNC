survival.rna <- function(gene_profile = NULL, clinData = NULL){
	#library("survival")
	#library("survminer")
	intersection<-intersect(colnames(gene_profile),clinData[,1])
	length(intersection)
	a=0;d=0
	for(i in 1:length(intersection)){
	a[i]=which(clinData[,1]==intersection[i])  ##clinData
	d[i]=which(colnames(gene_profile)==intersection[i])  ##profileData
	}
	clinData<-clinData[a,]
	rownames(clinData)<-c(1:nrow(clinData))
	gene_profile<-gene_profile[,d]
	ab<-c()
	cox_gene<-c()
	p_gene<-c()
	cox<-c()
	p<-c()
	#cox_three<-c()
	#p_three<-c()
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
	#p_gene
	colnames(p_gene) <- c("genes","PValue")
	p_gene_type <- merge(p_gene,fRNC::IDsymbol[,c("genes","type")], by="genes", all.x=TRUE)
	p_gene_type <- p_gene_type[!duplicated(p_gene_type),]
	rownames(p_gene_type) <- p_gene_type$genes
	p_gene_type <- p_gene_type[order(p_gene_type$PValue),]
	surivallist=list(rna_p = p_gene_type, algorithm = "univariable Cox")
	#write.table(p_mir,file = file_miR, sep ="\t", row.names = F,col.names = F, quote = F)
	#write.table(p_gene, file = file_rna, sep ="\t", row.names = F,col.names = F, quote = F)
	return(surivallist)
}




