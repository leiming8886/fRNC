#' @title Performe univariate Cox regression analysis and plot Kaplan-Meier curve of the RNC module
#'
#' @description univariate Cox regression analysis using survival with survival data and expression data,and Patients with higher and lower than the median risk score of the he 
#' dysregulated RNC are classified into different groups. Kaplan-Meier survival analysis was used to assess the clinical significance between the comparison groups
#' @param gene_profile the expression value matrix, in which the row name is gene id and the column name is sample id
#' @param clinData the survival data, in which the column name is sample id, Survival(time) and Status(0,1)
#' @param genes the gene in the RNC, The name of mature miRNA is miRbase ID, and name of lncRNA and RBP is Ensemble gene ID
#' @param filename the output figure name in Kaplan-Meier curve

#' @return contain the figure with the pdf format in Kaplan-Meier curve
#' @export


survival.km <- function(gene_profile = NULL, clinData = NULL,genes= NULL,filename="temp"){
	#library("survival")
	#library("survminer")
	#gene_profile = case_exp
	#clinData = clinData
	#genes = module_genes
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
	
	normalizeMatrix <- function (mRNA_matrix, d = 1) 
	{
	#mRNA_matrix<-data_t[,1:variableL]
    norm_mRNA_matrix <- matrix(nrow = dim(mRNA_matrix)[1], ncol = dim(mRNA_matrix)[2], 
        data = NA)
    rownames(norm_mRNA_matrix) <- rownames(mRNA_matrix)
    colnames(norm_mRNA_matrix) <- colnames(mRNA_matrix)
    mu <- apply(mRNA_matrix, d, mean)
    sigma <- apply(mRNA_matrix, d, sd)
    if (d == 1) {#行
        for (i in 1:dim(mRNA_matrix)[d]) {
            norm_mRNA_matrix[i, ] <- (mRNA_matrix[i, ] - mu[i])/sigma[i]
        }
    }
    if (d == 2) {#列
        for (i in 1:dim(mRNA_matrix)[d]) {
            norm_mRNA_matrix[, i] <- (mRNA_matrix[, i] - mu[i])/sigma[i]
        }
    }
    return(norm_mRNA_matrix)
	}
	#0替换0.00001
	gene_profile[gene_profile == 0]<-0.00001
	gene_profile_log<-log2(gene_profile)
	gene_profile_log_zscore <- normalizeMatrix(gene_profile_log,1)
	
	
	gene_profile_log_zscore_contain_genes <- gene_profile_log_zscore[genes,]
	for ( i in genes){
		test1 <- list(time=clinData$Survival,
				status=clinData$Status,
				x=as.numeric(gene_profile_log_zscore_contain_genes[i,])
		)
	t1<-coxph(survival::Surv(time, status) ~ x , test1)
	cox_gene<-c(cox_gene,summary(t1)[[7]][1])
	}
	gene_profile_log_zscore_contain_genes_plus_r <- cox_gene%*%gene_profile_log_zscore_contain_genes
	
	fenlei<-c()
	#identical(as.character(colnames(gene_profile_log_zscore_contain_genes_plus_r)),as.character(clinData[,1]))
	m<-median(gene_profile_log_zscore_contain_genes_plus_r)
	for(j in 1:length(gene_profile_log_zscore_contain_genes_plus_r)){
		if(gene_profile_log_zscore_contain_genes_plus_r[j]>=m){fenlei[j]<-1}else{fenlei[j]<-0}
	}
	y <- Surv(clinData$Survival,clinData$Status==1)
	kmfit2 <- survfit(y~as.numeric(fenlei))
#
	diff<-survdiff(Surv(clinData$Survival, clinData$Status) ~ as.numeric(fenlei))
	Figer<-diff
	getPvalue<-function(x=Figer,digits = max(options()$digits - 4, 3)){### 函数5＿  得到p倿 
	if (is.matrix(x$obs)){
    otmp <- apply(x$obs,1,sum)
    etmp <- apply(x$exp,1,sum)
	}else{
    otmp <- x$obs
    etmp <- x$exp
	}
	df <- (sum(1*(etmp>0))) -1
	pvalue<-format(signif(1-pchisq(x$chisq, df),digits))
	return(pvalue);
	}

	p<-getPvalue(Figer);
  filename<-paste(filename,".pdf",sep="")
pdf(file=filename);
#plot(kmfit2, mgp=c(3.5,1,0),lty = 2:3,xlab="Survival Time (days)",ylab="Survival #Probability",main=y,col=c("#ff9605","#067a04"),cex=2,cex.lab=2.2,cex.main=1.8,cex.axis=1.8,las=1,lwd=2)
mainlabel<-"Survival curves"                  #标题的名秿
par(bg="white", mai=c(1,1.2,1,1.2),family="serif");
plot(kmfit2, mgp=c(3.5,1,0.1),lty =c(1,1),xlab="Time (months)",ylab="Probability",main=mainlabel,mark.time=TRUE,col=c("black","red"),cex=1.2,cex.lab=1.5,cex.main=1.2,cex.axis=1.2,las=1,lwd=1)
text(30,0.17,"Expression",cex=1); 
legend(10,0.15,c("Low rish", "High rish"),lty =c(1,1),col=c("black","red"),bty="n")
#plot(kmfit2)

cFile<-as.matrix(clinData)
Survalue<-as.double(cFile[!is.na(cFile[,c("Survival")]),c("Survival")])
px<-max(Survalue)
if((px>100)&(px<1000)) px=px-30
if(px>1000) px=px-2000
if(px<100)  px=px-3
text(px-5,0.93,paste("Log Rank P=",p,sep=""),cex=1); 
dev.off()

#pdf(file="ROC.pdf")
#modelroc <- roc(factor(clinData$Survival)~as.numeric(gene_profile_log_zscore_contain_genes_plus_r))

#plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), grid.col=c("green", "red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
#dev.off()
#library(survivalROC)
#nobs <- NROW(clinData)
#Mayo= survivalROC(Stime=clinData$Survival,##生存时间
#                     status=clinData$Status,## 终止事件    
#                     marker = as.numeric(gene_profile_log_zscore_contain_genes_plus_r), ## marker value    
 #                    predict.time = 72,## 预测时间截点
#                     span = 0.25*nobs^(-0.20))##span,NNE法的namda
#str(Mayo)## list结构
## List of 6
##  $ cut.values  : num [1:313] -Inf 4.58 4.9 4.93 4.93 ...
##  $ TP          : num [1:313] 1 0.997 0.995 0.993 0.99 ...
##  $ FP          : num [1:313] 1 0.997 0.994 0.99 0.987 ...
##  $ predict.time: num 365
##  $ Survival    : num 0.929
##  $ AUC         : num 0.931
## 绘图
#plot(Mayo$FP, Mayo$TP, ## x=FP,y=TP
#     type="l",col="red", ##线条设置
 #    xlim=c(0,1), ylim=c(0,1),   
 #    xlab=paste( "FP", "\n", "AUC = ",round(Mayo$AUC,3)), ##连接
 #    ylab="TP",
 #    main="Mayoscore 4, Method = NNE \n  Year = 1")## \n换行符
#abline(0,1,col="gray",lty=2)##线条颜色
}
