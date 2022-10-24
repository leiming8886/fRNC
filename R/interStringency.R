#' @title Extract interactions according to stringency and interaction type
#'
#' @description interactions were extracted according to stringency and interaction type in the database of ENCORI

#' @param type	a character string indicating which interaction type is to be choosed
#' . One of "Protein" (RBP-circRNA, RBP-lncRNA, miRNA-circRNA, miRNA-lncRNA, RBP-miRNA, RBP-RBP), "transcription (miRNA-circRNA,miRNA-lncRNA,miRNA-RBP))": can be abbreviated

#' @param spec	a character string indicating which species is to be choosed
#' . One of "hg", "mm": can be abbreviated

#' @param stringency	a character string indicating which interaction stringency is to be choosed
#' . One of "low" ( number of supported experiments > = 1 or combined_score >= 150 or score >=0 ), "medium ( > = 2 or combined_score >= 400 or score >=5)","high ( > = 3 or combined_score >= 700 or score >=10)","strict ( > = 5 or combined_score >= 900, or score >=20)"
#'
#' @return Satisfactory interaction matrix. It contain five colnames:"node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name"

#' @examples
#' \dontrun{
#'  interac <- interStringency(type = "Protein",spec ="hg", stringency = "strict")
#'  interac <- interac[,c("node_gene_ID","target_gene_ID")]
#' }
#'
#' @export
interStringency <- function(type = c("Protein","transcription"),spec = c("hg","mm"),stringency = c("low","medium","high","strict") ){
  #data(dataM2C)
  #data(dataM2L)
  #data(dataM2R)
  #data(dataR2C)
  #data(dataR2L)
  spec <- match.arg(spec)
  specID <- NULL
  switch( spec,
          "hg"={
				specID <- "hg"
          },
          "mm"={
            specID <- "mm"
          },

          stop("Enter something that switches me!")
  )
  type <- match.arg(type)
  interaction <- NULL
  R2R <- NULL
  R2M <- NULL
  switch( type,
          "Protein"={
			if( specID == "hg"){
			load(system.file("extdata/", "hg_r2c.Rdata", package = "fRNC", mustWork = TRUE))#dataR2C
			load(system.file("extdata/", "hg_r2r.Rdata", package = "fRNC", mustWork = TRUE))#dataR2R
			load(system.file("extdata/", "hg_r2l.Rdata", package = "fRNC", mustWork = TRUE))#dataR2L
			load(system.file("extdata/", "hg_r2m.Rdata", package = "fRNC", mustWork = TRUE))#dataR2M
			#load(system.file("extdata/", "hg_m2r.Rdata", package = "fRNC", mustWork = TRUE))#dataM2R
			load(system.file("extdata/", "hg_m2l.Rdata", package = "fRNC", mustWork = TRUE))#dataM2L
			load(system.file("extdata/", "hg_m2c.Rdata", package = "fRNC", mustWork = TRUE))#dataM2C
            interaction <-rbind(dataR2C,dataR2L)
			R2R <- dataR2R
			R2M <- dataR2M
            #interaction <-rbind(interaction,dataM2R)
            interaction <-rbind(interaction,dataM2L)
			interaction <-rbind(interaction,dataM2C)
			}
			if( specID == "mm"){
			load(system.file("extdata/", "mm_r2c.Rdata", package = "fRNC", mustWork = TRUE))#dataR2C
			load(system.file("extdata/", "mm_r2r.Rdata", package = "fRNC", mustWork = TRUE))#dataR2R
			load(system.file("extdata/", "mm_r2l.Rdata", package = "fRNC", mustWork = TRUE))#dataR2L
			load(system.file("extdata/", "mm_r2m.Rdata", package = "fRNC", mustWork = TRUE))#dataR2M
			#load(system.file("extdata/", "hg_m2r.Rdata", package = "fRNC", mustWork = TRUE))#dataM2R
			load(system.file("extdata/", "mm_m2l.Rdata", package = "fRNC", mustWork = TRUE))#dataM2L
			load(system.file("extdata/", "mm_m2c.Rdata", package = "fRNC", mustWork = TRUE))#dataM2C
            interaction <-rbind(dataR2C,dataR2L)
			R2R <- dataR2R
			R2M <- dataR2M
            #interaction <-rbind(interaction,dataM2R)
            interaction <-rbind(interaction,dataM2L)
			interaction <-rbind(interaction,dataM2C)
			}
          },
          "transcription"={
			if( specID == "hg"){
			load(system.file("extdata/", "hg_m2r.Rdata", package = "fRNC", mustWork = TRUE))#dataM2R
			load(system.file("extdata/", "hg_m2l.Rdata", package = "fRNC", mustWork = TRUE))#dataM2L
			load(system.file("extdata/", "hg_m2c.Rdata", package = "fRNC", mustWork = TRUE))#dataM2C
            #interaction <-rbind(dataR2C,dataR2L)
			#R2R <- dataR2R
            interaction <-rbind(dataM2L,dataM2R)
            #interaction <-rbind(interaction,dataM2L)
			interaction <-rbind(interaction,dataM2C)
			}
			if( specID == "mm"){
			load(system.file("extdata/", "mm_m2r.Rdata", package = "fRNC", mustWork = TRUE))#dataM2R
			load(system.file("extdata/", "mm_m2l.Rdata", package = "fRNC", mustWork = TRUE))#dataM2L
			load(system.file("extdata/", "mm_m2c.Rdata", package = "fRNC", mustWork = TRUE))#dataM2C
            #interaction <-rbind(dataR2C,dataR2L)
			#R2R <- dataR2R
            interaction <-rbind(dataM2L,dataM2R)
            #interaction <-rbind(interaction,dataM2L)
			interaction <-rbind(interaction,dataM2C)
			}
          },
          stop("Enter something that switches me!")
  )


	if ( ! "clipExpNum" %in% colnames(interaction)){
		stop("clipExpNum must be in the col names!\n")
	}
	
	
	stringency <- match.arg(stringency)
	switch( stringency,
        "low"={
            interaction1 <- interaction[which(interaction$clipExpNum>=1),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
			if (! is.null(R2R)){
				interaction2 <- R2R[which(R2R$combined_score>=150),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				interaction1 <-rbind(interaction1,interaction2)
			
			}
			if (! is.null(R2M)){
				interaction3 <- R2M[which(R2M$score>=0),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				interaction1 <-rbind(interaction1,interaction3)
			}
        },
        "medium"={ 
				interaction1 <- interaction[which(interaction$clipExpNum>=2),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				if (! is.null(R2R)){
					interaction2 <- R2R[which(R2R$combined_score>=400),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
					interaction1 <-rbind(interaction1,interaction2)
				}
				if (! is.null(R2M)){
					interaction3 <- R2M[which(R2M$score>=5),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
					interaction1 <-rbind(interaction1,interaction3)
				}
        },
		"high"={ 
				interaction1 <- interaction[which(interaction$clipExpNum>=3),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				if (! is.null(R2R)){
					interaction2 <- R2R[which(R2R$combined_score>=700),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
					interaction1 <-rbind(interaction1,interaction2)
			
				}
				if (! is.null(R2M)){
					interaction3 <- R2M[which(R2M$score>=10),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
					interaction1 <-rbind(interaction1,interaction3)
				}
        },
		"strict"={
				interaction1 <- interaction[which(interaction$clipExpNum >= 5),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				if (! is.null(R2R)){
					interaction2 <- R2R[which(R2R$combined_score >=900),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
					interaction1 <-rbind(interaction1,interaction2)
			
				}
				if (! is.null(R2M)){
					interaction3 <- R2M[which(R2M$score>=20),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
					interaction1 <-rbind(interaction1,interaction3)
				}
        },
        stop("Enter something that switches me!")
    )

	return(interaction1)

}
