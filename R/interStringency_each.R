interStringency_each <- function(type = c("M2L","M2C","M2R","R2L","R2M","R2C","R2R"),spec = c("hg","mm"),stringency = c("low","medium","high","strict") ){
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

  switch( type,
        "M2L"={
			if( specID == "hg"){
			dataM2L <- read.table(system.file("extdata/hg/", "M2L.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataM2L
			}
			if( specID == "mm"){
			dataM2L <- read.table(system.file("extdata/mm/", "M2L.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataM2L
			}
			}
		,	
		"M2C"={
			if( specID == "hg"){
			dataM2C <- read.table(system.file("extdata/hg/", "M2C.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataM2C
			}
			if( specID == "mm"){
			dataM2C <- read.table(system.file("extdata/mm/", "M2C.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataM2C
			}
			}
		,	
		"M2R"={
			if( specID == "hg"){
			dataM2R <- read.table(system.file("extdata/hg/", "M2R.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataM2R
			}
			if( specID == "mm"){
			dataM2R <- read.table(system.file("extdata/mm/", "M2R.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataM2R
			}
			}
		,
		"R2L"={
			if( specID == "hg"){
			dataR2L <- read.table(system.file("extdata/hg/", "R2L.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2L
			}
			if( specID == "mm"){
			dataR2L <- read.table(system.file("extdata/mm/", "R2L.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2L
			}
			}	
		,	
		"R2M"={
			if( specID == "hg"){
			dataR2M <- read.table(system.file("extdata/hg/", "R2M.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2M
			}
			if( specID == "mm"){
			dataR2M <- read.table(system.file("extdata/mm/", "R2M.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2M
			}
			}
		,	
		"R2C"={
			if( specID == "hg"){
			dataR2C <- read.table(system.file("extdata/hg/", "R2C.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2C
			}
			if( specID == "mm"){
			dataR2C <- read.table(system.file("extdata/mm/", "R2C.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2C
			}
			}
		,	
		"R2R"={
			if( specID == "hg"){
			dataR2R <- read.table(system.file("extdata/hg/", "R2R.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2R
			}
			if( specID == "mm"){
			dataR2R <- read.table(system.file("extdata/mm/", "R2R.txt", package = "fRNC", mustWork = TRUE),header = T,sep="\t")
            interaction <-dataR2R
			}
			}
		,
		stop("Enter something that switches me!")
	)
	stringency <- match.arg(stringency)
	switch( stringency,
        "low"={
			score <- colnames(interaction)[6]
			if(score == "clipExpNum"){
				interaction1 <- interaction[which(interaction$clipExpNum>=1),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "combined_score"){
				interaction1 <- interaction[which(interaction$combined_score>=150),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "score"){
				interaction1 <- interaction[which(interaction$score>=0),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
        },
        "medium"={ 
			score <- colnames(interaction)[6]
			if(score == "clipExpNum"){
				interaction1 <- interaction[which(interaction$clipExpNum>=2),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "combined_score"){
				interaction1 <- interaction[which(interaction$combined_score>=400),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "score"){
				interaction1 <- interaction[which(interaction$score>=5),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
        },
		"high"={ 
			score <- colnames(interaction)[6]
			if(score == "clipExpNum"){
				interaction1 <- interaction[which(interaction$clipExpNum>=3),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "combined_score"){
				interaction1 <- interaction[which(interaction$combined_score>=700),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "score"){
				interaction1 <- interaction[which(interaction$score>=10),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
        },
		"strict"={
			score <- colnames(interaction)[6]
			if(score == "clipExpNum"){
				interaction1 <- interaction[which(interaction$clipExpNum>=5),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "combined_score"){
				interaction1 <- interaction[which(interaction$combined_score>=900),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}
			else if(score == "score"){
				interaction1 <- interaction[which(interaction$score>=20),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
				}		
        },
        stop("Enter something that switches me!")
    )

	return(interaction1)

}
