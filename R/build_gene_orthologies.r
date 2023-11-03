## Orthology functions

#' Build orthology table
#' 
#' This is the main function to call for building orthologies between any set of species
#'
#' @param taxIDs A named vector of all species with associated taxonomy IDs (names are common names to show up in columns, and values are NCBI taxa IDs)
#' @param primaryTaxID A numeric vector of the species to compare against. This MUST be included in the taxIDs. Default is 9606 (human)
#' @param geneOrthologs The gene_orthologs.gz file updated daily. A local file in the same format can be provided, but otherwise leave the default link as is.
#' @param addEnsemblID Should EnsemblIDs be added for each geneID?
#' @param addGeneInfo Gene info to include with each gene for the primary taxonomy (NULL for none). For species other than mouse and human, this can be very slow.  Currently this is only compatible with mammalial species.  Sensible defaults are included, but all options include: "Symbol", "LocusTag", "Synonyms", "dbXrefs", "chromosome", "map_location", "description", "type_of_gene", "Symbol_from_nomenclature_authority", "Full_name_from_nomenclature_authority", "Nomenclature_status", "Other_designations", "Modification_date", "Feature_type"
#' @param outputFilePrefix Prefix for outputted file name (or NULL to not output file). File will be outputFilePrefix_[date].csv.
#' @param returnTable Should the table be returned?
#'
#' @return Table of gene orthology information including EntrezID, Ensembl ID, and gene symbol for ever available species.
#'
#' @import data.table
#'
#' @export

build_orthology_table <- function(
  taxIDs           = setNames(c(9606,9483,10090),c("human","marmoset","mouse")),
  primaryTaxID     = 9606,
  geneOrthologs    = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",  
  addEnsemblID     = TRUE,
  addGeneInfo      = c("Symbol","description","type_of_gene"),
  outputFilePrefix = "ortholog_table",                                            
  returnTable      = TRUE)
{
  ## Function prep and checks
  library(data.table)
  if(!is.element(primaryTaxID,taxIDs)) stop(paste("Please include",primaryTaxID,"in the taxIDs vector."))
  if(length(taxIDs)<2) stop("Please include at least two species in taxIDs")
  w = which(taxIDs==primaryTaxID)
  taxIDs <- taxIDs[c(w,setdiff(1:length(taxIDs),w))]

  ## Download and read the current orthologs table 
  orthologs <- fread(geneOrthologs)
  orthologs <- orthologs[,c(1,2,4,5)]  # Omit unnecessary "relationship column
  colnames(orthologs) <- c("tax_id_1","GeneID_1","tax_id_2","GeneID_2")  # Edit column names
  
  ## Subset to only include rows where both gene IDs are in an inputted tax ID
  orthologs <- orthologs[is.element(orthologs$tax_id_1,taxIDs)&is.element(orthologs$tax_id_2,taxIDs),]
  
  ## Put the primary ortholog species first
  orths <- NULL
  for (tax2 in setdiff(taxIDs,primaryTaxID)){
    orthTmp <- homologeneOneToOne(primaryTaxID,tax2,orthologs)
    if(!is.na(as.character(orthTmp[1])[1]))
      orths <- rbind(orths,orthTmp)
  }
  genes <- sort(unique(orths$GeneID_1))

  ## Build the table
  orthologTable <- as.data.frame(matrix(NA,nrow=length(genes),ncol=length(taxIDs)))
  colnames(orthologTable) <- names(taxIDs)
  orthologTable[,1] <- genes
  for (i in 2:length(taxIDs)){
    orthTmp <- orths[orths$tax_id_2==taxIDs[i],]
    orthologTable[,i] <- orthTmp[match(orthologTable[,1],orthTmp$GeneID_1),"GeneID_2"]
  }
  colnames(orthologTable) <- paste0(colnames(orthologTable),"_geneid")

  ## Add Ensembl information
  if (addEnsemblID){
    convertEns = fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz")
    EnsemblTable = orthologTable
    for(i in 1:length(taxIDs)){
      convertTmp       <- convertEns[convertEns$`#tax_id`==taxIDs[i],2:3]
      EnsemblTable[,i] <- convertTmp[match(orthologTable[,i],convertTmp$GeneID),]$"Ensembl_gene_identifier"
    }
    colnames(EnsemblTable) <- gsub("geneid","EnsemblID",colnames(EnsemblTable))
    orthologTable <- cbind(orthologTable,EnsemblTable)
  }
  
  ## Add Gene Info, if not set to NULL
  if(!is.null(addGeneInfo)){
    if(primaryTaxID==9606){
      geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz")
    } else if(primaryTaxID==10090){
      geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz")
    } else {
      geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz")
      geneInfo <- geneInfo[geneInfo$`#tax_id`==primaryTaxID,]
    }
    
    addGeneInfo = intersect(addGeneInfo,colnames(geneInfo))
    if(length(addGeneInfo)>0){
      geneInfo <- as.data.frame(geneInfo)[match(orthologTable[,1],geneInfo$GeneID),addGeneInfo]
      colnames(geneInfo) <- paste(names(taxIDs)[1],colnames(geneInfo),sep="_")
      orthologTable <- cbind(orthologTable,geneInfo)
    }
  }
  
  ## Reorder columns
  cn <- colnames(orthologTable)
  cn <- c(cn[grepl(names(taxIDs)[1],cn)],sort(cn[!grepl(names(taxIDs)[1],cn)]))
  orthologTable <- orthologTable[,cn]
  
  ## Output the table if requested
  if(!is.null(outputFilePrefix)){
    date_suffix = gsub("-","",Sys.Date())
    fout <- paste0(outputFilePrefix,"_",date_suffix,".csv")
    fwrite(orthologTable,fout,na="NA")
  } 

  ## Return the table if requested
  if(returnTable)
    orthologTable
} 


#' Internal function for finding one-to-one orthologs between a pair of species
#' 
#' @param tax1 the primary taxonomy ID
#' @param tax2 the secondary taxonomy ID
#' @param orthologs An orthology table of gene IDs and taxonomy IDs
#'
#' @return a potentially reordered and filtered taxonomy table
#'
homologeneOneToOne <- function(tax1, tax2, orthologs){
  # Subset to gene ids from desired species
  ortholog2 <- subset(orthologs, (tax_id_1 == tax1 & tax_id_2 == tax2) | (tax_id_1 == tax2 & tax_id_2 == tax1))
  
  # Return NA if there are no homologs
  if(dim(ortholog2)[1]==0) return(NA)
  
  # Put all of tax1 values in column 1 if needed
  kp <- ortholog2$tax_id_1==tax2
  if(length(kp)>0)
    ortholog2[kp,] = ortholog2[kp,c(3,4,1,2)]
  
  # Delete any orthologs that don't have one-to-one matches
  kp1 <- table(as.character(ortholog2$GeneID_1))
  kp1 <- is.element(as.character(ortholog2$GeneID_1),names(kp1)[kp1>1])
  kp2 <- table(as.character(ortholog2$GeneID_2))
  kp2 <- is.element(as.character(ortholog2$GeneID_2),names(kp1)[kp2>1])
  kp  <- kp1|kp2 
  if(sum(kp)>0)
    ortholog2 <- ortholog2[!kp,]
  
  # Return the table
  ortholog2
}

