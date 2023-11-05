## Orthology functions

#' Build orthology table
#' 
#' This is the main function to call for building orthologies between any set of species
#'
#' @param taxIDs A named vector of all species with associated taxonomy IDs (names are common names to show up in columns, and values are NCBI taxa IDs)
#' @param primaryTaxID A numeric vector of one or more species to compare against. Default is 9606 (human). These all MUST be included in the taxIDs. If more than one are provided, any additional species are only used to fill in "NA" slots in the initial table, in the order provided.  Possible options include 7955 (zebrafish), 9031 (chicken), 9606 (human), 9615 (dog), 9685 (cat), 9823 (pig), 9913 (cow), 10090 (mouse), and 10116 (rat), since (as of 5 November 2023), these are the only species listed as primary species in the NCBI orthology table.
#' @param geneOrthologs The gene_orthologs.gz file updated daily. A local file in the same format can be provided, but otherwise leave the default link as is.
#' @param addEnsemblID Should EnsemblIDs be added for each geneID?
#' @param addGeneInfoFirst Gene info to include with each gene for the first primary species. To skip gene info, set to NULL.  Sensible defaults are included, but all options include: "Symbol", "LocusTag", "Synonyms", "dbXrefs", "chromosome", "map_location", "description", "type_of_gene", "Symbol_from_nomenclature_authority", "Full_name_from_nomenclature_authority", "Nomenclature_status", "Other_designations", "Modification_date", "Feature_type"
#' @param addGeneInfoAll Gene info to include with each gene for every other species (default is "Symbol"). Same options as above are allowed. Ignored if addGeneInfoFirst is set to NULL.
#' @param includeNonMammalianSpecies Are any non-mammalian species used (default is FALSE). This is only used to determine which gene info file to download. The mammalian file is ~190MB, while the all species file in ~1GB.
#' @param verbose Report status (default FALSE). Note that fread may report status of large file downloads regardless of TRUE/FALSE setting.
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
  addGeneInfoFirst = c("Symbol","type_of_gene","description"),
  addGeneInfoAll   = c("Symbol"),
  includeNonMammalianSpecies = FALSE,
  verbose          = FALSE,
  outputFilePrefix = "ortholog_table",   
  returnTable      = FALSE)
{
  ## Function prep and checks
  library(data.table)
  primaryTaxID <- intersect(primaryTaxID,taxIDs)
  if(length(primaryTaxID)==0) stop(paste("Please include",primaryTaxID,"in the taxIDs vector."))
  if(length(taxIDs)<2) stop("Please include at least two species in taxIDs")
  w = which(taxIDs==primaryTaxID[1])
  taxIDs <- taxIDs[c(w,setdiff(1:length(taxIDs),w))]
  additionalPrimaryTaxIDs = NULL
  if(length(primaryTaxID)>1) {
    additionalPrimaryTaxIDs = primaryTaxID[2:length(primaryTaxID)]
    primaryTaxID = primaryTaxID[1]
  }

  ## Download and read the current orthologs table 
  if(verbose) print("Downloading and subsetting NCBI orthology table...")
  orthologs <- fread(geneOrthologs)
  orthologs <- orthologs[,c(1,2,4,5)]  # Omit unnecessary "relationship column
  colnames(orthologs) <- c("tax_id_1","GeneID_1","tax_id_2","GeneID_2")  # Edit column names
  
  ## Subset to only include rows where both gene IDs are in an inputted tax ID
  orthologs <- orthologs[is.element(orthologs$tax_id_1,taxIDs)&is.element(orthologs$tax_id_2,taxIDs),]
  
  ## Report missing taxIDs
  missingIDs <- setdiff(taxIDs,orthologs$tax_id_2)
  print("The following species are not included in NCBI's gene ortholog table and will be omitted:")
  print(paste(names(taxIDs)[is.element(taxIDs,missingIDs)]))
  taxIDs <- taxIDs[!is.element(taxIDs,missingIDs)]
  
  ## Put the primary ortholog species first
  if(verbose) print("Building the output orthology table...")
  if(verbose) print("...putting first primary taxonomy")
  orths <- NULL
  for (tax2 in setdiff(taxIDs,primaryTaxID)){
    orthTmp <- homologeneOneToOne(primaryTaxID,tax2,orthologs)
    if(!is.na(as.character(orthTmp[1])[1]))
      orths <- rbind(orths,orthTmp)
  }
  genes <- sort(unique(orths$GeneID_1))

  ## Build the table
  if(verbose) print("...building the table")
  orthologTable <- as.data.frame(matrix(NA,nrow=length(genes),ncol=length(taxIDs)))
  colnames(orthologTable) <- names(taxIDs)
  orthologTable[,1] <- genes
  for (i in 2:length(taxIDs)){
    orthTmp <- orths[orths$tax_id_2==taxIDs[i],]
    orthologTable[,i] <- orthTmp[match(orthologTable[,1],orthTmp$GeneID_1),"GeneID_2"]
  }
  colnames(orthologTable) <- paste0(colnames(orthologTable),"_geneid")

  ## Address additional primary taxonomies, if any
  if(!is.null(additionalPrimaryTaxIDs)){
    if(verbose) print("...addressing additional primary taxonomies")
    for (i in 1:length(additionalPrimaryTaxIDs)){
      orths <- NULL
      for (tax2 in setdiff(taxIDs,additionalPrimaryTaxIDs[i])){
        orthTmp <- homologeneOneToOne(additionalPrimaryTaxIDs[i],tax2,orthologs)
        if(!is.na(as.character(orthTmp[1])[1]))
          orths <- rbind(orths,orthTmp)
      }
      column <- which(taxIDs==additionalPrimaryTaxIDs[i])
      # Add new gene rows
      newGenes <- sort(unique(setdiff(orths$GeneID_1,orthologTable[,column])))
      if(length(newGenes)>0){
        newRows <- orthologTable[1:length(newGenes),,drop=FALSE]*NA
        newRows[,column] <- newGenes
        colnames(newRows) <- colnames(orthologTable)
        orthologTable <- rbind(orthologTable,newRows)
      }
      # Replace any NAs, one species at a time
      for (j in setdiff(1:length(taxIDs),column)){
        orthTmp <- orths[orths$tax_id_2==taxIDs[j],c(2,4)]
        if(length(orthTmp)>0){
          orthTmp = orthTmp[match(orthologTable[,column],orthTmp$GeneID_1),]
          kp <- is.na(orthologTable[,j])
          orthologTable[kp,j] <- orthTmp[kp,2]
          print(paste("...added",sum(!is.na(orthTmp[kp,2])),"orthologs for",names(taxIDs)[column],"and",names(taxIDs)[j]))
        }
      }
    }
  }
  
  
  ## Add Ensembl information
  if (addEnsemblID){
    if(verbose) print("Adding Ensembl IDs...")
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
  if(!is.null(addGeneInfoFirst)){
    if(verbose) print("Adding requested gene info...")
    #if(primaryTaxID==9606){  # If we go back to reporting for only primary species, we can speed up the code by putting this back
    #  geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz")
    #} else if(primaryTaxID==10090){
    #  geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz")
    #} else {
    #  geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz")
    #  geneInfo <- geneInfo[geneInfo$`#tax_id`==primaryTaxID,]
    #}
    if(includeNonMammalianSpecies){
      if(verbose) print("...non-mammalian species indicated, so this step may be very slow...")
      geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz") # a 1GB file!
    } else {
      geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz")
    }
    
    addGeneInfoFirst = intersect(addGeneInfoFirst,colnames(geneInfo))
    if(length(addGeneInfoFirst)>0){
      for (i in 1){
        geneInfoTmp <- geneInfo[geneInfo$`#tax_id`==taxIDs[i],]
        geneInfoTmp <- as.data.frame(geneInfoTmp)[match(orthologTable[,i],geneInfoTmp$GeneID),addGeneInfoFirst,drop=FALSE]
        colnames(geneInfoTmp) <- paste(names(taxIDs)[i],colnames(geneInfoTmp),sep="_")
        orthologTable <- cbind(orthologTable,geneInfoTmp)
      }
    }
    addGeneInfoAll   = intersect(addGeneInfoAll,colnames(geneInfo))
    if(length(addGeneInfoAll)>0){
      for (i in 2:length(taxIDs)){
        geneInfoTmp <- geneInfo[geneInfo$`#tax_id`==taxIDs[i],]
        geneInfoTmp <- as.data.frame(geneInfoTmp)[match(orthologTable[,i],geneInfoTmp$GeneID),addGeneInfoAll,drop=FALSE]
        colnames(geneInfoTmp) <- paste(names(taxIDs)[i],colnames(geneInfoTmp),sep="_")
        orthologTable <- cbind(orthologTable,geneInfoTmp)
      }
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
    if(verbose) print(paste("Outputting file to",fout))
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

