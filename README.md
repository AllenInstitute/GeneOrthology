# GeneOntology

## Overview

GeneOntology is an R package and associated csv file for creating snapshots of the NCBI Gene ontology table (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz) for use with **[MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells)** and other purposes.  In addition, Ensembl IDs and other gene info are pulled from other tables from the [NCBI Gene page](https://www.ncbi.nlm.nih.gov/gene) FTP server.  These NCBI files are all updated daily and so the date is appended to all outputted files in the format YYYYMMDD.  Unless otherwise noted, all ortholog tables provided in csv format are anchored to human genes.

## Downloading precomputed tables

Click on this button on the right side of the screen after you click on one of the links below to download to precomputed orthology table:

![image](https://github.com/AllenInstitute/GeneOrthology/assets/25486679/3d176b70-70f1-4a09-b5d4-741b4ea714e3)

**Snapshots created on 13 November 2023.**  These files contain human gene symbols (and other info), with Ensembl IDs and NCBI gene IDs for every species. 
* Gene conversions between mouse, human, marmoset, and macaque (rhesus): **[mouse_human_marmoset_macaque_orthologs_20231113.csv](https://github.com/AllenInstitute/GeneOrthology/blob/main/csv/mouse_human_marmoset_macaque_orthologs_20231113.csv)**
* Gene conversions between 27 mammalian species, anchored to all available mammalian species (see below): **[mammalian_orthologs_20231113](https://github.com/AllenInstitute/GeneOrthology/blob/main/csv/mammalian_orthologs_20231113.csv).** (Note that we include dog as a stand-in for coyote and Vaquita as a stand-in for porpoise.)

## Using the R package

Install the package:
```
install.packages("remotes", repos='http://cran.us.r-project.org')
remotes::install_github("AllenInstitute/GeneOrthology")
```

Run the code for a simple example to create the mouse/human/marmoset/macaque orthology table above. This is the most common use case.  
```
library(GeneOrthology)
taxIDs <- setNames(c(9606,10090,9483,9544),
                   c("human","mouse","marmoset","rhesus_macaque"))
build_orthology_table(taxIDs = taxIDs,  primaryTaxID = 9606, 
                      outputFilePrefix="mouse_human_marmoset_macaque_orthologs")
```

Here is a bit more complicated example showing how to create the 27 mammal (+zebrafish) table above, which is anchored in multiple species and includes non-mammals.  
```
library(GeneOrthology)
taxIDs <- setNames(c(9669, 246437, 10116, 13616, 27679, 
                     9823, 9361, 9986, 60711, 9598, 
                     30608, 10181, 37293, 9545, 9544, 
                     10090, 30611, 9685, 9595, 9606, 
                     9483, 42100, 9515, 9614, 9407, 
                     9555, 9999, 9615, 7955),
                   c("Ferret", "Chinese.Treeshrew", "Rat", "Opossum", "Squirrel.monkey", 
                     "Pig", "Armadillo.Nine.banded", "Rabbit", "African.green.monkey", "Chimpanzee", 
                     "Mouse.lemur", "Naked.mole.rat", "Owl.monkey", "Macaque.pig.tailed", "Macaque.rhesus", 
                     "Mouse", "Galago", "Cat", "Gorilla", "Human", 
                     "Marmoset", "Vaquita", "Tufted.capuchin", "Coyote", "Egyptian.fruit.bat", 
                     "Olive.baboon", "Squirrel.arctic.ground","dog","zebrafish"))
build_orthology_table(taxIDs = taxIDs, primaryTaxID = c(9606,10090,10116,9615,9685,9823,9913,7955),  
                      outputFilePrefix="mammalian_orthologs",verbose=TRUE,
                      includeNonMammalianSpecies = TRUE)  # To include zebrafish gene symbols, but much slower.
# A few of these species do not have human homologies computed and are omitted from the output.
# It's also worth noting that anchoring to all of these additional species ONLY adds a total of
#   ~100 orthology pairs and probably is not necessary. 
```

## Mammalian species list

#### Current list

This list includes all mammals currently studied at the Allen Institute for Brain Science (as of 8 November 2023), as well as related species with NCBI orthologs to human:
|English Name|Species|NCBI TaxID|
| ----- | ----- | ----- |
|African green monkey|Chlorocebus sabaeus|60711|
|Armadillo (Nine-banded)|Dasypus novemcinctus|9361|
|Cat|Felis catus|9685|
|Chimpanzee|Pan troglodytes|9598|
|Chinese tree shrew|Tupaia chinensis|246437|
|Coyote|Canis latrans|9614|
|Dog|Canis lupus familiaris|9615|
|Egyptian fruit bat|Rousettus aegyptiacus|9407|
|Ferret|Mustela putorius furo|9669|
|Galago|Otolemur garnettii|30611|
|Gorilla|Gorilla gorilla gorilla|9595|
|Harbor porpoise|Phocoena phocoena|9742|
|Human|Homo sapiens|9606|
|Macaque (pig-tailed)|Macaca nemestrina|9545|
|Macaque (rhesus)|Macaca mulatta|9544|
|Marmoset|Callithrix jacchus|9483|
|Mouse lemur|Microcebus murinus|30608|
|Mouse|Mus musculus|10090|
|Naked mole rat|Heterocephalus glaber|10181|
|Olive baboon|Papio anubis|9555|
|Opossum|Monodelphis domestica|13616|
|Owl monkey|Aotus nancymaae|37293|
|Pig|Sus scrofa|9823|
|Rabbit|Oryctolagus cuniculus|9986|
|Rat|Rattus norvegicus|10116|
|Squirrel (Arctic ground)|Urocitellus parryii|9999|
|Squirrel monkey|Saimiri boliviensis|27679|
|Tufted capuchin|Sapajus apella|9515|
|Vaquita|Phocoena sinus|42100|

Available speciesfrom this list are included in the downloadable csv file.

#### Finding additional species with orthologs

Any species supported by NCBI can be included in the build_orthology_table function (mammalian or otherwise).  Currently the vast majority of orthologs are matched against either human or zebrafish.  Taxonomy IDs for other species can be found on the **[NCBI taxonomy website](https://www.ncbi.nlm.nih.gov/taxonomy)** manually.  

The **[taxonomizr R library](https://github.com/sherrillmix/taxonomizr/)** also provides a convenient wrapper for this information in R.  This script shows how to get additional information about taxa with NCBI orthologs and to search for species of interest.

```
## Install and load the taxonomizr
if(!require("taxonomizr", quietly = TRUE)) 
  install.packages("taxonomizr")
library(taxonomizr)
library(data.table)

## Find all available NCBI taxids with orthologs
orthologs <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz")
taxids    <- sort(unique(orthologs$Other_tax_id))

## Get the taxonomy lineage from taxonomizr
database  <- prepareDatabase(getAccessions=FALSE)
taxa      <- getTaxonomy(taxids,database)
print(head(taxa))

## Get common names from taxonomizr (for human and mouse)
getCommon(c(9606,10090),database)

## See which species with NCBI orthologs are a type of "squirrel monkey"
commons  <- getCommon(taxids,database)
isInNCBI <- unlist(lapply(commons,function(x) sum(grepl("squirrel monkey",x[,1]))))>0
getCommon(taxids[isInNCBI],database)
```

## Contributions and updates

#### License

The license for this package is available on Github at: https://github.com/AllenInstitute/GeneOrthology/blob/master/LICENSE

#### Level of Support

We do not anticipate updates to this tool, so long as NCBI retains current format of ontology table.  That said, we encourage submission of issues.

#### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/GeneOrthology/blob/master/CONTRIBUTION

#### Comments, issues, or suggestions?

Please make direct pull requests, submit issues, or contact [Jeremy Miller](mailto:jeremym@alleninstitute.org) with any input.
