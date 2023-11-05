# GeneOntology

## Overview

GeneOntology is an R package and associated csv file for creating snapshots of the NCBI Gene ontology table (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz) for use with **[MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells)** and other purposes.  In addition, Ensembl IDs and other gene info are pulled from other tables from the [NCBI Gene page](https://www.ncbi.nlm.nih.gov/gene) FTP server.  These NCBI files are all updated daily and so the date is appended to all outputted files in the format YYYYMMDD.  Unless otherwise noted, all ortholog tables provided in csv format are anchored to human genes.

## Downloading precomputed tables

Click on this button on the right side of the screen after you click on one of the links below to download to precomputed orthology table:
![image](https://github.com/AllenInstitute/GeneOrthology/assets/25486679/3d176b70-70f1-4a09-b5d4-741b4ea714e3)

**Snapshots created on 3 November 2023.**  These files contain human gene symbols (and other info), with Ensembl IDs and NCBI gene IDs for every species. 
* Gene conversions between mouse and human: **[mouse_human_orthologs_20231103.csv](https://github.com/AllenInstitute/GeneOrthology/blob/main/csv/mouse_human_orthologs_20231103.csv)**
* Gene conversions between 27 mammalian species (see below): **[mammalian_orthologs_20231103](https://github.com/AllenInstitute/GeneOrthology/blob/main/csv/mammalian_orthologs_20231103.csv).** (Note that a few of these species do not have human homologies computed and therefore have NA in every row.)

## Using the R package

Install the package:
```
install.packages("remotes", repos='http://cran.us.r-project.org')
remotes::install_github("AllenInstitute/GeneOrthology", build_vignettes = TRUE)
```

Run the code
```
# Load the library
library(GeneOrthology)

# Create the mouse/human orthology table
build_orthology_table(taxIDs = setNames(c(9606,10090),c("human","mouse")), 
                      primaryTaxID = 9606, outputFilePrefix="mouse_human_orthologs",returnTable=FALSE)

# Create the 27 mammal table (anchored in human)
taxIDs <- setNames(c(9669, 37347, 10116, 13616, 39432, 
                     9823, 9361, 9986, 60711, 9598, 
                     30608, 10181, 37293, 9545, 9544, 
                     10090, 30611, 9685, 9595, 9606, 
                     9483, 9742, 9516, 9614, 9407, 
                     9555, 9999),
                     c("Ferret", "Treeshrew", "Rat", "Opossum", "Squirrel.monkey", 
                     "Pig", "Armadillo.Nine.banded", "Rabbit", "African.green.monkey", "Chimpanzee", 
                     "Mouse.lemur", "Naked.mole.rat", "Owl.monkey", "Macaque.pig-tailed", "Macaque.rhesus", 
                     "Mouse", "Galago", "Cat", "Gorilla", "Human", 
                     "Marmoset", "Harbor.porpoise", "Capuchin", "Coyote", "Egyptian.fruit.bat", 
                     "Olive.baboon", "Squirrel.Arctic ground."))
build_orthology_table(taxIDs = taxIDs, primaryTaxID = 9606, outputFilePrefix="mammalian_orthologs",returnTable=FALSE)
# Note that a few of these species do not have human homologies computed and therefore have NA in every row.
```

## Mammalian species list

This list includes all mammals currently studied at the Allen Institute for Brain Science (as of 3 November 2023):
|English Name|Species|NCBI TaxID|
| ----- | ----- | ----- |
|African green monkey|Chlorocebus sabaeus|60711|
|Armadillo (Nine-banded)|Dasypus novemcinctus|9361|
|Capuchin|Sapajus apella|9516|
|Cat|Felis catus|9685|
|Chimpanzee|Pan troglodytes|9598|
|Coyote|Canis latrans|9614|
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
|Squirrel monkey|Saimiri boliviensis boliviensis|39432|
|Treeshrew|Tupaia belangeri|37347|

This species list is included in the downloadable file, but **any species supported by NCBI can be entered into this function**.  Taxonomy IDs for other species can be found on NCBI: https://www.ncbi.nlm.nih.gov/taxonomy.  The taxonomizr R library provides a convenient wrapper for this information (https://github.com/sherrillmix/taxonomizr/).  Adding gene information is only currently supported for mammals, but please let me know if you'd like this to be extended to other species.

## Contributions and updates

#### License

The license for this package is available on Github at: https://github.com/AllenInstitute/GeneOrthology/blob/master/LICENSE

#### Level of Support

We do not anticipate updates to this tool, so long as NCBI retains current format of ontology table.  That said, we encourage submission of issues.

#### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/GeneOrthology/blob/master/CONTRIBUTION

#### Comments, issues, or suggestions?

Please make direct pull requests, submit issues, or contact [Jeremy Miller](mailto:jeremym@alleninstitute.org) with any input.
