# GeneOntology

### Overview

GeneOntology is an R package and associated csv file for creating snapshots of the NCBI Gene ontology table (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz) for use with [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells) and other purposes.

### How to use the package

#### Downloading precomputed tables

Snapshots created 3 November 2023
* Gene conversions between mouse and human: mouse_human_orthologs_11022023.csv
* Gene conversions between 26 mammalian species (see below): mammalian_orhologs_11022023

#### Using the R package

Install the package:
```
install.packages("remotes", repos='http://cran.us.r-project.org')
remotes::install_github("AllenInstitute/GeneOrthology", build_vignettes = TRUE)
```

Open R and follow the vignette:
```
library(GeneOrthology)
browse_vignettes("GeneOrthology")
```

The vignette can also be viewed in HTML format at this link.

### Mammalian species list

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




### License

The license for this package is available on Github at: https://github.com/AllenInstitute/allen_institute_nomenclature/blob/master/LICENSE

### Level of Support

We do not anticipate updates to this tool, so long as NCBI retains current format of ontology table.

### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/allen_institute_nomenclature/blob/master/CONTRIBUTION
