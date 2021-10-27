source("./scripts/ABC-createRefDB.R")

myDB <- dbAddProtein(myDB, jsonlite::fromJSON("./myScripts/MBP1_CUTOL.json"))
myDB <- dbAddTaxonomy(myDB, jsonlite::fromJSON("./myScripts/CUTOLtaxonomy.json"))
