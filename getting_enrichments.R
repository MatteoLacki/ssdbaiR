library(readr)
library(data.table)
# library(httr)

source('string_api.R')

ncbi = 318829

D = data.table(read_delim("string_test_results.csv", 
               ";", escape_double = FALSE, trim_ws = TRUE))

ID = map_identifiers(D$uniprot, ncbi, limit=1)
View(ID)
dim(ID) == length(D$uniprot)

ENRICH = get_functional_enrichment(ID$stringId, ncbi)

saveNetworkPlot(ID$stringId, ncbi, 'test.png')
# saveNetworkPlot(ID$stringId, ncbi, 'test.png', query_limit=20)

saveNetworkPlot(D$uniprot, ncbi, 'test2.png')
