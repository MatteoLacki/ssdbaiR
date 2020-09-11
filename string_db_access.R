# install.packages('httr')
library(httr)


get_string_db_data_ = function(x, method, format='tsv', ...){
  http = file.path("https://string-db.org/api", format, method, fsep='/')
  body = list(identifiers=paste(x, collapse='\r'), ...)
  data = POST(http, body=body)
  Sys.sleep(1)
  data
}


getNetworkPlot = function(identifiers, path, highres=FALSE, query_limit=2000, ...){
  img_type = tools::file_ext(path)
  if(!(img_type %in% c('svg','png'))) stop('Unsupported type of image.')
  if(img_type == 'png'){
    if(highres) img_type = 'highres_image'
    else img_type = 'image'
  }
  if(length(identifiers) > query_limit) stop('Too many proteins for an API query. Consider subsets.')
  DATA = get_string_db_data_(identifiers, method='network', format=img_type, ...)
  if(length(DATA$content) > 0) writeBin(DATA$content, path) else print("empty result!")
  invisible()
}


#` Generate a sequence of N assignments into a minimal number of groups of similar size below K.
split_below = function(N, K) sort(suppressWarnings(rep.int(0, N) + 1:(N %/% K + 1)))

#` Download data from the STRING-DB server.
#`
#` Detail on \url{https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers}
get_string_db_data = function(identifiers, method, format='tsv', query_limit=2000, ...)
  lapply(
    split(identifiers, 
          split_below(length(identifiers), query_limit)),
    get_string_db_data_,
    method=method,
    format=format,
    ...)


content2table = function(post_data) fread(readBin(post_data$content, 'character'))

#` Download interaction partners from STRING-DB for a given set of identifiers.
#`
#` @params identifiers Protein identifiers.
#` @params species Which species to query?
#` @params query_limit The limit on the individual query on STRING-DB. As of April 2020, it's 2000 proteins.
#` @params ... other parameters to the API, see \url{https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers}
get_interaction_partners = function(identifiers, species, query_limit=2000, ...){
  DATA = get_string_db_data(identifiers,
                            method='interaction_partners',
                            format='tsv', 
                            query_limit=query_limit,
                            species=species,
                            ...)
  unique(rbindlist(lapply(DATA, content2table)))
}


get_string_ids = function(species, accessions, entries, query_limit=2000, ...){
  acc_ent = c(accessions, entries)
  acc_ent_dt = data.table(accession=accessions, entry=entries)
  
  # there is no limit on how many strings we query here!
  combined_query = content2table(get_string_db_data_(acc_ent,
                                                     'get_string_ids',
                                                     'tsv', 
                                                     # echo_query=1,
                                                     # limit=1,
                                                     species=species,
                                                     ...))
  
  combined_query[,origin:=ifelse(queryItem %in% accessions, 'A', 'E')]
  
  # the code below is typical R inferno:
  # the problem is that we find a protein either by accession or by entry
  # and both possibilities are valid
  # and an accession that is tied to the same entry can point to a different stringId
  # and there can be several stringIds for the same accession or the same entry.
  # Welcome to the Wild Wild West of Proteomics, hahahaha!
  
  # A | A & E
  res = combined_query[origin == 'A']
  setnames(res, 'queryItem','accession')
  res = merge(res, acc_ent_dt, all.y=F)

  # check if we have all:
  out = combined_query[, .(what=paste0(sort(origin), collapse='')), .(stringId)]
 
  # E \ A & E
  E_only = combined_query[origin == 'E' & !(queryItem %in% res$entry)]
  setnames(E_only, 'queryItem','entry')
  E_only = merge(E_only, acc_ent_dt, all.y=F)
  res = rbind(res, E_only)
  
  missing = setdiff(paste(accessions, entries, sep='.'), paste(res$accession, res$entry, sep='.'))
  list(found=res, missing=missing, origin_stats = table(out$what))
}

# identifiers = Z1@accession
FD_050 = FD_maxOcc_short[folderNo=='2019-050']
res = get_string_ids(species=10090, FD_050$accession, FD_050$entry)

stringsId = unique(res$found$stringId)
length(stringsId)
N0 = get_interaction_partners(stringsId, species=10090)
# N1  = get_interaction_partners(stringsId, species=10090, add_nodes=1)
# N10 = get_interaction_partners(stringsId, species=10090, add_nodes=10)

Zubidubi = res$found[FD_050[,.(accession, LTtest_pval, Ftest_pval_FDR)], on='accession']
Zubidubi[,origin:=NULL]
Zubidubi = Zubidubi[,.(stringId, LTtest_pval, Ftest_pval_FDR)]
Zubidubi[,stringId:=str_sub(stringId, 7)]

N0_pval = merge(N0, Zubidubi, by.x ='stringId_A', by.y='stringId')
# N0 = merge(N0, Zubidubi, by.x ='stringId_B', by.y='stringId')

degrees = N0[,.N,stringId_A]
Zubidubi = merge(Zubidubi, degrees, by.x = 'stringId', by.y='stringId_A')

plot( Zubidubi[LTtest_pval < .05]$LTtest_pval, log2(Zubidubi[LTtest_pval < .05]$N))
plot( Zubidubi$LTtest_pval, log2(Zubidubi$N))
plot(-log2(Zubidubi$LTtest_pval), log2(Zubidubi$N))
plot(-log2(Zubidubi$Ftest_pval_FDR), log2(Zubidubi$N))

ggplot(Zubidubi, aes(-log2(Ftest_pval_FDR), log2(N))) +
  geom_point() +
  geom_smooth()

ggplot(Zubidubi, aes(Ftest_pval_FDR, N)) +
  geom_point() +
  geom_smooth()

ggplot(Zubidubi, aes(LTtest_pval, log2(N))) +
  geom_point() +
  geom_smooth()


# accessions = Z1[pval_fdr < .05]$accession
# all_accessions = identifiers = Z1$accession
# 
# Z5perc = Z1[pval_fdr < .05]
# acc5perc = split(Z5perc$accession, Z5perc$folderNo)
# 
# acc = split(Z1$accession, Z1$folderNo)
# GG = lapply(acc, get_interaction_partners, species=10090)
# 
# getNetworkPlot(accessions,
#                'res/LTtest_pval_5perc/net.png',
#                species=10090,
#                hide_disconnected_nodes=1,
#                block_structure_pics_in_bubbles=1)
# 
# # getting multiple figures.
# library(purrr)
# 
# map2(acc5perc, 
#      file.path('res/LTtest_pval_5perc', paste0(names(acc), '.png')),
#      getNetworkPlot, 
#      species=10090,
#      hide_disconnected_nodes=1,
#      block_structure_pics_in_bubbles=1)
# 
# Z5perc[,.N,folderNo]
# 
# ncbi_id = 10090
# G = get_interaction_partners(all_accessions, species=10090)
# 
# get_interaction_partners(all_accessions, species=10090)
# 
# 
# GGG = graph_from_edgelist(as.matrix(G[,.(stringId_A, stringId_B)]))
# plot(GGG)
# write(accessions, file='res/LTtest_pval_5perc/all_prots.tsv')
# 
