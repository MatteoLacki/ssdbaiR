# install.packages("data.table")
# install.packages("httr")

#' General method to connect to STRING-DB and retrieve output.
#'
#' @param identifiers Character array of protein identifiers.
#' @param format Change to 'json' or 'xml' if 'tsv' does not work.
#' @return POST object: the content field contains binary outputs.
get_string_db_data_ = function(x, method, format='tsv', ...){
  http = file.path("https://string-db.org/api", format, method, fsep='/')
  body = list(identifiers=paste(x, collapse='\r'), ...)
  data = httr::POST(http, body=body)
  Sys.sleep(1)
  return(data)
}

readers = list(tsv=data.table::fread)

#' Get tabular data from the post.
get_tabular_data = function(x, method, format='tsv', ...){
  X = get_string_db_data_(x, method, format, ...)
  binary = readBin(X$content, "character")
  X = readers[[format]](binary)
  return(data.table::data.table(X))
}

## If data.tables does not work, or string-db will stop tsv service:
# install.packages('jsonlite')
# install.packages('xml2')
# readers = list(tsv=data.table::fread,
#                json=jsonlite::fromJSON,
#                xml=xml2::read_html)


#' Map identifiers to string-db format.
#'
#' @param identifiers Character array of protein identifiers.
#' @param species_ncbi_id The NCBI identification of the organism.
#' @param format Change to 'json' or 'xml' if 'tsv' does not work.
#' @param echo_query 1 - identifiers are appended to the output table. 0 - they ain't.
#' @param ... Other settings listed in https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
#' @return data.table with translated entries. Some might be MISSING or there might be MULTIPLE VALID TRANSLATIONS. So think what about what you do.
map_identifiers = function(identifiers,
                           species_ncbi_id,
                           format='tsv',
                           echo_query=1,
                           ...)
  get_tabular_data(identifiers,
                   "get_string_ids",
                   format,
                   species=species_ncbi_id,
                   echo_query=echo_query,
                   ...)

#' Get functional enrichment.
#'
#' @param identifiers Character array of protein identifiers or stringIDs.
#' @param species_ncbi_id The NCBI identification of the organism.
#' @param format Change to 'json' or 'xml' if 'tsv' does not work.
#' @param ... Other settings listed in https://string-db.org/cgi/help.pl?subpage=api%23getting-functional-enrichment
#' @return data.table with translated entries.
get_functional_enrichment = function(identifiers,
                                     species_ncbi_id,
                                     format='tsv',
                                     caller_identity="TenzerLab",                                              ...)
  get_tabular_data(identifiers,
                   "enrichment",
                   format,
                   species=species_ncbi_id,
                   caller_identity="TenzerLab",
                   ...)


#' Save a network plot showing functional dependencies.
#'
#' @param identifiers Character array of protein identifiers or stringIDs.
#' @param species_ncbi_id The NCBI identification of the organism.
#' @param path Where to save the image. It must have a 'png' or 'svg' extension.
#' @param highres Download high resolution image?
#' @param ... Other settings listed in https://string-db.org/cgi/help.pl?subpage=api%23getting-functional-enrichment
#' @return data.table with translated entries.
saveNetworkPlot = function(identifiers,
                           species_ncbi_id,
                           path,
                           highres=FALSE, ...){
  img_type = tools::file_ext(path)
  if(!(img_type %in% c('svg','png'))) stop('Unsupported type of image.')
  if(img_type == 'png'){
    if(highres) img_type = 'highres_image'
    else img_type = 'image'
  }
  DATA = get_string_db_data_(identifiers,
                             method='network',
                             format=img_type,
                             species=species_ncbi_id,
                             ...)
  if(length(DATA$content) > 0) writeBin(DATA$content, path) else print("empty result!")
  return(invisible())
}

