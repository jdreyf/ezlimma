#' Transform SMPDB protein and/or metabolite pathway data frames into a pathway object
#'
#' Transform SMPDB protein and/or metabolite pathway data frames, read in as CSV files from http://smpdb.ca/downloads,
#' into a pathway object that matches Pathway Commons, using Gene.Name to identify proteins and ChEBI IDs for metabolites.
#' Format of SMPDB downloads has changed, so now it has numerous semi-redundant fatty acid pathways, so this will need to
#' be modified if want to stay up to date with SMPDB.
#'
#' @param smpdb.prot Data frame of proteins per pathway from SMPDB.
#' @param smpdb.met Data frame of metabolites per pathway from SMPDB.
#' @param exclude.pwy.subj Vector of subjects (or types) of pathways to exclude

SMPDB2G <- function(smpdb.prot, smpdb.met, exclude.pwy.subj=NA){
  stopifnot(c('Pathway.Subject', 'Gene.Name') %in% colnames(smpdb.prot),
            c('Pathway.Type', 'ChEBI.ID') %in% colnames(smpdb.met))

  if (!is.na(exclude.pwy.subj[1])){
    smpdb.prot <- smpdb.prot[!(smpdb.prot$Pathway.Subject %in% exclude.pwy.subj),]
    smpdb.met <- smpdb.met[!(smpdb.met$Pathway.Type %in% exclude.pwy.subj),]
  }
  
  #remove analytes w/ no ID
  smpdb.prot <- smpdb.prot[smpdb.prot$Gene.Name!="",]
  smpdb.met <- smpdb.met[!is.na(smpdb.met$ChEBI.ID),]
  
  smpdb.prot$feat.name <- smpdb.prot$Gene.Name
  smpdb.met$feat.name <- paste0('CHEBI:', smpdb.met$ChEBI.ID)

  smpdb.tab <- rbind(smpdb.prot[,c('Pathway.Name', 'feat.name')], smpdb.met[,c('Pathway.Name', 'feat.name')])

  smpdb.tab <- smpdb.tab[!is.na(smpdb.tab$feat.name) & smpdb.tab$feat.name!="",]

  all.feats <- sort(unique(smpdb.tab$feat.name))
  all.pwys <- sort(unique(smpdb.tab$Pathway.Name))
  
  # need for loop to add names
  G <- list()
  for (pwy.u in unique(smpdb.tab$Pathway.Name)){
    G[[pwy.u]] <- list(name = pwy.u, description=NA, genes=smpdb.tab[smpdb.tab$Pathway.Name == pwy.u, "feat.name"])
  }
  return(G)
}
