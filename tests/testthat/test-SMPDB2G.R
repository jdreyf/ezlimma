
smpdb.met <- data.frame(Pathway.Name=c("Alanine Metabolism", "Aspartate Metabolism", "Aspartate Metabolism",
                        "Refsum Disease"), Pathway.Type=c("Metabolic", "Metabolic", "Metabolic", "Disease"),
                        ChEBI.ID=c("15956", "15366", "16027", "15351"), stringsAsFactors = FALSE)
smpdb.prot <- data.frame(Pathway.Name=c("Alanine Metabolism", "Aspartate Metabolism", "Aspartate Metabolism",
                        "Primary Hyperoxaluria Type"), Pathway.Subject=c("Metabolic", "Metabolic", "Metabolic",
                        "Disease"), Gene.Name=c("AARS", "ABAT", "ASL", "AARS"), stringsAsFactors = FALSE)

# example subset
g <- SMPDB2G(smpdb.prot=smpdb.prot, smpdb.met=smpdb.met, exclude.pwy.subj="Disease")
expect_equal(length(g[[1]]$genes), 2)
expect_equal(length(g[[2]]$genes), 4)
expect_true("AARS" %in% g[[1]]$genes)
expect_true("CHEBI:15956" %in% g[[1]]$genes)