library(covr)
library(limma) #for EList-class
library(testthat)

# example from limma::contrasts.fit
set.seed(0)
M <- matrix(rnorm(100*6, sd=0.3), nrow=100, ncol=6)
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
design <- cbind(First3Arrays=c(1,1,1,0,0,0), Last3Arrays=c(0,0,0,1,1,1))
grp <- rep(c("First3", "Last3"), each=3)
M[1, 1:3] <- M[1, 1:3] + 2

el <- new("EList")
el$E <- M
ww <- matrix(rexp(n=nrow(M)*ncol(M)), ncol=ncol(M), nrow=nrow(M))
el$weights <- ww

phenotype <- pheno.v <- rnorm(ncol(M))
names(phenotype) <- names(pheno.v) <- colnames(M)
pheno.v[1:3] <- pheno.v[1:3]-1
pheno2 <- pheno.v
pheno2[1] <- NA
covar <- rnorm(length(pheno.v))

grp2 <- batch2design(grp)

contr.v <- c(First3="First3", Last3="Last3", Last3vsFirst3="Last3-First3")
eztt <- limma_contrasts(object=M, grp = grp, contrast.v = contr.v)
lc <- limma_cor(object=M, phenotype = pheno.v)

G <- list(list(name="pwy1", description=NA, genes=paste0("gene", 1:10)),
          list(name="pwy2", description=NA, genes=paste0("gene", 11:20)),
          list(name="pwy3", description=NA, genes=paste0("gene", 21:30)))

rcn.f <- roast_contrasts(object=M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry")
set.seed(0)
rcn.m <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast")
rcr.f <- roast_cor(M, G=G, feat.tab=eztt, pheno=pheno.v, fun="fry")
set.seed(0)
rcr.m <- roast_cor(M, G=G, feat.tab=eztt, pheno=pheno.v, fun="mroast")

fl <- lapply(G, FUN=function(x) x$genes)
names(fl) <- lapply(G, FUN=function(x) x$name)

fit <- limma::eBayes(lmFit(M))
fit2 <- limma::eBayes(lmFit(M, design = design))
