context("p-value utils")

tab <- matrix(cbind(p=(1:10)/10, fdr=(1:10)/10), ncol=2, dimnames=list(letters[1:10], c("p", "fdr")))
tab2 <- cbind(foo.p=(1:9)/9, bar.p=(9:1)/9)
# ensure grep for pv.cols not confused by "pheno.cor"
tab3 <- data.frame(pheno.cor=(-2:3)/3, pheno.p=(0:5)/5, x.logFC=1:6, x.p=c(-1, (0:4)/5), 
                   Direction=rep(c("Up", "Down"), each=3), p=(6:1)/6, stringsAsFactors = FALSE)
tab3[tab3[, "pheno.cor"] == 0, "pheno.p"] <- 1
tab3.ss <- tab4 <- tab3[-1,]
tab4[2, "Direction"] <- "less"
tab3.sss <- tab3[-(1:2),]

# stat.cols="logFC|slope|cor|Direction"; p.cols="p|PValue"; alternative="greater"; nperm=NULL

# test vs Hui's fcn
test_that("combine_pvals", {
  # comb_pv
  cp <- as.numeric(c("1", "0.626422003782904", "0.593346207732144", "0.581518654612405", "0.578313236039981", 
                     "0.581518654612405", "0.593346207732144", "0.626422003782904", "1"))
  expect_equal(names(combine_pvalues(tab)), rownames(tab))
  expect_equal(combine_pvalues(tab2), cp)
  # NULL = 0 rows
  expect_error(combine_pvalues(NULL))
  expect_error(combine_pvalues(tab3))
  expect_warning(cp1 <- combine_pvalues(tab3.ss, alternative = "Up"))
  expect_warning(cp2 <- combine_pvalues(tab3.ss, alternative = "Down"))
})

test_that("grep_cols", {
  expect_error(grep_cols(tab3, p.cols = "p|PVal"))
  expect_silent(tmp <- grep_cols(tab3.ss, p.cols = "p|PVal"))
  expect_equal(tmp, c("pheno.p", "x.p", "p"))
  expect_error(grep_cols(tab3, p.cols = ".p"))
  expect_error(grep_cols(matrix(1:9, nrow=3), p.cols = "p"))
})

# test vs & roast fcns
test_that("2 --> 1 tailed", {
  expect_error(two2one_tailed(tab4))
  tto1 <- two2one_tailed(tab3.ss, alternative = "greater")
  expect_equal(tto1[1, "pheno.p"], 1-tab3.ss[1, "pheno.p"]/2)
  # sign = 0
  expect_equal(tto1[2, "pheno.p"], tab3.ss[2, "pheno.p"]/2)
  expect_equal(tto1[5, "pheno.p"], tab3.ss[5, "pheno.p"]/2)
  expect_equal(tto1[1, "x.p"], 0)
  # Direction
  expect_equal(tto1[1, "p"], tab3.ss[1, "p"]/2)
  expect_equal(tto1[5, "p"], 1-tab3.ss[5, "p"]/2)
  
  
})