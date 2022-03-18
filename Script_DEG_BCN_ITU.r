rm(list=ls()); gc(); options(digits=4)
setwd("G:/Mi unidad/NASIR/Barcelona/DevelopedData2022/Resultados HCC NASIR firmas/data_code/")

#libraries
require(reshape); library(ggplot2); library(ggpubr); library(dplyr)

# load dataset
load("data/hcc_data(3).Rdata")
ls()

clin <- dat_tumor
# clin <- clin[clin$Benefit3!="ED",]
# clin$Benefit3 <- factor(clin$Benefit3, levels=c("Remission", "Stable disease", "Progression"))
table(clin$Benefit2)
clin$Benefit2[clin$Benefit2=="Remission/Stable disease"] <- "Remission.Stabilization"
clin$Benefit2 <- factor(clin$Benefit2,
                        levels=c("Remission.Stabilization", "Progression"))

clin$Benefit2[clin$Patient %in% c("1009","6011")] = "Progression"

ex <- ex_counts_pc
summary(as.numeric(as.matrix(ex)))


### -------------------------------------------
### DEG analysis
### -------------------------------------------

### perform Differential expression
library(edgeR)

## Design matrix
clin <- clin[order(clin$id_ex),]
ex <- ex[, colnames(ex) %in% clin$id_ex]
ex <- ex[, order(colnames(ex))]
identical(as.character(clin$id_ex), as.character(colnames(ex)))


# DGE list object
library(edgeR)
d0 <- DGEList(counts=ex, genes=rownames(ex), group=as.factor(clin$Benefit2))
d0


### -------------------------------------------
### Filter low expression
### -------------------------------------------


# log counts per million reads (CPM)
cpm <- cpm(d0, log = TRUE)

# filter low expressed genes
keep.exprs <- filterByExpr(d0)
d <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d)
cpm_filt <- cpm(d, log = TRUE)



### -------------------------------------------
### Normalize before limma
### -------------------------------------------

# calculate normalization factors
d <- calcNormFactors(d, method="TMM")
cpm_norm <- cpm(d, log=TRUE)

dff = data.frame(cpm_norm)
dff2 = data.frame(Genenames = rownames(dff), dff)

ttt = as.numeric(dff2[rownames(dff2) == "RIMS2", ])[-1]
boxplot(ttt ~ clin$Benefit2, main = "RIMS2")


### -------------------------------------------
### Differential expression
### -------------------------------------------

# estimate dispersion
y <- estimateDisp(d)
sqrt(y$common.dispersion) # biological coefficient of variation

et <- exactTest(y)
topall = topTags(et, n = nrow(et$table))
summary(decideTests(et))

tt = topall$table
tt2 = data.frame(Gen = rownames(tt), tt)

openxlsx::write.xlsx(tt2, file = "NASIR_DEG_bcn_checked.xlsx", overwrite = T)




