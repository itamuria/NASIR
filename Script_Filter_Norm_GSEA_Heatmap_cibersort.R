setwd("G:/Mi unidad/NASIR/NASIR_counts/")
setwd("/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/")
dir()

# http://yulab-smu.top/clusterProfiler-book/chapter8.html
# https://biit.cs.ut.ee/gprofiler/page/apis (frogatzeko)

load("202101_counts_subread_together.RData")

se <- as.matrix(temp[,-1])
rownames(se) <- temp$GeneID

se2 <- se[,14:43]

sample <- data.frame(colnames(se2), rep("T",30))
names(sample) <- c("Sample", "Type")

main <- openxlsx::read.xlsx("../MAIN_FOLDER_NASIR/MAIN_TABLE_v007_nueva_nomenclatura.xlsx",2)
main <- main[!is.na(main$ID_1),]

# Functions ----------------------------------------------------------------

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

# deg_sh2 <- function(se2, main, contraste = "PD_IHdivEH_01", val_padj = 0.1)
# {

library("DESeq2")
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(officer)
library(EnhancedVolcano)

val_contr <- NULL
contraste = "BENEFIT2_39"
val_padj = 0.1

sample <- data.frame(main[,c("ID_1","TUMOR_RNA_T_R_10", contraste)])

sample[,contraste] <- factor(sample[,contraste])
sele1 <- which(!is.na(sample[,contraste]))
length(sele1)
sample$sele1 <- 888
sample$sele1[sele1] <- "SI"
sele2 <- which(!is.na(sample$TUMOR_RNA_T_R_10))
length(sele2)
sample$sele2 <- 888
sample$sele2[sele2] <- "SI"
sele <- ifelse(sample$sele1 == "SI" & sample$sele2 == "SI", TRUE, FALSE)
table(sele)

condition001 <- !is.na(sample$TUMOR_RNA_T_R_10)
tum_id <- sample$ID_1[condition001]
tum_rna <- sample$TUMOR_RNA_T_R_10[condition001]
groups <- sample[condition001, contraste]

se <- se2[,which(colnames(se2) %in% tum_rna)]
contr <- sample$TUMOR_RNA_T_R_10[sele]
se3 <- se[,which(colnames(se) %in% contr)]
colnames(se3) <- tum_id


# Ejemplos para tener en cuenta qué filtro utilizar -----------------------

library(clusterProfiler)
library(GeneStructureTools)

# Para ver como se distribuyen los genes en los datos crudos
check_contrast_gene <- function(gene_to_change = "CHP2",
                                se = se, 
                                tum_id, groups, masdecounts)
{
  gene_ens <- unname(mapIds(org.Hs.eg.db,
                            keys=gene_to_change,
                            column="ENSEMBL",
                            keytype="SYMBOL",
                            multiVals="first"))
  
  dd = data.frame(a = colnames(se), b = tum_id, c = groups)
  
  # se3[1:5,1:5]
  
  sele_counts <- se3[removeVersion(rownames(se3)) == gene_ens,]
  dd2 <- data.frame(counts = sele_counts, groups, patient = names(sele_counts))
  
  library(ggplot2)
  n0 <- nrow(dd2[dd2$groups == 0,])
  n1 <- nrow(dd2[dd2$groups == 1,])
  n00 <- nrow(dd2[dd2$groups == 0 & dd2$counts > masdecounts,])
  n01 <- nrow(dd2[dd2$groups == 1 & dd2$counts > masdecounts,])
  n00b <- ifelse(n00 == 0, 0, n00)
  n01b <- ifelse(n01 == 0, 0, n01)
  
  s0 <- sum(dd2$counts[dd2$groups == 0])
  s1 <- sum(dd2$counts[dd2$groups == 1])
  
  m0 <- mean(dd2$counts[dd2$groups == 0], na.rm = T)
  m1 <- mean(dd2$counts[dd2$groups == 1], na.rm = T)
  
  dd2 <- dd2[dd2$counts > masdecounts,]
  
  plot_title <- paste0(gene_to_change, " - ", contraste, " cutoff: ",masdecounts, " - Gr0:", n00b, "/", n0, "(",s0,"-",m0,") - Grp1:", n01b, "/", n1, "(", s1, "-", m1, ")")
  
  gg <- ggplot(dd2, aes(fill=groups, y=counts, x=patient)) + theme_classic() + 
    geom_bar(position="dodge", stat="identity") + coord_flip() + ggtitle(plot_title)
  
  return(list(dd2, gg))
}


masc = 5

check_contrast_gene(gene_to_change = "CHP2",
                    se = se, 
                    tum_id, groups, masdecounts = masc)


check_contrast_gene(gene_to_change = "XIST",
                    se = se, 
                    tum_id, groups, masdecounts = masc)

check_contrast_gene(gene_to_change = "CYP51A1",
                    se = se, 
                    tum_id, groups, masdecounts = masc)

check_contrast_gene(gene_to_change = "TM4SF20",
                    se = se, 
                    tum_id, groups, masdecounts = masc)


# Idea de cuantos genes quedan --------------------------------------------
dim(se)
countak <- 25
plot_title = paste0("Min counts: ", countak)
count1 <- apply(se, 1, function(h) sum(h > countak))
# sum(count1)
count1_df <- as.data.frame(table(count1))
ggplot(count1_df, aes(y=Freq, x=count1)) + theme_classic() + 
  geom_bar(position="dodge", stat="identity") + ggtitle(plot_title)  

values_count = NULL
for(g in 0:500)
{
  print(g)
  count1 <- apply(se, 1, function(h) sum(h > g))
  count1_df <- as.data.frame(table(count1))
  values_count <- c(values_count, 
                    count1_df$Freq[count1_df$count1 == 0],
                    count1_df$Freq[count1_df$count1 == 30])
}

df_counts = data.frame(matrix(values_count, byrow = T, ncol = 2))
df_counts$id <- 1:nrow(df_counts)
df_counts2 = reshape2::melt(df_counts,"id")

ggplot(df_counts2, aes(x=id, y=value, group=variable)) +
  geom_line(aes(color=variable)) + xlim(0,50)

# Filtro ------------------------------------------------------------------

# Vamos a filtrar:
# en la mitad de los samples de alguno de los grupos más de 5 counts
# es decir si tenemos 10 zeros y 20 unos, debe haber más de 7 con más de 5 counts y/o más de 12 con más de 5 counts

mas_de_5_por_grupo <- function(se3 = se3, grupos = sample$BENEFIT2_39[sele])
{
  g0 <- se3[,grupos == 0]
  g1 <- se3[,grupos == 1]
  
  coun0 <- apply(g0, 1, function(h) sum(h > 5))
  g0_fin <- g0[coun0 > 5,]
  
  coun1 <- apply(g1, 1, function(h) sum(h > 5))
  g1_fin <- g1[coun1 > 5,]
  
  namess <- unique(c(rownames(g0_fin), rownames(g1_fin)))
  return(namess)
}


# Funcion completa --------------------------------------------------------




deg_sh2 <- function(se2, main, contraste = "BENEFIT2_39", 
                    val_padj = 0.1, 
                    wd = "/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/Results_sensibility",
                    version1 = "004")
{
  
  library("DESeq2")
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(officer)
  library(EnhancedVolcano)
  
  val_contr <- NULL
  # contraste = "BENEFIT2_39"
  # val_padj = 0.1
  powerpointfile <- paste0(wd, "/",contraste,"_", version1,"_PowerPoint_Results.pptx")
  
  sample <- data.frame(main[,c("ID_1","TUMOR_RNA_T_R_10", contraste)])
  
  sample[,contraste] <- factor(sample[,contraste])
  sele1 <- which(!is.na(sample[,contraste]))
  length(sele1)
  sample$sele1 <- 888
  sample$sele1[sele1] <- "SI"
  sele2 <- which(!is.na(sample$TUMOR_RNA_T_R_10))
  length(sele2)
  sample$sele2 <- 888
  sample$sele2[sele2] <- "SI"
  sele <- ifelse(sample$sele1 == "SI" & sample$sele2 == "SI", TRUE, FALSE)
  table(sele)
  
  condition001 <- !is.na(sample$TUMOR_RNA_T_R_10)
  tum_id <- sample$ID_1[condition001]
  tum_rna <- sample$TUMOR_RNA_T_R_10[condition001]
  groups <- sample[condition001, contraste]
  
  se <- se2[,which(colnames(se2) %in% tum_rna)]
  contr <- sample$TUMOR_RNA_T_R_10[sele]
  se3 <- se[,which(colnames(se) %in% contr)]
  colnames(se3) <- tum_id  
  
  
  # Con los que nos quedamos
  remove_names <- mas_de_5_por_grupo(se3 = se3, grupos = sample[sele,contraste])
  
  # se3_froga <- se3[apply(se3, 1, sum) > 5,]
  se3_froga <- se3[rownames(se3) %in% remove_names,]  
  # se3_froga <- se3
  
  # Solo tumores ------------------------------------------------------------
  
  var_formula <- as.formula(paste0("~ ", contraste))
  dds <- DESeqDataSetFromMatrix(countData = se3_froga, colData = sample[sele,], design = var_formula)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  # save(dds, file = "dds_contraste_benefit2.RData")
  # save(sample, file = "sample_contraste_benefit2.RData")
  # save(sele2, file = "sele_contraste_benefit2.RData")
  
  boxplot(log(se3+1))
  boxplot(log(se3_froga+1))
  summary(normalized_counts)
  boxplot_logcounts <- boxplot(log(normalized_counts+1))
  
  # write.table(normalized_counts, file=paste0("SandraIbon/",contraste,"_normalized_counts2.txt"), 
  #             sep="\t", quote=F, col.names=NA)
  
  # Diff expression
  dds <- DESeq(dds)
  res <- results(dds)
  # head(results(dds, tidy=TRUE)) #let's look at the results table
  summary(res)
  
  res <- res[order(res$pvalue),]
  head(res)
  res2 = data.frame(res)
  res2_temp = data.frame(res)
  
  na_padj <- which(is.na(res2$padj))
  
  for(h in na_padj)
  {
    res2[h,"padj"] <- res2[h-1,"padj"]
  }
  
  table(res2$padj < 0.1)
  table(res2$pvalue < 0.001)
  
  
  # jpeg(paste0("Images/", contraste, "examples_RNAseq.jpg"))
  # par(mfrow=c(2,3))
  plotcoung1 <- plotCounts(dds, gene=rownames(res)[1], intgroup=contraste, col = 2, pch = 16)
  plotcoung2 <- plotCounts(dds, gene=rownames(res)[2], intgroup=contraste, col = 2, pch = 16)
  plotcoung3 <- plotCounts(dds, gene=rownames(res)[3], intgroup=contraste, col = 2, pch = 16)
  plotcoung4 <- plotCounts(dds, gene=rownames(res)[4], intgroup=contraste, col = 2, pch = 16)
  plotcoung5 <- plotCounts(dds, gene=rownames(res)[5], intgroup=contraste, col = 2, pch = 16)
  plotcoung6 <- plotCounts(dds, gene=rownames(res)[6], intgroup=contraste, col = 2, pch = 16)
  # par(mfrow=c(1,1))
  # dev.off()
  
  genesymbol <- read.table("/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/gene_symbol.txt", sep ="\t")
  
  res_frame <- data.frame(trans = rownames(res2), data.frame(res2))
  
  symbol <- mapIds(org.Hs.eg.db,
                   keys=gsub("\\..*","",rownames(res_frame)),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")
  symbol <- unname(ifelse(is.na(symbol), names(symbol), symbol))
  res_frame2 <- data.frame(sym = data.frame(symbol), res_frame)
  
  # plot
  
  # boxplot(normalized_counts[rownames(normalized_counts) == "ENSG00000001630.16"] ~ sample[sele,"BENEFIT2_39"], 
  #         ylab= "Gene expression", xlab= "Benefit2", main = "CYP51A1")
  # 
  # boxplot(normalized_counts[rownames(normalized_counts) == "ENSG00000005073.5"] ~ sample[sele,"BENEFIT2_39"], 
  #         ylab= "Gene expression", xlab= "Benefit2", main = "HOXA11")
  
  # openxlsx::write.xlsx(res_frame2, file = paste0(wd, "/",contraste,"_", version,"_DEG_with_Full_.xlsx"), overwrite = T)
  
  dtemp <- data.frame(normalized_counts[rownames(normalized_counts) == "ENSG00000005073.5"], sample[sele,"BENEFIT2_39"])
  
  
  # Volcano -----------------------------------
  mm <- -log10(min(res_frame2$padj, na.rm = TRUE))
  
  library(EnhancedVolcano)
  
  table(res_frame2$padj < 0.1)
  table(res_frame2$pvalue < 0.001)
  
  row0 =nrow(res_frame2[res_frame2$padj < 0.1 & res_frame2$log2FoldChange < -1,])
  row1 =nrow(res_frame2[res_frame2$padj < 0.1 & res_frame2$log2FoldChange > 1,])
  
  title_plot = paste0("Benefit B(1) vs B(0) ----- B(0) ", row0, " -- B(1) ", row1, " - padj")
  
  
  volc1 <- EnhancedVolcano(res_frame2,
                           lab = res_frame2$symbol,
                           x = 'log2FoldChange',
                           y = 'padj',
                           # xlim = c(-8, 8),
                           ylim = c(0,mm),
                           title = title_plot,
                           pCutoff = 1e-1,
                           FCcutoff = 1,
                           pointSize = 3.0,
                           labSize = 3.0)
  
  
  res_val3 <- res_frame2[res_frame2$padj < val_padj,]
  res_val3 <- res_val3[!is.na(res_val3$symbol),]
  
  res_val3_pvalue <- res_frame2[res_frame2$pvalue < 0.05,]
  
  # sensibility analysis of pvalue and qvalue
  save_qp_values <- NULL
  cutofss <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
  
  for(qp in cutofss)
  {
    sen_qpvalue <- qp
    qq <- nrow(res_frame2[res_frame2$padj < sen_qpvalue,])
    pp <- nrow(res_frame2[res_frame2$pvalue < sen_qpvalue,])
    save_qp_values <- c(save_qp_values, pp, qq)
  }
  df_qp <- data.frame(cutofss, matrix(save_qp_values, byrow = T, ncol = 2))
  names(df_qp) <- c("Cutoff", "Count pvalue", "Count padj")
  
  # export
  openxlsx::write.xlsx(df_qp, file = paste0(wd, "/",contraste,"_", version1,"_DEG_summary.xlsx"), overwrite = T)
  openxlsx::write.xlsx(res_val3, file = paste0(wd, "/",contraste,"_", version1,"_DEG_with_Filt_v002.xlsx"), overwrite = T)
  openxlsx::write.xlsx(res_frame2, file = paste0(wd, "/",contraste,"_", version1,"_DEG_with_Filt_v002_sinfiltrar.xlsx"), overwrite = T)
  
  
  # ClusterProfiler
  
  res_frame2b <- res_frame2[order(res_frame2$log2FoldChange, decreasing = TRUE),]
  d <- res_frame2b[,c("symbol","log2FoldChange")]
  geneList <- d[,2]
  names(geneList) <- as.character(d[,1])
  # geneList <- sort(geneList, decreasing = TRUE)
  
  
  # geneList2 <- data.frame(id = names(geneList), FC = geneList)
  # geneList2b <- geneList2[abs(geneList2$FC) > 2,]
  # geneList2b$id <- as.character(paste0(substr(geneList2b$id,1,1), tolower(substr(geneList2b$id,2,nchar(as.character(geneList2b$id))))))
  # gene <- geneList2b$id
  
  
  
  #   # Wikipath ------------------------------------------------------------
  
  
  
  library(magrittr)
  library(clusterProfiler)
  
  # data(geneList, package="DOSE")
  # gene <- names(geneList)[abs(geneList) > 2]
  
  # genesymbol
  
  # library(clusterProfiler)
  # bitr(gsub("\\..*","",rownames(resframe8)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  # Warning message:
  # In bitr(gsub("\\..*", "", rownames(resframe8)), fromType = "ENSEMBL",  :
  #           20.52% of input gene IDs are fail to map...
  
  
  resframe8 <- res_frame[res_frame$pvalue < 0.05,]
  # resframe8 <- resframe8[!is.na(resframe8$trans),]
  entrezid <- mapIds(org.Hs.eg.db,
                     keys=gsub("\\..*","",rownames(resframe8)),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
  entrezid2 <- unname(entrezid)
  table(is.na(entrezid))
  prop.table(table(is.na(entrezid)))
  dborra <- data.frame(entrezid2, resframe8$trans)
  entrezid2 <- entrezid2[!is.na(entrezid2)]
  
  
  # https://wikipathways-data.wmcloud.org/current/gmt/
  wpgmtfile <- "/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/wikipath/wikipathways-20210510-gmt-Homo_sapiens.gmt"
  wp2gene <- read.gmt(wpgmtfile)
  ont <- "term"
  wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  
  ewp <- enricher(entrezid2, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  head(ewp)
  ewpdf <- data.frame(ewp)
  # openxlsx::write.xlsx(ewpdf, file = paste0(wd, "/",contraste,"_", version1,"_enricher001.xlsx"), overwrite = T)
  
  # GSEA
  
  # res_val3 <- res_val3[order(res_val3$log2FoldChange, decreasing = TRUE),]
  
  
  resframeGSEA <- res_frame[order(res_frame$log2FoldChange, decreasing = TRUE),]
  resframeGSEA <- resframeGSEA[!is.na(resframeGSEA$trans),]
  dGSEA <- resframeGSEA$trans
  
  entrezidGSEA <- mapIds(org.Hs.eg.db,
                         keys=gsub("\\..*","",dGSEA),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  entrezid2GSEA <- unname(entrezidGSEA)
  
  geneListGSEA <- resframeGSEA$log2FoldChange
  names(geneListGSEA) <- entrezid2GSEA
  
  entrezid2GSEA_sinNa <- !is.na(names(geneListGSEA))
  geneListGSEA_sinNA <- geneListGSEA[entrezid2GSEA_sinNa]
  geneListGSEA_sinNA2 <- geneListGSEA_sinNA[!is.na(geneListGSEA_sinNA)]
  
  ewp2 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
  head(ewp2)
  
  
  
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggupset)
  
  ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
  head(ewp)
  
  if(nrow(ewp) > 0)
  {
    # library(enrichplot)
    heatmap1 <- heatplot(ewp, foldChange=geneListGSEA_sinNA2)
    # emapplot(ewp2)
    upset1 <- upsetplot(ewp)
  }
  
  
  # ridgeplot(ewp)
  gs1 <- gseaplot(ewp2, geneSetID = 1, title = ewp2$Description[1])
  gs2 <- gseaplot2(ewp2, geneSetID = 1, title = ewp2$Description[1])
  
  # gseaplot2(ewp2, geneSetID = 1:3, pvalue_table = TRUE,
  #           color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
  
  openxlsx::write.xlsx(data.frame(ewp), file = paste0(wd, "/",contraste,"_", version1,"_Wikipath_hipergeometrico2.xlsx"), overwrite = T)
  openxlsx::write.xlsx(data.frame(ewp2), file = paste0(wd, "/",contraste,"_", version1,"_Wikipath_gsea2.xlsx"), overwrite = T)
  
  # ploting
  
  library(enrichplot)
  barpl1 <- barplot(ewp, showCategory=4)
  dot1 <- dotplot(ewp, showCategory=30) + ggtitle("Wikipath for Progresion")
  
  dot2 <- dotplot(ewp2, showCategory=30) + ggtitle("Wikipath for Progresion")
  
  
  ### Plot heatmap
  
  df2 <- data.frame(ewp2)
  df2_genes <- df2$core_enrichment[2]
  
  sel_gen <- unlist(strsplit(df2_genes, "/"))
  
  normalized_counts2 <- counts(dds, normalized=TRUE)
  
  sample[,contraste] <- as.character(sample[,contraste])
  TTPlar <- ifelse(is.na(sample[,contraste]), "NA", sample[,contraste])
  group_df <- data.frame(TTPlar[sele2])
  rownames(group_df) <- sample$ID_1[sele2]
  names(group_df) <- contraste
  
  symbol <- mapIds(org.Hs.eg.db,
                   keys=gsub("\\..*","",rownames(normalized_counts2)),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")
  
  symbol <- unname(ifelse(is.na(symbol), names(symbol), symbol))
  
  heatdf <- normalized_counts2[symbol %in% sel_gen,]
  
  symbol_df <- data.frame(symbol)
  # rownames(symbol_df) <- rownames(heatdf)
  symbol2 <- symbol[symbol %in% sel_gen]
  
  phet1 <- pheatmap(heatdf, scale = "row", labels_col = sample$ID[sele2], cluster_cols = FALSE,
                    annotation_col = group_df, labels_row = symbol2, fontsize_row = 5) 
  
  phet2 <- pheatmap(heatdf, scale = "row", labels_col = sample$ID[sele2], cluster_cols = TRUE,
                    annotation_col = group_df, labels_row = symbol2, fontsize_row = 10)
  # # Cell marker -----------------------------------------------------------
  
  
  
  cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
  cell_markers
  
  y <- enricher(entrezid2, TERM2GENE=cell_markers, minGSSize=1)
  
  y <- setReadable(y, org.Hs.eg.db, keyType = "ENTREZID")
  DT::datatable(as.data.frame(y))
  
  openxlsx::write.xlsx(  data.frame(y), file = paste0(wd, "/",contraste,"_", version1,"_CellMarker.xlsx"), overwrite = T)
  
  # MSigDb
  
  
  library(msigdbr)
  msigdbr_show_species()
  
  m_df <- msigdbr(species = "Homo sapiens")
  head(m_df, 2) %>% as.data.frame
  
  m_df_desc <- unique(m_df[,c("gs_name","gs_description")])
  
  # m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  #   dplyr::select(gs_name, entrez_gene)
  # head(m_t2g)
  # 
  # em <- enricher(entrezid2, TERM2GENE=m_t2g)
  # em2 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = m_t2g)
  # head(em)
  
  library(org.Hs.eg.db)
  # em <- setReadable(em, org.Hs.eg.db, keyType = "ENTREZID")
  # em2 <- setReadable(em2, org.Hs.eg.db, keyType = "ENTREZID")
  # head(em)
  # head(em2)
  
  # openxlsx::write.xlsx(data.frame(em), file = paste0(wd, "/",contraste,"_", version1,"_MSigDb_hipergeometrico.xlsx"), overwrite = T)
  # openxlsx::write.xlsx(data.frame(em2), file = paste0(wd, "/",contraste,"_", version1,"_MSigDb_gsea.xlsx"), overwrite = T)
  # 
  # C-s
  # m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  #   dplyr::select(gs_name, entrez_gene)
  # head(m_t2g)
  # 
  # em3 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = m_t2g)
  # em3 <- setReadable(em3, org.Hs.eg.db, keyType = "ENTREZID")
  # d9 <- data.frame(em3)
  # gseaplot(em3, geneSetID = 7, title = em3$Description[7])
  
  for(call in c(paste0("C",1:7),"H"))
  {
    print(call)
    m_t2g <- NULL
    m_t2g <- msigdbr(species = "Homo sapiens", category = call) %>% 
      dplyr::select(gs_name, entrez_gene)
    em3 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = m_t2g)
    em3 <- setReadable(em3, org.Hs.eg.db, keyType = "ENTREZID")
    em <- enricher(entrezid2, TERM2GENE=m_t2g)
    em <- setReadable(em, org.Hs.eg.db, keyType = "ENTREZID")
    
    # head(em3)
    m_df_desc[m_df_desc$gs_name == "HALLMARK_BILE_ACID_METABOLISM",]
    
    em3b <- merge(em3, m_df_desc, by.x = "ID", by.y = "gs_name")
    emb <- merge(em, m_df_desc, by.x = "ID", by.y = "gs_name")
    
    openxlsx::write.xlsx(data.frame(em3b), file = paste0(wd, "/",contraste,"_", version1,"_MSigDb_gsea_",call,".xlsx"), overwrite = T)
    openxlsx::write.xlsx(data.frame(emb), file = paste0(wd, "/",contraste,"_", version1,"_MSigDb_enrich_",call,".xlsx"), overwrite = T)
    print(paste0(call," - ",dim(data.frame(em3b))[1], " counts in gsea"))
    print(paste0(call," - ",dim(data.frame(emb))[1], " counts in enrich"))
  }
  
  # Disease
  library(DOSE)
  x <- enrichDO(gene          = entrezid2,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                universe      = names(geneListGSEA_sinNA2),
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
  # head(x)
  x <- setReadable(x, org.Hs.eg.db, keyType = "ENTREZID")
  
  # gene2 <- names(geneListGSEA_sinNA2)[abs(geneListGSEA_sinNA2) < 3]
  # ncg <- enrichNCG(gene2)
  # head(ncg)
  # openxlsx::write.xlsx(data.frame(ncg), file = paste0(wd, "/",contraste,"_", version1,"_ncg_.xlsx"), overwrite = T)
  
  dgn <- enrichDGN(entrezid2)
  head(dgn)
  dgn <- setReadable(dgn, org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(dgn), file = paste0(wd, "/",contraste,"_", version1,"_enrich_dgn_.xlsx"), overwrite = T)
  
  library(DOSE)
  
  y <- gseDO(geneListGSEA_sinNA2,
             nPerm         = 100000,
             minGSSize     = 120,
             pvalueCutoff  = 0.2,
             pAdjustMethod = "BH",
             verbose       = FALSE)
  head(y, 3)
  y <- setReadable(y, org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(y), file = paste0(wd, "/",contraste,"_", version1,"_gseDO_gsea.xlsx"), overwrite = T)
  
  
  # gene2 <- names(geneListGSEA_sinNA2)[abs(geneListGSEA_sinNA2) > 2]
  ncg1 <- enrichNCG(entrezid2)
  
  
  ncg2 <- gseNCG(geneListGSEA_sinNA2,
                 nPerm         = 10000,
                 minGSSize     = 120,
                 pvalueCutoff  = 0.2,
                 pAdjustMethod = "BH",
                 verbose       = FALSE)
  ncg1 <- setReadable(ncg1, 'org.Hs.eg.db')
  ncg2 <- setReadable(ncg2, 'org.Hs.eg.db')
  head(ncg1, 3)
  openxlsx::write.xlsx(data.frame(ncg1), file = paste0(wd, "/",contraste,"_", version1,"_ncg_enrich.xlsx"), overwrite = T)
  openxlsx::write.xlsx(data.frame(ncg2), file = paste0(wd, "/",contraste,"_", version1,"_ncg_gsea.xlsx"), overwrite = T)
  
  # p1 <- dotplot(ncg1, showCategory=30) + ggtitle("dotplot for Enrich")
  # p2 <- dotplot(ncg2, showCategory=30) + ggtitle("dotplot for GSEA")
  # plot_grid(p1, p2, ncol=2)
  
  
  
  
  
  dgn <- gseDGN(geneListGSEA_sinNA2,
                nPerm         = 10000,
                minGSSize     = 120,
                pvalueCutoff  = 0.2,
                pAdjustMethod = "BH",
                verbose       = FALSE)
  dgn <- setReadable(dgn, 'org.Hs.eg.db')
  head(dgn, 3)
  openxlsx::write.xlsx(data.frame(dgn), file = paste0(wd, "/",contraste,"_", version1,"_dgn_gsea.xlsx"), overwrite = T)
  
  
  # Gene Ontology Analysis
  
  library(clusterProfiler)
  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  gene.df <- bitr(entrezid2, fromType = "ENTREZID",
                  toType = c("ENSEMBL", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  head(gene.df)
  
  ggo <- groupGO(gene     = entrezid2,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 2,
                 readable = TRUE)
  
  head(ggo)
  
  openxlsx::write.xlsx(data.frame(ggo), file = paste0(wd, "/",contraste,"_", version1,"_ggo.xlsx"), overwrite = T)
  
  
  ego <- enrichGO(gene          = entrezid2,
                  universe      = names(geneListGSEA_sinNA2),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "none",
                  readable      = TRUE)
  head(ego)
  
  # ez det ulertzen ezberdintasuna
  
  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "none")
  
  head(ego2)
  
  ego <- setReadable(ego, OrgDb = org.Hs.eg.db)
  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
  
  openxlsx::write.xlsx(data.frame(ego), file = paste0(wd, "/",contraste,"_", version1,"_ego.xlsx"), overwrite = T)
  openxlsx::write.xlsx(data.frame(ego2), file = paste0(wd, "/",contraste,"_", version1,"_ego2.xlsx"), overwrite = T)
  
  
  ego3 <- gseGO(geneList     = geneListGSEA_sinNA2,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                nPerm        = 1000,
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
  
  ego3 <- setReadable(ego3, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(ego3), file = paste0(wd, "/",contraste,"_", version1,"_ego3.xlsx"), overwrite = T)
  
  
  # KEGG
  
  # data(geneList, package="DOSE")
  # gene <- names(geneList)[abs(geneList) > 2]
  
  kk <- enrichKEGG(gene         = entrezid2,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  
  kk <- setReadable(ego3, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(kk), file = paste0(wd, "/",contraste,"_", version1,"_KEGG_enrich.xlsx"), overwrite = T)
  
  kk2 <- gseKEGG(geneList     = geneListGSEA_sinNA2,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  head(kk2)
  kk2 <- setReadable(kk2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(kk2), file = paste0(wd, "/",contraste,"_", version1,"_KEGG_gsea.xlsx"), overwrite = T)
  
  mkk <- enrichMKEGG(gene = entrezid2,
                     organism = 'hsa')
  mkk <- setReadable(mkk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(mkk), file = paste0(wd, "/",contraste,"_", version1,"_KEGG_mkk.xlsx"), overwrite = T)
  
  mkk2 <- gseMKEGG(geneList = geneListGSEA_sinNA2,
                   organism = 'hsa')
  mkk2 <- setReadable(mkk2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(mkk2), file = paste0(wd, "/",contraste,"_", version1,"_KEGG_mkk2_gsea.xlsx"), overwrite = T)
  
  # REACTOME
  
  # library(GOSemSim)
  # hsGO <- godata('org.Hs.eg.db', ont="MF")
  # 
  # library(AnnotationHub)
  # hub <- AnnotationHub()
  # q <- query(hub, "Cricetulus")
  # id <- q$ah_id[length(q)]
  # Cgriseus <- hub[[id]]
  # 
  # goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Jiang")
  # goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
  # 
  # go1 = c("GO:0004022","GO:0004024","GO:0004174")
  # go2 = c("GO:0009055","GO:0005515")
  # mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
  # 
  # mgoSim(go1, go2, semData=hsGO, measure="Wang", combine="BMA")
  # 
  # geneSim("241", "251", semData=hsGO, measure="Wang", combine="BMA")
  # 
  # hsGO2 <- godata('org.Hs.eg.db', keytype = "ENTREZID", ont="MF", computeIC=FALSE)
  # # genes <- c("CDC45", "MCM10", "CDC20", "NMU", "MMP1")
  # 
  # mgeneSim(entrezid2[1:4], semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
  
  
  library(ReactomePA)
  
  x <- enrichPathway(gene=entrezid2, pvalueCutoff = 0.05, readable=TRUE)
  head(x)
  x <- setReadable(x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(x), file = paste0(wd, "/",contraste,"_", version1,"_ReactomePA_hiper.xlsx"), overwrite = T)
  
  y <- gsePathway(geneListGSEA_sinNA2, 
                  pvalueCutoff = 0.2,
                  pAdjustMethod = "BH", 
                  verbose = FALSE)
  head(y)
  y <- setReadable(y, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(y), file = paste0(wd, "/",contraste,"_", version1,"_ReactomePA_gsea.xlsx"), overwrite = T)
  
  # library(graphite)
  # viewPathway(" The citric acid (TCA) cycle and respiratory electron transport", 
  #             readable = TRUE, 
  #             foldChange = geneListGSEA_sinNA2)
  
  
  # library("pathview")
  # hsa04080 <- pathview(gene.data  = geneListGSEA_sinNA2,
  #                      pathway.id = "hsa04080",
  #                      species    = "hsa",
  #                      limit      = list(gene=max(abs(geneListGSEA_sinNA2)), cpd=1))
  
  # ------------------------------------------------------------------------
  
  
  
  # Heatmap from pathways ---------------------------------------------------
  
  
  ### Plot heatmap
  
  t_gsea1a <- openxlsx::read.xlsx(paste0(wd, "/",contraste,"_", version1,"_Wikipath_hipergeometrico2.xlsx"), 1)
  t_gsea1b <- openxlsx::read.xlsx(paste0(wd, "/",contraste,"_", version1,"_Wikipath_gsea2.xlsx"), 1)
  
  
  # df2 <- data.frame(ewp2)
  df2_genes <- t_gsea1b$core_enrichment[2]
  
  sel_gen <- unlist(strsplit(df2_genes, "/"))
  
  normalized_counts2 <- counts(dds, normalized=TRUE)
  
  sample[,contraste] <- as.character(sample[,contraste])
  TTPlar <- ifelse(is.na(sample[,contraste]), "NA", sample[,contraste])
  group_df <- data.frame(TTPlar[sele2])
  # rownames(group_df) <- sample$ID_1[sele2]
  rownames(group_df) <- sample$ID_1[sele2]
  names(group_df) <- contraste
  
  symbol <- mapIds(org.Hs.eg.db,
                   keys=gsub("\\..*","",rownames(normalized_counts2)),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")
  
  symbol <- unname(ifelse(is.na(symbol), names(symbol), symbol))
  
  heatdf <- normalized_counts2[symbol %in% sel_gen,]
  
  symbol_df <- data.frame(symbol)
  # rownames(symbol_df) <- rownames(heatdf)
  symbol2 <- symbol[symbol %in% sel_gen]
  
  ph3 <- pheatmap(heatdf, scale = "row", labels_col = sample$ID[sele2], cluster_cols = FALSE,
                  annotation_col = group_df, labels_row = symbol2, fontsize_row = 5) 
  
  ph4 <- pheatmap(heatdf, scale = "row", labels_col = sample$ID[sele2], cluster_cols = TRUE,
                  annotation_col = group_df, labels_row = symbol2, fontsize_row = 10, cellheight = 15)
  
  save_pheatmap_pdf(ph3, filename = paste0(wd, "/",contraste,"_", version1,"_hipergeom_heatmap.pdf"), 
                    width=7, height=7)
  
  
  save_pheatmap_pdf(ph4, filename = paste0(wd, "/",contraste,"_", version1,"_gsea_heatmap.pdf"), 
                    width=7, height=7)
  
  # Volcano -----------------------------------
  mm <- -log10(min(t_gsea1b$pvalue, na.rm = TRUE))
  
  volcn2 <- EnhancedVolcano(t_gsea1b,
                            lab = t_gsea1b$ID,
                            x = 'NES',
                            y = 'pvalue',
                            # xlim = c(-8, 8),
                            ylim = c(0,mm),
                            # title = paste0(contraste, '_with_padj01'),
                            pCutoff = 5e-3,
                            FCcutoff = 1.5,
                            pointSize = 3.0,
                            labSize = 3.0,xlab = "NES")
  
  
  
  # Power point -------------------------------------------------------------
  library(officer)
  library("grid")
  library("ggplotify")
  
  doc <- read_pptx()
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = boxplot_logcounts, location = ph_location_fullsize())
  # doc <- ph_with(doc, value = "Volcano plot", location = ph_location_type(type = "title"))
  
  
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = plotcoung1, location = ph_location_fullsize())
  # 
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = plotcoung2, location = ph_location_fullsize())
  # 
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = plotcoung3, location = ph_location_fullsize())
  # 
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = plotcoung4, location = ph_location_fullsize())
  # 
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = plotcoung5, location = ph_location_fullsize())
  # 
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = plotcoung6, location = ph_location_fullsize())
  
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = volc1, location = ph_location_fullsize())
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = heatmap1, location = ph_location_fullsize())
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = upset1, location = ph_location_fullsize())
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = gs1, location = ph_location_fullsize())
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = gs2, location = ph_location_fullsize())
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = barpl1, location = ph_location_fullsize())
  
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = dot1, location = ph_location_fullsize())
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = dot2, location = ph_location_fullsize())
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = phet1[[4]], location = ph_location_fullsize())
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = phet2, location = ph_location_fullsize())
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = ph3, location = ph_location_fullsize())
  # doc <- add_slide(doc)
  # doc <- ph_with(doc, value = ph4, location = ph_location_fullsize())
  doc <- add_slide(doc)
  doc <- ph_with(doc, value = volcn2, location = ph_location_fullsize())
  
  print(doc, target = powerpointfile)
}



deg_sh2(se2, main, contraste = "BENEFIT2_39", 
        val_padj = 0.1, 
        wd = "/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/Results_sensibility",
        version1 = "004")


deg_sh2(se2, main, contraste = "Progresion_45", 
        val_padj = 0.1, 
        wd = "/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/Results_sensibility",
        version1 = "004")

otros <- c("RESPUESTA_43","Etiology_32","PVI_29","AFPMzy400_35")
otros <- c("PVI_29")

for(oo in otros)
{
  deg_sh2(se2, main, contraste = oo, 
          val_padj = 0.1, 
          wd = "/Volumes/GoogleDrive/My Drive/NASIR/NASIR_counts/Results_sensibility",
          version1 = "004")
}



# Heatmap general ---------------------------------------------------------

library("DESeq2")
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(officer)
library(EnhancedVolcano)

val_contr <- NULL
contraste = "BENEFIT2_39"
val_padj = 0.1
sample <- data.frame(main[,c("ID_1","TUMOR_RNA_T_R_10", contraste)])

sample[,contraste] <- factor(sample[,contraste])
sele1 <- which(!is.na(sample[,contraste]))
sample$sele1 <- 888
sample$sele1[sele1] <- "SI"
sele2 <- which(!is.na(sample$TUMOR_RNA_T_R_10))
sample$sele2 <- 888
sample$sele2[sele2] <- "SI"
sele <- ifelse(sample$sele1 == "SI" & sample$sele2 == "SI", TRUE, FALSE)

condition001 <- !is.na(sample$TUMOR_RNA_T_R_10)
tum_id <- sample$ID_1[condition001]
tum_rna <- sample$TUMOR_RNA_T_R_10[condition001]
groups <- sample[condition001, contraste]

se <- se2[,which(colnames(se2) %in% tum_rna)]
contr <- sample$TUMOR_RNA_T_R_10[sele]
se3 <- se[,which(colnames(se) %in% contr)]
colnames(se3) <- tum_id  


# Con los que nos quedamos
remove_names <- mas_de_5_por_grupo(se3 = se3, grupos = sample[sele,contraste])
se3_froga <- se3[rownames(se3) %in% remove_names,]  

# Solo tumores ------------------------------------------------------------

var_formula <- as.formula(paste0("~ ", contraste))
dds <- DESeqDataSetFromMatrix(countData = se3_froga, colData = sample[sele,], design = var_formula)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

dds <- DESeq(dds)
res <- results(dds)
# head(results(dds, tidy=TRUE)) #let's look at the results table
summary(res)

res <- res[order(res$pvalue),]
head(res)
res2 = data.frame(res)
res2_temp = data.frame(res)

na_padj <- which(is.na(res2$padj))

for(h in na_padj)
{
  res2[h,"padj"] <- res2[h-1,"padj"]
}

deg_subset <- res2[res2$padj < 0.1,]

library("RColorBrewer")
# vsd <- vst(dds, blind=FALSE)
# sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(normalized_counts))
sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$sample)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# DEG subset
normalized_counts_deg <- normalized_counts[rownames(normalized_counts) %in% rownames(deg_subset),]


# Agregando SYMBOL --------------------------------------------------------

entrezid <- mapIds(org.Hs.eg.db,
                   keys=gsub("\\..*","",rownames(normalized_counts)),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")  
table(is.na(entrezid))
normalized_counts_symbol <- data.frame(Symbol = entrezid, normalized_counts)

dim(normalized_counts_symbol)
dim(lm22)
normalized_counts_symbol$Symbol[normalized_counts_symbol$Symbol %in% lm22$V1]  
lm22$V1[!lm22$V1 %in% normalized_counts_symbol$Symbol]  

names(normalized_counts_symbol)
View(normalized_counts_symbol)

normalized_counts_symbol2 <- normalized_counts_symbol[!is.na(normalized_counts_symbol$Symbol),]
write_tsv(normalized_counts_symbol2, file = "NASIR_featureCount_norm.txt")

# which(rownames(normalized_counts_symbol) == "ENSG00000113520")
# temp <- normalized_counts_symbol
# temp$ensg <- rownames(temp)



# Covariates --------------------------------------------------------------


covr = data.frame(
  stringsAsFactors = FALSE,
  ID_1 = c(1001L,
           1008L,1009L,1010L,2001L,3001L,3002L,3003L,
           4002L,4003L,4004L,4005L,5001L,5002L,5003L,
           5006L,5007L,6007L,6011L,6012L,6013L,6014L,
           6015L,6016L,6017L,6018L,7001L,8003L,
           8004L,8005L),
  ID_PACIENTE_2 = c(1.001,
                    1.008,1.009,1.01,2.001,3.001,3.002,3.003,
                    4.002,4.003,4.004,4.005,5.001,5.002,5.003,
                    5.006,5.007,6.007,6.011,6.012,6.013,6.014,
                    6.015,6.016,6.017,6.018,7.001,8.003,
                    8.004,8.005),
  PVI_29 = c(0L,0L,0L,
             1L,0L,1L,0L,0L,0L,0L,0L,0L,1L,0L,1L,
             1L,0L,0L,0L,0L,0L,0L,1L,0L,0L,0L,
             1L,0L,0L,0L),
  TACE_30 = c(0L,0L,0L,
              0L,0L,0L,1L,0L,0L,1L,0L,1L,0L,0L,0L,
              0L,1L,0L,1L,1L,0L,0L,0L,1L,0L,0L,
              0L,0L,0L,0L),
  Sorafenib_31 = c(0L,0L,0L,
                   0L,0L,0L,0L,0L,1L,0L,0L,0L,0L,0L,0L,
                   0L,1L,0L,0L,0L,0L,1L,0L,0L,0L,0L,
                   0L,0L,0L,0L),
  Etiology_32 = c(1L,2L,1L,
                  1L,1L,1L,1L,2L,1L,1L,1L,1L,2L,1L,2L,
                  1L,2L,1L,1L,1L,1L,1L,1L,1L,1L,1L,
                  2L,2L,1L,2L),
  AgeMay65_33 = c(1L,1L,1L,
                  1L,0L,1L,1L,0L,1L,1L,1L,1L,0L,1L,0L,
                  0L,0L,1L,0L,0L,0L,0L,0L,0L,0L,0L,
                  0L,0L,1L,0L),
  AFPMzy400_35 = c(1L,0L,0L,
                   1L,1L,1L,1L,0L,1L,0L,0L,0L,0L,0L,1L,
                   0L,0L,0L,1L,0L,0L,0L,1L,0L,0L,1L,
                   0L,0L,1L,0L),
  MAAuptake_36 = c(1L,2L,2L,
                   1L,1L,2L,2L,1L,1L,1L,2L,1L,1L,1L,1L,
                   1L,2L,1L,1L,2L,1L,1L,1L,1L,1L,1L,
                   1L,2L,1L,2L),
  BENEFIT3_38 = c(2L,1L,1L,
                  3L,3L,3L,3L,2L,2L,3L,1L,2L,3L,3L,3L,
                  1L,2L,2L,1L,3L,1L,1L,1L,3L,1L,2L,
                  1L,2L,1L,1L),
  BENEFIT2_39 = c(1L,1L,1L,
                  0L,0L,0L,0L,1L,1L,0L,1L,1L,0L,0L,0L,
                  1L,1L,1L,1L,0L,1L,1L,1L,0L,1L,1L,
                  1L,1L,1L,1L),
  BOR_41 = c(3L,1L,2L,
             3L,2L,3L,3L,3L,3L,4L,2L,3L,3L,3L,4L,
             1L,3L,2L,2L,3L,2L,2L,1L,3L,1L,3L,
             1L,3L,2L,2L),
  RESPUESTA_43 = c(0L,1L,1L,
                   0L,1L,0L,0L,0L,0L,0L,1L,0L,0L,0L,0L,
                   1L,0L,1L,1L,0L,1L,1L,1L,0L,1L,0L,
                   1L,0L,1L,1L),
  CONTROL_44 = c(1L,1L,1L,
                 1L,1L,1L,1L,1L,1L,0L,1L,1L,1L,1L,0L,
                 1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,
                 1L,1L,1L,1L),
  Progresion_45 = c(1L,1L,1L,
                    1L,1L,1L,0L,0L,0L,1L,0L,1L,1L,1L,1L,
                    0L,1L,0L,1L,1L,0L,0L,0L,1L,1L,0L,
                    1L,1L,1L,1L),
  PatronPD_47 = c("NIH",
                  "IHG","IHG","IHG","EH","EH",NA,NA,NA,"EH",
                  NA,"NIH","IHG","NIH","EH",NA,"NIH",NA,
                  "EH","NIH",NA,NA,NA,"EH","IHG",NA,"NIH",
                  "IHG","NIH","EH"),
  PatronPD_123_47b = c(1L,3L,3L,
                       3L,2L,2L,NA,NA,NA,2L,NA,1L,3L,1L,2L,
                       NA,1L,NA,2L,1L,NA,NA,NA,2L,3L,NA,
                       1L,3L,1L,2L),
  PatronPD_en2_48 = c("IH","IH",
                      "IH","IH","EH","EH",NA,NA,NA,"EH",NA,
                      "IH","IH","IH","EH",NA,"IH",NA,"EH",
                      "IH",NA,NA,NA,"EH","IH",NA,"IH","EH","IH",
                      "EH"),
  PatronPD_en2_12_48b = c(1L,1L,1L,
                          1L,2L,2L,NA,NA,NA,2L,NA,1L,1L,1L,2L,
                          NA,1L,NA,2L,1L,NA,NA,NA,2L,1L,NA,
                          1L,2L,1L,2L),
  PFS_49 = c(1L,1L,1L,
             1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,
             0L,1L,0L,1L,1L,0L,0L,0L,1L,1L,0L,
             1L,1L,1L,1L),
  Muerte_51 = c(1L,1L,1L,
                1L,1L,0L,1L,1L,1L,1L,1L,0L,1L,0L,1L,
                0L,0L,0L,1L,1L,0L,0L,0L,1L,1L,0L,
                0L,1L,0L,0L),
  Resection_53 = c(0L,0L,0L,
                   0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                   0L,0L,0L,0L,0L,1L,0L,1L,0L,0L,1L,
                   0L,0L,0L,0L),
  Nivolumabpost_54 = c(0L,0L,0L,
                       0L,0L,0L,0L,1L,1L,0L,1L,0L,0L,0L,0L,
                       0L,0L,1L,1L,0L,0L,1L,0L,0L,0L,0L,
                       0L,0L,0L,0L),
  TKIpost_55 = c(0L,0L,0L,
                 0L,1L,1L,0L,0L,0L,NA,0L,NA,0L,1L,1L,
                 0L,0L,0L,0L,1L,0L,0L,0L,1L,0L,0L,
                 0L,0L,NA,NA),
  CRPR_SDPD_61 = c("SDPD",
                   "CRPR","CRPR","SDPD","CRPR","SDPD","SDPD",
                   "SDPD","SDPD","SDPD","CRPR","SDPD","SDPD",
                   "SDPD","SDPD","CRPR","SDPD","SDPD","SDPD",
                   "SDPD","CRPR","CRPR","CRPR","SDPD","CRPR",
                   "SDPD","CRPR","CRPR","CRPR","CRPR"),
  CRPR_SDPD_bin_61b = c(0L,1L,1L,
                        0L,1L,0L,0L,0L,0L,0L,1L,0L,0L,0L,0L,
                        1L,0L,0L,0L,0L,1L,1L,1L,0L,1L,0L,
                        1L,1L,1L,1L)
)

annotation_col = data.frame(
  GeneClass = factor(c(1001L, 1008L, 1009L, 1010L, 2001L, 3001L, 3002L, 3003L, 
                       4002L, 4003L, 4004L, 4005L, 5001L, 5002L, 5003L, 5006L, 
                       5007L, 6007L, 6011L, 6012L, 6013L, 6014L, 6015L, 6016L, 
                       6017L, 6018L, 7001L, 8003L, 8004L, 8005L)
  ))

rowizen = c(1001L, 1008L, 1009L, 1010L, 2001L, 3001L, 3002L, 3003L, 
            4002L, 4003L, 4004L, 4005L, 5001L, 5002L, 5003L, 5006L, 
            5007L, 6007L, 6011L, 6012L, 6013L, 6014L, 6015L, 6016L, 
            6017L, 6018L, 7001L, 8003L, 8004L, 8005L)

annotation_row = covr[,c(3,6,8,11,13,15)]

for(h in 1:ncol(annotation_row))
{
  annotation_row[,h] <- factor(annotation_row[,h])
}

rownames(annotation_row) = rownames(sampleDistMatrix)

# creat colours for each group
newCols <- colorRampPalette(grDevices::rainbow(2))
set.seed(2012)

PVI_29_col <- newCols(2)
names(PVI_29_col) <- unique(c("0","1"))

Etiology_32_col <- newCols(2)
names(Etiology_32_col) <- unique(c("1","2"))

AFPMzy400_35_col <- newCols(2)
names(AFPMzy400_35_col) <- unique(c("0","1"))

BENEFIT2_39_col <- newCols(2)
names(BENEFIT2_39_col) <- unique(c("0","1"))

RESPUESTA_43_col <- newCols(2)
names(RESPUESTA_43_col) <- unique(c("0","1"))

Progresion_45_col <- newCols(2)
names(Progresion_45_col) <- unique(c("0","1"))

annoCol <- list(PVI_29 = PVI_29_col, 
                Etiology_32 = Etiology_32_col, 
                AFPMzy400_35 = AFPMzy400_35_col, 
                BENEFIT2_39 = BENEFIT2_39_col, 
                RESPUESTA_43 = RESPUESTA_43_col, 
                Progresion_45 = Progresion_45_col
)




# # Diapo 5, cluster 1 ----------------------------------------------------



pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         # col=colors,
         annotation_row = annotation_row,
         labels_row = rowizen,
         labels_col = rowizen,
         annotation_colors = annoCol)

# Diapo 6, cluster with DEG

pheatmap(normalized_counts_deg,
         scale = "row",
         annotation_col = annotation_row,
         show_rownames =F,
         annotation_colors = annoCol
         
         # clustering_distance_rows=sampleDists2,
         # clustering_distance_cols=sampleDists2,
         # col=colors,
         # annotation_row = annotation_row,
         # labels_row = rowizen,
         # labels_col = rowizen,
         # annotation_colors = annoCol
)



# CIBERSORT ---------------------------------------------------------------

setwd("/Volumes/GoogleDrive/My Drive/PabloSarobe/Immunodeconv/")

library(Rserve)
library('e1071')
library('parallel')
library('colorRamps')  

Rserve(args="--no-save")

source("/Volumes/GoogleDrive/My Drive/PabloSarobe/Immunodeconv/CIBERSORT.R")

lm22 <- read.table("LM22.txt", sep = "\t")

# results <- CIBERSORT("LM22.txt", "20210121_symbol_fpkm_nasir.txt", perm=100, QN=TRUE)
results <- CIBERSORT("LM22.txt", "NASIR_featureCount_norm.txt", perm=100, QN=TRUE)

# save(results, file = "20211130_CIBERSORT_Results.RData")

load("20211130_CIBERSORT_Results.RData")    


results2 <- data.frame(results[,-c(23:25)])
results2$pat <- rownames(results2)
results2b <- reshape2::melt(results2, id = "pat")

library(ggplot2)

ggplot(aes(x = pat, y = value, fill=variable), data=results2b) + 
  geom_bar(stat="summary")


# Heatmap
sapply(results, class)
cib_mean = apply(results[,-c(23:25)], 2, mean)
cib_mean[cib_mean == 0]
results3 <- data.frame(results[,-which(colnames(results) %in% c("T cells gamma delta","Eosinophils"))])
results3 <- results3[,-c(21:23)]
results4 <- t(results3)
colnames(results4) <- gsub("X","",colnames(results4))

library(pheatmap)
pheatmap(results4,
         scale = "row",
         annotation_col = annotation_row,
         show_rownames = T,
         annotation_colors = annoCol)

# Tabla
results5_table <- data.frame(results3, annotation_row)
check_names <- data.frame(rownames(results3), rownames(annotation_row))

ciber_vector <- NULL
n <- 0 

for(h in colnames(annotation_row))
{
  print(h)
  for(x in colnames(results3))
  {
    print(x)
    n <- n + 1 
    print(n)
    
    temp <- data.frame(results3[,x], annotation_row[,h])
    names(temp) <- c(x, h)
    
    tt <- t.test(temp[,1] ~ temp[,2])
    tt$p.value
    ciber_vector <- c(ciber_vector, h, x, tt$p.value)
  }
  
}

ciber_table <- data.frame(matrix(ciber_vector, byrow = T, ncol = 3))
names(ciber_table) <- c("Contraste","Cibersort","p-valor")
head(ciber_table)

ciber_table2 <- reshape2::dcast(ciber_table, Cibersort ~ Contraste)
openxlsx::write.xlsx(ciber_table2, file = "20211201_Cibersort_pvalores_contraste.xlsx")

sapply(ciber_table2, class)
for(b in 2:ncol(ciber_table2))
{
  ciber_table2[,b] <- round(as.numeric(ciber_table2[,b]),3)
}

library(data.table)

library(dplyr)

library(formattable)

library(tidyr)

#Set a few color variables to make our table more visually appealing

customGreen0 = "#DeF7E9"

customGreen = "#71CA97"

customRed = "#ff7f7f"

formattable(ciber_table2)


library(formattable)

# df_try <- structure(list(Normal = structure(2:1, .Label = c("p-values", 
#                                                             "stistic"), class = "factor"), Jan = c(0.93069466, 0.05123532
#                                                             ), Feb = c(0.90404849, 0.01056474)), class = "data.frame", row.names = c("1", 
#                                                                                                                                      "2"))

first_col_formatter <- formatter("span", 
                                 style = ~ style(color = "grey",
                                                 font.weight = "bold"))

improvement_formatter <- formatter("span", 
                                   style = x ~ style(
                                     font.weight = "bold", 
                                     color = ifelse(x < 0.1, "red", "black")))

# improvement2_formatter <- formatter("span", 
#                                     style = x ~ style(
#                                       font.weight = "bold", 
#                                       color = ifelse(x > 0.05, "red", "black")))

formattable(df_try,
            align =c("l","c","c"),
            list(Normal = first_col_formatter,
                 area(row = 1, col = -1) ~ improvement_formatter,
                 area(row = 2, col = -1) ~ improvement_formatter))

summary(ciber_table2)

# formattable(ciber_table2, align =c("l","c","c","c","c", "c", "c"), list(
#   `Cibersort` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
#   `AFPMzy400_35`= color_tile(customRed, customGreen0),
#   `BENEFIT2_39`= color_tile(customRed, customGreen0),
#   `Etiology_32`= color_tile(customRed, customGreen0),
#   `Progresion_45`= color_tile(customRed, customGreen0),
#   `PVI_29`= color_tile(customRed, customGreen0),
#   `RESPUESTA_43`= color_tile(customRed, customGreen0)
# ))

formattable(ciber_table2, align =c("l","c","c","c","c", "c", "c"), list(
  `Cibersort` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  area(row = 1:20, col = 2) ~ improvement_formatter,
  area(row = 1:20, col = 3) ~ improvement_formatter,
  area(row = 1:20, col = 4) ~ improvement_formatter,
  area(row = 1:20, col = 5) ~ improvement_formatter,
  area(row = 1:20, col = 6) ~ improvement_formatter,
  area(row = 1:20, col = 7) ~ improvement_formatter
))


t.test(results3$Mast.cells.resting ~ annotation_row$Progresion_45)
t.test(results3$NK.cells.activated ~ annotation_row$Progresion_45)
t.test(results3$Monocytes ~ annotation_row$RESPUESTA_43)
t.test(results3$T.cells.CD8 ~ annotation_row$PVI_29)


# GSVA --------------------------------------------------------------------




# Firmas ------------------------------------------------------------------


library(gtools)

ann6 <- normalized_counts_symbol2

h <- 0
ann7 <- sapply(ann6[,-1], function(x) 
{
  h <<- h + 1
  print(paste0(h, " - ", names(ann6)[h]))
  decil <- quantcut(x, seq(0,1,by=0.1) )
  as.numeric(decil)
})

ann7b <- data.frame(Symbol = ann6$Symbol, ann7)


# Log10 -------------------------------------------------------------------

ann8_log <- log10(ann7b[,-1])
ann8_log <- data.frame(Symbol = ann7b$Symbol, ann8_log)

save(ann8_log, file = "20211201_ann8_log.RData")

# Groups mean -----------------------------------------------------------

#• Load groups of genes

gengr <- openxlsx::read.xlsx("/Volumes/GoogleDrive/My Drive/NASIR/SignaturesSandra/Gene-signature.xlsx",3)

# Check if we have all the interesting genes

gengr$Genes[!gengr$Genes %in% ann8_log$Symbol]

un_gr <- unique(gengr$Group)

un_gen <- unique(gengr$Genes)
un_gen <- un_gen[which(un_gen != "NGK7")]

ann8_log$Symbol <- gsub("-",".",ann8_log$Symbol )
un_gen <- gsub("-",".",un_gen)

grep("HLA",ann8_log$Pat, value = T)

un_gen <- c(un_gen, "CD4", "CD3D", "CD3E","CD3G", "CD247", "CD8A","CD8B", "CD8B2")
# un_gen <- c("CD4", "CD3D", "CD3E","CD3G", "CD247", "CD8A","CD8B", "CD8B2")


library(tidyverse)

t8 <- ann8_log[ann8_log$Symbol %in% un_gen,]
t8 <- t8[-which(rownames(t8) == 14668), ]
rownames(t8) <- t8$Symbol

ann9_fil <- data.frame(t(t8[,-1]))
# rownames(ann8_log) = ann8_log$Symbol
# ann8_log = t(ann8_log[,-1])
# ann9_fil <- data.frame(ann8_log[,colnames(ann8_log) %in% un_gen])

# openxlsx::write.xlsx(ann9_fil, file = "20210215_Sandra_CDxx.xlsx")

library(dplyr)

# temp <- data.frame(colnames(ann8_log),colnames(ann8_log))


names(ann9_fil)
grep("HLA",names(ann9_fil), value = T)

ann9_fil %>%
  mutate(Inflamatory = CD274+CD8A+LAG3+STAT1,
         Cytolytic32 = GZMA+PRF1,
         Gajewski33 = CCL2+CCL3+CCL4+CD8A+CXCL10+CXCL9+GZMK+HLA.DMA+HLA.DMB+HLA.DOA+HLA.DOB+ICOS+IRF1,
         InterferonGamaSig = CXCL10+CXCL9+HLA.DRA+IDO1+IFNG+STAT1,
         AntigenPresenting = CMKLR1+HLA.DQA1+HLA.DRB1+PSMB10,
         InterferonGamaBiology = CCL5+CD27+CXCL9+CXCR6+IDO1+STAT1,
         TcellExhaustion = CD276+CD8A+LAG3+PDCD1LG2+TIGIT,
         T_nK_sig = HLA.E,
         RibasGeneInterferon = CCR5+CXCL10+CXCL11+CXCL9+GZMA+HLA.DRA+IDO1+IFNG+PRF1+STAT1,
         Inflamatory_mean = (CD274+CD8A+LAG3+STAT1)/4,
         Cytolytic32_mean = (GZMA+PRF1)/2,
         Gajewski33_mean = (CCL2+CCL3+CCL4+CD8A+CXCL10+CXCL9+GZMK+HLA.DMA+HLA.DMB+HLA.DOA+HLA.DOB+ICOS+IRF1)/13,
         InterferonGamaSig_mean = (CXCL10+CXCL9+HLA.DRA+IDO1+IFNG+STAT1)/6,
         AntigenPresenting_mean = (CMKLR1+HLA.DQA1+HLA.DRB1+PSMB10)/4,
         InterferonGamaBiology_mean = (CCL5+CD27+CXCL9+CXCR6+IDO1+STAT1)/6,
         TcellExhaustion_mean = (CD276+CD8A+LAG3+PDCD1LG2+TIGIT)/5,
         T_nK_sig_mean = (HLA.E)/1,
         RibasGeneInterferon_mean = (CCR5+CXCL10+CXCL11+CXCL9+GZMA+HLA.DRA+IDO1+IFNG+PRF1+STAT1)/10) -> ann10_fin


# Checking t.test with log-cuantiles values

ciber_vector <- NULL
n <- 0 

# t9 <- data.frame(t(t8[,-1]))
# colnames(t9) <- t8$Symbol

for(h in colnames(ann10_fin))
{
  print(h)
  for(x in colnames(annotation_row))
  {
    print(x)
    n <- n + 1 
    print(n)
    
    temp <- data.frame(ann10_fin[,h], annotation_row[,x])
    names(temp) <- c( h, x)
    
    if(var(temp[,1]) != 0)
    {
      tt <- t.test(temp[,1] ~ temp[,2])
      tt$p.value
      ciber_vector <- c(ciber_vector, h, x, tt$p.value)
    }
    
    
  }
  
}

ciber_table <- data.frame(matrix(ciber_vector, byrow = T, ncol = 3))
names(ciber_table) <- c("Contraste","Firma","p-valor")
head(ciber_table)
table(ciber_table$Contraste)
ciber_table$`p-valor` <- as.numeric(ciber_table$`p-valor`)
ciber_table <- unique(ciber_table)

ciber_table2 <- reshape2::dcast(ciber_table, Firma ~ Contraste)
openxlsx::write.xlsx(ciber_table2, file = "20211201_FIRMAS_pvalores_contraste.xlsx")

# Significancia

ciber_table3 <- data.frame(colnames(ciber_table2[,-1]), t(ciber_table2[,-1]))
colnames(ciber_table3)  <- c("Genes",ciber_table2$Firma)

formattable(ciber_table3, align =c("l","c","c","c","c", "c", "c"), list(
  `Cibersort` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  
  area(row = 1:53, col = 2) ~ improvement_formatter,
  area(row = 1:53, col = 3) ~ improvement_formatter,
  area(row = 1:53, col = 4) ~ improvement_formatter,
  area(row = 1:53, col = 5) ~ improvement_formatter,
  area(row = 1:53, col = 6) ~ improvement_formatter,
  area(row = 1:53, col = 7) ~ improvement_formatter
  
))

ann11 <- data.frame(t(ann10_fin[,-grep("mean",colnames(ann10_fin))]))

ann12 <- ann11[apply(ann11, 1, var) > 0,]
colnames(ann12) <- gsub("X","", colnames(ann12))

pheatmap(ann12,
         scale = "row",
         annotation_col = annotation_row,
         show_rownames =T,
         annotation_colors = annoCol)



summary(ann10_fin)

