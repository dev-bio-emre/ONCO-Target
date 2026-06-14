# ==============================================================================
# SCRIPT: clusterProfiler ile GO ve KEGG Yolak (Pathway) Analizi
# ==============================================================================

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("Kütüphaneler yükleniyor...\n")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

report_path <- snakemake@input[["report"]]
out_go <- snakemake@output[["go_plot"]]
out_kegg <- snakemake@output[["kegg_plot"]]

alpha_val <- as.numeric(snakemake@params[["alpha"]])
lfc_val <- as.numeric(snakemake@params[["lfc_threshold"]])

cat("1. DESeq2 verileri okunuyor ve filtreleniyor...\n")
res_df <- read.table(report_path, header=TRUE, row.names=1, sep="\t")

# Sadece anlamlı değişen genleri alıyoruz
sig_genes <- res_df %>%
  filter(!is.na(padj) & padj < alpha_val & abs(log2FoldChange) > lfc_val)

# Ensembl ID'lerin versiyon takılarını temizle (.16 gibi)
ensembl_ids <- gsub("\\..*", "", rownames(sig_genes))

cat("2. Ensembl ID'ler Entrez ID formatına çevriliyor (KEGG için zorunlu)...\n")
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = ensembl_ids,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# NA olanları (eşleşmeyenleri) çıkar
entrez_ids <- na.omit(entrez_ids)

cat("3. Gene Ontology (GO) Analizi yapılıyor...\n")
go_results <- enrichGO(gene = entrez_ids,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP", # Biological Process (Biyolojik Süreç)
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)

if (!is.null(go_results) && nrow(go_results) > 0) {
  p_go <- dotplot(go_results, showCategory=15) + ggtitle("Gene Ontology (BP) Enrichment")
  ggsave(out_go, plot=p_go, width=10, height=8, dpi=300)
} else {
  cat("Anlamlı GO yolağı bulunamadı, boş grafik oluşturuluyor...\n")
  file.create(out_go)
}

cat("4. KEGG Yolak Analizi yapılıyor...\n")
kegg_results <- enrichKEGG(gene = entrez_ids,
                           organism = 'hsa', # Homo sapiens
                           pvalueCutoff = 0.05)

if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
  p_kegg <- dotplot(kegg_results, showCategory=15) + ggtitle("KEGG Pathway Enrichment")
  ggsave(out_kegg, plot=p_kegg, width=10, height=8, dpi=300)
} else {
  cat("Anlamlı KEGG yolağı bulunamadı, boş grafik oluşturuluyor...\n")
  file.create(out_kegg)
}

cat("BÖLÜM 3 TAMAMLANDI: Zenginleştirme analizi başarıyla bitti!\n")