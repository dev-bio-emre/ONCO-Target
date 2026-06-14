# ==============================================================================
# SCRIPT: Ensembl ID -> Gene Symbol Dönüşümü ve Top Hedeflerin Seçimi
# ==============================================================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
})

cat("1. DESeq2 Klinik Raporu okunuyor...\n")
report_path <- "results/reports/PAAD_clinical_significance_report.txt"
res_df <- read.table(report_path, header=TRUE, row.names=1, sep="\t")

# Biyoinformatik Hilesi: TCGA verilerindeki Ensembl ID'ler versiyon numarası içerir 
# (Örn: ENSG00000141510.16). Veritabanıyla eşleşmesi için noktadan sonrasını (.16) atıyoruz.
cat("2. Ensembl ID versiyonları temizleniyor...\n")
ensembl_ids <- gsub("\\..*", "", rownames(res_df))

cat("3. Veritabanından (org.Hs.eg.db) Gen Sembolleri çekiliyor...\n")
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Yeni kolonları veri çerçevesine ekle
res_df$Gene_Symbol <- gene_symbols
res_df$Ensembl_ID <- ensembl_ids

cat("4. En kritik potansiyel ilaç hedefleri filtreleniyor...\n")
# "Significant" kolonu yerine doğrudan p-değeri ve Fold Change eşiklerini kullanıyoruz
top_targets <- res_df %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  arrange(padj) %>%
  select(Gene_Symbol, Ensembl_ID, log2FoldChange, padj)

# Liste görünümünü temizlemek için NA (Sembolü bulunamayan pseudogenler) olanları filtreleyebiliriz
top_targets <- top_targets[!is.na(top_targets$Gene_Symbol), ]

cat("5. İlk 20 Hedef kaydediliyor...\n")
write.csv(head(top_targets, 20), "results/reports/Top20_PAAD_Targets.csv", row.names=FALSE)

cat("==============================================================\n")
cat("İŞLEM TAMAM! En iyi adaylar 'Top20_PAAD_Targets.csv' dosyasına yazıldı.\n")
cat("==============================================================\n")