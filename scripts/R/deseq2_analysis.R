# ==============================================================================
# SCRIPT: DESeq2 ile Tumor vs Normal Karşılaştırması ve Volcano Plot
# ==============================================================================

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("Kütüphaneler yükleniyor...\n")
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

# Dosya yolları ve parametreler
counts_path <- snakemake@input[["counts"]]
meta_path <- snakemake@input[["metadata"]]
out_report <- snakemake@output[["report"]]
out_plot <- snakemake@output[["plot"]]

alpha_val <- as.numeric(snakemake@params[["alpha"]])
lfc_val <- as.numeric(snakemake@params[["lfc_threshold"]])
min_count <- as.numeric(snakemake@params[["min_count"]])

cat("Veriler okunuyor...\n")
counts_df <- read.csv(counts_path, row.names=1, check.names=FALSE)
meta_df <- read.csv(meta_path, row.names=1, check.names=FALSE)

# Rownaems(meta) ve Colnames(counts) eşleşmesini garanti altına al
common_samples <- intersect(rownames(meta_df), colnames(counts_df))
counts_df <- counts_df[, common_samples]
meta_df <- meta_df[common_samples, ]

# Sınıflandırma faktörünü belirle ('definition' kolonu)
# Not: Boşluklu isimler sorun yaratmasın diye faktörleştiriyoruz
meta_df$definition <- as.factor(meta_df$definition)

cat("DESeq2 objesi oluşturuluyor...\n")
# Dizayn formülü: ~ definition (Yani ifadesi doku türüne göre değişen genleri bul)
dds <- DESeqDataSetFromMatrix(countData = counts_df, 
                              colData = meta_df, 
                              design = ~ definition)

# Düşük okunan (gürültü) genleri filtrele
keep <- rowSums(counts(dds)) >= min_count
dds <- dds[keep,]

cat("DESeq2 istatistiksel modeli çalıştırılıyor (Bu biraz zaman alabilir)...\n")
dds <- DESeq(dds)

# Sonuçları çıkar
res <- results(dds)
res_df <- as.data.frame(res) %>% 
          arrange(padj) # En anlamlıdan anlamsıza doğru sırala

cat("Klinik rapor kaydediliyor...\n")
write.table(res_df, file=out_report, sep="\t", quote=FALSE, row.names=TRUE)

cat("Volcano Plot çiziliyor...\n")
# Anlamlılık etiketlemesi (Hem P-değeri hem de Fold Change eşiğini geçenler)
res_df$Significant <- ifelse(!is.na(res_df$padj) & 
                             res_df$padj < alpha_val & 
                             abs(res_df$log2FoldChange) > lfc_val, "Yes", "No")

# ggplot2 ile klasik bir Volcano Plot tasarımı
p <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=Significant)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal() +
  labs(title="TCGA-PAAD: Primary Solid Tumor vs Solid Tissue Normal", 
       x="Log2 Fold Change", 
       y="-Log10 (Adjusted P-value)") +
  geom_vline(xintercept=c(-lfc_val, lfc_val), col="blue", linetype="dashed") +
  geom_hline(yintercept=-log10(alpha_val), col="blue", linetype="dashed")

ggsave(out_plot, plot=p, width=8, height=6, dpi=300)

cat("BÖLÜM 2 TAMAMLANDI: DESeq2 Analizi başarıyla bitti!\n")