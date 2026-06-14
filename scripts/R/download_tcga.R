# ==============================================================================
# SCRIPT: TCGAbiolinks ile TCGA RNA-Seq Verisi İndirme ve Matris Oluşturma
# ==============================================================================

# Loglama işlemini başlat (Hataları ve mesajları log dosyasına yazdırır)
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("TCGAbiolinks kütüphanesi yükleniyor...\n")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# Snakemake parametrelerini değişkene al
project_id <- snakemake@params[["project"]]
out_counts <- snakemake@output[["counts"]]
out_meta   <- snakemake@output[["metadata"]]

cat(sprintf("GDC sunucusuna sorgu atılıyor: %s\n", project_id))

# GDC API Sorgusu
query <- GDCquery(
  project = project_id,
  data.category = snakemake@params[["category"]],
  data.type = snakemake@params[["type"]],
  workflow.type = snakemake@params[["workflow"]]
)

cat("Veriler indiriliyor (Bu işlem boyut ve internet hızına göre sürebilir)...\n")
GDCdownload(query)

cat("Sayım (Count) matrisi ve Klinik veri hazırlanıyor...\n")
se <- GDCprepare(query)

# Ham okuma (Unstranded count) sayımlarını matris olarak çek
counts_matrix <- assay(se, "unstranded")

# Klinik ve fenotipik verileri (Metadata) çek
meta_data <- colData(se)

cat("Veriler yerel diske kaydediliyor...\n")
write.csv(counts_matrix, file = out_counts, quote = FALSE)

# CSV formatına uymayan 'list' tipli kolonları tespit edip düz metne çeviriyoruz
meta_df <- as.data.frame(meta_data)
meta_df[] <- lapply(meta_df, function(x) {
  if (is.list(x)) {
    sapply(x, paste, collapse = ", ")
  } else {
    x
  }
})

write.csv(meta_df, file = out_meta, quote = TRUE)

cat("BÖLÜM 1 TAMAMLANDI: TCGA verisi başarıyla indirildi ve temizlendi!\n")