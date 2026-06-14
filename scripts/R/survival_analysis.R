# ==============================================================================
# SCRIPT: Top 20 Hedefin Sağkalım Taraması ve En İyi Prognostik Genin Çizimi
# ==============================================================================

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("Klinik paketler yükleniyor...\n")
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
})

meta_path <- snakemake@input[["metadata"]]
counts_path <- snakemake@input[["counts"]]
targets_path <- snakemake@input[["top_targets"]]
out_plot <- snakemake@output[["plot"]]

cat("Veriler okunuyor...\n")
meta_df <- read.csv(meta_path, row.names=1, check.names=FALSE)
counts_df <- read.csv(counts_path, row.names=1, check.names=FALSE)
targets_df <- read.csv(targets_path)

# HİLE BURADA: Orijinal matrisi bozmuyoruz, ayrı bir arama vektörü oluşturuyoruz
clean_counts_ids <- gsub("\\..*", "", rownames(counts_df))

# Tümörlü hastaları filtrele ve sağkalım günlerini hazırla
tumor_meta <- meta_df[meta_df$definition == "Primary solid Tumor", ]
tumor_meta$time <- ifelse(tumor_meta$vital_status == "Dead", 
                          as.numeric(tumor_meta$days_to_death), 
                          as.numeric(tumor_meta$days_to_last_follow_up))
tumor_meta$status <- ifelse(tumor_meta$vital_status == "Dead", 1, 0)
tumor_meta <- tumor_meta[!is.na(tumor_meta$time) & tumor_meta$time > 0, ]

common_samples <- intersect(rownames(tumor_meta), colnames(counts_df))
tumor_meta <- tumor_meta[common_samples, ]

cat("20 Potansiyel Hedef Sağkalım (Log-Rank) Testinden Geçiriliyor...\n")
surv_results <- data.frame(Gene_Symbol=character(), Ensembl_ID=character(), P_Value=numeric(), stringsAsFactors=FALSE)

# HER BİR GEN İÇİN OTOMATİK TARAMA DÖNGÜSÜ
for (i in 1:nrow(targets_df)) {
  gene_id <- targets_df$Ensembl_ID[i]
  gene_sym <- targets_df$Gene_Symbol[i]
  
  # Arama vektöründe genin sırasını bul (İlk eşleşeni alır, kopya hatasını aşarız)
  match_idx <- match(gene_id, clean_counts_ids)
  
  # Eğer gen sayım matrisinde yoksa döngüyü atla
  if (is.na(match_idx)) next
  
  # Satır ismini değil, bulduğumuz sıra numarasını (match_idx) kullanarak veriyi çekiyoruz!
  gene_expr <- log2(as.numeric(counts_df[match_idx, common_samples]) + 1)
  temp_meta <- tumor_meta
  temp_meta$gene_expr <- gene_expr
  
  median_expr <- median(temp_meta$gene_expr, na.rm=TRUE)
  temp_meta$Expr_Group <- ifelse(temp_meta$gene_expr > median_expr, "High", "Low")
  temp_meta$Expr_Group <- factor(temp_meta$Expr_Group, levels=c("Low", "High"))
  
  surv_obj <- Surv(time = temp_meta$time, event = temp_meta$status)
  fit <- survfit(surv_obj ~ Expr_Group, data = temp_meta)
  
  pval <- surv_pvalue(fit, data = temp_meta)$pval
  surv_results <- rbind(surv_results, data.frame(Gene_Symbol=gene_sym, Ensembl_ID=gene_id, P_Value=pval))
}

# P-değerine göre sırala (En anlamlı olan en üste)
surv_results <- surv_results %>% arrange(P_Value)
write.csv(surv_results, "results/reports/Prognostic_Survival_Results.csv", row.names=FALSE)

cat(sprintf("GERÇEK KATİL BULUNDU: %s (P-Value: %f). Grafik çiziliyor...\n", surv_results$Gene_Symbol[1], surv_results$P_Value[1]))

# EN İYİ GEN İÇİN GRAFİK ÇİZİMİ (Listenin en üstündeki)
best_gene_id <- surv_results$Ensembl_ID[1]
best_gene_sym <- surv_results$Gene_Symbol[1]
best_p <- surv_results$P_Value[1]
best_idx <- match(best_gene_id, clean_counts_ids)

gene_expr <- log2(as.numeric(counts_df[best_idx, common_samples]) + 1)
tumor_meta$gene_expr <- gene_expr
median_expr <- median(tumor_meta$gene_expr, na.rm=TRUE)
tumor_meta$Expr_Group <- ifelse(tumor_meta$gene_expr > median_expr, "High", "Low")
tumor_meta$Expr_Group <- factor(tumor_meta$Expr_Group, levels=c("Low", "High"))

surv_obj <- Surv(time = tumor_meta$time, event = tumor_meta$status)
fit <- survfit(surv_obj ~ Expr_Group, data = tumor_meta)

plot_title <- paste0("True Prognostic Driver: ", best_gene_sym, " (p = ", signif(best_p, 3), ")")

p <- ggsurvplot(fit,
                data = tumor_meta,
                pval = TRUE,
                conf.int = FALSE,
                risk.table = TRUE,
                risk.table.col = "strata",
                palette = c("blue", "red"),
                title = plot_title,
                xlab = "Time (Days)",
                ylab = "Overall Survival Probability",
                legend.title = "Expression",
                legend.labs = c("Low", "High"))

png(out_plot, width = 800, height = 600, res = 120)
print(p)
dev.off()

cat("FAZ 1 FİNAL: Otomatik tarama tamamlandı!\n")