library(WGCNA) library(ggplot2) library(reshape2)
file1 <- "m_MTBLS406_21237_metabolite_profiling_mass_spectrometry_v2_maf.tsv" 
file2 <- "m_MTBLS980_21833__metabolite_profiling_mass_spectrometry_v2_maf.tsv"
data1_raw <- read.delim(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
data2_raw <- read.delim(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("MTBLS406データの次元:", dim(data1_raw), "\n") 
cat("MTBLS980データの次元:", dim(data2_raw), "\n")
preprocess_metabolomics_improved <- function(data, study_name) { cat("\n========================================\n") cat("=== ", study_name, "の改善版前処理開始 ===\n") cat("========================================\n\n")
# ステップ1: 空のヘッダー行の削除
if (data[1, "metabolite_identification"] == "") { data <- data[-1, ] cat("✓ 空のヘッダー行を削除しました\n") }
# ステップ2: サンプル列の抽出
metabolite_col <- "metabolite_identification" sample_cols <- grep("^PD\\.", colnames(data), value = TRUE) 
cat("✓ サンプル列数:", length(sample_cols), "\n")
sample_data <- data[, sample_cols, drop = FALSE]
# 数値変換
for (col in sample_cols) { sample_data[[col]] <- as.numeric(sample_data[[col]]) }
# ステップ3: 有効な代謝物の選択
metabolite_names <- data[[metabolite_col]] 
valid_rows <- metabolite_names != "" & !is.na(metabolite_names) 
sample_data <- sample_data[valid_rows, , drop = FALSE] 
metabolite_names <- metabolite_names[valid_rows] 
rownames(sample_data) <- make.unique(as.character(metabolite_names))
cat("✓ 有効な代謝物数:", nrow(sample_data), "\\n\\n")
# ステップ4: 転置（サンプル × 代謝物）
expr_data <- t(sample_data) cat("--- ステップ4: 転置完了 ---\n") 
cat(" 次元:", nrow(expr_data), "サンプル ×", ncol(expr_data), "代謝物\n\n")
# ステップ5: 欠損値の処理
cat("--- ステップ5: 欠損値処理 ---\n") 
na_proportion <- colSums(is.na(expr_data)) / nrow(expr_data) 
cat(" 欠損値を含む代謝物数:", sum(na_proportion > 0), "\n") 
cat(" 50%以上欠損の代謝物数:", sum(na_proportion >= 0.5), "\n")
# 50%以上NAの列を除去
expr_data <- expr_data[, na_proportion < 0.5]
# 残りのNAを最小値の半分で補完
for (i in 1:ncol(expr_data)) { if (any(is.na(expr_data[, i]))) { min_val <- min(expr_data[, i], na.rm = TRUE) expr_data[is.na(expr_data[, i]), i] <- min_val / 2 } } 
cat("✓ 欠損値処理完了\n") 
cat(" 処理後:", nrow(expr_data), "サンプル ×", ncol(expr_data), "代謝物\n\n")
# ステップ6: 分散ゼロの列を除去
cat("--- ステップ6: 分散チェック ---\n") 
vars <- apply(expr_data, 2, var) 
zero_var_count <- sum(vars == 0) 
cat(" 分散ゼロの代謝物数:", zero_var_count, "\n") 
expr_data <- expr_data[, vars > 0] 
cat("✓ 分散ゼロの列を除去\n\n")
# ステップ7: ログ変換
cat("--- ステップ7: ログ変換 ---\n") cat(" 変換前の値の範囲:\n") cat(" 最小値:", min(expr_data, na.rm = TRUE), "\n") cat(" 最大値:", max(expr_data, na.rm = TRUE), "\n") cat(" 中央値:", median(expr_data, na.rm = TRUE), "\n")
# 負の値やゼロがある場合の対処
min_positive <- min(expr_data[expr_data > 0], na.rm = TRUE) if (min(expr_data, na.rm = TRUE) <= 0) { cat(" 警告: 負の値またはゼロが検出されました\n") cat(" 最小の正の値:", min_positive, "を使用してシフトします\n") expr_data <- expr_data + abs(min(expr_data, na.rm = TRUE)) + min_positive }
# log2変換
expr_data <- log2(expr_data + 1) 
cat("✓ log2(x+1)変換完了\n") cat(" 変換後の値の範囲:\n") cat(" 最小値:", round(min(expr_data, na.rm = TRUE), 3), "\n") cat(" 最大値:", round(max(expr_data, na.rm = TRUE), 3), "\n") cat(" 中央値:", round(median(expr_data, na.rm = TRUE), 3), "\n\n")
# ステップ8: サンプル間正規化（中央値正規化）
cat("--- ステップ8: サンプル間正規化 ---\n") 
sample_medians <- apply(expr_data, 1, median, na.rm = TRUE) cat(" サンプル中央値の範囲:", round(min(sample_medians), 3), "-", round(max(sample_medians), 3), "\n")
global_median <- median(sample_medians)
scaling_factors <- global_median / sample_medians
expr_data <- expr_data * scaling_factors

cat("✓ 中央値正規化完了\\n")
cat("  正規化後のサンプル中央値範囲:", 
    round(min(apply(expr_data, 1, median)), 3), "-",
    round(max(apply(expr_data, 1, median)), 3), "\\n\\n")
# ステップ9: アウトライヤー検出（視覚的確認）
cat("--- ステップ9: アウトライヤー検出 ---\n")
# サンプル間距離の計算
sampleTree <- hclust(dist(expr_data), method = "average")
# デンドログラムの保存
png(filename = paste0("Sample_Clustering_", study_name, ".png"), width = 12, height = 6, units = "in", res = 300) par(mar = c(2, 5, 2, 1)) plot(sampleTree, main = paste(study_name, "- Sample Clustering to Detect Outliers"), sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# カットオフラインの追加
heights <- sampleTree$height cutoff_height <- mean(heights) + 2 * sd(heights) abline(h = cutoff_height, col = "red", lwd = 2, lty = 2) text(nrow(expr_data)/2, cutoff_height, labels = paste("Cutoff =", round(cutoff_height, 1)), pos = 3, col = "red", cex = 1.2) dev.off()
cat("✓ サンプルクラスタリング図を保存:", 
    paste0("Sample_Clustering_", study_name, ".png\\n"))
# アウトライヤーの自動検出
clust <- cutreeStatic(sampleTree, cutHeight = cutoff_height, minSize = 10) outlier_samples <- (clust == 0)
if (sum(outlier_samples) > 0) {
    cat("  警告: ", sum(outlier_samples), "個のアウトライヤーサンプルを検出\\n")
    cat("  アウトライヤー:", rownames(expr_data)[outlier_samples], "\\n")
    cat("  注意: これらのサンプルは除去されていません\\n")
    cat("       デンドログラムを確認して手動で判断してください\\n")
} else {
    cat("✓ 明確なアウトライヤーは検出されませんでした\\n")
}
cat("\\n")
# ステップ10: 主成分分析によるサンプル品質チェック
cat("--- ステップ10: PCA分析 ---\n") pca_result <- prcomp(expr_data, scale. = FALSE) variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
png(filename = paste0("PCA_Plot_", study_name, ".png"),
    width = 10, height = 8, units = "in", res = 300)
par(mar = c(5, 5, 4, 2))
plot(pca_result$x[, 1], pca_result$x[, 2],
     xlab = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
     main = paste(study_name, "- PCA Plot"),
     pch = 19, cex = 1.5, col = "steelblue")
text(pca_result$x[, 1], pca_result$x[, 2], 
     labels = rownames(expr_data), 
     cex = 0.6, pos = 3)
dev.off()

cat("✓ PCAプロットを保存:", paste0("PCA_Plot_", study_name, ".png\\n"))
cat("  PC1の寄与率:", round(variance_explained[1], 1), "%\\n")
cat("  PC2の寄与率:", round(variance_explained[2], 1), "%\\n\\n")
# ステップ11: WGCNAの品質チェック
cat("--- ステップ11: WGCNA品質チェック ---\n") gsg <- goodSamplesGenes(expr_data, verbose = 3)
if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
        cat("  除去する代謝物:", sum(!gsg$goodGenes), "\\n")
    }
    if (sum(!gsg$goodSamples) > 0) {
        cat("  除去するサンプル:", sum(!gsg$goodSamples), "\\n")
    }
    expr_data <- expr_data[gsg$goodSamples, gsg$goodGenes]
} else {
    cat("✓ すべてのサンプルと代謝物が品質基準を満たしています\\n")
}
cat("\\n")
# 最終サマリー
cat("========================================\n") cat("=== 前処理完了サマリー ===\n") cat("========================================\n") cat("最終データ次元:", nrow(expr_data), "サンプル ×", ncol(expr_data), "代謝物\n") cat("データ値の範囲:", round(min(expr_data), 3), "-", round(max(expr_data), 3), "\n") cat("データ値の中央値:", round(median(expr_data), 3), "\n") cat("データ値の平均:", round(mean(expr_data), 3), "\n") cat("データ値の標準偏差:", round(sd(as.vector(expr_data)), 3), "\n") cat("========================================\n\n")
return(expr_data)
}
expr1 <- preprocess_metabolomics_improved(data1_raw, "MTBLS406")
expr2 <- preprocess_metabolomics_improved(data2_raw, "MTBLS980")
