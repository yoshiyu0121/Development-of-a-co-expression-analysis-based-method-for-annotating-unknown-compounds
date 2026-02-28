# モジュール間の保存性を定量化する関数（Jaccard係数 + p値）
quantify_module_preservation <- function(module_df1, module_df2, study1_name = "MTBLS406", study2_name = "MTBLS980") {
cat("\\n=== モジュール保存性の定量化（Jaccard係数 + p値） ===\\n\\n")
# グレーモジュールを除く
modules1 <- module_df1[module_df1$Module != "grey", ] modules2 <- module_df2[module_df2$Module != "grey", ]
# 全代謝物のユニバース（共通代謝物のみ）
all_metabolites1 <- unique(modules1$Metabolite) all_metabolites2 <- unique(modules2$Metabolite) universe <- intersect(all_metabolites1, all_metabolites2) universe_size <- length(universe)
cat("共通代謝物の総数（ユニバース）:", universe_size, "\\n")
# ユニークなモジュールを取得
unique_modules1 <- unique(modules1$Module) unique_modules2 <- unique(modules2$Module)
cat(study1_name, "のモジュール数:", length(unique_modules1), "\\n")
cat(study2_name, "のモジュール数:", length(unique_modules2), "\\n\\n")
# 全ペアの組み合わせで計算
results_list <- list()
for (mod1 in unique_modules1) {
    metabolites1_all <- modules1$Metabolite[modules1$Module == mod1]
# ユニバース内の代謝物のみ
metabolites1 <- intersect(metabolites1_all, universe)
    for (mod2 in unique_modules2) {
        metabolites2_all <- modules2$Metabolite[modules2$Module == mod2]
# ユニバース内の代謝物のみ
metabolites2 <- intersect(metabolites2_all, universe)
# Jaccard係数の計算
intersection <- length(intersect(metabolites1, metabolites2)) union <- length(union(metabolites1, metabolites2)) jaccard <- if(union > 0) intersection / union else 0
# ハイパージオメトリック検定用の分割表
# a = 両方のモジュールに含まれる
# b = mod1のみに含まれる
# c = mod2のみに含まれる
# d = どちらにも含まれない
a <- intersection
b <- length(metabolites1) - a
c <- length(metabolites2) - a
d <- universe_size - a - b - c
# Fisher's exact testの実行
contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
# p値の計算
if (universe_size > 0 && (a + b) > 0 && (a + c) > 0) { fisher_result <- tryCatch({ fisher.test(contingency_table, alternative = "greater") }, error = function(e) { list(p.value = 1) }) p_value <- fisher_result$p.value } else { p_value <- 1 }
# 結果を保存
results_list[[length(results_list) + 1]] <- data.frame( Module_Study1 = mod1, Module_Study2 = mod2, Size_Study1 = length(metabolites1), Size_Study2 = length(metabolites2), Common_Metabolites = intersection, Total_Union = union, Jaccard_Coefficient = round(jaccard, 4), P_value = format(p_value, scientific = TRUE, digits = 4), P_value_numeric = p_value, Significant = ifelse(p_value < 0.05, "Yes", "No"), stringsAsFactors = FALSE ) } }
# データフレームに結合
preservation_scores <- do.call(rbind, results_list)
# 列名を研究名に変更
colnames(preservation_scores)[1:2] <- c( paste0("Module_", study1_name), paste0("Module_", study2_name) )
# Jaccard係数でソート（降順）
preservation_scores <- preservation_scores[order(-preservation_scores$Jaccard_Coefficient), ]
# P値の数値列を削除（表示用のフォーマット済み列のみ残す）
preservation_scores_output <- preservation_scores[, -9]
# 全体を保存
write.csv(preservation_scores_output, "Module_Preservation_Jaccard_Pvalue_Full.csv", row.names = FALSE)
cat("全ペアの保存性スコアを保存: Module_Preservation_Jaccard_Pvalue_Full.csv\\n")
# 有意なペアのみ抽出（p < 0.05 または Jaccard > 0.1）
significant_pairs <- preservation_scores[ preservation_scores$P_value_numeric < 0.05 | preservation_scores$Jaccard_Coefficient > 0.1, ] significant_pairs_output <- significant_pairs[, -9]
if (nrow(significant_pairs) > 0) {
    write.csv(significant_pairs_output,
             "Module_Preservation_Jaccard_Pvalue_Significant.csv",
             row.names = FALSE)
    cat("有意なペア（p < 0.05 または Jaccard > 0.1）を保存: Module_Preservation_Jaccard_Pvalue_Significant.csv\\n")
    cat("有意なペア数:", nrow(significant_pairs), "\\n\\n")
} else {
    cat("有意なペアはありませんでした\\n\\n")
}
# サマリー統計
cat("=== 保存性スコアの統計 ===\n") cat("平均 Jaccard係数:", round(mean(preservation_scores$Jaccard_Coefficient), 4), "\n") cat("中央値 Jaccard係数:", round(median(preservation_scores$Jaccard_Coefficient), 4), "\n") cat("最大 Jaccard係数:", round(max(preservation_scores$Jaccard_Coefficient), 4), "\n") cat("最小 Jaccard係数:", round(min(preservation_scores$Jaccard_Coefficient), 4), "\n\n")
cat("p < 0.05 のペア数:", sum(preservation_scores$P_value_numeric < 0.05), "\\n")
cat("p < 0.01 のペア数:", sum(preservation_scores$P_value_numeric < 0.01), "\\n")
cat("p < 0.001 のペア数:", sum(preservation_scores$P_value_numeric < 0.001), "\\n\\n")
# 上位10ペアを表示（Jaccard係数順）
cat("保存性が最も高い上位10ペア:\n") print(head(preservation_scores[, c(1, 2, 5, 7, 8, 10)], 10))
# 統計的に最も有意なペア（p値順）
preservation_scores_pvalue_sorted <- preservation_scores[ order(preservation_scores$P_value_numeric), ] cat("\n統計的に最も有意なペア（上位10）:\n") print(head(preservation_scores_pvalue_sorted[, c(1, 2, 5, 7, 8, 10)], 10))
# マトリックス形式で保存（Jaccard係数）
jaccard_matrix <- matrix(0, nrow = length(unique_modules1), ncol = length(unique_modules2), dimnames = list(unique_modules1, unique_modules2))
# p値のマトリックス
pvalue_matrix <- matrix(1, nrow = length(unique_modules1), ncol = length(unique_modules2), dimnames = list(unique_modules1, unique_modules2))
for (i in 1:nrow(preservation_scores)) {
    mod1 <- as.character(preservation_scores[i, 1])
    mod2 <- as.character(preservation_scores[i, 2])
    jaccard_matrix[mod1, mod2] <- preservation_scores$Jaccard_Coefficient[i]
    pvalue_matrix[mod1, mod2] <- preservation_scores$P_value_numeric[i]
}
# Jaccardマトリックス保存
jaccard_matrix_df <- as.data.frame(jaccard_matrix) jaccard_matrix_df <- cbind(Module = rownames(jaccard_matrix_df), jaccard_matrix_df) write.csv(jaccard_matrix_df, "Module_Preservation_Jaccard_Matrix.csv", row.names = FALSE)
# p値マトリックス保存
pvalue_matrix_df <- as.data.frame(pvalue_matrix) pvalue_matrix_df <- cbind(Module = rownames(pvalue_matrix_df), pvalue_matrix_df) write.csv(pvalue_matrix_df, "Module_Preservation_Pvalue_Matrix.csv", row.names = FALSE)
cat("\\nマトリックス形式を保存:\\n")
cat("  - Module_Preservation_Jaccard_Matrix.csv\\n")
cat("  - Module_Preservation_Pvalue_Matrix.csv\\n")
# 各研究のベストマッチを抽出（Jaccard係数基準）
best_matches1 <- data.frame() for (mod1 in unique_modules1) { subset_data <- preservation_scores[preservation_scores[[1]] == mod1, ] best_match <- subset_data[which.max(subset_data$Jaccard_Coefficient), ] best_matches1 <- rbind(best_matches1, best_match) }
best_matches2 <- data.frame()
for (mod2 in unique_modules2) {
    subset_data <- preservation_scores[preservation_scores[[2]] == mod2, ]
    best_match <- subset_data[which.max(subset_data$Jaccard_Coefficient), ]
    best_matches2 <- rbind(best_matches2, best_match)
}

best_matches1_output <- best_matches1[, -9]
best_matches2_output <- best_matches2[, -9]

write.csv(best_matches1_output,
         paste0("Best_Matches_from_", study1_name, ".csv"),
         row.names = FALSE)
write.csv(best_matches2_output,
         paste0("Best_Matches_from_", study2_name, ".csv"),
         row.names = FALSE)

cat("\\n各モジュールのベストマッチを保存:\\n")
cat("  - Best_Matches_from_", study1_name, ".csv\\n")
cat("  - Best_Matches_from_", study2_name, ".csv\\n")

return(list(
    full_scores = preservation_scores,
    significant_pairs = significant_pairs,
    jaccard_matrix = jaccard_matrix,
    pvalue_matrix = pvalue_matrix,
    best_matches_study1 = best_matches1,
    best_matches_study2 = best_matches2
))
}
module_df1 <- read.csv("Module_Metabolites_MTBLS406.csv", stringsAsFactors = FALSE) module_df2 <- read.csv("Module_Metabolites_MTBLS980.csv", stringsAsFactors = FALSE) jaccard_results <- quantify_module_preservation(module_df1, module_df2) cat("\n========================================") cat("\n解析完了！以下のファイルが作成されました：") cat("\n========================================\n") cat("1. Module_Preservation_Jaccard_Pvalue_Full.csv - 全ペアのJaccard係数とp値\n") cat("2. Module_Preservation_Jaccard_Pvalue_Significant.csv - 有意なペアのみ\n") cat("3. Module_Preservation_Jaccard_Matrix.csv - Jaccardマトリックス\n") cat("4. Module_Preservation_Pvalue_Matrix.csv - p値マトリックス\n") cat("5. Best_Matches_from_MTBLS406.csv - 各モジュールの最良マッチ\n") cat("6. Best_Matches_from_MTBLS980.csv - 各モジュールの最良マッチ\n") cat("\n指標の解釈:\n") cat(" Jaccard係数: 1.0 = 完全一致, 0.0 = 共通要素なし\n") cat(" p値: < 0.05 = 統計的に有意な重複\n") cat("========================================\n")
