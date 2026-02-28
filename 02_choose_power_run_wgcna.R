choose_soft_threshold <- function(expr_data, study_name) { cat("\n=== ", study_name, "のソフト閾値選択 ===\n") 
	powers <- c(1:10, seq(12, 20, by = 2))
	sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5) 
	par(mfrow = c(1, 2)) 
	plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", main = paste(study_name, ": Scale independence")) 
	abline(h = 0.85, col = "red") text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = 0.9, col = "red") 
	plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", main = paste(study_name, ": Mean connectivity")) 
	text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red") 
	recommended_power <- sft$fitIndices[sft$fitIndices[, 2] > 0.85, 1][1] if (is.na(recommended_power)) { recommended_power <- sft$fitIndices[which.max(sft$fitIndices[, 2]), 1] }
	 cat("推奨されるソフト閾値パワー:", recommended_power, "\n") return(list(power = recommended_power, sft = sft)) }
	 # 各データセットでソフト閾値を選択
	 sft1 <- choose_soft_threshold(expr1, "MTBLS406")
	 sft2 <- choose_soft_threshold(expr2, "MTBLS980")
	 
	 run_wgcna <- function(expr_data, power, study_name) { cat("\n=== ", study_name, "のWGCNA解析実行 ===\n")
# WGCNAのcor関数を明示的に使用
cor <- WGCNA::cor
# ネットワーク構築とモジュール検出
net <- blockwiseModules( expr_data, power = power, TOMType = "unsigned", minModuleSize = 5,
# 最小モジュールサイズ（代謝物数が少ない場合は調整）
reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = paste0("TOM_", study_name), corType = "pearson",
# 相関タイプを明示的に指定
maxBlockSize = 5000,
# ブロックサイズを指定
verbose = 3 )
# モジュールカラーの割り当て
moduleLabels <- net$colors moduleColors <- labels2colors(net$colors)
new_df <- data.frame(Metabolite = colnames(expr_data), Module = moduleColors) write.csv(new_df, paste0("Module_Metabolites_", study_name, ".csv"), row.names = FALSE) 
cat("✓ 最新のモジュール情報を保存:", paste0("Module_Metabolites_", study_name, ".csv\n"))
# モジュール数の表示
n_modules <- length(unique(moduleColors)) - 1
# グレー（モジュール未割り当て）を除く
cat("検出されたモジュール数:", n_modules, "\n") cat("モジュールカラー:", unique(moduleColors), "\n")
# モジュールごとの代謝物数
module_sizes <- table(moduleColors) print(module_sizes)
return(list(
    net = net,
    moduleLabels = moduleLabels,
    moduleColors = moduleColors,
    MEs = net$MEs
))
}
wgcna1 <- run_wgcna(expr1, sft1$power, "MTBLS406")
wgcna2 <- run_wgcna(expr2, sft2$power, "MTBLS980")

plot_dendrogram <- function(expr_data, wgcna_result, study_name) { cat("\n=== ", study_name, "のデンドログラム作成 ===\n")
# 階層的クラスタリングの計算
dissTOM <- 1 - TOMsimilarityFromExpr(expr_data, power = sft1$power) geneTree <- hclust(as.dist(dissTOM), method = "average")
# プロット
png(filename = paste0("Dendrogram_", study_name, ".png"), width = 12, height = 6, units = "in", res = 300)
plotDendroAndColors(
    geneTree,
    wgcna_result$moduleColors,
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste(study_name, "- Dendrogram and Module Colors")
)

dev.off()
cat("デンドログラムを保存:", paste0("Dendrogram_", study_name, ".png\\n"))
}
# 各データセットでデンドログラム作成
plot_dendrogram(expr1, wgcna1, "MTBLS406")
plot_dendrogram(expr2, wgcna2, "MTBLS980")