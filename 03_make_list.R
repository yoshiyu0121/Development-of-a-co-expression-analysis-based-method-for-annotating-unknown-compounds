create_expanded_module_summary <- function(module_df, study_name) { cat("\n=== ", study_name, "の展開型モジュールサマリー作成 ===\n")
# グレーモジュールを除く
modules_no_grey <- module_df[module_df$Module != "grey", ]
# モジュールごとの代謝物数を計算
module_counts <- as.data.frame(table(modules_no_grey$Module)) colnames(module_counts) <- c("Module", "Total_Count")
# 元のデータにカウント情報を追加
expanded_df <- merge(modules_no_grey, module_counts, by = "Module")
# モジュールでソート
expanded_df <- expanded_df[order(expanded_df$Module, expanded_df$Metabolite), ]
# 列の順序を整理
expanded_df <- expanded_df[, c("Module", "Total_Count", "Metabolite")]
# CSVとして保存
write.csv(expanded_df, file = paste0("Module_Summary_Expanded_", study_name, ".csv"), row.names = FALSE)
cat("展開型サマリーを保存:", paste0("Module_Summary_Expanded_", study_name, ".csv\\n"))
# サマリー統計を表示
cat("\nモジュール統計:\n") print(module_counts)
cat("\\n合計モジュール数:", nrow(module_counts), "\\n")
cat("合計代謝物数:", sum(module_counts$Total_Count), "\\n")

return(expanded_df)
}
# 両方のデータセットで実行
module_df1 <- read.csv("Module_Metabolites_MTBLS406.csv") 
module_df2 <- read.csv("Module_Metabolites_MTBLS980.csv")
expanded1 <- create_expanded_module_summary(module_df1, "MTBLS406")
expanded2 <- create_expanded_module_summary(module_df2, "MTBLS980")

#全代謝物名のリスト化
all_metabolites1 <- module_df1$Metabolite 
all_metabolites2 <- module_df2$Metabolite
#集合演算で共通と固有を定義
common_metabolites <- intersect(all_metabolites1, all_metabolites2) 
unique_to_study1 <- setdiff(all_metabolites1, all_metabolites2) 
unique_to_study2 <- setdiff(all_metabolites2, all_metabolites1)
#共通代謝物の詳細比較データの作成
common_detailed <- data.frame( Metabolite = common_metabolites, Module_Study1 = module_df1$Module[match(common_metabolites, module_df1$Metabolite)], Module_Study2 = module_df2$Module[match(common_metabolites, module_df2$Metabolite)], stringsAsFactors = FALSE )
#同じ色名かのグラフ作成
common_detailed$Same_Module <- (common_detailed$Module_Study1 == common_detailed$Module_Study2)
sule_count <- sum(common_detailed$Same_Module, na.rm = TRUE)
if (length(unique_to_study1) > 0) { unique_df1 <- module_df1[module_df1$Metabolite %in% unique_to_study1, c("Metabolite", "Module")] unique_df1 <- unique_df1[order(unique_df1$Module, unique_df1$Metabolite), ]
    write.csv(unique_df1,
             paste0("Unique_to_", study1_name, ".csv"),
             row.names = FALSE)
    
    cat(study1_name, "固有の代謝物リストを保存: Unique_to_", study1_name, ".csv\\n")
}

if (length(unique_to_study2) > 0) {
    unique_df2 <- module_df2[module_df2$Metabolite %in% unique_to_study2, 
                             c("Metabolite", "Module")]
    unique_df2 <- unique_df2[order(unique_df2$Module, unique_df2$Metabolite), ]
    
    write.csv(unique_df2,
             paste0("Unique_to_", study2_name, ".csv"),
             row.names = FALSE)
    
    cat(study2_name, "固有の代謝物リストを保存: Unique_to_", study2_name, ".csv\\n\\n")
}
summary_table <- data.frame( Category = c("Total_Metabolites", "Common_Metabolites", paste0("Unique_to_", study1_name), paste0("Unique_to_", study2_name)), Study1_Count = c(length(all_metabolites1), length(common_metabolites), length(unique_to_study1), 0), Study2_Count = c(length(all_metabolites2), length(common_metabolites), 0, length(unique_to_study2)), stringsAsFactors = FALSE )
colnames(summary_table) <- c("Category", study1_name, study2_name)

write.csv(summary_table,
         "Metabolite_Comparison_Summary.csv",
         row.names = FALSE)

cat("比較サマリーを保存: Metabolite_Comparison_Summary.csv\\n")
png("Metabolite_Overlap_Venn.png", width = 10, height = 8, units = "in", res = 300)
par(mar = c(2, 2, 3, 2)) plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-1.5, 2),
xlab = "", ylab = "", axes = FALSE, asp = 1,
main = "Metabolite Overlap Between Studies")
theta <- seq(0, 2*pi, length.out = 100) x1 <- cos(theta) - 0.5 y1 <- sin(theta) x2 <- cos(theta) + 0.5 y2 <- sin(theta) polygon(x1, y1, col = rgb(1, 0, 0, 0.3), border = "red", lwd = 2) polygon(x2, y2, col = rgb(0, 0, 1, 0.3), border = "blue", lwd = 2)
text(-0.8, 0, paste0(study1_name, "\nOnly\n", length(unique_to_study1)), cex = 1.2, font = 2) text(0.8, 0, paste0(study2_name, "\nOnly\n", length(unique_to_study2)), cex = 1.2, font = 2) text(0, 0, paste0("Common\n", length(common_metabolites)), cex = 1.5, font = 2) text(-1.5, -1.2, study1_name, col = "red", cex = 1.1, font = 2) text(1.5, -1.2, study2_name, col = "blue", cex = 1.1, font = 2)
dev.off() quartz_off_screen 2
cat("✓ ベン図を保存しました: Metabolite_Overlap_Venn.png\n")
