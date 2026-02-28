library(glasso)
data_scaled_980 <- scale(data_for_glasso_980)
n_980 <- nrow(data_scaled_980)
p_980 <- ncol(data_scaled_980)
cat("===================================\n")
cat("MTBLS980 データ\n")
cat("===================================\n")
cat("サンプル数:", n_980, "\n")
cat("メタボライト数:", p_980, "\n\n")
S_980 <- cor(data_scaled_980)
lambda_980 <- sqrt(log(p_980) / n_980)
cat("lambda =", lambda_980, "\n")
cat("Graphical Lassoを実行中...\n")
glasso_result_980 <- glasso(S_980, rho = lambda_980)
precision_980 <- glasso_result_980$wi
D_980 <- diag(1 / sqrt(diag(precision_980)))
partial_cor_980 <- -D_980 %% precision_980 %% D_980
diag(partial_cor_980) <- 1
rownames(partial_cor_980) <- colnames(partial_cor_980) <-
colnames(data_for_glasso_980)
cat("Graphical Lasso完了\n")
cat("ゼロでないエッジ数:",
sum(abs(partial_cor_980[upper.tri(partial_cor_980)]) > 1e-10), "\n\n")
cat("Debiasingステップを実行中...\n")
Theta_hat_980 <- precision_980
Theta_inv_980 <- solve(Theta_hat_980)
residual_980 <- S_980 - Theta_inv_980
T_hat_980 <- Theta_hat_980 - Theta_hat_980 %% Theta_hat_980 %%
residual_980
rownames(T_hat_980) <- colnames(T_hat_980) <-
colnames(data_for_glasso_980)
cat("Debiasing完了\n\n")
cat("p値を計算中...\n")
sigma_hat_980 <- sqrt(diag(Theta_hat_980 %*% Theta_hat_980) / n_980)
sigma_matrix_980 <- outer(sigma_hat_980, sigma_hat_980, "*")
z_scores_980 <- T_hat_980 / sigma_matrix_980
p_values_980 <- 2 * (1 - pnorm(abs(z_scores_980)))
diag(p_values_980) <- 1
rownames(p_values_980) <- colnames(p_values_980) <-
colnames(data_for_glasso_980)
cat("p値計算完了\n\n")
D_debiased_980 <- diag(1 / sqrt(abs(diag(T_hat_980))))
partial_cor_debiased_980 <- -D_debiased_980 %% T_hat_980 %%
D_debiased_980
diag(partial_cor_debiased_980) <- 1
rownames(partial_cor_debiased_980) <- colnames(partial_cor_debiased_980)
<- colnames(data_for_glasso_980)
idx_980 <- which(upper.tri(partial_cor_debiased_980), arr.ind = TRUE)
dspc_edge_list_980 <- data.frame(
Node1 = rownames(partial_cor_debiased_980)[idx_980[, 1]],
Node2 = colnames(partial_cor_debiased_980)[idx_980[, 2]],
Partial_Correlation = partial_cor_debiased_980[idx_980],
P_value = p_values_980[idx_980],
Abs_Partial_Correlation = abs(partial_cor_debiased_980[idx_980]),
stringsAsFactors = FALSE
)
cat("多重検定補正を実行中...\n")
dspc_edge_list_980$Adjusted_P_value <-
p.adjust(dspc_edge_list_980$P_value, method = "BH")
dspc_significant_980 <-
dspc_edge_list_980[dspc_edge_list_980$Adjusted_P_value < 0.05, ]
dspc_significant_980 <-
dspc_significant_980[order(dspc_significant_980$Adjusted_P_value), ]
rownames(dspc_significant_980) <- NULL
cat("===================================\n")
cat("MTBLS980 DSPC結果\n")
cat("===================================\n")
cat("全エッジ数:", nrow(dspc_edge_list_980), "\n")
cat("有意なエッジ数 (adjusted p < 0.05):", nrow(dspc_significant_980),
"\n")
if(nrow(dspc_significant_980) > 0) {
cat("正の相関:", sum(dspc_significant_980$Partial_Correlation > 0),
"\n")
cat("負の相関:", sum(dspc_significant_980$Partial_Correlation < 0),
"\n\n")
cat("上位10の有意なエッジ:\n")
print(head(dspc_significant_980[, c("Node1", "Node2",
"Partial_Correlation",
"P_value", "Adjusted_P_value")],
10))
}
