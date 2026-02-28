# Development-of-a-co-expression-analysis-based-method-for-annotating-unknown-compounds
# 非標的メタボロミクスにおけるネットワーク解析プログラム (DSPC & WGCNA)

本リポジトリは、学士論文 「非標的メタボロミクスにおけるネットワーク解析を用いた未知化合物のアノテーション手法の比較検証」で使用した解析プログラム一式を公開するものです。

## 概要
本研究では、共発現解析手法（WGCNA）を行い、モジュールの保存性を検証しました。
さらに、最適な共発現手法の検討として共発現解析手法（WGCNA）、GraphicalLASSO、偏相関解析手法（DSPC）の比較検証を行いました。特にDSPCについては、先行研究（Basu et al., 2018）のアルゴリズムに基づき、R言語を用いて独自に実装を行っています。

## ディレクトリ構成
解析のステップ順にファイルを分割しています。

1. **01_preprocessing.R** - データの読み込み、欠損値処理、`scale()` による標準化。
2. **02_choose_power_run_wgcna.R** - ソフト閾値の選定、モジュール抽出。
3. **03_make_list.R** -代謝物リストの作成。
4. **04_Jaccard_index.R** -およびJaccard係数による再現性評価。
5. **05_GraphicalLASSO_DSPC.R** -Graphical Lassoによる推定と、Debiasing処理（バイアス補正）による不偏推定量の算出、およびp値の計算。

## 使用言語・主なパッケージ
- **Language**: R (version [4.5.0 等])
- **Packages**: 
  - `glasso`: スパース精度行列の算出に使用
  - `WGCNA`: 共発現ネットワーク解析に使用
  - `ggplot2`: 図表の可視化に使用
  - `reshape2`: データ整形に使用

## 独自実装のポイント (DSPC)
`05__GraphicalLASSO_DSPC.R` では、Basu et al. (2018) の提唱する理論に基づき、以下のステップを実装しています。
- 正則化パラメータ $\lambda = \sqrt{\log(p)/n}$ の設定
- 行列演算によるDebiasingステップ
- 標準正規分布に基づくp値の算出と多重検定補正（FDR制御）

## 免責事項・ライセンス
本プログラムは学術研究を目的として作成されたものです。
ライセンス：MIT License
