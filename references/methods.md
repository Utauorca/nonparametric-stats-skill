# 無母數檢定速查

每個方法列出：虛無假設、檢定統計量、ties 處理、效果量、scipy 函式。
涵蓋本專案 slides.qmd（底本 Hollander, Wolfe & Chicken 2015）的 6 個方法；
完整假設、公式與 R code 見 `textbook-methods.md`。

## Wilcoxon 符號等級檢定 (signed-rank) — 單樣本 / 配對位置
- H0：配對差值分布對稱中心為 0。**假設差值分布對稱**。
- 統計量：T⁺ = 正差的等級和（對 |Z| 排名後取正差的等級和）。
- Ties：差值=0 多半剔除並縮減 n；|Z| 相同取平均等級（檢定變近似，變異數縮小）。
- 效果量：matched-pairs rank-biserial；位置估計用 Hodges-Lehmann（配對差中位數）。
- scipy：`scipy.stats.wilcoxon`。

## Mann-Whitney U / Wilcoxon rank-sum — 兩獨立樣本位置
- H0：兩組來自同分布（位置移動模型 Δ=0）。**要解讀成中位數差需兩組形狀相近**，否則只能說隨機優勢 P(X<Y)。
- 統計量：W = 第二組 Y 的等級和；等價 Mann-Whitney 計數 U，W = U + n(n+1)/2。
- Ties：等級取平均，建議用常態近似校正。
- 效果量：rank-biserial r = 1 - 2U/(n1·n2)；位置差用 Hodges-Lehmann。
- scipy：`scipy.stats.mannwhitneyu`。

## Ansari-Bradley — 兩獨立樣本散布 (scale)
- H0：尺度比 γ² = 1（兩組散布相同）。**額外要求兩母體中位數相等**，否則結果會混進位置差異。
- 統計量：合併排序後由兩端往中間給分（最小與最大=1、次小與次大=2…），C = Y 的分數和。較集中的樣本落中央得高分，較分散落兩端得低分。
- Ties：取平均分數，檢定變近似（變異數公式隨 N 奇偶不同）。
- 效果量：尺度比估計 γ² = var(X)/var(Y)。
- scipy：`scipy.stats.ansari`。

## Kolmogorov-Smirnov 雙樣本 (Smirnov) — 兩分布一般差異
- H0：F(t)=G(t) 對所有 t（偵測任何差異，天生雙尾）。
- 統計量：D = max_t |F_m(t)−G_n(t)|（兩經驗分布的最大垂直落差）；精確檢定用整數化的 J = (mn/d)·D，d=gcd(m,n)。
- Ties：連續母體假設下理論上無 ties；實務有 ties 時 p 值為近似。
- 效果量：D 本身即落差大小（0~1）可當效果量。
- scipy：`scipy.stats.ks_2samp`。

## Kruskal-Wallis — k 組獨立位置
- H0：k 組來自同分布（位置模型 τ 全相等）。
- 統計量：H，**合併全部 N 個觀測排名**後衡量各組平均等級偏離 (N+1)/2 的程度，近似 χ²(k−1)。k=2 時退化為雙尾 rank-sum。
- Ties：取平均等級後除以收縮因子校正 H′。
- 事後（顯著後）：Dunn's test + Bonferroni / Holm 校正（`scikit-posthocs`）。
- 效果量：epsilon² = (H − k + 1)/(n − k) 或 eta²_H。
- scipy：`scipy.stats.kruskal`。

## Friedman — k 組相依 / 重複測量 (RCBD)
- H0：k 個 treatment 效果相同。
- 統計量：S，**在每個 block 內各自排名** 1…k（消掉 block 效果）後衡量各 treatment 平均等級偏離 (k+1)/2 的程度，近似 χ²(k−1)。k=2 時退化為符號檢定。
  - 與 Kruskal-Wallis 的關鍵差異：KW 聯合排名；Friedman 在每個 block 內排名。
- Ties：block 內取平均等級，用 tie-corrected S。
- 事後（顯著後）：Nemenyi test（`scikit-posthocs`）。
- 效果量：Kendall's W（一致性係數）。
- scipy：`scipy.stats.friedmanchisquare`。
