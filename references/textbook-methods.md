# 無母數方法詳解（教材版）

> 本檔內容萃取自本專案的 `slides.qmd`，底本為
> Hollander, Wolfe & Chicken (2015), *Nonparametric Statistical Methods* (3rd ed.)。
> 回答時可引用此處的假設、統計量、ties 處理與判斷規則作為依據。

---

## 1. 單樣本 / 配對位置問題 — Wilcoxon 符號等級檢定 (signed-rank)

**資料與假設.** 對 $n$ 個受試者取配對 $(X_i, Y_i)$，化為差值 $Z_i = Y_i - X_i$。
- A1：$Z_1,\ldots,Z_n$ 互相獨立。
- A2：每個 $Z_i$ 來自連續母體，且**對共同中位數 $\theta$ 對稱**：$F_i(\theta+t)+F_i(\theta-t)=1$。
  （$\theta$ = treatment effect；純單樣本時 $Z_i$ 即原始觀測值，檢定中位數是否等於 $\theta_0$。）

**假設檢定.** $H_0:\theta=0$；對立：$H_1:\theta>0$、$H_2:\theta<0$、$H_3:\theta\neq 0$。

**統計量.** 對 $|Z_i|$ 由小到大排名得 $R_i$，令 $\psi_i=1$ 若 $Z_i>0$ 否則 $0$：
$$T^{+}=\sum_{i=1}^{n}\psi_i R_i \quad(\text{正差的等級和,兼含符號與大小}).$$

**Ties.** 差值為 0 → 剔除並縮減 $n$；$|Z|$ 相同 → 取平均等級（檢定變近似）。
有 ties 時 $E_0(T^+)$ 不變，變異數縮小：
$$\operatorname{Var}_0(T^{+})=\tfrac{1}{24}\big[n(n+1)(2n+1)-\tfrac12\sum_j t_j(t_j-1)(t_j+1)\big].$$

**判斷（精確）.** $t_\alpha$ 為上 $\alpha$ 臨界值（R：`qsignrank`/`psignrank`）。
| 對立 | 拒絕 $H_0$ 若 |
|---|---|
| $H_1:\theta>0$ | $T^{+}\ge t_\alpha$ |
| $H_2:\theta<0$ | $T^{+}\le \tfrac{n(n+1)}{2}-t_\alpha$ |
| $H_3:\theta\neq0$ | $T^{+}\ge t_{\alpha/2}$ 或 $\le \tfrac{n(n+1)}{2}-t_{\alpha/2}$ |

**大樣本近似.** $E_0(T^+)=\tfrac{n(n+1)}{4}$，$\operatorname{Var}_0(T^+)=\tfrac{n(n+1)(2n+1)}{24}$，
$T^{*}=\dfrac{T^{+}-E_0}{\sqrt{\operatorname{Var}_0}}\xrightarrow{d}N(0,1)$。

**R.** `wilcox.test(post, pre, paired = TRUE, alternative = "less")`
（等價：`wilcox.test(post - pre, alternative = "less")`）。

---

## 2. 兩樣本位置問題 — Wilcoxon 等級和 / Mann–Whitney U

**資料與假設.** $X_1,\ldots,X_m$、$Y_1,\ldots,Y_n$，各自 i.i.d.、兩組互相獨立、連續母體。
位置移動模型 $G(t)=F(t-\Delta)$：兩分布**形狀相同、只差平移 $\Delta$**。

**假設檢定.** $H_0:\Delta=0$（即 $F=G$）；對立 $H_1:\Delta>0$、$H_2:\Delta<0$、$H_3:\Delta\neq0$。

**統計量.** 合併排名；$S_j$ = $Y_j$ 的聯合排名。
$$W=\sum_{j=1}^{n}S_j \quad(\text{第二組 }Y\text{ 的等級和}).$$
等價形式：Mann–Whitney 計數 $U=\sum_i\sum_j \mathbb{1}(X_i<Y_j)$，與 $W=U+\tfrac{n(n+1)}{2}$。

**Ties.** 取平均等級，檢定變近似；$E_0(W)$ 不變，變異數縮小：
$$\operatorname{Var}_0(W)=\tfrac{mn}{12}\big[m+n+1-\tfrac{\sum_j(t_j-1)t_j(t_j+1)}{(m+n)(m+n-1)}\big].$$

**判斷（精確）.** $w_\alpha$ 上 $\alpha$ 臨界值（R：`qwilcox`）。
| 對立 | 拒絕 $H_0$ 若 |
|---|---|
| $H_1:\Delta>0$ | $W\ge w_\alpha$ |
| $H_2:\Delta<0$ | $W\le n(m+n+1)-w_\alpha$ |
| $H_3:\Delta\neq0$ | $W\ge w_{\alpha/2}$ 或 $\le n(m+n+1)-w_{\alpha/2}$ |

**大樣本近似.** $E_0(W)=\tfrac{n(m+n+1)}{2}$，$\operatorname{Var}_0(W)=\tfrac{mn(m+n+1)}{12}$，標準化後 $\to N(0,1)$。

**R.** `wilcox.test(treatment, control)`（雙尾，小樣本自動用精確分布）。

**重要解讀界線.** 要把結果講成「中位數差異」需兩組**形狀相近**；否則只能解讀為隨機優勢 $P(X<Y)$。

---

## 3. 兩樣本散布問題 — Ansari–Bradley 檢定（scale）

**資料與假設.** 同兩樣本 i.i.d. + 互相獨立 + 連續；**額外要求 A3：兩母體中位數相等** $\theta_1=\theta_2$
（讓檢定只反映散布差異）。位置-尺度模型 $F(t)=H(\tfrac{t-\theta_1}{\eta_1})$、$G(t)=H(\tfrac{t-\theta_2}{\eta_2})$。

**假設檢定.** 尺度比 $\gamma=\eta_1/\eta_2$（若變異數存在，$\gamma^2=\operatorname{var}(X)/\operatorname{var}(Y)$）。
$H_0:\gamma^2=1$；對立 $H_1:\gamma^2>1$、$H_2:\gamma^2<1$、$H_3:\gamma^2\neq1$。

**統計量.** 合併排序後**由兩端往中間給分**：最小與最大給 1、次小與次大給 2…
令 $R_j$ 為 $Y_j$ 的分數，$C=\sum_{j=1}^n R_j$。
（概念：較集中的樣本落在中央＝高分，較分散的落在兩端＝低分。）

**Ties.** 取平均分數，檢定變近似（變異數公式隨 $N$ 奇偶不同，見 slides）。

**判斷（精確）.** $c_\alpha$ 上 $\alpha$ 臨界值。大樣本近似的 $E_0(C)$、$\operatorname{Var}_0(C)$ **依 $N$ 奇偶而異**：
- $N$ 偶：$E_0(C)=\tfrac{n(N+2)}{4}$，$\operatorname{Var}_0(C)=\tfrac{mn(N+2)(N-2)}{48(N-1)}$。
- $N$ 奇：$E_0(C)=\tfrac{n(N+1)^2}{4N}$，$\operatorname{Var}_0(C)=\tfrac{mn(N+1)(3+N^2)}{48N^2}$。

**R.** `ansari.test(y, x)`（`alternative` 可指定單尾）。

---

## 4. 兩分布一般差異 — Kolmogorov–Smirnov 雙樣本檢定 (Smirnov)

**假設.** 兩獨立樣本、連續母體，分布函數 $F$、$G$。

**假設檢定.** $H_0:F(t)=G(t)\ \forall t$ vs $H_1:F(t)\neq G(t)$ 至少一個 $t$（偵測**任何**差異，天生雙尾）。

**統計量.** 經驗分布 $F_m,G_n$，取最大垂直落差：
$$J=\tfrac{mn}{d}\max_t|F_m(t)-G_n(t)|,\quad d=\gcd(m,n).$$
最大值只需在 $N=m+n$ 個資料點檢查。

**判斷.** $J\ge j_\alpha$ 時拒絕。大樣本：以 $D=\max_t|F_m-G_n|$，$\sqrt{mn/N}\,D$ 比對 Kolmogorov 臨界值。

**R.** `ks.test(x, y)`。

---

## 5. 單因子配置（$k$ 組獨立）— Kruskal–Wallis

**假設.** $N=\sum n_j$ 個觀測分 $k$ 組；互相獨立；各組 i.i.d. 連續；位置模型 $F_j(t)=F(t-\tau_j)$。
（即 rank 版 one-way ANOVA：$X_{ij}=\theta+\tau_j+e_{ij}$。）

**假設檢定.** $H_0:\tau_1=\cdots=\tau_k$ vs $H_1$ 至少兩者不等。

**統計量.** **合併全部 $N$ 個觀測排名**，$R_j$ = 第 $j$ 組等級和，$\bar R_{\cdot j}=R_j/n_j$：
$$H=\tfrac{12}{N(N+1)}\sum_j n_j\big(\bar R_{\cdot j}-\tfrac{N+1}{2}\big)^2=\tfrac{12}{N(N+1)}\sum_j \tfrac{R_j^2}{n_j}-3(N+1).$$
$H$ 大則拒絕（$k=2$ 時退化成雙尾 Wilcoxon rank-sum）。

**Ties.** 取平均等級後除以收縮因子校正：
$$H'=\dfrac{H}{1-\dfrac{\sum_j(t_j^3-t_j)}{N^3-N}}.$$

**判斷.** 精確：$H\ge h_\alpha$（R：`cKW`）。大樣本：$H\xrightarrow{d}\chi^2_{k-1}$，$H\ge\chi^2_{k-1,\alpha}$ 拒絕。

**事後（顯著後）.** Dunn's test + 多重比較校正（Bonferroni/Holm）。

**R.** `kruskal.test(value ~ group)`。

---

## 6. 雙因子配置（$k$ 組相依 / RCBD）— Friedman 檢定

**假設.** 隨機完全區集設計：$n$ 個 block、每 block 有 $k$ 個 treatment 各一觀測。
互相獨立、連續；加法模型 $F_{ij}(u)=F(u-\beta_i-\tau_j)$（$\beta_i$ = block 效果，$\tau_j$ = treatment 效果）。
（即 rank 版 two-way / repeated-measures ANOVA。）

**假設檢定.** $H_0:\tau_1=\cdots=\tau_k$ vs $H_1$ 至少兩者不等。

**統計量.** **在每個 block 內各自排名** $1\ldots k$（這樣消掉 block 效果 $\beta_i$），$R_j$ = treatment $j$ 的等級和：
$$S=\tfrac{12n}{k(k+1)}\sum_j\big(\bar R_{\cdot j}-\tfrac{k+1}{2}\big)^2=\tfrac{12}{nk(k+1)}\sum_j R_j^2-3n(k+1).$$
$S$ 大則拒絕（$k=2$ 時退化成符號檢定）。

> **與 Kruskal–Wallis 的關鍵差異**：KW 把全部觀測**聯合排名**；Friedman **在每個 block 內排名**。

**Ties.** block 內取平均等級，並用 tie-corrected $S$（分母扣 $\tfrac{1}{k-1}\sum_i\sum_j(t_{ij}^3-t_{ij})$）。

**判斷.** 精確：$S\ge s_\alpha$。大樣本：$S\xrightarrow{d}\chi^2_{k-1}$。

**事後（顯著後）.** Nemenyi test。

**R.** `friedman.test(times)`（`times` 為 block×treatment 矩陣）。
