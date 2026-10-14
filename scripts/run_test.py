"""無母數檢定計算腳本。

用 scipy 實際計算常見無母數檢定的統計量、p 值與效果量，
避免模型心算或憑空給數字。資料從 CSV 讀入 (long 格式: 一欄值、一欄組別)。

範例:
    # 兩獨立樣本 Mann-Whitney U
    python run_test.py mannwhitney --csv data.csv --value score --group grp

    # 配對 Wilcoxon signed-rank (寬格式: 兩欄為前後測)
    python run_test.py wilcoxon --csv data.csv --col-a pre --col-b post

    # 多組 Kruskal-Wallis
    python run_test.py kruskal --csv data.csv --value score --group grp

    # EDA 視覺化 (檢定前務必先畫；圖型 auto 依參數判斷)
    python run_test.py plot --csv data.csv --value score --group grp   # 組間 → 箱型圖
    python run_test.py plot --csv data.csv --col-a pre --col-b post    # 配對 → 前後 + 差值
    python run_test.py plot --csv data.csv --x height --y weight       # 相關 → 散布圖

    # Phase 0 常態/對稱檢查：直方圖 + 常態曲線 + Q-Q 圖
    python run_test.py plot --kind hist --csv data.csv --value score --group grp
    python run_test.py plot --kind hist --csv data.csv --col-a pre --col-b post  # 配對差值
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats


def _groups(df, value, group):
    """把 long 格式拆成各組的數值陣列 list。"""
    return [g[value].to_numpy() for _, g in df.groupby(group)]


def mannwhitney(df, args):
    a, b = _groups(df, args.value, args.group)
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    n1, n2 = len(a), len(b)
    rank_biserial = 1 - (2 * res.statistic) / (n1 * n2)
    return {
        "test": "Mann-Whitney U",
        "U": float(res.statistic),
        "p_value": float(res.pvalue),
        "rank_biserial_r": float(rank_biserial),
        "n": [n1, n2],
        "note": "要解讀成中位數差，需兩組分布形狀相近；否則僅為隨機優勢。",
    }


def wilcoxon(df, args):
    a = df[args.col_a].to_numpy()
    b = df[args.col_b].to_numpy()
    res = stats.wilcoxon(a, b)
    diff = a - b
    return {
        "test": "Wilcoxon signed-rank",
        "W": float(res.statistic),
        "p_value": float(res.pvalue),
        "median_diff": float(np.median(diff)),
        "n": int(len(diff)),
        "note": "假設配對差值分布對稱。",
    }


def kruskal(df, args):
    groups = _groups(df, args.value, args.group)
    res = stats.kruskal(*groups)
    n = sum(len(g) for g in groups)
    k = len(groups)
    epsilon_sq = (res.statistic - k + 1) / (n - k)
    return {
        "test": "Kruskal-Wallis",
        "H": float(res.statistic),
        "p_value": float(res.pvalue),
        "epsilon_squared": float(epsilon_sq),
        "k_groups": k,
        "note": "顯著後做 Dunn's test 事後比較並校正多重檢定。",
    }


def plot(df, args):
    """依參數自動選圖型做 EDA，存成 PNG。檢定前務必先畫圖看分布/假設。"""
    try:
        import matplotlib
    except ImportError as e:
        raise SystemExit("未安裝 matplotlib，請執行：pip install matplotlib") from e
    matplotlib.use("Agg")  # 不開視窗，直接存檔
    # 盡量讓中文標題可顯示（Windows 11 內建 Microsoft JhengHei）；找不到則自動退回
    matplotlib.rcParams["font.sans-serif"] = [
        "Microsoft JhengHei", "Microsoft YaHei", "PingFang TC", "SimHei", "DejaVu Sans",
    ]
    matplotlib.rcParams["axes.unicode_minus"] = False
    import matplotlib.pyplot as plt

    kind = args.kind
    if kind == "auto":
        if args.col_a and args.col_b:
            kind = "paired"
        elif args.x and args.y:
            kind = "scatter"
        elif args.value and args.group:
            kind = "box"
        else:
            raise SystemExit(
                "plot 需提供下列之一：--value+--group（組間）、"
                "--col-a+--col-b（配對）、--x+--y（散布）。"
            )

    out = args.out or "nonparametric_eda.png"

    if kind == "box":
        groups = list(df.groupby(args.group))
        labels = [str(k) for k, _ in groups]
        data = [g[args.value].to_numpy() for _, g in groups]
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.boxplot(data, showmeans=True)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        rng = np.random.default_rng(0)
        for i, arr in enumerate(data, start=1):  # 疊上抖動散點，看形狀與離群
            jitter = (rng.random(len(arr)) - 0.5) * 0.15
            ax.scatter(np.full(len(arr), i) + jitter, arr, alpha=0.6, s=20, color="tab:blue")
        ax.set_xlabel(args.group)
        ax.set_ylabel(args.value)
        ax.set_title("組間分布（箱型圖 + 散點）— 檢查各組形狀是否相近")
        note = "形狀相近才可把 Mann-Whitney / Kruskal-Wallis 解讀成中位數差；否則只談隨機優勢。"

    elif kind == "paired":
        a = df[args.col_a].to_numpy()
        b = df[args.col_b].to_numpy()
        diff = b - a
        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        for ai, bi in zip(a, b):  # 左：每個受試者的前後連線
            axes[0].plot([0, 1], [ai, bi], color="gray", alpha=0.6, marker="o")
        axes[0].set_xticks([0, 1])
        axes[0].set_xticklabels([args.col_a, args.col_b])
        axes[0].set_title("配對變化（before → after）")
        axes[1].hist(diff, bins="auto", color="tab:blue", alpha=0.8, edgecolor="white")
        axes[1].axvline(0, color="red", linestyle="--", label="0")
        axes[1].axvline(np.median(diff), color="green", label=f"中位數={np.median(diff):.3g}")
        axes[1].set_title(f"差值分布（{args.col_b} - {args.col_a}）- 檢查對稱性")
        axes[1].legend()
        note = "右圖差值若明顯不對稱，Wilcoxon signed-rank 的對稱假設可疑，需謹慎解讀。"

    elif kind == "scatter":
        x = df[args.x].to_numpy()
        y = df[args.y].to_numpy()
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(x, y, alpha=0.7, color="tab:blue")
        ax.set_xlabel(args.x)
        ax.set_ylabel(args.y)
        ax.set_title("散布圖 - 檢視是否單調關係（Spearman / Kendall）")
        note = "散布圖主要配合相關分析；看關係是否單調，再決定用 Spearman 或 Kendall。"

    elif kind == "hist":
        # Phase 0 常態/對稱檢查：每個樣本/組 畫 直方圖+常態曲線 與 Q-Q 圖
        if args.col_a and args.col_b:
            diff = df[args.col_b].to_numpy() - df[args.col_a].to_numpy()
            series = [(f"{args.col_b}-{args.col_a} 差值", diff)]
        elif args.value and args.group:
            series = [(str(k), g[args.value].to_numpy()) for k, g in df.groupby(args.group)]
        elif args.value:
            series = [(args.value, df[args.value].to_numpy())]
        else:
            raise SystemExit(
                "hist 需提供：--value（可加 --group 分組）或 --col-a+--col-b（配對差值）。"
            )

        nrow = len(series)
        fig, axes = plt.subplots(nrow, 2, figsize=(10, 3.3 * nrow), squeeze=False)
        for r, (label, arr) in enumerate(series):
            arr = arr[~np.isnan(arr)]
            mu, sigma = float(np.mean(arr)), float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
            # 左：直方圖 + 常態密度曲線
            ax_h = axes[r][0]
            ax_h.hist(arr, bins="auto", density=True, color="tab:blue", alpha=0.7, edgecolor="white")
            if sigma > 0:
                xs = np.linspace(arr.min(), arr.max(), 200)
                ax_h.plot(xs, stats.norm.pdf(xs, mu, sigma), "r-", lw=2, label="常態曲線")
            ax_h.axvline(mu, color="green", linestyle="--", label=f"平均={mu:.3g}")
            ax_h.axvline(float(np.median(arr)), color="orange", linestyle=":", label=f"中位數={np.median(arr):.3g}")
            ax_h.set_title(f"{label}：直方圖 + 常態曲線（看鐘型/對稱）")
            ax_h.legend(fontsize=8)
            # 右：Q-Q 圖（對常態）
            ax_q = axes[r][1]
            stats.probplot(arr, dist="norm", plot=ax_q)
            ax_q.set_title(f"{label}：Q-Q 圖（點貼直線=近常態）")
        note = (
            "直方圖近鐘型且 Q-Q 點貼直線 → 近常態，n 夠時考慮母數法；偏態/厚尾/離群 → 無母數。"
            "平均與中位數差很多代表偏態；配對差值不對稱時 Wilcoxon signed-rank 假設可疑。"
        )

    else:
        raise SystemExit(f"未知圖型：{kind}")

    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return {"test": "EDA 視覺化", "kind": kind, "output": out, "note": note}


DISPATCH = {
    "mannwhitney": mannwhitney,
    "wilcoxon": wilcoxon,
    "kruskal": kruskal,
    "plot": plot,
}


def main():
    p = argparse.ArgumentParser(description="無母數檢定計算")
    p.add_argument("test", choices=DISPATCH.keys())
    p.add_argument("--csv", required=True)
    p.add_argument("--value", help="數值欄名 (long 格式)")
    p.add_argument("--group", help="組別欄名 (long 格式)")
    p.add_argument("--col-a", help="配對第一欄 (wilcoxon / plot paired)")
    p.add_argument("--col-b", help="配對第二欄 (wilcoxon / plot paired)")
    p.add_argument("--x", help="散布圖 X 欄 (plot scatter)")
    p.add_argument("--y", help="散布圖 Y 欄 (plot scatter)")
    p.add_argument("--out", help="圖檔輸出路徑 (plot)，預設 nonparametric_eda.png")
    p.add_argument(
        "--kind",
        choices=["auto", "box", "paired", "scatter", "hist"],
        default="auto",
        help="plot 圖型，預設 auto 依參數自動判斷；hist=直方圖+常態曲線+Q-Q（Phase 0 常態/對稱檢查）",
    )
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    result = DISPATCH[args.test](df, args)

    print(f"\n=== {result['test']} ===")
    for key, val in result.items():
        if key == "test":
            continue
        print(f"{key:>18}: {val}")
    print()


if __name__ == "__main__":
    sys.exit(main())
