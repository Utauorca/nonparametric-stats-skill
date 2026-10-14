# nonparametric-stats-skill

A Claude Agent Skill that acts as a **nonparametric / distribution-free statistics
consultant**. It helps students and researchers pick the right rank-based test,
explains the assumptions and limits, runs the computation with real code, and
reports effect sizes — not just a p-value.

The statistical content is aligned with this project's course slides
(`slides.qmd`), based on **Hollander, Wolfe & Chicken (2015),
*Nonparametric Statistical Methods* (3rd ed.)**.

## Scope — 6 methods

| Problem | Method |
|---|---|
| One-sample / paired location | Wilcoxon signed-rank |
| Two independent samples (location) | Mann–Whitney U / Wilcoxon rank-sum |
| Two independent samples (dispersion) | Ansari–Bradley |
| General difference between two distributions | Kolmogorov–Smirnov (two-sample) |
| k independent groups | Kruskal–Wallis (post-hoc: Dunn) |
| k dependent groups (repeated measures) | Friedman (post-hoc: Nemenyi) |

Out of scope: nominal-data tests, correlation, goodness-of-fit, randomness tests.

## File structure

```text
nonparametric-stats-skill/
├── SKILL.md                    # main file: agent behavior & answering workflow
├── references/
│   ├── methods.md              # quick reference: H0, statistic, ties, effect size, scipy fn
│   ├── textbook-methods.md     # textbook detail: full assumptions, formulas, R code
│   └── decision-tree.md        # method-selection tree (Phase 0 + main tree)
├── scripts/
│   └── run_test.py             # actual computation (scipy) + plotting (matplotlib)
├── decision-tree-phase0.png    # rendered Phase 0 flowchart (used by the slides)
├── decision-tree-methods.png   # rendered method-selection tree
├── skill-overview.qmd          # Quarto Beamer slides introducing the skill
└── skill-overview.pdf          # rendered slides
```

The skill uses **progressive disclosure**: `SKILL.md` is lightweight and always
loaded; detailed knowledge lives under `references/` and is read on demand;
numeric work is delegated to `scripts/run_test.py`.

## How the agent behaves

- After data is imported it **asks before plotting** (even a scatter plot) and
  **never runs a test without explicit consent** — one step at a time.
- Workflow: clarify key points → describe data & offer a scatter plot →
  recommend a method (ask before running) → assumptions & details →
  effect sizes → parametric vs. nonparametric trade-off.

## scripts/run_test.py

Requires `numpy`, `pandas`, `scipy`, and `matplotlib`.

```bash
# Tests (returns statistic, p-value, effect size)
python scripts/run_test.py mannwhitney --csv data.csv --value score --group grp
python scripts/run_test.py wilcoxon    --csv data.csv --col-a pre --col-b post
python scripts/run_test.py kruskal     --csv data.csv --value score --group grp

# EDA plots (saved as PNG)
python scripts/run_test.py plot --kind scatter --csv data.csv --x height --y weight
python scripts/run_test.py plot --kind box     --csv data.csv --value score --group grp
python scripts/run_test.py plot --kind paired  --csv data.csv --col-a pre --col-b post
python scripts/run_test.py plot --kind hist    --csv data.csv --value score --group grp
```

## Slides

`skill-overview.qmd` is a Quarto Beamer deck introducing the skill. Render with:

```bash
quarto render skill-overview.qmd --to beamer
```

(Requires Quarto + a XeLaTeX engine with a CJK font.)

## Installation as a Claude skill

Copy this folder into your skills directory, e.g.:

```text
.claude/skills/nonparametric-stats/
```

## Reference

Hollander, M., Wolfe, D. A., & Chicken, E. (2015).
*Nonparametric Statistical Methods* (3rd ed.). Wiley.
