# EMP-Web 平台产出（Week 6，16S）

本目录是 EasyMultiProfiler-Web v9.0.4 上实际运行的产出，平台操作由学生本人完成。
**这些是平台自己的输出**；R 侧的独立分析在上一级 `results/` 中，两者的对照见报告 3.4 节。

平台同步到 GitHub 时 alpha、降维、富集三项结果为 4 字节空文件，因此平台的这些产出以本目录为准。

---

## 文件说明

| 文件 | 内容 | 样本覆盖 | 来源 |
|---|---|---|---|
| `16S_level_7_EMP_EMPT.rds` | 导入后的平台实验对象 | 132（含 2 个 Group 为 NA 的样本） | 界面下载 |
| `prep-preview.csv` | 过滤后 assay 的预览 | 132 | 界面下载 |
| `alpha-result.csv` | 7 种 alpha 指数 | 132 | 分析页下载 |
| `dim-result.csv` | PCoA（Bray-Curtis）坐标 | 132 | 分析页下载 |
| `diff-result.csv` | **IBS_before vs IBS_after** 差异分析，128 个特征 | IBS 72 | 分析页下载 |
| `20260925-003258.zip` | 一键运行结果包（勾选「仅保留 2 组」之后） | **仅 IBS 72** | 一键页下载 |
| `20260925-004748_RunAll_全部132样本.zip` | 一键运行结果包（勾选之前）：alpha、PCoA/PCA/NMDS、top15 柱状图、top40 热图 | **全部 132** | 取自平台会话目录 `bundles/` |
| `UC_diff-result.csv` | **UC_before vs UC_after** 差异分析，120 个特征 | UC 58 | 取自 GitHub 上的平台同步文件 |

两个结果包的差异分析表：`003258` 中有 IBS 的 top40 差异表；`004748` 在 4 个组上运行，
差异分析步骤未生成表（两两比较需要恰好两组）。

---

## 解读前必须知道的三点

1. **log2FC 为正 = 在参照组（`*_before`）更高**，与常规习惯相反。已用原始数据核验。
   判断哪组更高请看 `sign_group` 列。
2. **「仅保留 2 组」会直接从实验对象中删除其余样本**，之后的一键运行和同步都只含那两组。
3. **平台保留了 `K_XYL_F_0009_03`、`K_XYL_F_0035_03`**，二者 Group 为 NA，
   参与 alpha 与排序计算，不进入两组比较。R 侧已删除。

## 对应的 GitHub 同步

| 运行 | 内容 |
|---|---|
| `EMP2026/Week_06/microbiome_16s/weekly/runs/2026-09-24T16-35-32-319Z-res6qh` | IBS，72 样本 |
| `EMP2026/Week_06/microbiome_16s/weekly/runs/2026-09-24T16-48-50-412Z-oqnxna` | UC，58 样本 |
