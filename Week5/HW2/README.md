# Week 5 Homework 2 — EasyMultiProfiler-Web RNA-seq 分析

**Zixuan Zeng (曾梓轩) · SUAT24000114**

数据：`EasyMultiProfiler-Web/tests/RNAseq_output.csv` + `RNAseq_mapping.csv`
24,394 基因 × 24 样本，6 组 × 4 重复：

```
DMSO · DMSO+LIPUS · T4400 · T4400+LIPUS · T3976 · T3976+LIPUS
```

---

## ⚠️ 本目录文件的来源必须分清

平台 Sync 时 `results/` 下三个文件被导成了**空文件**（4 字节，内容为 `""`）：

```
EMP2026/Week_05/.../results/RNAseq_output_diff_analysis.csv  ->  ""
EMP2026/Week_05/.../results/RNAseq_output_dimension.csv      ->  ""
EMP2026/Week_05/.../results/RNAseq_output_enrichment.csv     ->  ""
```

因此本目录分两部分存放，**来源不同，不可混为一谈**：

### `platform_outputs/` — 平台自身的输出

| 文件 | 内容 | 来源 |
|---|---|---|
| `tx-analysis-result.csv` | GSEA 结果，109 条 KEGG 通路（Mouse） | EMP-Web 界面导出 |
| `barplot.png` | Top20 丰度概况 | EMP-Web 界面导出 |
| `prep-preview.csv` | 预处理预览（29 基因 × 24 样本原始 counts） | EMP-Web 界面导出 |

过滤后的表达矩阵（14,245 × 24，log2）已由平台 Sync 至
`EMP2026/Week_05/transcriptomics/weekly/runs/2026-09-23T05-39-08-324Z-mvq13l/data/`，
此处不重复存放。

### `r_rebuild/` — 用 R 重建的缺失产物

**这些不是平台的输出。** 它们是用**同一份源数据**在本机 R 中独立重跑得到的，
用于补上平台导出失败的 PCA、差异表与火山图。统计实现、默认参数与版本
与平台不同，数值不保证与平台一致。

| 文件 | 内容 |
|---|---|
| `hw2_rebuild_analysis.R` | 重建脚本，含输入断言与满秩检查 |
| `hw2_pca_all24.png` | 全部 24 样本的 PCA（VST 转换） |
| `hw2_volcano_T4400_vs_DMSO.png` | 火山图 |
| `hw2_deseq2_results_T4400_vs_DMSO.csv` | 完整差异结果表（含不显著基因） |
| `hw2_direction_check.csv` | 方向核验：收缩 LFC 符号 vs 原始 counts 组均值 |
| `hw2_pca_coordinates.csv` | PCA 坐标 |
| `hw2_rebuild_summary.csv` | 参数与结果摘要 |
| `hw2_session_info.txt` | 运行环境 |

**特征集的处理**：没有去猜平台的过滤规则（试过 11 种常见规则均不匹配 14,245——
被丢弃的基因里存在 rowSum 143 且 24 样本全非零者，说明不是简单计数阈值）。
改为**直接读取平台导出 assay 的基因名**作为特征集，再回到原始整数 counts 上分析，
从而保证特征集与平台完全一致。

**输入用的是原始整数 counts，不是 log 矩阵**——DESeq2 要求原始计数。

---

## 重建结果

**PCA（全部 24 样本）**：PC1 56.4%，PC2 13.9%。
六组未形成清晰分离，T4400 组在 PC1 上跨度很大（约 −2 到 +53），**组内变异明显**。
这一点在解读差异基因时必须一并考虑。

**差异分析 T4400 vs DMSO**（参考水平 DMSO，8 样本，模型矩阵满秩 2/2）：

| 项 | 值 |
|---|---:|
| 阈值 | `padj < 0.05` 且 `\|log2FC\| >= 1` |
| 上调（T4400 更高） | **116** |
| 下调（T4400 更低） | **12** |
| 合计显著 | **128** |
| `padj` 为 NA（独立过滤所致） | 832 |
| 仅过 padj、未过效应量 | 2,068 |
| 方向核验 | top20 中 **20/20** 与原始组均值一致 |

`resultsNames(dds)` = `Intercept` | `Group_T4400_vs_DMSO`（系数名经显式核对，未假设）。

---

## 报告

`作业二_EMP-Web分析报告.docx` —— 七节完整报告，第 4.3 节两点补充解读由学生本人撰写。

## 尚未完成

- [ ] 平台侧重新 Sync，使 `EMP2026/Week_05/.../results/` 三个文件不再为空。
      已同步三次（05:01 / 05:39 / 06:06），三次均为 4 字节空文件；
      本目录 `r_rebuild/` 已用 R 补齐这些产物。

## 未做的分析（可选扩展）

本次只跑了 T4400 vs DMSO 一个对比。数据本身是 3 药 × 2 超声的析因设计，
若要回答超声的作用需另跑 `DMSO+LIPUS vs DMSO`；若要判断药物与超声是否存在交互，
需在模型中引入交互项（`~ Drug * LIPUS`），平台界面未提供该选项，但 R 中可行。
