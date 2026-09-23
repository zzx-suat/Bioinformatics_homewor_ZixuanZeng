# Week 5 — Transcriptomics: RNA-seq and Differential Expression with DESeq2

**Zixuan Zeng (曾梓轩) · SUAT24000114**
Bioinformatics: From Multi-Omics Data to Discovery · Dr. Liwei Xie · SUAT · Fall 2026

本目录是 **Homework 1**（本地 R + DESeq2 分析）的提交物。
**Homework 2**（EasyMultiProfiler-Web 平台分析）按作业说明经平台 Sync 提交，
产物在本仓库 `EMP2026/Week_05/transcriptomics/weekly/` 下。

---

## Homework 1 — 八个必交项

| 文件 | 内容 |
|---|---|
| `week5_deseq2_analysis.R` | 完整分析脚本，含输入完整性断言与模型矩阵满秩检查 |
| `week5_deseq2_results.csv` | apeglm 收缩后的完整结果表（含不显著基因） |
| `week5_pca.png` | VST 转换后的 PCA |
| `week5_de_plot.png` | 差异表达图（标注了坐标轴与阈值） |
| `week5_interpretation.md` | 150 词解读 |
| `week5_AI_verification_log.md` | AI 使用与独立核验记录 |
| `week5_deseq2_object.rds` | 拟合后的 DESeq2 对象 |
| `session_info.txt` | 运行环境 |

## 附带材料（非必交，作为可复现性与独立判断的证据）

| 文件 | 内容 |
|---|---|
| `学生独立判断记录_J1J2J3.md` | **在看到任何 AI 代码之前**写下的三项设计判断原文；后续更正另起段落，未覆盖原文 |
| `run_log.txt` | 本机实际运行日志 |
| `run_week5.R` | 运行入口 |
| `design_and_checks.tsv` | 设计与检查项记录 |
| `week5_direction_check.csv` | 方向核验：top20 收缩 LFC 符号 vs 原始 counts 组均值 |
| `week5_top20_features.csv` | Top20 特征 |

---

## 分析要点

- 输入 1,000 基因 × 12 样本，全整数；设计完全平衡（3 批次 × 2 条件 × 2 重复）
- 设计式 `~ batch + condition`，`control` 为参考水平；模型矩阵满秩 4/4
- `resultsNames(dds)` = Intercept / batch_B_vs_A / batch_C_vs_A / **condition_treated_vs_control**
- 过滤规则：`>= 10 counts 于 >= 3 样本` → 保留 **989 / 1000**
- 显著（`padj < 0.05` 且 `|log2FC| >= 1`）：**上调 36、下调 24、合计 60**；`padj` 为 NA：**0**
- 另有 23 个基因通过 padj 阈值但未过效应量阈值 —— 两个门槛筛掉的是不同类型的弱证据
- 方向核验：top20 收缩 LFC 符号与原始 counts 组均值 **20/20 一致**

运行环境：R 4.6.1 (2026-06-24 ucrt)，Windows 11 x64，Bioconductor 3.23。

## 数据性质的说明

Homework 1 的输入是**合成数据**（基因名为 `Gene0001` 形式，biotype 单一），
因此解读中不作机制性断言，这一点已写入 `week5_interpretation.md` 的局限性部分。
