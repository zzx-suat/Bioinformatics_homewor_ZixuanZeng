# Week 5 — AI Verification Log

**Course:** Bioinformatics: From Multi-Omics Data to Discovery — Dr. Liwei Xie, SUAT
**Student:** Zixuan Zeng (曾梓轩) · SUAT24000114
**Date:** 2026-09-23
**Run environment:** R 4.6.1 (2026-06-24 ucrt), Windows 11 x64, Bioconductor 3.23 — 见 `session_info.txt`

（格式沿用 Week 4 的 `Week4_AI_Analysis.md`：AI-assisted workflow → Verification → Conclusion）

## What this file is, and what it is not

本文件只记录 Week 5 作业中**使用 AI 的那一部分**：提示词、AI 给出的东西、被拒绝或更正之处、
我独立完成的核查，以及一条尚未解决的不确定性。分析方法与结论写在 `week5_interpretation.md`，
代码在 `week5_deseq2_analysis.R`。

**顺序声明**：J1/J2/J3 三项设计判断是我在看到任何 AI 生成的代码**之前**写下的，
原文逐字保存在 `学生独立判断记录_J1J2J3.md`，未被后续修改覆盖。

---

## Task assigned to AI — 分析脚本与火山图代码

### My independent reasoning, recorded before any AI output

在请求代码之前我先自己确定了三件事（原文见 `学生独立判断记录_J1J2J3.md`）：

- **J1** 设计公式必须含 batch，用于扣除批次带来的系统性差异；本次每批 control/treated 各 2 个，
  设计平衡，batch 与 condition 不混杂，两个效应可分别估计
- **J2** 过滤阈值取「>=10 counts 于 >=3 个样本」，3 介于最小实验单元 2 与每组样本数 6 之间
- **J3** 火山图：x 轴用 apeglm **收缩后**的 log2FC，y 轴用 **padj**，
  水平线 `-log10(0.05) = 1.3010`，垂直线 `log2FC = ±1`，显著须两条件同时满足

### AI-assisted workflow

AI 编写了 `week5_deseq2_analysis.R` 全脚本（含第 8b 节火山图代码）与一键运行器 `一键运行.bat`。
**AI 未能在其自身环境中安装 DESeq2，因此脚本从未由 AI 运行过**；全部运行与验证由我在本机完成。

### Verification — 我做的五项独立核查

| # | 核查项 | 怎么核的 | 实测结果 | 我的判断 |
|---|---|---|---|---|
| 1 | 火山图是否按我事先指定的三个量绘制 | 对照 `week5_de_plot.png` 与代码第 8b 节 | x 轴 shrunken log2FC，y 轴 -log10(padj)，横线标注 y=1.301，竖线 ±1，显著取交集 | 是 |
| 2 | 系数名是否先查再用，而非凭猜 | 看 `run_log.txt` 与代码第 6 节的 `stopifnot` | `resultsNames(dds)` 打印为 Intercept / batch_B_vs_A / batch_C_vs_A / **condition_treated_vs_control**；断言通过 | 先查再用 |
| 3 | 系数方向是否正确（会不会反） | 用 `week5_direction_check.csv` 比对 top20 基因的收缩 LFC 符号与原始 counts 组均值 | **20 / 20 全部一致**。例：Gene0035 LFC=1.737，control 均值 47.3 → treated 均值 166.8 | 正确，没有反 |
| 4 | 收缩是否改变了 padj | 代码第 7 节 `all.equal(res$padj, resLFC$padj)` 断言 | 运行未报错，padj 未被收缩改变 | 没有改变 |
| 5 | padj 为 NA 的基因如何处理 | 统计结果表中 padj 为 NA 的行数 | **0 个**，本次不存在「把 NA 误当不显著」的风险 | 本次没有NA的基因 |

补充记录：另有 **23 个基因 padj < 0.05 但 |log2FC| < 1**，被效应量门槛排除；
反向（|log2FC| >= 1 但 padj >= 0.05）为 **0 个**。

### Errors and corrections — AI 出错或被我更正的地方

1. **我自己的初始判断被更正（J2）**：我最初把过滤规则理解成「删掉 3 个样本」。
   核对代码后确认它过滤的是**基因**而非样本（`rowSums(counts >= 10) >= 3`），12 个样本全程保留。
   原文与更正都保留在 `学生独立判断记录_J1J2J3.md`，未覆盖原文。
2. **AI 对评分标准的判断出错并被推翻**：AI 最初依据阅读材料里的 100 分评分表，
   声称按作业说明书交付会丢 30 分。核对 PPT 第 49 页后（老师红框圈出的作业文件只有两个 `.md` 说明书）
   该判断被推翻并更正。
3. **AI 无法在其运行环境中安装 DESeq2**（其环境为 Windows，conda 无 win-64 构建，
   Bioconductor 官方源在其网络中不可达），因此脚本中依赖 DESeq2 的部分**未经 AI 自身运行验证**。
   实际运行验证由我在本机完成，见 `run_log.txt` 与 `session_info.txt`。
4. **AI 另用 pydeseq2（DESeq2 的独立 Python 重实现）做了交叉核对**，其结果与我本机 R DESeq2
   **完全一致**（同为 36 上调 / 24 下调 / 合计 60，padj 为 NA 均为 0，方向 20/20）。
   交叉核对结果存放于 `crosscheck_pydeseq2_NOT_submission/`，**不作为提交物**。

### Conclusion the AI supports

在上述五项核查通过的前提下，AI 生成的代码支持这一结论：在 `~ batch + condition` 设计下，
treated 相对 control 有 36 个基因上调、24 个下调（padj < 0.05 且 |log2FC| >= 1）。
它**不能**支持任何机制性结论——数据为合成数据，基因标识为占位符，且相关不等于因果。

---

## One unresolved uncertainty

- 合成数据的离散度结构是否与真实 RNA-seq 一致，这直接影响这套流程迁移到真实数据时的表现

---

## Appendix

### A1. Prompts used

提示词当时以口头/对话形式给出，以下为事后依对话记录的重述（非逐字原文）：

> 写 Week 5 的分析脚本和火山图代码：x 轴用 apeglm 收缩后的 log2FC，y 轴用 -log10(padj)，
阈值线画在 ±1 和 -log10(0.05)，显著要求两个条件同时满足

### A2. Thresholds and their justification

- `padj < 0.05`：BH 校正后的 FDR 阈值。1000 个基因各做一次检验，即使全无差异也约有 50 个
  原始 p < 0.05，必须校正
- `|log2FC| >= 1`：效应量阈值，即至少 2 倍变化，排除统计显著但幅度无实际意义的基因
- 两者取**交集**，分别控制「可信度」与「幅度」。本次数据中有 23 个基因只满足前者

### A3. Reproducing this analysis

双击 `一键运行.bat`（自动定位 Rscript、按需安装 DESeq2/apeglm、运行分析、输出到 `run_log.txt`），
或在 R 中 `source("week5_deseq2_analysis.R")`。输入为 `../Homework/for_student/` 下的
count matrix 与 metadata，**原始输入文件未被修改**。

### A4. Declaration of AI use

AI 参与了：分析脚本与火山图代码的编写、一键运行器的编写、数据完整性的初步核查、
以及 pydeseq2 交叉核对。
AI 未参与：J1/J2/J3 三项设计判断（在 AI 输出之前完成）、上表五项核查的执行、
脚本的实际运行、以及「未解决的不确定性」的选择。
