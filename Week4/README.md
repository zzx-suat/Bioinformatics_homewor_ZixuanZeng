# Week 4 — Genomics: Sequencing, Variant Interpretation, and AI-Assisted Analysis

**Zixuan Zeng (曾梓轩) · SUAT24000114**
Bioinformatics: From Multi-Omics Data to Discovery · Dr. Liwei Xie · SUAT · Fall 2026

---

## 提交结构：两份文档，两个声音

作业要求每题包含四部分：用 AI 之前的推理 → AI 辅助流程 → 核验 → 最终结论。
本次提交把 **「我的判断」** 和 **「AI 的贡献」** 拆成两份独立文档，避免两者混在一起
看不出谁是谁。

| 文件 | 内容 | 对应作业要求 |
|---|---|---|
| **`Week4_手写部分_曾梓轩.docx`** | 四道题的 *Reasoning before AI*，本人手写填写 | 每题第 1 部分 |
| **`Week4_AI分析_曾梓轩.docx`** | AI 辅助流程、核验、最终结论（Word 版） | 每题第 2–4 部分 |
| `docs/Week4_AI_Analysis.md` | 同上（Markdown 版，内容一致） | 同上 |

手写文档在撰写时不包含任何 AI 结论，四道题的先验判断均在阅读 AI 分析之前完成。
AI 文档中每题开头均标注了对应手写部分的题号（`see worksheet, Q1 ①②③`）。

---

## 四张必需插图

| 图 | 对应题目 |
|---|---|
| `figures/Q1_workflow.png` | Q1 实验策略流程图（按证据层级分级：观察 → 能力 → 必需性） |
| `figures/Q2_fastqc_audit.png` | Q2 对 demo FASTQ 的独立审计（实测，非引用） |
| `figures/Q3_locus_chain.png` | Q3 整合位点图（可及性 → 染色质状态 → 甲基化 → 3D 接触 → 表达 → 干预） |
| `figures/Q4_prioritization.png` | Q4 过滤级联 + 为什么 ClinVar 不能作为首个过滤条件 |

---

## 可复现

```bash
cd Week4
"D:\R-4.6.1\bin\Rscript.exe" R/q4_variant_prioritization.R   # Q4 分析 + 图
py R/q2_fastq_audit.py                                        # Q2 FASTQ 审计 + 图
py R/q1_q3_figures.py                                         # Q1 / Q3 图
py R/q4_verify_databases.py                                   # 实时查询 Ensembl + ClinVar
py R/build_worksheet_docx.py                                  # 生成手写工作簿（空白版）
py R/build_docx.py                                            # 由 md 生成 AI 分析 Word 版
```

环境：R 4.6.1（readr, dplyr, ggplot2, ggrepel, patchwork）、Python 3.12
（matplotlib 3.11.1, python-docx）。会话信息见 `results/q4_sessionInfo.txt`。

`results/` 下所有文件均由上述脚本生成，可逐字节复现（随机数已固定种子）。

---

## 核验说明

Q4 的核验通过实时查询两个权威资源完成，结果存于 `results/q4_verification.tsv`：

- **Ensembl REST**（GRCh38.p14 与 GRCh37 两个端点）—— 基因重叠与参考碱基比对
- **NCBI ClinVar E-utilities** —— 所引 VCV 编号的真实记录

核验发现记录在 AI 分析文档 Q4 的 Verification 一节。数据文件 `data/variants_q4.tsv`
头部已声明为合成教学数据（*"Positions are illustrative"*），相关不一致属教学数据性质，
非需要上报的错误；核验的价值在于展示流程本身。
