# Week 6 — 16S 微生物组分析

**Zixuan Zeng (曾梓轩) · SUAT24000114**
Bioinformatics: From Multi-Omics Data to Discovery · Dr. Liwei Xie · SUAT · Fall 2026

数据：`EasyMultiProfiler-Web/tests/16S_level-7.csv` + `16S_mapping.csv`

---

## 主文档

**`Week6_分析报告_曾梓轩.docx`** —— 被评文档，所有必交元素均在其中：

| 章节 | 内容 |
|---|---|
| 一 | 数据与设计：三处数据不一致的处置、配对结构的识别 |
| 二 | 方法与选择理由 |
| 三 | 结果：alpha 多样性、beta 多样性与 PERMANOVA、差异丰度；**3.4 EMP-Web 平台结果与 R 对照** |
| **四** | **科学假设（本人撰写）** |
| 五 | 局限 |
| 六 | AI 使用声明与审计表 |
| 附录 A | 全部分析代码原文 |
| 附录 B | 完整运行日志 |
| 附录 C | sessionInfo() |

`docs/Week6_分析报告.md` 为同一报告的 Markdown 源文件。

---

## 核心发现

数据是 **64 位受试者 × before/after 的配对设计**，不是四个独立组。

同一份数据、同一个 PERMANOVA，只改置换方式：

| 疾病 | 无限制置换 p | 受试者内受限置换 p |
|---|---:|---:|
| IBS | 0.277 | **0.018** |
| UC | 0.276 | **0.022** |

受试者身份单独解释群落方差 R² = 0.613。无限制置换会跨受试者打乱标签，
把个体间差异灌进零分布，从而漏掉真实的时间点效应。

---

## 目录

```
Week6/
├── Week6_分析报告_曾梓轩.docx
├── docs/     报告 md 源文件、EMP 操作单
├── R/        00 环境检查 · 01 清洗 · 02 alpha · 03 beta · 04 差异丰度 · 05 平台对照 · run_all · build_docx
├── figures/  02_alpha.png · 03_beta.png · 04_paired_diff.png · 05_platform_vs_R.png
└── results/  清洗日志、各步日志与结果表、完整运行日志、sessionInfo
    └── EMPresult/  EMP-Web 平台的原始产出（文件说明见其中 README.md）
```

## 复现

```bash
cd Week6
"D:\R-4.6.1\bin\Rscript.exe" R/run_all.R    # 01-05，约 5 秒；05 需要 results/EMPresult/
py R/build_docx.py                           # 重新生成 docx
```

随机数已固定种子（`set.seed(1)`），`results/` 可逐字节复现。

---

## EMP-Web 平台提交

作业要求经 EMP-web 提交最终结果。平台操作由本人完成，两次有效同步：

| 同步运行 | 内容 |
|---|---|
| `EMP2026/Week_06/microbiome_16s/weekly/runs/2026-09-24T16-35-32-319Z-res6qh` | IBS before vs after，72 样本 |
| `EMP2026/Week_06/microbiome_16s/weekly/runs/2026-09-24T16-48-50-412Z-oqnxna` | UC before vs after，58 样本（LATEST） |

同步文件中 alpha、降维、富集为空（平台只导出实验对象里的结果，分析页结果另存），
因此平台的全部实际产出另存于 `results/EMPresult/`，与 R 的对照见报告 3.4 节：
平台（不配对）与 R（配对）在名义显著的特征上方向 100% 一致；
平台自己算的 UC Shannon 用配对检验 p=0.020，按不配对处理 p=0.214。

操作步骤见 `docs/EMP操作单.md`。
