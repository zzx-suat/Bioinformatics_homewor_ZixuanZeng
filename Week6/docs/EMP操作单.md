# Week 6 — EMP-Web 操作单（16S）

作业要求原文：*"Use 16S data files in EMP tests folder and EMP-web to run analysis.
You need to complete whole analysis procedures and submit the final result via the EMP-web system."*

**提交通道是 EMP Sync，必须在平台上跑。** 本目录的 R 分析是并行的方法学交叉核对，不替代平台操作。

---

## 上传的两个文件

```
EasyMultiProfiler-Web/tests/16S_level-7.csv     丰度表（470 taxa × 132 列）
EasyMultiProfiler-Web/tests/16S_mapping.csv     元数据（130 行；SampleID, Group, Group_sub）
```

> ⚠ 丰度表比元数据多两列（`K_XYL_F_0009_03`、`K_XYL_F_0035_03`，无元数据）。
> 平台如何处理这两列不确定——**导入后第一件事就是核对样本数**，看它保留了 130 还是 132。
> 这个数字要记下来写进报告，因为它决定后续所有统计的分母。

---

## 分析步骤

| # | 步骤 | 要点 |
|---|---|---|
| 1 | 导入两个 CSV | 核对样本数（130 还是 132）与特征数（470） |
| 2 | 过滤 | 记下平台用的规则与保留特征数。R 侧用的是"流行度 ≥10%"，保留 136 |
| 3 | 抽平 / 标准化 | 深度相差 12.5 倍（1624~20283）。若平台提供 rarefaction，记录抽平深度 |
| 4 | Alpha 多样性 | 至少出 Shannon 与 Observed；分组选 `Group`（四组）或 `Group_sub`（八组） |
| 5 | Beta 多样性 / 降维 | 若可选距离，优先 Aitchison 或 CLR 后的欧氏；否则 Bray-Curtis 并在报告中说明局限 |
| 6 | 差异分析 | 分组选择见下方"方向"一节 |
| 7 | 富集 / 功能预测 | 若平台有 PICRUSt 之类的功能预测，可做，但结论要标为"预测非实测" |
| 8 | **Export → Sync** | 见下方"同步"一节 |

---

## 三个必须当心的地方

### ① 分组变量选哪个

元数据有两列可选：

```
Group      : IBS_before / IBS_after / UC_before / UC_after        （4 组）
Group_sub  : ..._poor / ..._great                                  （8 组）
```

**平台大概率只能做两两独立组比较，做不了配对。** 这不是你的操作问题，是界面能力所限。
因此：在平台上做 `IBS_before` vs `IBS_after` 是可以的，但报告里必须写明
**这个比较把同一受试者的两个样本当成了独立样本**，而 R 侧的配对分析给出了正确版本
（受限置换 p=0.018/0.022，无限制置换 p≈0.28）。

这个对比本身就是报告里很有价值的一段——不是缺陷，是你发现了平台的方法学边界。

### ② 对比方向

Week 5 在这里翻过车：REFERENCE 与 TEST 设反了，基因列表不变但所有 log2FC 符号颠倒。

**这次先定死**：REFERENCE = `*_before`，TEST = `*_after`。
即正值表示"after 更高"。导出的 CSV 里通常不含方向元数据，所以**在界面上截图存证**。

### ③ ORGANISM 下拉会复位

Week 5 记录过：安装完物种注释包后页面重渲染，ORGANISM 自动跳回 Human。
16S 数据是人体肠道菌群，若平台要求选物种，**每次点运行前重新确认一次**。

---

## 同步（最容易出问题的一步）

Week 5 同步了 **五次**，`results/` 下三个文件**每次都是 4 字节空文件**：

```
RNAseq_output_diff_analysis.csv  ->  ""
RNAseq_output_dimension.csv      ->  ""
RNAseq_output_enrichment.csv     ->  ""
```

推测原因：EMP 预处理是"叠加"模式，log/过滤会产生新快照；同步导出的是**当前 session 绑定的那个 experiment**，若分析结果挂在另一个快照上，导出就是空的。

**这次的做法**：

1. 点 Sync **之前**，先在界面上确认当前选中的 experiment **能看到**差异分析和降维的结果
2. 看得到再同步；看不到就切换快照
3. 同步后**立刻核对**，别只看成功提示：

```bash
git -C "D:/BioLession in total/Bioinformatics/week1" fetch origin main
git -C "D:/BioLession in total/Bioinformatics/week1" ls-tree -r -l origin/main -- "EMP2026/Week_06/"
```

看每个文件的字节数。**4 字节就是空的**，成功提示不作数。

> 参见记忆条目：EMP 的同步和 git push 是两条独立通道，平台弹成功提示不等于文件到了 GitHub。

---

## 平台跑完之后要记录的数字

这些要填进报告，用来和 R 侧结果对照：

- [ ] 导入后的样本数与特征数
- [ ] 过滤规则与保留特征数
- [ ] 抽平深度（若做了）
- [ ] Alpha 多样性：各组的中位数
- [ ] Beta：用的什么距离、前两轴解释比例
- [ ] 差异分析：REFERENCE/TEST 方向、显著特征数、用的什么校正方法
- [ ] Sync 的 commit 链接
