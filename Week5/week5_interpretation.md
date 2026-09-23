# Week 5 Interpretation

Student: Zixuan Zeng (曾梓轩) · Course: Bioinformatics: From Multi-Omics Data to Discovery
Comparison: treated versus control (reference level = control), design `~ batch + condition`
Word count: 150

This analysis compares treated versus control in a 1,000-gene count matrix using DESeq2,
design `~ batch + condition`, control as reference. On the PCA of VST counts, PC1 (24%
variance) separated the conditions completely, while PC2 (9%) showed no clean batch grouping
and no outlying sample. After filtering (>=10 counts in >=3 samples; 989 of 1,000 retained)
and apeglm shrinkage, 36 genes were higher and 24 lower in treated samples at padj < 0.05
and |log2FC| >= 1, with no undefined adjusted p values. A further 23 genes passed the
adjusted-p but not the effect-size threshold, showing the two criteria exclude different
weak evidence. The main limitation is that these are simulated counts with placeholder
identifiers and one biotype, so no mechanistic claim is possible. AI wrote the script to
thresholds I set beforehand but could not install DESeq2, so I ran it and verified direction
against raw group means.

<!-- 学生原始中文草稿（保留备查，提交前可删）：
1. QC：PC1 占 24%，把 12 个样本按 condition 完全分开无重叠；PC2 占 9%，批次 C 有三个样本偏上但 TC1 不在，没形成干净的批次分层。没有离群样本。
2. 双阈值：火山图上水平线以上、夹在两条竖线之间的那些灰点就是这 23 个——这正好说明为什么两个门槛都要有。
3. AI 分工原稿：「AI 负责帮我完成 R 语言代码和跑通这个代码」。
   更正：AI 写了代码，但它无法在自己的环境中安装 DESeq2，代码从未由 AI 运行过；
   实际运行与验证由学生在本机 R 4.6.1 完成（见 run_log.txt / session_info.txt）。
-->
