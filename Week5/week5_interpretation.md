# Week 5 Interpretation

Student: Zixuan Zeng (曾梓轩) · Course: Bioinformatics: From Multi-Omics Data to Discovery
Comparison: treated versus control (reference level = control), design `~ batch + condition`
Word count: 148

This analysis compares treated versus control in a 1,000-gene count matrix using DESeq2,
design `~ batch + condition`, control as reference. On the PCA of VST counts, PC1 (24%
variance) separated the conditions completely, while PC2 (9%) showed no clean batch grouping
and no outlying sample. After filtering (>=10 counts in >=3 samples; 989 of 1,000 retained)
and apeglm shrinkage, 36 genes were higher and 24 lower in treated samples at padj < 0.05
and |log2FC| >= 1, with no undefined adjusted p values. This suggests the treatment mainly
activates gene expression, and that its effect is larger than the batch effect. The main
limitation is that these are simulated counts with placeholder identifiers and one biotype,
so no mechanistic claim is possible. AI wrote the script to thresholds I set beforehand but
could not install DESeq2, so I ran it and verified direction against raw group means.
