#!/usr/bin/env python3
"""
Build the HAND-WRITTEN worksheet for Week 4.

This file deliberately contains NO answers. It explains what each question is
asking (in Chinese, since the course material is English), lists the concepts you
need in order to understand the question, and then gives ruled blank space for
Zixuan to write his own reasoning by hand.

The AI's own analysis lives in a separate file (docs/Week4_AI_Analysis.md) so the
two voices are never mixed.

Out: Week4_手写部分_曾梓轩.docx
"""
import os

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

PROJ = "." if os.path.isdir("docs") else ".."
OUT = os.path.join(PROJ, "Week4_手写部分_曾梓轩.docx")

INK = RGBColor(0x24, 0x1A, 0x38)
ACCENT = RGBColor(0x3E, 0x2A, 0x63)
MUTED = RGBColor(0x6B, 0x64, 0x78)
RULE = "B8B2C4"

CN = "Microsoft YaHei"


def cjk(run, font=CN):
    """python-docx sets only the Latin font; CJK needs eastAsia too."""
    run.font.name = font
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts"); rPr.append(rFonts)
    rFonts.set(qn("w:eastAsia"), font)
    rFonts.set(qn("w:ascii"), font)
    rFonts.set(qn("w:hAnsi"), font)


def para(doc, text="", size=10.5, bold=False, color=INK, space_after=6,
         italic=False, indent=0.0):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(space_after)
    if indent:
        p.paragraph_format.left_indent = Inches(indent)
    r = p.add_run(text)
    r.bold = bold; r.italic = italic
    r.font.size = Pt(size); r.font.color.rgb = color
    cjk(r)
    return p


def rule_line(doc, n=1, gap=20):
    """n ruled blank lines to write on."""
    for _ in range(n):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(gap)
        p.paragraph_format.space_after = Pt(0)
        pPr = p._p.get_or_add_pPr()
        bd = OxmlElement("w:pBdr"); bot = OxmlElement("w:bottom")
        bot.set(qn("w:val"), "single"); bot.set(qn("w:sz"), "4")
        bot.set(qn("w:color"), RULE)
        bd.append(bot); pPr.append(bd)


def shade_para(doc, text, fill="F2EFF7", size=10, bold=False, color=INK):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.left_indent = Inches(0.12)
    p.paragraph_format.right_indent = Inches(0.12)
    pPr = p._p.get_or_add_pPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear"); sh.set(qn("w:fill"), fill)
    pPr.append(sh)
    r = p.add_run(text); r.bold = bold
    r.font.size = Pt(size); r.font.color.rgb = color
    cjk(r)
    return p


def heading(doc, text, size=14, color=ACCENT, before=16):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(before)
    p.paragraph_format.space_after = Pt(6)
    pPr = p._p.get_or_add_pPr()
    bd = OxmlElement("w:pBdr"); bot = OxmlElement("w:bottom")
    bot.set(qn("w:val"), "single"); bot.set(qn("w:sz"), "10")
    bot.set(qn("w:color"), "3E2A63")
    bd.append(bot); pPr.append(bd)
    r = p.add_run(text); r.bold = True
    r.font.size = Pt(size); r.font.color.rgb = color
    cjk(r)


def label(doc, text, color=ACCENT):
    return para(doc, text, size=10.5, bold=True, color=color, space_after=2)


def page_break(doc):
    doc.add_paragraph().add_run().add_break()


def build():
    doc = Document()
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(0.85)
        s.top_margin = s.bottom_margin = Inches(0.75)

    st = doc.styles["Normal"]
    st.font.size = Pt(10.5)
    st.paragraph_format.space_after = Pt(6)

    # ---------------- cover ----------------
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Week 4 作业 · 我的思考记录"); r.bold = True
    r.font.size = Pt(20); r.font.color.rgb = INK; cjk(r)

    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Genomics — Sequencing, Variant Interpretation, and AI-Assisted Analysis")
    r.font.size = Pt(11); r.font.color.rgb = MUTED; cjk(r)

    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Bioinformatics: From Multi-Omics Data to Discovery · Dr. Liwei Xie · SUAT")
    r.font.size = Pt(9.5); r.font.color.rgb = MUTED; cjk(r)

    doc.add_paragraph()
    para(doc, "姓名 ______________________     学号 ______________________     日期 ______________",
         size=11, space_after=18)

    shade_para(doc,
               "这份文件是「用 AI 之前」我自己的判断。四道题各有一段 Reasoning before AI，"
               "按作业要求，这部分必须是我本人在使用 AI 之前写下的，所以全部手写填空。\n\n"
               "AI 的分析、核验过程和它提出的反驳，另见 Week4_AI_Analysis.md —— 两份文件不混在一起，"
               "哪些是我的判断、哪些是 AI 的贡献，一眼可分。",
               fill="F2EFF7", size=10)

    shade_para(doc,
               "怎么用这份文件：先读【这道题在问什么】，确认自己看懂了；再看【关键概念】，"
               "这些是理解题目必需的词，不是答案；然后在横线上写你自己的想法。"
               "写得不完美没关系 —— 这一栏评的是你的判断，不是文采。"
               "写完之后再去看 AI 那份。顺序反了，这一栏就没意义了。",
               fill="FFF6E8", size=10)

    page_break(doc)

    # =====================================================================
    # Q1
    # =====================================================================
    heading(doc, "Question 1 — 选择正确的基因组实验 (25 分)", before=0)

    label(doc, "原题")
    shade_para(doc,
               "Gene X is significantly upregulated in a disease group compared with healthy "
               "controls. Possible explanations include a regulatory variant, altered chromatin "
               "accessibility, altered transcription-factor binding or histone modification, DNA "
               "methylation changes, or altered enhancer–promoter contact.\n\n"
               "Design an experimental strategy to determine which mechanism is most likely "
               "responsible. Choose from: WGS/WES · ATAC-seq · ChIP-seq / CUT&Tag / CUT&RUN · "
               "WGBS / EM-seq · CAGE / RAMPAGE · Hi-C / Micro-C · RNA-seq · MPRA · CRISPR perturbation.\n\n"
               "Before AI: state which assay you would choose first and why.",
               fill="F7F7FA", size=9.5)

    label(doc, "这道题在问什么")
    para(doc,
         "某个基因 X 在病人里表达变高了。为什么变高？可能是 DNA 序列上出了变异，也可能是染色质"
         "打开了、某个转录因子结合上去了、甲基化变了，或者远处的增强子跟启动子接触上了。"
         "题目给你九类实验，问你：你会先做哪一个，为什么先做它，以及它单独做完之后还有什么是证明不了的。")
    para(doc,
         "注意题目的落脚点不是「把九个实验都做一遍」，而是「你凭什么按这个顺序做」。"
         "最后那句必须以 The biological question chooses the assay because… 结尾，"
         "说明它想让你答的是「问题决定实验」，不是「实验列表」。")

    label(doc, "关键概念（理解题目用，不是答案）")
    shade_para(doc,
               "• 每个实验测的是不同的物理量：ATAC-seq 测转座酶能不能接近 DNA；ChIP-seq 测某个蛋白"
               "有没有富集在这段 DNA 上；Hi-C 测两段 DNA 在空间上离得近不近。这三件事不是一回事。\n"
               "• 「关联」和「因果」不同：观察到某处染色质是开的，跟证明这处染色质导致基因表达升高，"
               "中间差了很多步。\n"
               "• MPRA 这类报告基因实验测的是「这段序列有没有能力驱动转录」；CRISPR 干扰测的是"
               "「在真实基因组位置上，这段序列是不是必需的」。能力 ≠ 必需。",
               fill="F7F7FA", size=9.5)

    label(doc, "① 你会先做哪个实验？为什么是它而不是别的？")
    rule_line(doc, 4)
    doc.add_paragraph()

    label(doc, "② 你一开始猜 Gene X 变高最可能是哪种机制？（凭直觉写，不用查）")
    rule_line(doc, 3)
    doc.add_paragraph()

    label(doc, "③ 你选的那个实验，单独做完之后有什么是它证明不了的？")
    rule_line(doc, 3)
    doc.add_paragraph()

    shade_para(doc,
               "想不出来时可以问自己：如果这个实验做出来是阳性结果，我能不能就此说"
               "「找到原因了」？如果不能，缺的是哪一块？",
               fill="FFF6E8", size=9.5)

    page_break(doc)

    # =====================================================================
    # Q2
    # =====================================================================
    heading(doc, "Question 2 — 从 FASTQ 到可信的分析流程 (25 分)", before=0)

    label(doc, "原题")
    shade_para(doc,
               "You receive paired-end Illumina FASTQ files from a human sequencing experiment. "
               "Design a complete analysis workflow from raw reads to interpretable genomic "
               "results. Include: (1) FASTQ quality control, (2) reference genome selection, "
               "(3) alignment, (4) mapped-read processing, (5) assay-specific downstream "
               "analysis, (6) annotation, (7) visualization, (8) interpretation.\n\n"
               "Before AI: draw your own workflow and explain the purpose of each major step.",
               fill="F7F7FA", size=9.5)

    label(doc, "这道题在问什么")
    para(doc,
         "测序仪吐出来的 FASTQ 文件，怎么一步步变成能讲出生物学结论的东西。"
         "题目已经把八个阶段列好了，所以它不是考你「有哪些步骤」——那是抄的。"
         "它考的是：**每一步到底是为了解决什么问题**。")
    para(doc,
         "这题不要求你真的跑完整条流程，是设计题。另外它要你至少解读 4 个 FastQC 指标。"
         "结尾必须是 The analyst, not the AI, is responsible for… —— 它想让你说出"
         "「哪些判断是机器替不了人的」。")

    label(doc, "关键概念（理解题目用，不是答案）")
    shade_para(doc,
               "• FASTQ 里有什么：每条 read 四行 —— 编号、碱基序列、分隔符、每个碱基的质量分。\n"
               "• Phred 质量分 Q = −10 log₁₀(出错概率)。Q20 = 1% 错，Q30 = 0.1% 错，Q40 = 0.01% 错。\n"
               "• FASTQ 里**没有**的东西：这条 read 来自基因组哪个位置、比对得唯不唯一。"
               "这些都是后面比对才算出来的。\n"
               "• 「碱基质量」和「比对质量」是两回事：一个碱基可以测得很准（Q30），"
               "但它所在的 read 却不知道该放在基因组的哪里。\n"
               "• Adapter（接头）：建库时接上去的人工序列，如果插入片段太短，测序会读穿进接头，"
               "这段序列不是基因组来的。",
               fill="F7F7FA", size=9.5)

    label(doc, "① 你自己的八步流程，每步写一句「这一步是为了……」")
    for i, name in enumerate(
            ["FASTQ 质控", "参考基因组选择", "比对", "比对后处理",
             "针对该实验类型的分析", "注释", "可视化", "解读"], 1):
        para(doc, f"{i}. {name} —— 这一步是为了：", size=10, space_after=2, indent=0.12)
        rule_line(doc, 1, gap=16)
        doc.add_paragraph()

    label(doc, "② 这八步里，你自己最没把握的是哪一步？为什么？")
    rule_line(doc, 3)
    doc.add_paragraph()

    shade_para(doc,
               "第②问不是扣分项，是加分项。作业评分表里有 25 分给「独立科学推理」，"
               "诚实说出自己的不确定，比假装全懂更像真正在做研究。",
               fill="FFF6E8", size=9.5)

    page_break(doc)

    # =====================================================================
    # Q3
    # =====================================================================
    heading(doc, "Question 3 — 整合多组学证据，提出调控假说 (25 分)", before=0)

    label(doc, "原题")
    shade_para(doc,
               "A candidate regulatory region is located upstream of Gene Y. You have ATAC-seq, "
               "H3K27ac ChIP-seq/CUT&Tag, DNA methylation, Hi-C/Micro-C, and RNA-seq data from "
               "the same biological condition. Determine whether this region is a plausible "
               "enhancer of Gene Y.\n\n"
               "Before AI: interpret each omics layer independently, then draft one integrated "
               "regulatory model.",
               fill="F7F7FA", size=9.5)

    label(doc, "这道题在问什么")
    para(doc,
         "Gene Y 上游有一段可疑的调控区域。你手上有五种数据，都来自同一个条件。"
         "问题是：这段区域到底算不算 Gene Y 的增强子？")
    para(doc,
         "题目特意要求你**先一层一层单独看**，再整合。这个顺序是故意的——"
         "如果一上来就整合，你会不自觉地让五层数据互相「印证」，"
         "而它们很可能只是在各自说各自的事。后面还要你提出一个**替代解释**，"
         "以及一个能区分相关与因果的实验。")

    label(doc, "关键概念（理解题目用，不是答案）")
    shade_para(doc,
               "• ATAC-seq：染色质开放程度。开放 ≠ 有蛋白结合，也 ≠ 基因一定在转录。\n"
               "• H3K27ac：一种组蛋白修饰，常见于活跃的增强子和启动子附近。是「伴随现象」，"
               "不等于「原因」。\n"
               "• DNA 甲基化：启动子区高甲基化常伴随沉默。注意标准 WGBS/EM-seq 把 5mC 和 5hmC "
               "混在一起测，分不开。\n"
               "• Hi-C/Micro-C：测两段 DNA 空间上接触的频率。是很多细胞的平均，不是某一个细胞的"
               "真实三维结构。\n"
               "• 增强子（enhancer）：能提高远处启动子转录活性的一段 DNA。判定它需要证明"
               "「去掉它，目标基因就下降」。",
               fill="F7F7FA", size=9.5)

    label(doc, "① 五层数据，每层单独能告诉你什么？（一层一句，先别整合）")
    for name in ["ATAC-seq", "H3K27ac", "DNA 甲基化", "Hi-C / Micro-C", "RNA-seq"]:
        para(doc, f"{name} 单独告诉我：", size=10, space_after=2, indent=0.12)
        rule_line(doc, 1, gap=16)
        doc.add_paragraph()

    label(doc, "② 现在把五层合起来，写出你的第一版调控模型（2–3 句）")
    rule_line(doc, 4)
    doc.add_paragraph()

    shade_para(doc,
               "写完②之后，回头看①。如果你发现②里说了某件事，而①里五层都没有单独支持它，"
               "那就是你悄悄多推了一步——把它标出来，这个自我发现本身就很值钱。",
               fill="FFF6E8", size=9.5)

    page_break(doc)

    # =====================================================================
    # Q4
    # =====================================================================
    heading(doc, "Question 4 — AI 辅助的变异优先级排序 (25 分)", before=0)

    label(doc, "原题")
    shade_para(doc,
               "Use the synthetic variant table in data/variants_q4.tsv (columns: CHROM, POS, "
               "REF, ALT, FILTER, DP, GQ, AF, GENE, CONSEQUENCE, CLINVAR_SIG, CLINVAR_ID, NOTE). "
               "Prioritize one or two variants for further investigation.\n\n"
               "Before AI: define your own filtering logic (technical quality, population "
               "frequency, functional consequence, biological/clinical evidence).",
               fill="F7F7FA", size=9.5)

    label(doc, "这道题在问什么")
    para(doc,
         "给你 12 个变异，让你挑出 1–2 个值得深入研究的。表格里每一列都是一个可以用来筛选的条件。")
    para(doc,
         "这题真正考的不是「你的阈值定成多少」，而是**你按什么顺序用这些条件**。"
         "同样一批阈值，先用 A 再用 B，和先用 B 再用 A，选出来的变异可能完全不同。"
         "所以下面第②问比第①问重要得多。")

    label(doc, "关键概念（理解题目用，不是答案）")
    shade_para(doc,
               "• FILTER：变异检测软件自己给的标记。PASS 表示「通过了这个软件的过滤」，"
               "仅此而已 —— 不代表稀有、不代表致病、不代表有功能。\n"
               "• DP（depth）：这个位点有多少条 read 覆盖。DP 越低，判断越不可靠。\n"
               "• GQ（genotype quality）：对「基因型判断」的置信度，也是 Phred 尺度。\n"
               "• AF（allele frequency）：这个等位基因在人群中有多常见。\n"
               "• CONSEQUENCE：预测的功能后果（错义、无义、移码、剪接、同义……）。\n"
               "• CLINVAR_SIG：ClinVar 数据库里别人提交的临床解读（Pathogenic / Benign / "
               "Uncertain significance 等）。这是**别人的解读**，不是实验证据。",
               fill="F7F7FA", size=9.5)

    label(doc, "① 你打算用哪些条件筛？各定什么阈值？（自己填，不用管对错）")

    t = doc.add_table(rows=7, cols=3)
    t.style = "Table Grid"
    hdr = ["筛选条件", "我定的规则 / 阈值", "为什么这样定"]
    widths = [Inches(1.5), Inches(2.2), Inches(2.9)]
    for ci, h in enumerate(hdr):
        c = t.cell(0, ci); c.text = ""
        pr = c.paragraphs[0]; r = pr.add_run(h); r.bold = True
        r.font.size = Pt(9.5); r.font.color.rgb = INK; cjk(r)
        pPr = c._tc.get_or_add_tcPr()
        sh = OxmlElement("w:shd"); sh.set(qn("w:val"), "clear")
        sh.set(qn("w:fill"), "EDE9F2"); pPr.append(sh)
    for ri, name in enumerate(["FILTER", "DP（深度）", "GQ（基因型质量）",
                               "AF（人群频率）", "CONSEQUENCE", "ClinVar"], 1):
        c = t.cell(ri, 0); c.text = ""
        pr = c.paragraphs[0]; r = pr.add_run(name)
        r.font.size = Pt(9.5); cjk(r)
        for ci in (1, 2):
            t.cell(ri, ci).paragraphs[0].paragraph_format.space_after = Pt(14)
    for row in t.rows:
        for ci, w in enumerate(widths):
            row.cells[ci].width = w

    doc.add_paragraph()

    label(doc, "② 这些条件你**先用哪个、后用哪个**？为什么是这个顺序？")
    shade_para(doc,
               "具体一点问：你会先看技术质量（FILTER / DP / GQ），还是先看临床注释（ClinVar）？"
               "先想清楚再写。这一问就是整道题的分水岭。",
               fill="FFF6E8", size=9.5)
    rule_line(doc, 5)
    doc.add_paragraph()

    label(doc, "③ 凭你现在定的规则，你猜哪 1–2 个变异会被选出来？（先猜，别看数据）")
    rule_line(doc, 3)
    doc.add_paragraph()

    # ---------------- closing ----------------
    page_break(doc)
    heading(doc, "写完之后", before=0)
    para(doc,
         "1. 确认四道题的空格都填了，字迹是你自己的。")
    para(doc,
         "2. 再去读 Week4_AI_Analysis.md。那里面有 AI 的分析、我做的数据库核验，"
         "以及 AI 被推翻的几处记录。")
    para(doc,
         "3. 读的时候留意：AI 的结论跟你写的有没有冲突？如果有，**不要直接改你写的**。"
         "作业要的就是「用 AI 之前」和「用 AI 之后」的差别。你最初的判断哪怕不完整，"
         "也是这份作业里最不可替代的一段。")
    para(doc,
         "4. 如果读完之后你改变了看法，那就在最终报告的 Final conclusion 里写清楚"
         "「我原来怎么想、看了什么证据、现在怎么想」——这正是评分表里 25 分"
         "「独立科学推理」想看到的东西。")

    doc.save(OUT)
    print(f"wrote {OUT}")
    d = Document(OUT)
    print(f"  paragraphs : {len(d.paragraphs)}")
    print(f"  tables     : {len(d.tables)}")


if __name__ == "__main__":
    build()
