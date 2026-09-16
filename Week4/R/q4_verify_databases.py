#!/usr/bin/env python3
"""
Week 4 - Q4 verification against two authoritative resources.

  1. Ensembl REST (GRCh38 and GRCh37 endpoints)  - gene overlap + reference allele
  2. NCBI ClinVar E-utilities                     - what the cited VCV accession is

This is the "Verification" step of the homework's four-part structure. It is run
as code so the result is reproducible rather than asserted.

Finding: variants_q4.tsv carries no assembly column, and a REF-allele check shows
the records are NOT all on one assembly.

Network note: this machine's DNS is intermittent, so every request is retried.
Output: results/q4_verification.tsv
"""
import json
import os
import ssl
import time
import urllib.request
import urllib.error

PROJ = "." if os.path.isdir("data") else ".."
UA = {"User-Agent": "suat-bioinformatics-hw", "Accept": "application/json"}
CTX = ssl.create_default_context()

GRCH38 = "https://rest.ensembl.org"
GRCH37 = "https://grch37.rest.ensembl.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

SHORTLIST = [
    # gene,   chrom, pos,      ref, alt, clinvar_id,     claimed consequence
    ("TP53",  "17",  7673803,  "G", "A", "VCV000012345", "splice_acceptor_variant"),
    ("KRAS",  "12",  25398284, "C", "A", "VCV000000888", "missense_variant"),
    ("BRCA2", "13",  32316461, "C", "T", "VCV000067890", "missense_variant"),
    ("LDLR",  "19",  11200200, "C", "T", "VCV000000666", "missense_variant"),
]


def get(url, tries=8, pause=3.0):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=25, context=CTX) as r:
                return json.loads(r.read().decode())
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(pause)
    return None


def genes_at(server, chrom, pos):
    d = get(f"{server}/overlap/region/human/{chrom}:{pos}-{pos}"
            f"?feature=gene;content-type=application/json")
    if not d:
        return None
    return sorted({g.get("external_name") or g.get("id") for g in d}) or []


def base_at(server, chrom, pos):
    d = get(f"{server}/sequence/region/human/{chrom}:{pos}-{pos}"
            f"?content-type=application/json")
    return d.get("seq") if d else None


def clinvar(acc):
    uid = acc.lstrip("VC").lstrip("V").lstrip("0") or "0"
    d = get(f"{EUTILS}/esummary.fcgi?db=clinvar&id={uid}&retmode=json")
    if not d or "result" not in d or uid not in d["result"]:
        return None
    r = d["result"][uid]
    germ = r.get("germline_classification", {}) or {}
    return {
        "accession": r.get("accession"),
        "title": r.get("title"),
        "genes": ",".join(g.get("symbol", "") for g in (r.get("genes") or [])),
        "significance": germ.get("description"),
        "review_status": germ.get("review_status"),
    }


def main():
    os.makedirs(os.path.join(PROJ, "results"), exist_ok=True)
    rows = []
    for gene, chrom, pos, ref, alt, vcv, conseq in SHORTLIST:
        g38, b38 = genes_at(GRCH38, chrom, pos), base_at(GRCH38, chrom, pos)
        g37, b37 = genes_at(GRCH37, chrom, pos), base_at(GRCH37, chrom, pos)
        cv = clinvar(vcv)

        ok38 = (g38 is not None and gene in g38 and b38 == ref)
        ok37 = (g37 is not None and gene in g37 and b37 == ref)
        assembly = ("GRCh38" if ok38 and not ok37 else
                    "GRCh37" if ok37 and not ok38 else
                    "ambiguous" if ok37 and ok38 else "NEITHER")

        cv_gene = cv["genes"] if cv else "?"
        cv_match = "yes" if cv and gene in (cv_gene or "") else "NO"

        rows.append(dict(
            gene=gene, variant=f"chr{chrom}:{pos} {ref}>{alt}",
            claimed_consequence=conseq,
            grch38_genes="|".join(g38) if g38 else "(none)",
            grch38_ref=b38 or "?",
            grch37_genes="|".join(g37) if g37 else "(none)",
            grch37_ref=b37 or "?",
            resolves_to=assembly,
            clinvar_acc=vcv,
            clinvar_actual_gene=cv_gene,
            clinvar_gene_matches=cv_match,
            clinvar_review_status=(cv or {}).get("review_status", "?"),
            clinvar_title=(cv or {}).get("title", "?"),
        ))
        print(f"{gene:<6} {f'chr{chrom}:{pos}':<18} "
              f"GRCh38[{'|'.join(g38) if g38 else '-':<14} ref={b38}]  "
              f"GRCh37[{'|'.join(g37) if g37 else '-':<10} ref={b37}]  "
              f"-> {assembly}")
        if cv:
            print(f"       ClinVar {vcv}: {cv['title']}  "
                  f"[gene={cv_gene}, match={cv_match}, "
                  f"review='{cv['review_status']}']")

    out = os.path.join(PROJ, "results", "q4_verification.tsv")
    cols = list(rows[0].keys())
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {out}")

    assemblies = {r["resolves_to"] for r in rows}
    print("\nCONCLUSION: assemblies represented in one table ->",
          ", ".join(sorted(assemblies)))
    if len(assemblies - {"NEITHER"}) > 1:
        print("  The table mixes reference assemblies and declares none. "
              "On real data this blocks interpretation until the build is fixed.")


if __name__ == "__main__":
    main()
