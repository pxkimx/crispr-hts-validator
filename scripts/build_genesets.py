"""Build the bundled gene-set files (GMT, gzip) from public sources with open licences.

  Gene Ontology  go-basic.obo + goa_human.gaf.gz + mgi.gaf.gz   (geneontology.org, CC BY 4.0)
  Reactome       ReactomePathways.gmt                            (reactome.org, CC BY 4.0)

Annotations marked NOT are left out; every other evidence code is kept. Each GO term collects the genes annotated to
it or to any of its descendants (is_a and part_of). Only sets with 10 to 500 genes are written.

    python scripts/build_genesets.py <download folder>
"""
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

MIN, MAX = 10, 500
OUT = Path(__file__).resolve().parent.parent / "crisprhts" / "data" / "genesets"


def read_obo(p):
    terms, cur = {}, None
    for ln in open(p, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln == "[Term]":
            cur = {"parents": set(), "obsolete": False}
        elif ln.startswith("[") and ln.endswith("]"):
            cur = None
        elif cur is not None:
            if ln.startswith("id: "):
                cur["id"] = ln[4:]
                terms[cur["id"]] = cur
            elif ln.startswith("name: "):
                cur["name"] = ln[6:]
            elif ln.startswith("namespace: "):
                cur["ns"] = ln[11:]
            elif ln.startswith("is_a: "):
                cur["parents"].add(ln[6:].split(" ")[0])
            elif ln.startswith("relationship: part_of "):
                cur["parents"].add(ln.split(" ")[2])
            elif ln.startswith("is_obsolete: true"):
                cur["obsolete"] = True
    return terms


def ancestors(terms):
    memo = {}

    def up(t):
        if t in memo:
            return memo[t]
        out = {t}
        for p in terms.get(t, {}).get("parents", ()):
            out |= up(p)
        memo[t] = out
        return out
    sys.setrecursionlimit(10000)
    return up


def read_gaf(p):
    ann = defaultdict(set)
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("!"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or "NOT" in f[3].split("|"):
                continue
            ann[f[2]].add(f[4])
    return ann


def write_gmt(path, sets, source):
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for name, desc, genes in sets:
            if MIN <= len(genes) <= MAX:
                out.write("\t".join([name, desc] + sorted(genes)) + "\n")
                n += 1
    print(f"{path.name}: {n} sets ({source})")


def main(dl):
    dl = Path(dl)
    terms = read_obo(dl / "go-basic.obo")
    version = next((ln.split(": ", 1)[1].strip() for ln in open(dl / "go-basic.obo") if ln.startswith("data-version")), "")
    up = ancestors(terms)
    ns_short = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
    for species, gaf in (("human", "goa_human.gaf.gz"), ("mouse", "mgi.gaf.gz")):
        ann = read_gaf(dl / gaf)
        by_term = defaultdict(set)
        for gene, gos in ann.items():
            for go in gos:
                for a in up(go):
                    by_term[a].add(gene)
        for ns, short in ns_short.items():
            sets = [(f"{t} {terms[t]['name']}", f"GO:{short}", g) for t, g in by_term.items()
                    if t in terms and terms[t].get("ns") == ns and not terms[t]["obsolete"]]
            write_gmt(OUT / f"GO_{short}_{species}.gmt.gz", sorted(sets), f"Gene Ontology {version}")
    sets = []
    for ln in open(dl / "ReactomePathways.gmt", encoding="utf-8"):
        f = ln.rstrip("\n").split("\t")
        if len(f) > 2 and f[1].startswith("R-HSA-"):
            sets.append((f"{f[1]} {f[0]}", "Reactome", {g for g in f[2:] if re.match(r"^[A-Z0-9][A-Za-z0-9.\-]*$", g)}))
    write_gmt(OUT / "Reactome_human.gmt.gz", sets, "Reactome")
    (OUT / "SOURCES.txt").write_text(
        "Gene sets bundled with CRISPR HTS Validator (built by scripts/build_genesets.py)\n\n"
        f"Gene Ontology ({version}): go-basic.obo, goa_human.gaf, mgi.gaf from geneontology.org. Licence: CC BY 4.0.\n"
        "  Cite: The Gene Ontology Consortium, Genetics 2023; Ashburner et al., Nat Genet 2000.\n"
        "Reactome (current release at build time): ReactomePathways.gmt from reactome.org. Licence: CC BY 4.0.\n"
        "  Cite: Milacic et al., Nucleic Acids Res 2024.\n"
        f"Sets with {MIN} to {MAX} genes. GO sets include genes annotated to descendant terms (is_a, part_of); "
        "NOT annotations are excluded.\n")


if __name__ == "__main__":
    main(sys.argv[1])
