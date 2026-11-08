"""Reading pooled-screen inputs: count tables, guide libraries, sample sheets and BAGEL screen-info files.

Everything is tolerant of small mess (header or no header, tab or comma, Windows line ends) and reports what it
skipped, so the user can see exactly what was read.
"""
from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ACGT = re.compile(r"^[ACGTNacgtn]+$")


@dataclass
class CountTable:
    guides: list[str]
    genes: list[str]
    samples: list[str]
    counts: np.ndarray                 # guides x samples, float
    fmt: str = "generic"               # "MAGeCK", "BAGEL" or "generic"
    source: str = ""
    notes: list[str] = field(default_factory=list)
    reads: dict | None = None          # sample -> total reads (FASTQ input only)
    mapped: dict | None = None         # sample -> mapped reads (FASTQ input only)

    def frame(self) -> pd.DataFrame:
        df = pd.DataFrame(self.counts, columns=self.samples)
        df.insert(0, "gene", self.genes)
        df.insert(0, "guide", self.guides)
        return df


def _text(src) -> str:
    if isinstance(src, (bytes, bytearray)):
        b = bytes(src)
    else:
        p = Path(src)
        b = p.read_bytes()
        if p.suffix.lower() == ".gz":
            import gzip
            b = gzip.decompress(b)
    if b[:3] == b"\xef\xbb\xbf":
        b = b[3:]
    return b.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")


def _delim(lines: list[str]) -> str:
    head = "\n".join(lines[:20])
    if head.count("\t") >= max(1, len(lines[:20]) // 2):
        return "\t"
    if head.count(",") >= max(1, len(lines[:20]) // 2):
        return ","
    return "\t" if "\t" in head else ("," if "," in head else None)


def _isnum(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_count_table(src, name: str = "") -> CountTable:
    """Read a MAGeCK- or BAGEL-style count table (guide, gene, sample columns...).

    Follows `mageck test`'s reader where it matters for agreement: the first line is a header when its sample
    fields are not numbers; a repeated guide ID keeps the first record; rows with a different number of fields
    are skipped and reported.
    """
    text = _text(src)
    lines = [ln for ln in text.split("\n") if ln.strip() != ""]
    if not lines:
        raise ValueError("The count table is empty.")
    sep = _delim(lines)
    if sep is None:
        raise ValueError("Could not find columns in the count table (expected tab- or comma-separated text).")
    rows = list(csv.reader(lines, delimiter=sep)) if sep == "," else [ln.split("\t") for ln in lines]
    rows = [[c.strip() for c in r] for r in rows]
    first = rows[0]
    if len(first) < 3:
        raise ValueError("A count table needs at least three columns: guide ID, gene and one sample.")
    has_header = not all(_isnum(x) for x in first[2:] if x != "")
    notes: list[str] = []
    if has_header:
        samples = [s or f"sample{i + 1}" for i, s in enumerate(first[2:])]
        body = rows[1:]
        h0 = first[0].upper()
    else:
        samples = [f"sample{i + 1}" for i in range(len(first) - 2)]
        body = rows
        h0 = ""
    if h0 in ("SEQUENCE", "GENE_CLONE", "SEQ", "CLONE"):
        fmt = "BAGEL"
    elif h0 in ("SGRNA", "SGRNA_ID", "GUIDE", "GUIDE_ID", "ID", "SGRNAID"):
        fmt = "MAGeCK"
    else:
        fmt = "generic"
    n = len(samples)
    seen: set[str] = set()
    guides, genes, vals = [], [], []
    bad_len = bad_num = dup = empty_gene = 0
    for i, r in enumerate(body):
        if len(r) - 2 != n:
            bad_len += 1
            continue
        g = r[0]
        if g in seen:
            dup += 1
            continue
        if r[1] == "":
            empty_gene += 1
            continue
        try:
            v = [float(x) for x in r[2:]]
        except ValueError:
            bad_num += 1
            continue
        seen.add(g)
        guides.append(g)
        genes.append(r[1])
        vals.append(v)
    if not guides:
        raise ValueError("No guide rows could be read from the count table.")
    if bad_len:
        notes.append(f"{bad_len} row(s) skipped: wrong number of columns.")
    if dup:
        notes.append(f"{dup} row(s) skipped: repeated guide ID (the first record was kept, as MAGeCK does).")
    if bad_num:
        notes.append(f"{bad_num} row(s) skipped: a count was not a number.")
    if empty_gene:
        notes.append(f"{empty_gene} row(s) skipped: empty gene name.")
    counts = np.asarray(vals, dtype=float)
    if (counts < 0).any():
        raise ValueError("The count table has negative counts.")
    if np.any(np.abs(counts - np.round(counts)) > 1e-9):
        notes.append("Counts are not whole numbers (normalised or averaged counts?). They are used as given.")
    return CountTable(guides, genes, samples, counts, fmt, name, notes)


@dataclass
class Library:
    ids: list[str]
    seqs: list[str]
    genes: list[str]
    chrom: list[str] | None = None
    pos: list[float] | None = None
    notes: list[str] = field(default_factory=list)


def read_library(src, name: str = "") -> Library:
    """Read a guide library: ID, sequence, gene (MAGeCK order), with or without a header, CSV or TSV.

    The sequence column is found by content (only A/C/G/T/N), so other column orders work too. Optional
    chromosome and position columns (named chr/chrom/chromosome and start/pos/position) are kept for
    copy-number correction. Like `mageck count`: repeated IDs and repeated sequences keep the first record.
    """
    text = _text(src)
    lines = [ln for ln in text.split("\n") if ln.strip() != ""]
    if not lines:
        raise ValueError("The library file is empty.")
    sep = "," if (name or str(src)).lower().endswith(".csv") else (_delim(lines) or "\t")
    rows = list(csv.reader(lines, delimiter=sep)) if sep == "," else [ln.split("\t") for ln in lines]
    rows = [[c.strip() for c in r] for r in rows if len(r) >= 3]
    if not rows:
        raise ValueError("A library file needs at least three columns: ID, sequence and gene.")
    sample = rows[1:200] if len(rows) > 1 else rows
    ncol = min(len(r) for r in sample)
    seq_col = None
    for c in range(ncol):
        if all(ACGT.match(r[c] or "x") for r in sample):
            seq_col = c
            break
    if seq_col is None:
        raise ValueError("Could not find a column of guide sequences (A, C, G, T) in the library file.")
    header = rows[0] if not ACGT.match(rows[0][seq_col] or "x") else None
    body = rows[1:] if header else rows
    others = [c for c in range(ncol) if c != seq_col]
    id_col, gene_col = others[0], others[1]
    chrom_col = pos_col = None
    if header:
        low = [h.lower() for h in header]
        for c, h in enumerate(low):
            if h in ("gene", "gene_id", "genesymbol", "gene_symbol", "symbol", "target", "target_gene"):
                gene_col = c
            if h in ("id", "sgrna", "sgrna_id", "guide", "guide_id", "name"):
                id_col = c
            if h in ("chr", "chrom", "chromosome", "chrm"):
                chrom_col = c
            if h in ("start", "pos", "position", "cut", "cut_position", "startpos"):
                pos_col = c
        if gene_col == id_col:
            gene_col = [c for c in others if c != id_col][0]
    ids, seqs, genes, chrom, pos = [], [], [], [], []
    seen_id: set[str] = set()
    seen_seq: set[str] = set()
    dup_id = dup_seq = 0
    for r in body:
        if len(r) <= max(seq_col, id_col, gene_col):
            continue
        i, s, g = r[id_col], r[seq_col].upper(), r[gene_col]
        if i in seen_id:
            dup_id += 1
            continue
        if s in seen_seq:
            dup_seq += 1
            continue
        seen_id.add(i)
        seen_seq.add(s)
        ids.append(i)
        seqs.append(s)
        genes.append(g)
        if chrom_col is not None and pos_col is not None and len(r) > max(chrom_col, pos_col):
            chrom.append(r[chrom_col])
            try:
                pos.append(float(r[pos_col]))
            except ValueError:
                pos.append(float("nan"))
    notes = []
    if dup_id:
        notes.append(f"{dup_id} guide(s) skipped: repeated ID.")
    if dup_seq:
        notes.append(f"{dup_seq} guide(s) skipped: same sequence as an earlier guide (as `mageck count` does).")
    has_pos = chrom_col is not None and pos_col is not None and len(pos) == len(ids)
    return Library(ids, seqs, genes, chrom if has_pos else None, pos if has_pos else None, notes)


# ------------------------------------------------------------------ samples and comparisons
REF_HINT = re.compile(r"(^|[^a-z])(t0|d0|day0|day_0|plasmid|pdna|initial|pre|baseline|library|lib|ctrl0)([^a-z]|$)",
                      re.I)


REP_TOKEN = re.compile(r"^(?:r|rep|replicate)?(?:[A-Za-z]|\d{1,2})$", re.I)


def guess_conditions(samples: list[str]) -> dict[str, str]:
    """Group replicates by dropping one short replicate token, when that leaves a name shared by other samples.

    T18_A, T18_B -> T18; T18_A_Starved, T18_B_Starved -> T18_Starved; HL60.initial stays itself.
    """
    split = {s: re.split(r"([._\- ])", s.strip()) for s in samples}
    cands: dict[str, list[str]] = {}
    for s, parts in split.items():
        opts = []
        for i in range(0, len(parts), 2):
            if REP_TOKEN.match(parts[i]) and len(parts) > 1:
                rest = parts[:i] + parts[i + 2:] if i + 1 < len(parts) else parts[:max(0, i - 1)]
                name = "".join(rest).strip("._- ")
                if name:
                    opts.append(name)
        if re.search(r"\d[A-Za-z]$", s.strip()):      # HAP1_T18A -> HAP1_T18
            opts.append(s.strip()[:-1])
        cands[s] = opts
    count: dict[str, int] = {}
    for opts in cands.values():
        for o in set(opts):
            count[o] = count.get(o, 0) + 1
    out = {}
    for s in samples:
        shared = [o for o in cands[s] if count.get(o, 0) >= 2]
        out[s] = shared[0] if shared else s
    return out


def guess_design(samples: list[str]) -> dict:
    """A first guess at the sample sheet; the user confirms or edits it in the app."""
    rows = []
    gc = guess_conditions(samples)
    for s in samples:
        cond = gc[s]
        rows.append({"sample": s, "condition": cond, "reference": bool(REF_HINT.search(s))})
    conds = {}
    for r in rows:
        conds.setdefault(r["condition"], []).append(r["sample"])
    for r in rows:
        r["replicate"] = conds[r["condition"]].index(r["sample"]) + 1
    refs = sorted({r["condition"] for r in rows if r["reference"]})
    comps = []
    if refs:
        def tokens(x):
            return {t.lower() for t in re.split(r"[._\- ]", x) if t}
        for c in conds:
            if c in refs:
                continue
            # pair with the reference that shares most name parts (HL60.final -> HL60.initial)
            ref = max(refs, key=lambda r: (len(tokens(r) & tokens(c)), -refs.index(r)))
            comps.append({"name": f"{c} vs {ref}", "treatment": conds[c], "control": conds[ref]})
    return {"samples": rows, "comparisons": comps}


def read_screen_info(src) -> list[dict]:
    """BAGEL screen-info file: name <tab> replicates (comma) <tab> controls (comma)."""
    comps = []
    for ln in _text(src).split("\n"):
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        f = [x.strip() for x in ln.split("\t") if x.strip() != ""]
        if len(f) < 3:
            f = ln.split()
        if len(f) < 3:
            continue
        comps.append({"name": f[0], "treatment": [x for x in f[1].split(",") if x],
                      "control": [x for x in f[2].split(",") if x]})
    if not comps:
        raise ValueError("No screens found in the screen-info file (expected: name, replicates, controls).")
    return comps


def read_gene_list(src) -> set[str]:
    out = set()
    for i, ln in enumerate(_text(src).split("\n")):
        f = ln.strip().split("\t")[0].split(",")[0].strip()
        if not f or (i == 0 and f.upper() in ("GENE", "GENES", "SYMBOL")):
            continue
        out.add(f)
    return out


def to_csv_bytes(df: pd.DataFrame, sep: str = ",") -> bytes:
    buf = io.StringIO()
    df.to_csv(buf, index=False, sep=sep)
    return buf.getvalue().encode()
