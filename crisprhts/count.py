"""Counting guide reads in FASTQ files against a guide library.

Default behaviour reproduces `mageck count` 0.5.9.5 (exact match): the guide's start in the read is found
automatically from the first 100,000 reads (the most common 5' trim lengths), and each read is matched at those
offsets, longest guide length first. Extras, all off unless chosen or needed:
  * one mismatch: a read that matches no guide exactly is assigned to the single guide one base away (reads one
    base away from two guides are left unassigned);
  * reverse-complement reads (when the guide is read on the other strand);
  * scaffold search: find a flanking sequence (default the SpCas9 scaffold start GTTTTAGA) and take the guide
    right before it, for reads where the guide does not sit at a fixed offset.
"""
from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .pooled_io import Library

COMP = str.maketrans("ACGTN", "TGCAN")
SCAFFOLD = "GTTTTAGA"
TEST_READS = 100_000


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def _open(path):
    p = str(path)
    if p.lower().endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")


def iter_reads(path, limit: int | None = None):
    with _open(path) as fh:
        n = 0
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()
                n += 1
                if limit and n >= limit:
                    return


@dataclass
class CountResult:
    sample: str
    file: str
    counts: np.ndarray
    reads: int
    mapped: int
    mode: str
    trims: list[int] = field(default_factory=list)
    reverse: bool = False
    mismatch_reads: int = 0
    notes: list[str] = field(default_factory=list)

    @property
    def mapped_fraction(self) -> float:
        return self.mapped / self.reads if self.reads else 0.0


class Matcher:
    def __init__(self, lib: Library):
        self.index = {s: i for i, s in enumerate(lib.seqs)}
        self.lengths = sorted({len(s) for s in lib.seqs}, reverse=True)
        self.minlen = min(self.lengths)
        self._halves = None
        self.seqs = lib.seqs

    # ---- auto trim, exactly as mageckcount_trim5_auto
    def trim_auto(self, path, reverse: bool = False) -> tuple[list[int], float, int]:
        trimkey: Counter = Counter()
        nread = 0
        idx, lens, minlen = self.index, self.lengths, self.minlen
        for r in iter_reads(path, TEST_READS):
            nread += 1
            if reverse:
                r = revcomp(r)
            for t in range(len(r) - minlen + 1):
                rt = r[t:]
                hit = False
                for L in lens:
                    if len(rt) < L:
                        continue
                    if rt[:L] in idx:
                        trimkey[t] += 1
                        hit = True
                        break
                if hit:
                    break
        mapped = sum(trimkey.values())
        frac = mapped / nread if nread else 0.0
        cands = []
        if mapped:
            ks = sorted(trimkey.items(), key=lambda z: z[1], reverse=True)
            fracs = []
            last = 0.01
            for t, c in ks:
                f = c / mapped
                if last > f * 3.0:
                    break
                cands.append(t)
                fracs.append(f)
                if sum(fracs) > 0.75:
                    break
                if cands and f < 0.05:
                    break
                last = f
        return cands, frac, nread

    def scaffold_test(self, path, flank: str, reverse: bool = False) -> float:
        idx, hit, n = self.index, 0, 0
        for r in iter_reads(path, 20_000):
            n += 1
            if reverse:
                r = revcomp(r)
            if self._scaffold_match(r, flank) is not None:
                hit += 1
        return hit / n if n else 0.0

    def _scaffold_match(self, r: str, flank: str):
        p = r.find(flank)
        while p >= 0:
            for L in self.lengths:
                if p >= L:
                    i = self.index.get(r[p - L:p])
                    if i is not None:
                        return i
            p = r.find(flank, p + 1)
        return None

    # ---- one mismatch (pigeonhole on the two halves, unique best hit only)
    def _build_halves(self):
        left: dict = {}
        right: dict = {}
        for i, s in enumerate(self.seqs):
            h = len(s) // 2
            left.setdefault((len(s), s[:h]), []).append(i)
            right.setdefault((len(s), s[h:]), []).append(i)
        self._halves = (left, right)

    def one_mismatch(self, frag: str, L: int):
        if self._halves is None:
            self._build_halves()
        left, right = self._halves
        h = L // 2
        cand = set(left.get((L, frag[:h]), ())) | set(right.get((L, frag[h:]), ()))
        hit = None
        for i in cand:
            s = self.seqs[i]
            d = 0
            for a, b in zip(s, frag):
                if a != b:
                    d += 1
                    if d > 1:
                        break
            if d == 1:
                if hit is not None and hit != i:
                    return None          # ambiguous
                hit = i
        return hit


def count_fastq(path, lib: Library, sample: str = "", mode: str = "auto", trim: int | None = None,
                mismatch: bool = False, flank: str = SCAFFOLD, progress=None, cancel=None) -> CountResult:
    """Count one FASTQ file. mode: 'auto' (MAGeCK offset detection, then reverse complement / scaffold search
    if the guide is not found), 'offset' (fixed `trim`), 'scaffold' (search the flank)."""
    m = Matcher(lib)
    notes: list[str] = []
    reverse = False
    trims: list[int] = []
    use_scaffold = mode == "scaffold"
    if mode == "offset":
        trims = [int(trim or 0)]
    elif mode == "auto":
        trims, frac, _ = m.trim_auto(path)
        if frac < 0.5:
            rtrims, rfrac, _ = m.trim_auto(path, reverse=True)
            if rfrac > max(frac * 1.5, 0.05):
                trims, frac, reverse = rtrims, rfrac, True
                notes.append("Guides were found on the reverse strand of the reads, so reads were reverse-complemented.")
        if frac < 0.5:
            sf = max(m.scaffold_test(path, flank), m.scaffold_test(path, flank, reverse=True))
            if sf > frac * 1.2 and sf > 0.05:
                use_scaffold = True
                reverse = m.scaffold_test(path, flank, reverse=True) > m.scaffold_test(path, flank)
                notes.append(f"The guide is not at a fixed position; it was located by the scaffold sequence {flank}.")
        if not trims and not use_scaffold:
            raise ValueError(f"{Path(path).name}: fewer than 0.1% of the first {TEST_READS:,} reads contain a library "
                             "guide. Check that the library file belongs to these reads.")
        if trims and not use_scaffold:
            notes.append("Guide starts at read position " + ", ".join(str(t + 1) for t in trims) +
                         " (found automatically from the first 100,000 reads, as mageck count does).")
    elif use_scaffold:
        reverse = m.scaffold_test(path, flank, reverse=True) > m.scaffold_test(path, flank)

    counts = np.zeros(len(lib.seqs), dtype=np.int64)
    idx, lens = m.index, m.lengths
    nread = mapped = mm = 0
    for r in iter_reads(path):
        nread += 1
        if nread % 200_000 == 0:
            if progress:
                progress(nread)
            if cancel and cancel():
                raise RuntimeError("Cancelled.")
        if reverse:
            r = revcomp(r)
        hit = None
        if use_scaffold:
            hit = m._scaffold_match(r, flank)
        else:
            for t in trims:
                rt = r[t:]
                for L in lens:
                    if len(rt) < L:
                        continue
                    hit = idx.get(rt[:L])
                    if hit is not None:
                        break
                if hit is not None:
                    break
            if hit is None and mismatch:
                for t in trims:
                    rt = r[t:]
                    for L in lens:
                        if len(rt) < L:
                            continue
                        frag = rt[:L]
                        if "N" in frag and frag.count("N") > 1:
                            continue
                        hit = m.one_mismatch(frag, L)
                        if hit is not None:
                            mm += 1
                            break
                    if hit is not None:
                        break
        if hit is not None:
            counts[hit] += 1
            mapped += 1
    if progress:
        progress(nread)
    if mismatch:
        notes.append(f"{mm:,} reads were assigned with one mismatch.")
    return CountResult(sample or Path(path).name, str(path), counts, nread, mapped,
                       "scaffold" if use_scaffold else ("offset" if mode == "offset" else "auto"),
                       trims, reverse, mm, notes)
