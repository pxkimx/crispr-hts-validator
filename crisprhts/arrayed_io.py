"""Reading arrayed (plate-based) screens.

Supported inputs, all tolerant of mess and reporting what they skipped:
  * cellHTS-style sets: a plate list (Filename, Plate, Replicate[, Channel]), one file per plate and replicate
    (plate name, well, value), a plate configuration (new wildcard format or the old one-row-per-well format),
    an annotation (Plate, Well, reagent ID, gene) and a screen log of flagged wells;
  * per-well tables (CSV/TSV/Excel-exported text) with plate and well (or row + column) columns and one or more
    value columns;
  * plate-shaped grids such as EnVision exports (a header row 01..24 followed by rows A..P), one or more blocks
    per file.
Everything ends up as one long table of wells: plate, replicate, channel, well, row, col, value, content,
reagent, gene, flagged.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .pooled_io import _text

ROWS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
FORMATS = {96: (8, 12), 384: (16, 24), 1536: (32, 48)}
WELL_RX = re.compile(r"^\s*([A-Za-z]{1,2})\s*0*(\d{1,2})\s*$")


def parse_well(w: str) -> tuple[int, int] | None:
    """'B03' -> (1, 2) as zero-based (row, col); 'AA01' rows for 1536-well plates."""
    m = WELL_RX.match(str(w))
    if not m:
        return None
    r = m.group(1).upper()
    if len(r) == 2 and (r[0] != "A" or r[1] > "F"):
        return None                      # two-letter rows exist only on 1536-well plates (AA-AF)
    row = ROWS.index(r) if len(r) == 1 else 26 + ROWS.index(r[1])
    col = int(m.group(2)) - 1
    if col < 0 or col > 47:
        return None
    return row, col


def well_name(row: int, col: int) -> str:
    r = ROWS[row] if row < 26 else "A" + ROWS[row - 26]
    return f"{r}{col + 1:02d}"


def plate_format(n_rows: int, n_cols: int) -> int:
    for wells, (r, c) in FORMATS.items():
        if n_rows <= r and n_cols <= c:
            return wells
    raise ValueError(f"Plates with {n_rows} rows and {n_cols} columns are not supported.")


@dataclass
class ArrayedScreen:
    wells: pd.DataFrame          # one row per plate x well x replicate x channel
    n_wells: int
    name: str = ""
    notes: list[str] = field(default_factory=list)
    description: dict = field(default_factory=dict)

    @property
    def dims(self) -> tuple[int, int]:
        return FORMATS[self.n_wells]

    def summary(self) -> dict:
        w = self.wells
        return {"plates": int(w["plate"].nunique()), "replicates": sorted(int(x) for x in w["replicate"].unique()),
                "channels": sorted(int(x) for x in w["channel"].unique()), "wells": self.n_wells,
                "measurements": int(w["value"].notna().sum()), "contents": w.drop_duplicates(["plate", "well"])["content"].value_counts().to_dict(),
                "genes": int(w.loc[w["content"] == "sample", "gene"].replace("", np.nan).dropna().nunique()),
                "flagged": int(w["flagged"].sum()), "notes": self.notes, "name": self.name,
                "description": self.description}


# --------------------------------------------------------------------------- small readers
def _rows(src, sep=None) -> list[list[str]]:
    text = _text(src)
    lines = [ln for ln in text.split("\n")]
    if sep is None:
        head = "\n".join(lines[:30])
        sep = "\t" if head.count("\t") >= head.count(",") else ","
    if sep == ",":
        return [[c.strip() for c in r] for r in csv.reader(lines)]
    return [[c.strip().strip('"') for c in ln.split("\t")] for ln in lines]


def read_plate_list(src) -> pd.DataFrame:
    rows = [r for r in _rows(src) if any(c for c in r)]
    head = [h.lower() for h in rows[0]]
    if "filename" not in head or "plate" not in head:
        raise ValueError("A plate list needs Filename and Plate columns.")
    df = pd.DataFrame([r[:len(head)] for r in rows[1:] if len(r) >= 2], columns=head[:max(len(r) for r in rows)] if False else None)
    df = pd.DataFrame([dict(zip(head, r)) for r in rows[1:] if len(r) >= 2])
    df["plate"] = df["plate"].astype(int)
    df["replicate"] = df.get("replicate", pd.Series(["1"] * len(df))).replace("", "1").astype(int)
    df["channel"] = df.get("channel", pd.Series(["1"] * len(df))).replace("", "1").astype(int)
    return df[["filename", "plate", "replicate", "channel"]]


def read_plate_file(src) -> tuple[list[tuple[str, float]], int]:
    """cellHTS plate file: plate name, well, value (tab). Returns [(well, value)], skipped lines."""
    out, bad = [], 0
    for r in _rows(src):
        if not any(r):
            continue
        if len(r) >= 3 and parse_well(r[1]):
            try:
                out.append((r[1], float(r[2]) if r[2] not in ("", "NA", "NaN") else float("nan")))
            except ValueError:
                bad += 1
        elif len(r) >= 2 and parse_well(r[0]):
            try:
                out.append((r[0], float(r[1])))
            except ValueError:
                bad += 1
        else:
            bad += 1
    return out, bad


def read_plate_conf(src, n_plates: int | None = None) -> tuple[callable, int | None, list[str]]:
    """Plate configuration. New cellHTS format: 'Wells: 384', 'Plates: 57', then Plate / Well / Content lines with
    '*' or regular expressions; later lines override earlier ones. Old format: Batch / Well / Content, one row per
    well. Returns a function (plate, well) -> content, the number of wells (if stated) and notes."""
    rows = [r for r in _rows(src, "\t") if any(r)]
    notes: list[str] = []
    n_wells = None
    rules = []
    i = 0
    while i < len(rows) and ":" in rows[i][0]:
        k, v = rows[i][0].split(":", 1)
        if k.strip().lower() == "wells":
            n_wells = int(v.strip() or rows[i][1])
        i += 1
    head = [h.lower() for h in rows[i]]
    body = rows[i + 1:]
    if "batch" in head:                     # old format: same layout on every plate (per batch)
        table = {r[1]: r[2] for r in body if len(r) >= 3}
        notes.append("Old-style plate configuration: the same layout was applied to every plate.")
        return (lambda p, w: table.get(w, "sample")), n_wells, notes
    for r in body:
        if len(r) < 3:
            continue
        pp, ww, content = r[0], r[1], r[2]
        prx = None if pp == "*" else re.compile("^(?:" + pp + ")$")
        wrx = None if ww == "*" else re.compile("^(?:" + ww + ")$")
        rules.append((prx, wrx, content))

    def content_of(plate, well):
        c = "sample"
        for prx, wrx, cc in rules:
            if (prx is None or prx.match(str(plate))) and (wrx is None or wrx.match(well)):
                c = cc
        return c
    return content_of, n_wells, notes


def read_annotation(src) -> tuple[pd.DataFrame, list[str]]:
    """Plate, Well, reagent ID, gene (more columns kept as extra). Ragged rows are read for the columns present."""
    rows = [r for r in _rows(src, "\t") if any(c for c in r)]
    head = rows[0]
    hl = [h.lower() for h in head]
    if "plate" not in hl or "well" not in hl:
        raise ValueError("An annotation file needs Plate and Well columns.")
    ip, iw = hl.index("plate"), hl.index("well")
    gene_i = next((hl.index(k) for k in ("geneid", "gene", "genesymbol", "symbol", "gene_symbol") if k in hl), None)
    reag_i = next((i for i, h in enumerate(hl) if i not in (ip, iw, gene_i)), None)
    out, ragged, bad = [], 0, 0
    for n, r in enumerate(rows[1:], start=2):
        if len(r) != len(head):
            ragged += 1
            r = (r + [""] * len(head))[:len(head)]
        try:
            p = int(r[ip])
        except ValueError:
            bad += 1
            continue
        if not parse_well(r[iw]):
            bad += 1
            continue
        g = r[gene_i] if gene_i is not None else ""
        a = r[reag_i] if reag_i is not None else ""
        out.append((p, r[iw].upper(), "" if a in ("NA", "") else a, "" if g in ("NA", "") else g))
    notes = []
    if ragged:
        notes.append(f"Annotation: {ragged} row(s) had a different number of columns than the header "
                     "(e.g. an extra tab); they were read for the columns present.")
    if bad:
        notes.append(f"Annotation: {bad} row(s) skipped (no plate number or well).")
    return pd.DataFrame(out, columns=["plate", "well", "reagent", "gene"]), notes


def read_screen_log(src) -> pd.DataFrame:
    rows = [r for r in _rows(src, "\t") if any(c for c in r)]
    hl = [h.lower() for h in rows[0]]
    recs = []
    for r in rows[1:]:
        d = dict(zip(hl, r))
        try:
            recs.append({"plate": int(d.get("plate")), "replicate": int(d.get("sample") or d.get("replicate") or 0),
                         "channel": int(d.get("channel") or 0), "well": d.get("well", "").upper(),
                         "comment": d.get("comment", "")})
        except (TypeError, ValueError):
            continue
    return pd.DataFrame(recs, columns=["plate", "replicate", "channel", "well", "comment"])


def read_description(src) -> dict:
    out = {}
    for ln in _text(src).split("\n"):
        if ":" in ln and not ln.startswith("["):
            k, v = ln.split(":", 1)
            if v.strip():
                out[k.strip()] = v.strip()
    return out


def read_grids(src) -> list[dict]:
    """Plate-shaped blocks: a header row of column numbers (01..24) then rows starting with A, B, ...
    Returns [{"title", "values" (rows x cols array)}]; the title is the nearest text line above the block."""
    rows = _rows(src)
    blocks = []
    i = 0
    last_title = ""
    while i < len(rows):
        r = rows[i]
        cells = [c for c in r]
        nums = [c for c in cells[1:] if c != ""]
        if cells and cells[0] == "" and nums and all(re.match(r"^\d{1,2}$", c) for c in nums) and int(nums[0]) == 1:
            ncol = len(nums)
            vals = []
            j = i + 1
            while j < len(rows) and rows[j] and re.match(r"^[A-Z]{1,2}$", rows[j][0] or ""):
                row = rows[j][1:ncol + 1]
                vals.append([float(x) if re.match(r"^-?[\d.]+(e-?\d+)?$", x or "", re.I) else np.nan for x in row]
                            + [np.nan] * (ncol - len(row)))
                j += 1
            if vals:
                blocks.append({"title": last_title, "values": np.array(vals, dtype=float)})
            i = j
            continue
        if cells and cells[0] and not re.match(r"^[\d,.\s]+$", cells[0]) and len([c for c in cells if c]) <= 2:
            last_title = cells[0]
        i += 1
    return blocks


# --------------------------------------------------------------------------- assembling a screen
def _finish(df: pd.DataFrame, n_wells: int | None, name: str, notes: list[str], description=None) -> ArrayedScreen:
    rc = df["well"].map(parse_well)
    df["row"] = [x[0] for x in rc]
    df["col"] = [x[1] for x in rc]
    df["well"] = [well_name(r, c) for r, c in zip(df["row"], df["col"])]
    if n_wells is None:
        n_wells = plate_format(int(df["row"].max()) + 1, int(df["col"].max()) + 1)
    for c, d in (("content", "sample"), ("reagent", ""), ("gene", ""), ("flagged", False), ("comment", "")):
        if c not in df:
            df[c] = d
    df["content"] = df["content"].fillna("sample").astype(str)
    df["gene"] = df["gene"].fillna("").astype(str)
    df["reagent"] = df["reagent"].fillna("").astype(str)
    df = df.sort_values(["channel", "replicate", "plate", "row", "col"]).reset_index(drop=True)
    return ArrayedScreen(df, n_wells, name, notes, description or {})


def load_cellhts(files: dict[str, Path], name: str = "") -> ArrayedScreen:
    """`files`: file name (case-insensitive) -> path, containing at least the plate list and the plate files."""
    low = {k.lower(): v for k, v in files.items()}
    pl_key = next((k for k in low if re.match(r"^platelist", k)), None)
    if not pl_key:
        raise ValueError("No plate list (Platelist.txt) among the files.")
    pl = read_plate_list(low[pl_key])
    notes: list[str] = []
    recs = []
    missing = []
    skipped = 0
    for _, r in pl.iterrows():
        p = low.get(r["filename"].lower())
        if p is None:
            missing.append(r["filename"])
            continue
        vals, bad = read_plate_file(p)
        skipped += bad
        for w, v in vals:
            recs.append((r["plate"], r["replicate"], r["channel"], w, v))
    if missing:
        notes.append(f"{len(missing)} plate file(s) named in the plate list were not provided: {', '.join(missing[:5])}"
                     + ("…" if len(missing) > 5 else ""))
    if skipped:
        notes.append(f"{skipped} line(s) in plate files were not well/value lines and were skipped.")
    if not recs:
        raise ValueError("None of the plate files listed in the plate list could be read.")
    df = pd.DataFrame(recs, columns=["plate", "replicate", "channel", "well", "value"])
    n_wells = None
    conf_key = next((k for k in low if re.match(r"^plateconf", k)), None)
    if conf_key:
        content_of, n_wells, cn = read_plate_conf(low[conf_key])
        notes += cn
        cache = {}
        df["content"] = [cache.setdefault((p, w), content_of(p, w)) for p, w in zip(df["plate"], df["well"])]
    else:
        notes.append("No plate configuration: every well is treated as a sample until controls are marked.")
    ann_key = next((k for k, v in low.items() if k not in (pl_key, conf_key) and _looks_like_annotation(v)), None)
    if ann_key:
        ann, an = read_annotation(low[ann_key])
        notes += an
        ann["well"] = [well_name(*parse_well(w)) for w in ann["well"]]
        df["well"] = [well_name(*parse_well(w)) for w in df["well"]]
        df = df.merge(ann, on=["plate", "well"], how="left")
    log_key = next((k for k in low if re.match(r"^screenlog", k)), None)
    if log_key:
        lg = read_screen_log(low[log_key])
        df["well"] = [well_name(*parse_well(w)) for w in df["well"]]
        lg["well"] = [well_name(*parse_well(w)) if parse_well(w) else w for w in lg["well"]]
        flag = set()
        for _, r in lg.iterrows():
            flag.add((r["plate"], r["replicate"], r["channel"], r["well"]))
        df["flagged"] = [((p, rp, ch, w) in flag) or ((p, rp, 0, w) in flag)
                         for p, rp, ch, w in zip(df["plate"], df["replicate"], df["channel"], df["well"])]
        df["comment"] = np.where(df["flagged"], "screen log", "")
        if len(lg):
            notes.append(f"Screen log: {int(df['flagged'].sum())} measurement(s) flagged and left out "
                         f"({', '.join(sorted(set(lg['comment'].astype(str))))[:80]}).")
    desc_key = next((k for k in low if re.match(r"^description", k)), None)
    desc = read_description(low[desc_key]) if desc_key else {}
    return _finish(df, n_wells, name or desc.get("Screen", ""), notes, desc)


def _looks_like_annotation(p: Path) -> bool:
    try:
        head = _text(p).split("\n", 1)[0].lower()
    except OSError:
        return False
    return "plate" in head and "well" in head and ("gene" in head or "id" in head) and "content" not in head and "comment" not in head


def load_well_table(src, name: str = "") -> ArrayedScreen:
    """A per-well table. Plate and well columns are found by name (or row + column); every other numeric column is
    a measurement (named like rep1/replicate2/R1 -> replicates, otherwise channels)."""
    rows = [r for r in _rows(src) if any(c for c in r)]
    head = rows[0]
    hl = [h.lower().strip() for h in head]
    df = pd.DataFrame([(r + [""] * len(head))[:len(head)] for r in rows[1:]], columns=head)
    notes = []

    def col(*names):
        return next((head[i] for i, h in enumerate(hl) if h in names), None)
    pc = col("plate", "plate_id", "plateid", "plate number", "plate_number", "barcode")
    wc = col("well", "well_id", "wellid", "position", "well name")
    rc, cc = col("row"), col("col", "column")
    if pc is None or (wc is None and (rc is None or cc is None)):
        raise ValueError("The table needs a plate column and a well column (or row and column columns).")
    if wc is None:
        df["__well"] = [f"{r.strip().upper() if not r.strip().isdigit() else ROWS[int(r) - 1]}{int(float(c)):02d}"
                        for r, c in zip(df[rc], df[cc])]
        wc = "__well"
    gene_c = col("gene", "geneid", "gene_symbol", "symbol", "target", "gene id")
    reag_c = col("reagent", "sirna", "sgrna", "guide", "id", "reagent_id", "hfaid", "oligo")
    cont_c = col("content", "type", "control", "well type", "welltype")
    plates = pd.factorize(df[pc])[0] + 1 if not pd.to_numeric(df[pc], errors="coerce").notna().all() else pd.to_numeric(df[pc]).astype(int)
    if not pd.to_numeric(df[pc], errors="coerce").notna().all():
        notes.append("Plate names were numbered in order of appearance.")
    meas = []
    for h in head:
        if h in (pc, wc, rc, cc, gene_c, reag_c, cont_c):
            continue
        v = pd.to_numeric(df[h], errors="coerce")
        if v.notna().mean() > 0.5:
            meas.append(h)
    if not meas:
        raise ValueError("No numeric measurement columns were found.")
    rep_like = [h for h in meas if re.search(r"(rep|replicate|^r)\s*_?\d+$", h, re.I)]
    recs = []
    bad = int((~df[wc].map(lambda w: parse_well(w) is not None)).sum())
    ok = df[wc].map(lambda w: parse_well(w) is not None)
    for k, h in enumerate(meas):
        rep, ch = (rep_like.index(h) + 1, 1) if h in rep_like else (1, k + 1)
        if rep_like and h not in rep_like:
            continue
        sub = pd.DataFrame({"plate": np.asarray(plates)[ok.to_numpy()], "replicate": rep, "channel": ch,
                            "well": df.loc[ok, wc].to_numpy(), "value": pd.to_numeric(df.loc[ok, h], errors="coerce").to_numpy()})
        if gene_c:
            sub["gene"] = df.loc[ok, gene_c].to_numpy()
        if reag_c:
            sub["reagent"] = df.loc[ok, reag_c].to_numpy()
        if cont_c:
            sub["content"] = [_norm_content(x) for x in df.loc[ok, cont_c]]
        recs.append(sub)
    if bad:
        notes.append(f"{bad} row(s) skipped: the well could not be read.")
    used = rep_like or meas
    notes.append("Measurements: " + ", ".join(f"{h} → {'replicate ' + str(rep_like.index(h) + 1) if h in rep_like else 'channel ' + str(meas.index(h) + 1)}" for h in used))
    return _finish(pd.concat(recs, ignore_index=True), None, name, notes)


def _norm_content(x: str) -> str:
    s = str(x).strip()
    l = s.lower()
    if l in ("", "sample", "samples", "library", "test", "compound", "na", "nan"):
        return "sample"
    if l in ("neg", "negative", "negative control", "neg ctrl", "negctrl", "nc", "ntc", "non-targeting", "nontargeting"):
        return "neg"
    if l in ("pos", "positive", "positive control", "pos ctrl", "posctrl", "pc"):
        return "pos"
    if l in ("empty", "blank", "none"):
        return "empty"
    return s


def load_grids(files: dict[str, Path], name: str = "", block: int | None = None) -> ArrayedScreen:
    """Plate-shaped exports. With a plate list, each file is one plate/replicate; without, every grid block in every
    file is a separate plate, numbered in order."""
    low = {k.lower(): v for k, v in files.items()}
    pl_key = next((k for k in low if re.match(r"^platelist", k)), None)
    notes = []
    recs = []
    order = []
    if pl_key:
        pl = read_plate_list(low[pl_key])
        order = [(r["filename"], r["plate"], r["replicate"], r["channel"]) for _, r in pl.iterrows()]
    else:
        data_files = [k for k in files if k.lower() not in (pl_key,) and not re.match(r"^(plateconf|screenlog|description)", k.lower())]
        order = [(f, None, 1, 1) for f in sorted(data_files)]
    plate_no = 0
    for fname, plate, rep, ch in order:
        p = low.get(fname.lower())
        if p is None:
            notes.append(f"{fname} (named in the plate list) was not provided.")
            continue
        blocks = read_grids(p)
        if not blocks:
            continue
        titles = [b["title"] for b in blocks]
        if plate is not None:
            pick = block if block is not None and block < len(blocks) else _default_block(titles)
            blocks = [blocks[pick]]
            if len(titles) > 1:
                notes.append(f"{fname}: {len(titles)} measurement blocks found; used “{titles[pick] or 'block ' + str(pick + 1)}”. "
                              f"Others: {'; '.join(t for i, t in enumerate(titles) if i != pick)[:160]}")
        for b in blocks:
            plate_no += 1
            pn = plate if plate is not None else plate_no
            v = b["values"]
            for i in range(v.shape[0]):
                for j in range(v.shape[1]):
                    recs.append((pn, rep, ch, well_name(i, j), v[i, j]))
    if not recs:
        raise ValueError("No plate-shaped data blocks (a row of column numbers followed by rows A, B, C…) were found.")
    df = pd.DataFrame(recs, columns=["plate", "replicate", "channel", "well", "value"])
    conf_key = next((k for k in low if re.match(r"^plateconf", k)), None)
    n_wells = None
    if conf_key:
        content_of, n_wells, cn = read_plate_conf(low[conf_key])
        notes += cn
        df["content"] = [content_of(p, w) for p, w in zip(df["plate"], df["well"])]
    return _finish(df, n_wells, name, notes)


def _default_block(titles: list[str]) -> int:
    """Prefer the instrument's main result over derived blocks (EnVision: 'Results for …' over 'Calculated results')."""
    for i, t in enumerate(titles):
        if t.lower().startswith("results for"):
            return i
    return 0


def detect_and_load(files: dict[str, Path], name: str = "") -> ArrayedScreen:
    """Choose the reader from the files given."""
    keys = [k.lower() for k in files]
    has_pl = any(re.match(r"^platelist", k) for k in keys)
    if has_pl:
        pl = read_plate_list(next(v for k, v in files.items() if re.match(r"^platelist", k.lower())))
        first = next((v for k, v in files.items() if k.lower() == pl["filename"].iloc[0].lower()), None)
        if first is not None and read_grids(first):
            return load_grids(files, name)
        return load_cellhts(files, name)
    data = [v for k, v in files.items()]
    if len(data) == 1:
        p = data[0]
        if read_grids(p):
            return load_grids(files, name)
        return load_well_table(p, name)
    if all(read_grids(p) for p in data):
        return load_grids(files, name)
    raise ValueError("Several files without a plate list: add Platelist.txt (Filename, Plate, Replicate), "
                     "or combine the measurements into one per-well table.")
