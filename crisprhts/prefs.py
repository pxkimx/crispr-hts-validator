"""Settings that belong to the person, not to a project: the quality thresholds they have chosen.

Kept in ~/Library/Application Support/CRISPR HTS Validator/preferences.json so a change survives closing the
app and loading other data. Only values whose names appear in the given defaults are stored, and only when they
differ from the default, so a later version's new thresholds are picked up instead of being frozen at today's.
Nothing here is sent anywhere; a damaged file is ignored rather than stopping the app.
"""
from __future__ import annotations

import json
import os
import tempfile
import threading
from pathlib import Path

HOME = Path.home() / "Library" / "Application Support" / "CRISPR HTS Validator"
_lock = threading.RLock()


def path() -> Path:
    """Read at every call, not at import, so a test run can point CH_PREFS elsewhere."""
    return Path(os.environ.get("CH_PREFS") or (HOME / "preferences.json"))


def _read() -> dict:
    try:
        with open(path(), encoding="utf-8") as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def _write(d: dict) -> None:
    f = path()
    f.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(f.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(d, fh, indent=1, sort_keys=True)
        os.replace(tmp, f)              # atomic: a crash mid-write cannot leave a half file
    except OSError:
        try:
            os.unlink(tmp)
        except OSError:
            pass


def get(section: str, defaults: dict) -> dict:
    """The saved values for this section that are still known thresholds and still differ from the default."""
    with _lock:
        saved = _read().get(section) or {}
    out = {}
    for k, v in saved.items():
        if k not in defaults:
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if f != float(defaults[k]):
            out[k] = f
    return out


def set(section: str, values: dict, defaults: dict) -> dict:
    """Save the values that differ from the default; drop the rest. Returns what is now saved."""
    keep = {}
    for k, v in (values or {}).items():
        if k not in defaults:
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if f != float(defaults[k]):
            keep[k] = f
    with _lock:
        d = _read()
        if keep:
            d[section] = keep
        else:
            d.pop(section, None)
        _write(d)
    return keep


def clear(section: str) -> None:
    with _lock:
        d = _read()
        if d.pop(section, None) is not None:
            _write(d)
