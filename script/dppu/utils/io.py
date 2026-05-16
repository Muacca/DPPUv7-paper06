"""
I/O Utilities
=============

Shared helpers for reading and writing CSV, JSON, and text files, plus
a UTC timestamp utility.  All functions handle directory creation and
enforce UTF-8 encoding to avoid platform-specific codec issues.

Usage::

    from dppu.utils.io import (
        now_utc_iso,
        read_text, write_text,
        read_json, write_json,
        read_csv, write_csv,
        upsert_csv_rows,
    )
"""

from __future__ import annotations

import csv
import json
import os
from datetime import datetime, timezone
from typing import Any, Iterable


def now_utc_iso() -> str:
    """Return a UTC ISO-8601 timestamp (second precision) for artifacts."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


# ── Text ──────────────────────────────────────────────────────────────────────

def read_text(path: str) -> str:
    """Read a UTF-8 text artifact."""
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def write_text(path: str, text: str) -> None:
    """Write a UTF-8 text artifact, creating parent directories as needed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


# ── JSON ──────────────────────────────────────────────────────────────────────

def read_json(path: str) -> dict[str, Any]:
    """Read a UTF-8 JSON artifact."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: str, payload: dict[str, Any]) -> None:
    """Write deterministic UTF-8 JSON (sorted keys, 2-space indent)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2, sort_keys=True)
        f.write("\n")


# ── CSV ───────────────────────────────────────────────────────────────────────

def read_csv(path: str, missing_ok: bool = False) -> list[dict[str, str]]:
    """
    Read a CSV file and return its rows as a list of dicts.

    Args:
        path:       Absolute or relative path to the CSV file.
        missing_ok: When True, return [] instead of raising FileNotFoundError
                    if the file does not exist.

    Returns:
        List of row dicts mapping column name to string value.
    """
    if missing_ok and not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(
    path: str,
    fieldnames: Iterable[str],
    rows: Iterable[dict[str, Any]],
) -> None:
    """
    Write rows to a CSV file, creating parent directories as needed.

    None values in rows are written as empty strings.  Extra fields in row
    dicts beyond *fieldnames* are silently ignored.

    Args:
        path:       Absolute or relative path for the output file.
        fieldnames: Column names in order.
        rows:       Iterable of row dicts.
    """
    fieldnames = list(fieldnames)
    rows = list(rows)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {k: "" if row.get(k) is None else row.get(k) for k in fieldnames}
            )


def upsert_csv_rows(
    path: str,
    fieldnames: list[str],
    rows: Iterable[dict[str, Any]],
    key_fields: tuple[str, ...],
) -> None:
    """
    Upsert rows into a CSV file using a stable composite key.

    Existing rows whose key does not appear in *rows* are preserved.
    Rows with a matching key are replaced.  Output is sorted by key.

    Args:
        path:        Path to the CSV file (created if absent).
        fieldnames:  Ordered column names.
        rows:        New or updated rows to merge in.
        key_fields:  Field names that form the unique row key.
    """
    existing = read_csv(path, missing_ok=True)
    by_key: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in existing:
        by_key[tuple(row.get(k, "") for k in key_fields)] = row
    for row in rows:
        by_key[tuple(str(row.get(k, "")) for k in key_fields)] = row
    ordered = sorted(
        by_key.values(),
        key=lambda r: tuple(str(r.get(k, "")) for k in key_fields),
    )
    write_csv(path, fieldnames, ordered)
