#!/usr/bin/env python3
"""Single data entry: extract the analysis tables from the Mendeley Data workbook.

Reads   data/Shenyang_metro_entrance_acoustic_environment_data.xlsx
Writes  intermediate/clean/<table>.csv  for the four tables the analysis reads

Each cell is written exactly as the workbook stores it (integers without a
decimal point, other numbers in their shortest round-trip form, blanks empty),
so every module re-infers column types from the CSV exactly as the working
pipeline did from its own clean tables.
"""
from __future__ import annotations

import csv
import datetime as dt
import sys
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / "Shenyang_metro_entrance_acoustic_environment_data.xlsx"
OUT = ROOT / "intermediate" / "clean"

# workbook sheet -> clean-table file name the analysis modules read
SHEETS = {
    "station_master": "station_master.csv",
    "station_hour_panel": "station_hour_panel.csv",
    "noise_metrics_station": "noise_metrics_station.csv",
    "landuse_name_map": "landuse_column_map.csv",
}


def cell_text(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, dt.datetime):
        return v.strftime("%Y-%m-%d") if (v.hour, v.minute, v.second) == (0, 0, 0) else v.isoformat(sep=" ")
    if isinstance(v, dt.date):
        return v.isoformat()
    return str(v)


def main() -> None:
    if not XLSX.exists():
        sys.exit(f"Workbook not found: {XLSX}\nDownload it from Mendeley Data (DOI 10.17632/7rwt35tchz) and place it in data/ — see data/README.md")
    OUT.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    for sheet, fname in SHEETS.items():
        ws = wb[sheet]
        rows = ws.iter_rows(values_only=True)
        header = [cell_text(c) for c in next(rows)]
        n = 0
        with open(OUT / fname, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(header)
            for r in rows:
                if all(c is None for c in r):
                    continue
                w.writerow([cell_text(c) for c in r])
                n += 1
        print(f"  {sheet:24s} -> intermediate/clean/{fname}  ({n} rows x {len(header)} cols)")


if __name__ == "__main__":
    main()
