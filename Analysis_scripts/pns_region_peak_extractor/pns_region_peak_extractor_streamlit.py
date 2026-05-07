#!/usr/bin/env python3
"""
PNS region peak extractor

CLI example:
    python pns_region_peak_extractor_streamlit.py \
        --data-dir /mnt/d/Snyder_bams/CH01_PNS/PNS/ \
        --bed input_regions.bed \
        --peak-flank-bp 1000 \
        --out-prefix CH01_test

Streamlit example:
    streamlit run pns_region_peak_extractor_streamlit.py

Required Python packages:
    pyBigWig pandas numpy streamlit
"""

from __future__ import annotations

import argparse
import bisect
import io
import os
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import pyBigWig

try:
    import streamlit as st
except ImportError:
    st = None


DEFAULT_COVERAGE_BW = "CH01_chrAll_PNS_mode167_lower120_upper180_coverage_meanNorm_x100.bw"
DEFAULT_PNS_BW = "CH01_chrAll_PNS_mode167_lower120_upper180_pns_meanPeakNorm_x100.bw"
DEFAULT_BREAKPOINT_BB = "CH01_chrAll_PNS_mode167_lower120_upper180_breakpoint_peaks_meanPeakNorm_x100.bb"
DEFAULT_NUCLEOSOME_BB = "CH01_chrAll_PNS_mode167_lower120_upper180_nucleosome_regions_meanPeakNorm_x100.bb"


@dataclass(frozen=True)
class BedRegion:
    chrom: str
    start: int
    end: int

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


@dataclass(frozen=True)
class PeakRecord:
    chrom: str
    start: int
    end: int
    center: int
    score: str


# -----------------------------
# BED parsing
# -----------------------------

def parse_bed_text(text: str) -> tuple[list[BedRegion], list[str]]:
    regions: list[BedRegion] = []
    skipped: list[str] = []

    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line.lower().startswith("track"):
            continue
        if line.lower().startswith("browser"):
            continue

        parts = line.split()
        if len(parts) < 3:
            skipped.append("line " + str(line_no) + ": fewer than 3 columns: " + raw)
            continue

        chrom = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            skipped.append("line " + str(line_no) + ": non-integer start/end: " + raw)
            continue

        if start < 0 or end <= start:
            skipped.append("line " + str(line_no) + ": invalid interval: " + raw)
            continue

        regions.append(BedRegion(chrom=chrom, start=start, end=end))

    return regions, skipped


# -----------------------------
# Chromosome handling
# -----------------------------

def resolve_chrom_name(chrom: str, chroms: dict[str, int]) -> Optional[str]:
    """Strict chromosome matching only."""
    if chrom in chroms:
        return chrom
    return None


# -----------------------------
# bigBed peak extraction using UCSC bigBedToBed
# -----------------------------

def parse_bigbedtobed_line(line: str) -> Optional[PeakRecord]:
    parts = line.rstrip().split("	")
    if len(parts) < 3:
        return None

    try:
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None

    score = "."
    if len(parts) >= 5 and parts[4] != "":
        score = parts[4]

    center = (start + end) // 2
    if len(parts) >= 8:
        try:
            center = int(parts[6])
        except ValueError:
            center = (start + end) // 2

    return PeakRecord(
        chrom=chrom,
        start=start,
        end=end,
        center=center,
        score=score,
    )


def query_bigbed_peaks_with_ucsc(
    bigbed_path: str,
    chrom: str,
    center: int,
    flank_bp: int,
    bigbedtobed_exe: str,
) -> tuple[str, list[PeakRecord]]:
    query_start = max(0, center - flank_bp)
    query_end = center + flank_bp

    cmd = [
        bigbedtobed_exe,
        bigbed_path,
        "stdout",
        "-chrom=" + chrom,
        "-start=" + str(query_start),
        "-end=" + str(query_end),
    ]

    result = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "bigBedToBed failed for "
            + chrom
            + ":"
            + str(query_start)
            + "-"
            + str(query_end)
            + chr(10)
            + result.stderr
        )

    peaks: list[PeakRecord] = []
    for line in result.stdout.splitlines():
        peak = parse_bigbedtobed_line(line)
        if peak is not None:
            peaks.append(peak)

    peaks.sort(key=lambda peak: peak.center)
    return chrom, peaks


def get_upstream_downstream(peaks: list[PeakRecord], center: int) -> tuple[Optional[PeakRecord], Optional[PeakRecord]]:
    if not peaks:
        return None, None

    centers = [peak.center for peak in peaks]
    insert_index = bisect.bisect_left(centers, center)

    upstream = None
    downstream = None

    if insert_index - 1 >= 0:
        upstream = peaks[insert_index - 1]
    if insert_index < len(peaks):
        downstream = peaks[insert_index]

    return upstream, downstream


def add_peak_fields(row: dict, prefix: str, upstream: Optional[PeakRecord], downstream: Optional[PeakRecord]) -> None:
    for label, peak in [("upstream", upstream), ("downstream", downstream)]:
        base = prefix + "_" + label
        if peak is None:
            row[base + "_center"] = "NA"
            row[base + "_score"] = "NA"
            row[base + "_distance_from_region_center"] = "NA"
        else:
            row[base + "_center"] = peak.center
            row[base + "_score"] = peak.score
            row[base + "_distance_from_region_center"] = peak.center - row["region_center"]


def make_flanking_peak_table(
    regions: list[BedRegion],
    nucleosome_bigbed_path: str,
    breakpoint_bigbed_path: str,
    flank_bp: int,
    bigbedtobed_exe: str = "bigBedToBed",
) -> pd.DataFrame:
    rows: list[dict] = []

    for index, region in enumerate(regions, start=1):
        row: dict = {
            "region_index": index,
            "chrom": region.chrom,
            "start": region.start,
            "end": region.end,
            "region_center": region.center,
            "peak_query_flank_bp": flank_bp,
        }

        nuc_chrom, nuc_peaks = query_bigbed_peaks_with_ucsc(
            nucleosome_bigbed_path,
            region.chrom,
            region.center,
            flank_bp,
            bigbedtobed_exe,
        )
        bp_chrom, bp_peaks = query_bigbed_peaks_with_ucsc(
            breakpoint_bigbed_path,
            region.chrom,
            region.center,
            flank_bp,
            bigbedtobed_exe,
        )

        row["nucleosome_query_chrom"] = nuc_chrom
        row["breakpoint_query_chrom"] = bp_chrom
        row["nucleosome_peaks_found_in_query"] = len(nuc_peaks)
        row["breakpoint_peaks_found_in_query"] = len(bp_peaks)

        nuc_up, nuc_down = get_upstream_downstream(nuc_peaks, region.center)
        bp_up, bp_down = get_upstream_downstream(bp_peaks, region.center)

        add_peak_fields(row, "nucleosome", nuc_up, nuc_down)
        add_peak_fields(row, "breakpoint", bp_up, bp_down)
        rows.append(row)

    return pd.DataFrame(rows)


# -----------------------------
# bigWig signal extraction using UCSC bigWigToBedGraph
# -----------------------------

def query_bigwig_bedgraph_with_ucsc(
    bigwig_path: str,
    chrom: str,
    start: int,
    end: int,
    bigwigtobedgraph_exe: str,
) -> list[tuple[int, int, float]]:
    cmd = [
        bigwigtobedgraph_exe,
        bigwig_path,
        "stdout",
        "-chrom=" + chrom,
        "-start=" + str(start),
        "-end=" + str(end),
    ]

    result = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "bigWigToBedGraph failed for "
            + chrom
            + ":"
            + str(start)
            + "-"
            + str(end)
            + chr(10)
            + result.stderr
        )

    intervals: list[tuple[int, int, float]] = []
    for line in result.stdout.splitlines():
        parts = line.rstrip().split("\t")
        if len(parts) < 4:
            continue
        try:
            interval_start = int(parts[1])
            interval_end = int(parts[2])
            value = float(parts[3])
        except ValueError:
            continue
        intervals.append((interval_start, interval_end, value))

    return intervals


def get_bigwig_values_for_region_ucsc(
    bigwig_path: str,
    region: BedRegion,
    missing_value: str,
    bigwigtobedgraph_exe: str,
) -> list[str]:
    region_width = region.end - region.start
    values: list[Optional[float]] = [None] * region_width

    intervals = query_bigwig_bedgraph_with_ucsc(
        bigwig_path=bigwig_path,
        chrom=region.chrom,
        start=region.start,
        end=region.end,
        bigwigtobedgraph_exe=bigwigtobedgraph_exe,
    )

    for interval_start, interval_end, value in intervals:
        fill_start = max(interval_start, region.start)
        fill_end = min(interval_end, region.end)
        for pos in range(fill_start, fill_end):
            values[pos - region.start] = value

    output_values: list[str] = []
    for value in values:
        if value is None:
            output_values.append(missing_value)
        else:
            output_values.append(format(float(value), ".8g"))

    return output_values


def make_signal_matrix(
    regions: list[BedRegion],
    bigwig_path: str,
    missing_value: str,
    bigwigtobedgraph_exe: str = "bigWigToBedGraph",
) -> pd.DataFrame:
    rows: list[list[str | int]] = []
    max_width = 0

    for region in regions:
        values = get_bigwig_values_for_region_ucsc(
            bigwig_path=bigwig_path,
            region=region,
            missing_value=missing_value,
            bigwigtobedgraph_exe=bigwigtobedgraph_exe,
        )
        if len(values) > max_width:
            max_width = len(values)
        rows.append([region.chrom, region.start, region.end] + values)

    columns = ["chrom", "start", "end"]
    for i in range(max_width):
        columns.append("pos_" + str(i))

    padded_rows = []
    for row in rows:
        if len(row) < len(columns):
            row = row + [missing_value] * (len(columns) - len(row))
        padded_rows.append(row)

    return pd.DataFrame(padded_rows, columns=columns)


# -----------------------------
# Path helpers
# -----------------------------

def build_paths(data_dir: str, coverage_bw: str, pns_bw: str, breakpoint_bb: str, nucleosome_bb: str) -> dict[str, str]:
    return {
        "coverage bigWig": os.path.join(data_dir, coverage_bw),
        "PNS bigWig": os.path.join(data_dir, pns_bw),
        "breakpoint peak bigBed": os.path.join(data_dir, breakpoint_bb),
        "nucleosome peak bigBed": os.path.join(data_dir, nucleosome_bb),
    }


def path_report(paths: dict[str, str]) -> str:
    lines = []
    for name, path in paths.items():
        lines.append(name + ": " + path)
    return chr(10).join(lines)


def missing_paths(paths: dict[str, str]) -> list[str]:
    missing = []
    for name, path in paths.items():
        if not os.path.exists(path):
            missing.append(name)
    return missing


def run_extraction(
    regions: list[BedRegion],
    paths: dict[str, str],
    flank_bp: int,
    missing_value: str,
    bigbedtobed_exe: str = "bigBedToBed",
    bigwigtobedgraph_exe: str = "bigWigToBedGraph",
):
    flanking_df = make_flanking_peak_table(
        regions=regions,
        nucleosome_bigbed_path=paths["nucleosome peak bigBed"],
        breakpoint_bigbed_path=paths["breakpoint peak bigBed"],
        flank_bp=flank_bp,
        bigbedtobed_exe=bigbedtobed_exe,
    )
    coverage_df = make_signal_matrix(regions, paths["coverage bigWig"], missing_value, bigwigtobedgraph_exe)
    pns_df = make_signal_matrix(regions, paths["PNS bigWig"], missing_value, bigwigtobedgraph_exe)
    return flanking_df, coverage_df, pns_df


# -----------------------------
# CLI
# -----------------------------

def parse_cli_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract upstream/downstream nucleosome and breakpoint peaks plus per-base bigWig signal.")
    parser.add_argument("--data-dir", default=".")
    parser.add_argument("--coverage-bw", default=DEFAULT_COVERAGE_BW)
    parser.add_argument("--pns-bw", default=DEFAULT_PNS_BW)
    parser.add_argument("--breakpoint-bb", default=DEFAULT_BREAKPOINT_BB)
    parser.add_argument("--nucleosome-bb", default=DEFAULT_NUCLEOSOME_BB)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--peak-flank-bp", type=int, default=1000)
    parser.add_argument("--missing-value", default="nan")
    parser.add_argument("--bigbedtobed", default="bigBedToBed")
    parser.add_argument("--bigwigtobedgraph", default="bigWigToBedGraph")
    parser.add_argument("--out-prefix", default="output")
    return parser.parse_args()


def run_cli() -> None:
    args = parse_cli_args()

    paths = build_paths(
        data_dir=args.data_dir,
        coverage_bw=args.coverage_bw,
        pns_bw=args.pns_bw,
        breakpoint_bb=args.breakpoint_bb,
        nucleosome_bb=args.nucleosome_bb,
    )
    paths["input BED"] = args.bed

    missing = missing_paths(paths)
    if missing:
        raise FileNotFoundError("Missing required file(s): " + ", ".join(missing) + chr(10) + path_report(paths))

    with open(args.bed, "r", encoding="utf-8") as handle:
        bed_text = handle.read()

    regions, skipped = parse_bed_text(bed_text)
    if not regions:
        raise ValueError("No valid BED regions found")

    flanking_df, coverage_df, pns_df = run_extraction(
        regions=regions,
        paths=paths,
        flank_bp=args.peak_flank_bp,
        missing_value=args.missing_value,
        bigbedtobed_exe=args.bigbedtobed,
        bigwigtobedgraph_exe=args.bigwigtobedgraph,
    )

    flanking_path = args.out_prefix + "_flanking.tsv"
    coverage_path = args.out_prefix + "_coverage.tsv"
    pns_path = args.out_prefix + "_pns.tsv"

    flanking_df.to_csv(flanking_path, sep="	", index=False)
    coverage_df.to_csv(coverage_path, sep="	", index=False)
    pns_df.to_csv(pns_path, sep="	", index=False)

    if skipped:
        skipped_path = args.out_prefix + "_skipped_lines.txt"
        with open(skipped_path, "w", encoding="utf-8") as handle:
            handle.write(chr(10).join(skipped))
            handle.write(chr(10))
        print("Wrote " + skipped_path)

    print("Wrote " + flanking_path)
    print("Wrote " + coverage_path)
    print("Wrote " + pns_path)


# -----------------------------
# Streamlit
# -----------------------------

def dataframe_to_tsv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(sep="	", index=False).encode("utf-8")


def make_output_zip(outputs: dict[str, pd.DataFrame], skipped: list[str]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for filename, df in outputs.items():
            zf.writestr(filename, df.to_csv(sep="	", index=False))
        if skipped:
            zf.writestr("skipped_lines.txt", chr(10).join(skipped) + chr(10))
    buffer.seek(0)
    return buffer.getvalue()


def default_data_dir() -> str:
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.getcwd()


def run_streamlit() -> None:
    if st is None:
        raise ImportError("Streamlit is not installed")

    st.set_page_config(page_title="PNS region extractor", layout="wide")
    st.title("PNS region extractor")
    st.write("Extract upstream/downstream nucleosome and breakpoint peak centres, peak scores, coverage, and PNS signal.")

    with st.sidebar:
        st.header("Data files")
        data_dir = st.text_input("Data directory", value=default_data_dir())
        coverage_bw = st.text_input("Coverage bigWig", value=DEFAULT_COVERAGE_BW)
        pns_bw = st.text_input("PNS bigWig", value=DEFAULT_PNS_BW)
        breakpoint_bb = st.text_input("Breakpoint peak bigBed", value=DEFAULT_BREAKPOINT_BB)
        nucleosome_bb = st.text_input("Nucleosome peak bigBed", value=DEFAULT_NUCLEOSOME_BB)

        paths = build_paths(data_dir, coverage_bw, pns_bw, breakpoint_bb, nucleosome_bb)
        st.caption("Current file paths")
        st.code(path_report(paths), language="text")

        st.header("Options")
        peak_flank_bp = st.number_input("Peak query flank size (bp)", min_value=100, max_value=100000, value=1000, step=100)
        missing_value = st.text_input("Missing bigWig value", value="nan")
        bigbedtobed_exe = st.text_input("bigBedToBed executable", value="bigBedToBed")
        bigwigtobedgr

    st.header("Input BED regions")
    uploaded_bed = st.file_uploader("Upload BED file", type=["bed", "txt", "tsv"])
    pasted_bed = st.text_area("Or paste BED regions", height=180, placeholder="chr20	100000	100200")

    run_button = st.button("Run extraction", type="primary")
    if not run_button:
        return

    missing = missing_paths(paths)
    if missing:
        st.error("Missing required file(s): " + ", ".join(missing))
        st.code(path_report(paths), language="text")
        return

    if uploaded_bed is not None:
        bed_text = uploaded_bed.getvalue().decode("utf-8")
    else:
        bed_text = pasted_bed

    regions, skipped = parse_bed_text(bed_text)
    if not regions:
        st.error("No valid BED regions found")
        if skipped:
            st.text(chr(10).join(skipped[:50]))
        return

    try:
        with st.spinner("Extracting..."):
            flanking_df, coverage_df, pns_df = run_extraction(regions, paths, int(peak_flank_bp), missing_value, bigbedtobed_exe)
    except Exception as exc:
        st.exception(exc)
        return

    st.success("Processed " + str(len(regions)) + " valid BED region(s)")

    if skipped:
        with st.expander("Skipped " + str(len(skipped)) + " input line(s)"):
            st.text(chr(10).join(skipped[:200]))

    st.subheader("Flanking peaks")
    st.dataframe(flanking_df, use_container_width=True)

    st.subheader("Coverage matrix preview")
    st.dataframe(coverage_df.head(20), use_container_width=True)

    st.subheader("PNS matrix preview")
    st.dataframe(pns_df.head(20), use_container_width=True)

    outputs = {
        "flanking_peaks.tsv": flanking_df,
        "coverage_matrix.tsv": coverage_df,
        "pns_matrix.tsv": pns_df,
    }

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.download_button("Download all outputs", data=make_output_zip(outputs, skipped), file_name="pns_region_extractor_outputs.zip", mime="application/zip")
    with col2:
        st.download_button("flanking_peaks.tsv", data=dataframe_to_tsv_bytes(flanking_df), file_name="flanking_peaks.tsv")
    with col3:
        st.download_button("coverage_matrix.tsv", data=dataframe_to_tsv_bytes(coverage_df), file_name="coverage_matrix.tsv")
    with col4:
        st.download_button("pns_matrix.tsv", data=dataframe_to_tsv_bytes(pns_df), file_name="pns_matrix.tsv")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_cli()
    else:
        run_streamlit()
