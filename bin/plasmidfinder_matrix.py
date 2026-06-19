#!/usr/bin/env python3
"""Build a samples x plasmids matrix from per-sample PlasmidFinder summaries.

Reads every *.plasmidfinder_summary.tsv in the working directory and writes:
  <prefix>.plasmidfinder_matrix.tsv  samples x plasmids, cells "identity% / coverage%"
  <prefix>.plasmidfinder_long.tsv    tidy one-row-per-hit table with computed coverage

Cells hold the best hit (highest identity) of each plasmid in each sample;
coverage is derived from PlasmidFinder's "Query / Template length" column.
Plasmids absent from a sample are left blank.
"""

import argparse
import glob
import re

import pandas as pd


def find_col(columns, want):
    for c in columns:
        if want in c.lower():
            return c
    return None


def coverage(value):
    nums = re.findall(r"[\d.]+", str(value))
    if len(nums) >= 2 and float(nums[1]) > 0:
        return round(100 * float(nums[0]) / float(nums[1]), 1)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="all_samples", help="Output file prefix")
    parser.add_argument(
        "--pattern",
        default="*.plasmidfinder_summary.tsv",
        help="Glob for input summary files",
    )
    args = parser.parse_args()

    frames = []
    for path in sorted(glob.glob(args.pattern)):
        try:
            frames.append(pd.read_csv(path, sep="\t", dtype=str))
        except pd.errors.EmptyDataError:
            continue

    long_cols = ["sample", "mag", "plasmid", "identity", "coverage", "accession", "note"]

    if not frames or all(f.empty for f in frames):
        # No plasmid hits anywhere: emit empty, well-formed outputs.
        pd.DataFrame(columns=["sample"]).to_csv(
            f"{args.prefix}.plasmidfinder_matrix.tsv", sep="\t", index=False
        )
        pd.DataFrame(columns=long_cols).to_csv(
            f"{args.prefix}.plasmidfinder_long.tsv", sep="\t", index=False
        )
        return

    df = pd.concat(frames, ignore_index=True)

    c_sample = find_col(df.columns, "sample")
    c_mag = find_col(df.columns, "mag")
    c_plasmid = find_col(df.columns, "plasmid")
    c_ident = find_col(df.columns, "identity")
    c_len = find_col(df.columns, "template length")
    c_acc = find_col(df.columns, "accession")
    c_note = find_col(df.columns, "note")

    df["identity"] = pd.to_numeric(df[c_ident], errors="coerce")
    df["coverage"] = df[c_len].apply(coverage) if c_len else None

    # Tidy long table (every hit, with computed coverage).
    long = pd.DataFrame(
        {
            "sample": df[c_sample],
            "mag": df[c_mag] if c_mag else "",
            "plasmid": df[c_plasmid],
            "identity": df["identity"],
            "coverage": df["coverage"],
            "accession": df[c_acc] if c_acc else "",
            "note": df[c_note] if c_note else "",
        }
    ).sort_values(["sample", "plasmid", "identity"], ascending=[True, True, False])
    long.to_csv(f"{args.prefix}.plasmidfinder_long.tsv", sep="\t", index=False)

    # Matrix: best hit (max identity) per sample x plasmid.
    best = df.sort_values("identity", ascending=False).drop_duplicates([c_sample, c_plasmid])
    best["cell"] = best.apply(
        lambda r: f"{r['identity']:.1f} / {r['coverage']:.1f}"
        if pd.notna(r["identity"]) and pd.notna(r["coverage"])
        else (f"{r['identity']:.1f}" if pd.notna(r["identity"]) else ""),
        axis=1,
    )
    matrix = best.pivot(index=c_sample, columns=c_plasmid, values="cell").fillna("")
    matrix.index.name = "sample"
    matrix.to_csv(f"{args.prefix}.plasmidfinder_matrix.tsv", sep="\t")


if __name__ == "__main__":
    main()
