"""Run the B constructor in parallel over a range of integer seeds."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from run_b import BOX_H, BOX_W, check_new_paths, csv_text, periodic_axis_delta, run_one, write_files

TARGET_S = 64.00918439502748


def signature(layout):
    """Heuristic grouping by rounded periodic distances and their two diameters."""
    pairs = []
    for i, (x, y, d) in enumerate(layout):
        for u, v, other_d in layout[:i]:
            distance = math.hypot(periodic_axis_delta(x - u, BOX_W),
                                  periodic_axis_delta(y - v, BOX_H))
            pairs.append((min(d, other_d), max(d, other_d), round(distance, 6)))
    return tuple(sorted(d for _, _, d in layout)), tuple(sorted(pairs))


def count_label(record):
    return " ".join(f"{d}x{count}" for d, count in record["counts"].items())


def save_sweep(records, out_dir):
    records.sort(key=lambda record: (-record["S_percent"], record["seed"]))
    groups = {}
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(["seed", "n_circles", "S_percent", "min_gap_nm", "counts", "equivalence_group"])
    json_lines = []
    outputs = {}
    hits = 0
    for record in records:
        key = signature(record["layout"])
        group = groups.setdefault(key, len(groups) + 1)
        writer.writerow([record["seed"], record["n"], f"{record['S_percent']:.9f}",
                         f"{record['min_gap_nm']:.9f}", count_label(record), group])
        json_lines.append(json.dumps({
            "seed": record["seed"], "n": record["n"], "S_percent": record["S_percent"],
            "counts": record["counts"], "min_gap_nm": record["min_gap_nm"],
            "equivalence_group": group, "layout_nm": record["layout"],
        }, allow_nan=False))
        if abs(record["S_percent"] - TARGET_S) < 1e-6:
            outputs[out_dir / f"hit_{record['seed']}.csv"] = csv_text(record["layout"])
            hits += 1
    outputs[out_dir / "summary.csv"] = stream.getvalue()
    outputs[out_dir / "layouts.jsonl"] = "\n".join(json_lines) + "\n"
    write_files(outputs)
    return hits


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example from the repository root:\n"
            "  python -B b_dmax3/sweep.py --start 200 --count 4 --workers 2 "
            "--out runs/b_dmax3/sweep\n\n"
            "Outputs: summary.csv, layouts.jsonl, and hit_<seed>.csv when S matches\n"
            "the reference area fraction. A match in S can have different coordinates."
        ),
    )
    parser.add_argument("--start", type=int, default=200, help="First seed (default: 200).")
    parser.add_argument("--count", type=int, default=24, help="Number of seeds (default: 24).")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 4),
                        help="Parallel processes (default: up to 8, based on available CPUs).")
    parser.add_argument("--out", type=Path, default=Path("runs/b_dmax3/sweep"),
                        help="Output directory (default: runs/b_dmax3/sweep).")
    args = parser.parse_args(argv)
    if args.count < 1 or args.workers < 1:
        parser.error("--count and --workers must be positive.")
    seeds = range(args.start, args.start + args.count)
    try:
        check_new_paths([args.out / "summary.csv", args.out / "layouts.jsonl"])
        check_new_paths(args.out / f"hit_{seed}.csv" for seed in seeds)
        print(f"Seeds {seeds[0]}..{seeds[-1]} | workers: {args.workers}", flush=True)
        started = time.perf_counter()
        records = []
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_one, seed): seed for seed in seeds}
            for done, future in enumerate(as_completed(futures), 1):
                seed = futures[future]
                try:
                    record = future.result()
                except Exception as exc:
                    print(f"[{done}/{len(seeds)}] Seed {seed} failed: {exc}", file=sys.stderr)
                    continue
                records.append(record)
                print(f"[{done}/{len(seeds)}] Seed {seed} | n={record['n']} | "
                      f"S={record['S_percent']:.6f}% | {count_label(record)}", flush=True)
        if not records:
            parser.exit(1, "Error: no seed completed successfully.\n")
        hits = save_sweep(records, args.out)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(f"Completed {len(records)}/{len(seeds)} seeds in {time.perf_counter() - started:.1f} s. "
          f"Reference S hits: {hits}.")
    print(f"Summary: {args.out / 'summary.csv'}")
    return int(len(records) != len(seeds))


if __name__ == "__main__":
    sys.exit(main())

