"""Reproduce, generate, or preview the B packing with diameters up to 3.0 nm."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import time
from contextlib import ExitStack
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE_CSV = HERE / "candidate_14.csv"
CONSTRUCTOR_FILE = HERE / "constructor_gen6.py"
REFERENCE_SEED = 17
BOX_W, BOX_H, GAP = 9.7, 11.1, 0.34
DIAMETERS = (1.0, 1.5, 2.0, 2.5, 3.0)
CSV_FIELDS = ("x_nm", "y_nm", "diameter_nm")

# Use the same geometry checks and SVG renderer as the root A command.
sys.path.insert(0, str(HERE.parent))
from run import draw_svg, inspect_layout


def constructor_hash():
    return hashlib.sha256(CONSTRUCTOR_FILE.read_bytes()).hexdigest()


def summarize(layout, seed):
    geometry = inspect_layout(layout, BOX_W, BOX_H, GAP, DIAMETERS)
    counts = {
        str(row["diameter_nm"]): row["count"]
        for row in geometry["diameter_counts"] if row["count"]
    }
    return {
        "seed": seed,
        "n": len(layout),
        "S_percent": geometry["nominal_area_fraction_percent"],
        "pore_area_nm2": geometry["circle_area_nm2"],
        "counts": counts,
        "min_gap_nm": geometry["minimum_periodic_gap_nm"],
        "layout": layout,
        "geometry": geometry,
    }


def run_one(seed):
    if not CONSTRUCTOR_FILE.is_file():
        raise FileNotFoundError("Missing constructor_gen6.py beside run_b.py.")
    import constructor_gen6 as constructor

    raw = constructor.run_experiment(
        random_seed=seed, box_w=BOX_W, box_h=BOX_H,
        diameters=list(DIAMETERS), safety_gap=GAP,
    )
    return summarize([(float(x) % BOX_W, float(y) % BOX_H, float(d)) for x, y, d in raw], seed)


def read_layout(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        for fields in (CSV_FIELDS, ("x", "y", "diameter")):
            if columns.issuperset(fields):
                return [tuple(float(row[field]) for field in fields) for row in reader]
    raise ValueError("Expected CSV columns x_nm,y_nm,diameter_nm or x,y,diameter (nm).")


def csv_text(layout):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(CSV_FIELDS)
    for x, y, diameter in layout:
        writer.writerow((f"{x:.10f}", f"{y:.10f}", f"{diameter:g}"))
    return stream.getvalue()


def check_new_paths(paths):
    for path in paths:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}. Choose a new --out directory.")


def write_files(outputs):
    created = []
    try:
        with ExitStack() as stack:
            handles = {}
            for path in outputs:
                path.parent.mkdir(parents=True, exist_ok=True)
                handles[path] = stack.enter_context(path.open("x", encoding="utf-8", newline=""))
                created.append(path)
            for path, content in outputs.items():
                handles[path].write(content)
    except OSError:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def output_paths(out_dir, name):
    return {ext: out_dir / f"{name}.{ext}" for ext in ("csv", "svg", "json")}


def write_outputs(record, paths, source):
    layout = record["layout"]
    payload = {key: value for key, value in record.items() if key != "layout"}
    payload.update(source=source, layout_nm=layout)
    write_files({
        paths["csv"]: csv_text(layout),
        paths["svg"]: draw_svg(layout, record["geometry"]),
        paths["json"]: json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
    })
    print(f"Circles: {record['n']} | S: {record['S_percent']:.9f}% | "
          f"Min periodic gap: {record['min_gap_nm']:.9f} nm")
    for ext, path in paths.items():
        print(f"{ext.upper()}: {path}")


def periodic_axis_delta(delta, length):
    return abs((delta + length / 2) % length - length / 2)


def verify():
    reference = summarize(read_layout(REFERENCE_CSV), REFERENCE_SEED)
    print(f"Reproducing seed {REFERENCE_SEED} ...", flush=True)
    started = time.perf_counter()
    generated = run_one(REFERENCE_SEED)
    tol = 1e-6  # nm
    unmatched = list(reference["layout"])
    for x, y, diameter in generated["layout"]:
        match = next((i for i, (u, v, d) in enumerate(unmatched)
                      if abs(diameter - d) <= tol
                      and periodic_axis_delta(x - u, BOX_W) <= tol
                      and periodic_axis_delta(y - v, BOX_H) <= tol), None)
        if match is None:
            raise ValueError("Seed 17 does not match the reference coordinates within 1e-6 nm.")
        unmatched.pop(match)
    if unmatched:
        raise ValueError("Seed 17 did not reproduce all reference circles.")
    print(f"PASS: {generated['n']} circles, S={generated['S_percent']:.9f}%, "
          f"min gap={generated['min_gap_nm']:.9f} nm. "
          f"Coordinates match within {tol:g} nm per axis ({time.perf_counter() - started:.1f} s).")
    return generated


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples from the repository root:\n"
            "  python -B b_dmax3/run_b.py --verify\n"
            "  python -B b_dmax3/run_b.py --seeds 17 --out runs/b_dmax3\n"
            "  python -B b_dmax3/run_b.py --input-csv b_dmax3/candidate_14.csv\n\n"
            "Fixed geometry: 9.7 x 11.1 nm, periodic in x and y, 0.34 nm minimum\n"
            "surface gap; allowed nominal diameters: 1, 1.5, 2, 2.5, and 3 nm."
        ),
    )
    parser.add_argument("--verify", action="store_true",
                        help="Regenerate seed 17 and compare it with candidate_14.csv.")
    parser.add_argument("--seeds", type=int, nargs="+", default=[],
                        help="Generate one CSV, SVG, and JSON bundle per seed.")
    parser.add_argument("--input-csv", type=Path,
                        help="Validate and plot an existing layout without regenerating it.")
    parser.add_argument("--out", type=Path, default=Path("runs/b_dmax3"),
                        help="Output directory (default: runs/b_dmax3). Existing files are protected.")
    args = parser.parse_args(argv)
    if not (args.verify or args.seeds or args.input_csv):
        parser.print_help()
        return 0
    if args.input_csv and (args.seeds or args.verify):
        parser.error("--input-csv is a separate mode from --verify and --seeds.")
    if len(set(args.seeds)) != len(args.seeds):
        parser.error("List each seed only once.")

    try:
        if args.input_csv:
            paths = output_paths(args.out, args.input_csv.stem)
            check_new_paths(paths.values())
            record = summarize(read_layout(args.input_csv), None)
            source = {"mode": "existing_csv", "filename": args.input_csv.name,
                      "sha256": hashlib.sha256(args.input_csv.read_bytes()).hexdigest()}
            write_outputs(record, paths, source)
            return 0

        seed_outputs = {seed: output_paths(args.out, f"seed_{seed}") for seed in args.seeds}
        for paths in seed_outputs.values():
            check_new_paths(paths.values())
        digest = constructor_hash()
        verified = verify() if args.verify else None
        for seed, paths in seed_outputs.items():
            if verified is not None and seed == REFERENCE_SEED:
                record = verified
            else:
                print(f"Generating seed {seed} ...", flush=True)
                record = run_one(seed)
            write_outputs(record, paths, {"mode": "generated", "constructor_sha256": digest})
    except (OSError, ValueError, TypeError, csv.Error) as exc:
        parser.exit(2, f"Error: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

