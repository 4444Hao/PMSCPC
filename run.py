"""Standalone CLI and SVG preview for the unchanged Shinka-evolved constructor.

Python standard library only. All geometry lengths are in nm.
"""
import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import sys
import time


TOL = 1e-7  # Accommodates rounding of the archived 8-decimal CSV.
COLORS = {
    1.0: "#46788E", 1.5: "#408C85", 2.0: "#6B8D50",
    2.5: "#A49B33", 3.0: "#D58C31", 3.5: "#A373AA",
}
FALLBACK = ["#46788E", "#408C85", "#6B8D50", "#A49B33", "#D58C31", "#A373AA"]


def read_layout(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"x_nm", "y_nm", "diameter_nm"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CSV must contain x_nm,y_nm,diameter_nm (nm units).")
        return [(float(row["x_nm"]), float(row["y_nm"]), float(row["diameter_nm"]))
                for row in reader]


def inspect_layout(layout, width, height, gap, diameters):
    if not layout:
        raise ValueError("Empty layout: no figure exported.")
    counts = {d: 0 for d in sorted(diameters)}
    pair_gap = math.inf
    pair = None
    self_gap = math.inf
    self_circle = None
    for i, (x, y, diameter) in enumerate(layout):
        if not all(math.isfinite(v) for v in (x, y, diameter)):
            raise ValueError(f"Circle {i + 1}: non-finite coordinate/diameter.")
        if not (0 <= x < width and 0 <= y < height):
            raise ValueError(f"Circle {i + 1}: center outside [0,W) x [0,H).")
        matches = [d for d in diameters if abs(d - diameter) <= TOL]
        if diameter <= 0 or not matches:
            raise ValueError(f"Circle {i + 1}: diameter not in --diameters.")
        counts[matches[0]] += 1
        candidate = min(width, height) - diameter
        if candidate < self_gap:
            self_gap, self_circle = candidate, i + 1
        for j in range(i):
            u, v, other_d = layout[j]
            dx, dy = abs(x - u), abs(y - v)
            distance = math.hypot(min(dx, width - dx), min(dy, height - dy))
            candidate = distance - (diameter + other_d) / 2
            if candidate < pair_gap:
                pair_gap, pair = candidate, [j + 1, i + 1]
    minimum = min(pair_gap, self_gap)
    if minimum < gap - TOL:
        raise ValueError(f"Periodic gap {minimum:.9g} nm < requested {gap:g} nm.")
    area = math.fsum(math.pi * (d / 2) ** 2 for _, _, d in layout)
    return {
        "circle_count": len(layout), "box_nm": [width, height],
        "allowed_diameters_nm": sorted(diameters),
        "used_diameter_kinds": sum(n > 0 for n in counts.values()),
        "diameter_counts": [{"diameter_nm": d, "count": n} for d, n in counts.items()],
        "circle_area_nm2": area, "box_area_nm2": width * height,
        "nominal_area_fraction_percent": 100 * area / (width * height),
        "required_gap_nm": gap,
        "minimum_intercircle_gap_nm": pair_gap if pair else None,
        "closest_intercircle_pair_ids": pair,
        "minimum_self_image_gap_nm": self_gap,
        "closest_self_image_circle_id": self_circle,
        "minimum_periodic_gap_nm": minimum,
        "geometry_tolerance_nm": TOL, "geometry_check": "PASS_WITH_STATED_TOLERANCE",
    }


def draw_svg(layout, report):
    width, height = report["box_nm"]
    scale = min(550 / width, 620 / height)
    left, top = 74, 88
    pw, ph = width * scale, height * scale
    sidebar = left + pw + 52
    rows = report["diameter_counts"]
    canvas_w = max(960, sidebar + 340)
    canvas_h = max(570, top + ph + 135, 480 + 28 * len(rows))
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{canvas_w:.0f}" '
           f'height="{canvas_h:.0f}" viewBox="0 0 {canvas_w:.0f} {canvas_h:.0f}" '
           'role="img" aria-labelledby="title description">',
           '<title id="title">Periodic multi-size circle packing</title>',
           '<desc id="description">True-scale periodic cell with nominal diameters, '
           'counts, area fraction and minimum periodic surface gap. Geometry only.</desc>',
           '<rect width="100%" height="100%" fill="white"/>',
           '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#253442}'
           '.muted{fill:#627382}.small{font-size:13px}.normal{font-size:16px}</style>',
           f'<defs><clipPath id="cell"><rect x="{left}" y="{top}" '
           f'width="{pw:.6f}" height="{ph:.6f}"/></clipPath></defs>']

    def text(x, y, value, size=16, anchor="start", extra=""):
        svg.append(f'<text x="{x:.4f}" y="{y:.4f}" font-size="{size}" '
                   f'text-anchor="{anchor}" {extra}>{html.escape(str(value))}</text>')

    text(left, 35, "Periodic multi-size circle packing", 24, extra='font-weight="bold"')
    text(left, 61, "Nominal geometry | lengths in nm | not a global-optimum certificate", 13)
    svg.append(f'<rect x="{left}" y="{top}" width="{pw:.6f}" height="{ph:.6f}" fill="#FAFCFD"/>')
    palette = {row["diameter_nm"]: COLORS.get(row["diameter_nm"], FALLBACK[i % len(FALLBACK)])
               for i, row in enumerate(rows)}
    svg.append('<g clip-path="url(#cell)">')
    for x, y, diameter in layout:
        key = min(palette, key=lambda d: abs(d - diameter))
        color = palette[key]
        for ix in (-1, 0, 1):
            for iy in (-1, 0, 1):
                cx, cy = x + ix * width, y + iy * height
                radius = diameter / 2
                if cx + radius <= 0 or cx - radius >= width or cy + radius <= 0 or cy - radius >= height:
                    continue
                svg.append(f'<circle cx="{left + cx * scale:.6f}" '
                           f'cy="{top + (height - cy) * scale:.6f}" r="{radius * scale:.6f}" '
                           f'fill="{color}" fill-opacity="0.18" stroke="{color}" stroke-width="1.7"/>')
    # Diameter labels refer to the unshifted primary centers; fragments are not extra circles.
    for x, y, diameter in layout:
        if diameter * scale >= 20:
            text(left + x * scale, top + (height - y) * scale + 4, f"{diameter:g}", 12, "middle")
    svg.append('</g>')
    svg.append(f'<rect x="{left}" y="{top}" width="{pw:.6f}" height="{ph:.6f}" '
               'fill="none" stroke="#253442" stroke-width="1.5"/>')
    for i in range(5):
        vx, vy = width * i / 4, height * i / 4
        px, py = left + vx * scale, top + (height - vy) * scale
        svg.append(f'<path d="M {px:.6f} {top+ph:.6f} v 5 M {left-5} {py:.6f} h 5" stroke="#253442"/>')
        text(px, top + ph + 24, f"{vx:.3g}", 13, "middle")
        text(left - 12, py + 4, f"{vy:.3g}", 13, "end")
    text(left + pw / 2, top + ph + 50, "x (nm)", 16, "middle")
    svg.append(f'<text transform="translate(22,{top+ph/2:.4f}) rotate(-90)" '
               'text-anchor="middle" font-size="16">y (nm)</text>')
    y = top + 8
    text(sidebar, y, "Layout summary", 20, extra='font-weight="bold"')
    lines = [f"Circles: {report['circle_count']}",
             f"Diameter kinds: {report['used_diameter_kinds']} used / {len(rows)} allowed",
             f"Cell: {width:g} x {height:g} nm",
             f"Area fraction S: {report['nominal_area_fraction_percent']:.4f}%",
             f"Total circle area: {report['circle_area_nm2']:.4f} nm^2",
             f"Required gap: {report['required_gap_nm']:g} nm",
             f"Min periodic gap: {report['minimum_periodic_gap_nm']:.6f} nm",
             "Geometry: PASS (tolerance 1e-7 nm)"]
    for j, line in enumerate(lines):
        text(sidebar, y + 34 + j * 28, line, 15)
    legend_y = y + 34 + len(lines) * 28 + 16
    text(sidebar, legend_y, "Diameter (nm) / count", 17, extra='font-weight="bold"')
    for j, row in enumerate(rows):
        ly = legend_y + 29 + j * 28
        color = palette[row["diameter_nm"]]
        svg.append(f'<circle cx="{sidebar+9:.4f}" cy="{ly-5:.4f}" r="7" '
                   f'fill="{color}" fill-opacity="0.3" stroke="{color}" stroke-width="1.5"/>')
        text(sidebar + 28, ly, f"{row['diameter_nm']:g} nm   x {row['count']}", 15)
    footer_y = max(top + ph + 83, legend_y + 40 + len(rows) * 28)
    text(left, footer_y, "Edge fragments are periodic images, not extra circles. Labels inside circles are diameters (nm).", 13)
    text(left, footer_y + 22, "S uses nominal disk areas; gaps include self-images. No atomic structure or MD performance is implied.", 13)
    svg.append('</svg>')
    return '\n'.join(svg) + '\n'


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Full-parameter generation example (PowerShell, from D:\\PMSCPC):\n"
            "  python -B run.py --seed 17 --box-w 9.7 --box-h 11.1 "
            "--diameters 1 1.5 2 2.5 3 3.5 --gap 0.34 --output runs/example_seed17\n"
            "\n"
            "All lengths are in nm. Existing outputs are never overwritten; "
            "choose a new --output prefix for another run."
        ),
    )
    parser.add_argument("--seed", type=int, default=17, help="Generator seed (default: 17).")
    parser.add_argument("--box-w", type=float, default=9.7, help="Periodic x length in nm.")
    parser.add_argument("--box-h", type=float, default=11.1, help="Periodic y length in nm.")
    parser.add_argument("--diameters", type=float, nargs="+", default=[1, 1.5, 2, 2.5, 3, 3.5],
                        help="Allowed nominal diameters in nm; classes need not all be used.")
    parser.add_argument("--gap", "--safety-gap", dest="gap", type=float, default=0.34,
                        help="Minimum nominal surface gap in nm.")
    parser.add_argument("--input-csv", type=Path, help="Render an existing CSV instead of searching.")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output prefix, e.g. runs/A_seed17 (not a filename extension).")
    return parser


def main_cli(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if not all(math.isfinite(v) and v > 0 for v in (args.box_w, args.box_h)):
            raise ValueError("Box lengths must be finite and positive.")
        if not math.isfinite(args.gap) or args.gap < 0:
            raise ValueError("Gap must be finite and nonnegative.")
        if any(not math.isfinite(d) or d <= 0 for d in args.diameters):
            raise ValueError("Diameters must be finite and positive.")
        if len(set(args.diameters)) != len(args.diameters):
            raise ValueError("Do not repeat diameter classes.")
        if any(d + args.gap > min(args.box_w, args.box_h) + TOL for d in args.diameters):
            raise ValueError("An allowed diameter + gap exceeds a cell length (self-image conflict).")
        paths = {ext: Path(str(args.output) + '.' + ext) for ext in ('csv', 'svg', 'json')}
        for path in paths.values():
            if path.exists():
                raise FileExistsError(f"Refusing to overwrite {path}; choose a new --output prefix.")
        started = time.perf_counter()
        if args.input_csv is None:
            import main as constructor
            layout = constructor.run_experiment(random_seed=args.seed, box_w=args.box_w,
                                                box_h=args.box_h, diameters=args.diameters,
                                                safety_gap=args.gap)
            source = {"mode": "generated", "seed": args.seed,
                      "constructor_sha256": hashlib.sha256(Path(constructor.__file__).read_bytes()).hexdigest()}
        else:
            layout = read_layout(args.input_csv)
            source = {"mode": "existing_csv", "input_csv": str(args.input_csv),
                      "input_sha256": hashlib.sha256(args.input_csv.read_bytes()).hexdigest(),
                      "seed": None}
        report = inspect_layout(layout, args.box_w, args.box_h, args.gap, args.diameters)
        report.update(source)
        report["elapsed_seconds"] = time.perf_counter() - started
        svg = draw_svg(layout, report)
        for path in paths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        # Reserve all files exclusively before writing any; never overwrite existing results.
        handles = {}
        try:
            for ext, path in paths.items():
                handles[ext] = path.open('x', encoding='utf-8', newline='')
        except BaseException:
            for ext, handle in handles.items():
                handle.close()
                paths[ext].unlink()  # Only newly created, still-empty reservations.
            raise
        try:
            writer = csv.writer(handles['csv'])
            writer.writerow(['x_nm', 'y_nm', 'diameter_nm'])
            writer.writerows(layout)
            handles['svg'].write(svg)
            json.dump(report, handles['json'], ensure_ascii=False, indent=2, allow_nan=False)
            handles['json'].write('\n')
        finally:
            for handle in handles.values():
                handle.close()
        print(f"Circles: {len(layout)} | S: {report['nominal_area_fraction_percent']:.6f}%")
        print('Diameter counts (nm): ' + ', '.join(f"{r['diameter_nm']:g} x {r['count']}" for r in report['diameter_counts']))
        print(f"Min periodic surface gap: {report['minimum_periodic_gap_nm']:.9f} nm | PASS (tol {TOL:g} nm)")
        for ext, path in paths.items():
            print(f"{ext.upper()}: {path}")
        return 0
    except (ValueError, OSError, csv.Error, OverflowError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == '__main__':
    sys.exit(main_cli())
