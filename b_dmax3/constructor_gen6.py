"""B task: nominal diameters <= 3.0 nm in the original periodic box.

Warm-start constructor: saved A best program, generation 26; source snapshots and fingerprints are in sources/.
Maximize nominal circular area fraction, not MD P/R. Fixed composition gives
constant S; positions determine feasibility and insertion opportunities.
The 0.34 nm gap is a design constraint, not proof of water-tightness.
Only the EVOLVE-BLOCK is mutable. The independent evaluator enforces bounds.
"""
import math
import random

# ==================== 任务常量（EVOLVE-BLOCK 之外，固定，LLM 不可变异）====================
DEFAULT_BOX_W = 9.7           # nm  (= 97 Å)
DEFAULT_BOX_H = 11.1          # nm  (= 111 Å)
DEFAULT_DIAMETERS = [1.0, 1.5, 2.0, 2.5, 3.0]
DEFAULT_SAFETY_GAP = 0.34     # nm

# EVOLVE-BLOCK-START
# ---------------------------------------------------------------------------
# Geometry helpers (periodic)
# ---------------------------------------------------------------------------

def _pdist2(x1, y1, x2, y2, W, H):
    dx = abs(x1 - x2)
    if dx > W - dx:
        dx = W - dx
    dy = abs(y1 - y2)
    if dy > H - dy:
        dy = H - dy
    return dx * dx + dy * dy


def _pdist(x1, y1, x2, y2, W, H):
    return math.sqrt(_pdist2(x1, y1, x2, y2, W, H))


def _area(layout):
    s = 0.0
    for (_x, _y, d) in layout:
        s += d * d
    return 0.25 * math.pi * s


def _fits(x, y, d, layout, W, H, gap):
    r = d * 0.5 + gap
    for (px, py, pd) in layout:
        need = r + pd * 0.5
        if _pdist2(x, y, px, py, W, H) < need * need - 1e-12:
            return False
    return True


# ---------------------------------------------------------------------------
# Maximal-empty-circle (MEC) machinery
# ---------------------------------------------------------------------------

def _slack(x, y, layout, W, H):
    if not layout:
        return float("inf")
    best = float("inf")
    for (px, py, pd) in layout:
        dx = abs(x - px)
        if dx > W - dx:
            dx = W - dx
        dy = abs(y - py)
        if dy > H - dy:
            dy = H - dy
        s = math.sqrt(dx * dx + dy * dy) - pd * 0.5
        if s < best:
            best = s
    return best


def _push_to_local_max(x, y, layout, W, H, iters=12, step0=None):
    if not layout:
        return x, y
    if step0 is None:
        step0 = 0.22 * min(W, H)
    step = step0
    cur_s = _slack(x, y, layout, W, H)
    for _it in range(iters):
        improved = False
        best_dx, best_dy, best_s = 0.0, 0.0, cur_s
        for k in range(8):
            ang = k * math.pi / 4.0
            nx = (x + step * math.cos(ang)) % W
            ny = (y + step * math.sin(ang)) % H
            s = _slack(nx, ny, layout, W, H)
            if s > best_s + 1e-12:
                best_s = s
                best_dx = nx - x
                best_dy = ny - y
                improved = True
        if improved:
            x = (x + best_dx) % W
            y = (y + best_dy) % H
            cur_s = best_s
        else:
            step *= 0.5
            if step < 5e-4:
                break
    return x, y


def _mec_sites(layout, W, H, rng, n_seeds, extra_corner=True):
    sites = []
    if not layout:
        nx = max(2, int(round(math.sqrt(n_seeds))))
        for i in range(nx):
            for j in range(nx):
                sites.append(((i + 0.5) * W / nx, (j + 0.5) * H / nx))
        return sites

    n = len(layout)
    pair_count = 0
    # (a) midpoints between nearby pairs
    for i in range(n):
        xi, yi, _di = layout[i]
        for j in range(i + 1, n):
            xj, yj, _dj = layout[j]
            dx = xi - xj
            if dx > W * 0.5:
                dx -= W
            elif dx < -W * 0.5:
                dx += W
            dy = yi - yj
            if dy > H * 0.5:
                dy -= H
            elif dy < -H * 0.5:
                dy += H
            if dx * dx + dy * dy < 32.0:
                sites.append(((xj + dx * 0.5) % W, (yj + dy * 0.5) % H))
                pair_count += 1
                if pair_count > n_seeds * 4:
                    break
        if pair_count > n_seeds * 4:
            break

    # (b) jittered circle-center seeds
    for (px, py, _pd) in layout:
        for _ in range(3):
            ox = rng.uniform(-0.28, 0.28) * W
            oy = rng.uniform(-0.28, 0.28) * H
            sites.append(((px + ox) % W, (py + oy) % H))

    # (c) uniform random
    for _ in range(n_seeds):
        sites.append((rng.uniform(0.0, W), rng.uniform(0.0, H)))

    # (d) fine grid seeds (more uniform than corners only)
    if extra_corner:
        g = 7
        for i in range(g):
            for j in range(g):
                sites.append((i * W / float(g), j * H / float(g)))

    pushed = []
    for (sx, sy) in sites:
        pushed.append(_push_to_local_max(sx, sy, layout, W, H, iters=12))

    uniq = []
    for (px, py) in pushed:
        keep = True
        for (qx, qy, _qd) in uniq:
            if _pdist2(px, py, qx, qy, W, H) < 4e-4:
                keep = False
                break
        if keep:
            uniq.append((px, py, _slack(px, py, layout, W, H)))
    uniq.sort(key=lambda t: -t[2])
    return [(u[0], u[1]) for u in uniq]


def _max_fit_diameter(x, y, layout, W, H, gap, diameters_desc):
    s = _slack(x, y, layout, W, H)
    dmax = 2.0 * (s - gap)
    for d in diameters_desc:
        if d <= dmax + 1e-9:
            if _fits(x, y, d, layout, W, H, gap):
                return d
    return None


# ---------------------------------------------------------------------------
# Triangular lattice backbones
# ---------------------------------------------------------------------------

def _lattice_fill(d, W, H, gap, rng, tries=10):
    p = d + gap
    nx = int(W / p)
    ny = int(H / p)
    if nx < 1 or ny < 1:
        return []
    row_h = p * math.sqrt(3.0) / 2.0
    best, best_n = [], 0
    for _ in range(tries):
        ox = rng.uniform(0.0, W)
        oy = rng.uniform(0.0, H)
        cand = []
        for j in range(ny + 1):
            y = (oy + j * row_h) % H
            off = 0.0 if j % 2 == 0 else p / 2.0
            for i in range(nx + 1):
                x = (ox + off + i * p) % W
                cand.append((x, y, d))
        kept = []
        for c in cand:
            if _fits(c[0], c[1], c[2], kept, W, H, gap):
                kept.append(c)
        if len(kept) > best_n:
            best_n, best = len(kept), kept
    return best


def _hex_fill(d, W, H, gap, rng, tries=12):
    p = d + gap
    if W < p - 1e-9 or H < p - 1e-9:
        return []
    lat_h = p * math.sqrt(3.0) * 0.5

    ncol_max = max(1, int(W / p + 1e-9))
    nrow_max = max(1, int(H / lat_h + 1e-9))

    col_opts = sorted(set([ncol_max, max(1, ncol_max - 1),
                           max(1, ncol_max + 1)]))
    row_opts = sorted(set([nrow_max, max(1, nrow_max - 1),
                           max(1, nrow_max + 1)]))

    best = []
    best_n = 0
    for ncol in col_opts:
        px = W / float(ncol)
        if px < p - 1e-9:
            continue
        for nrow in row_opts:
            if nrow >= 2:
                py = H / float(nrow)
                if py < lat_h - 1e-9:
                    continue
            else:
                py = H
            for _ in range(tries):
                ox = rng.uniform(0.0, px)
                oy = rng.uniform(0.0, py)
                cand = []
                if nrow >= 2:
                    for j in range(nrow):
                        y = (oy + j * py) % H
                        off = 0.0 if (j & 1) == 0 else px * 0.5
                        for i in range(ncol):
                            x = (ox + off + i * px) % W
                            cand.append((x, y, d))
                else:
                    y = oy % H
                    for i in range(ncol):
                        x = (ox + i * px) % W
                        cand.append((x, y, d))
                kept = []
                for c in cand:
                    if _fits(c[0], c[1], c[2], kept, W, H, gap):
                        kept.append(c)
                if len(kept) > best_n:
                    best_n = len(kept)
                    best = kept
    return best


def _hex_pair_backbone(d_big, d_small, W, H, gap, rng, tries=5):
    """Hex lattice of d_big + d_small at triangular interstitial sites."""
    p = d_big + gap
    if W < p - 1e-9 or H < p - 1e-9:
        return []
    lat_h = p * math.sqrt(3.0) * 0.5

    ncol_max = max(1, int(W / p + 1e-9))
    nrow_max = max(1, int(H / lat_h + 1e-9))

    col_opts = sorted(set([ncol_max, max(1, ncol_max - 1)]))
    row_opts = sorted(set([nrow_max, max(1, nrow_max - 1),
                           max(1, nrow_max - 2)]))

    best = []
    best_a = -1.0
    for ncol in col_opts:
        px = W / float(ncol)
        if px < p - 1e-9:
            continue
        for nrow in row_opts:
            if nrow >= 2:
                py = H / float(nrow)
                if py < lat_h - 1e-9:
                    continue
            else:
                py = H
            for _ in range(tries):
                ox = rng.uniform(0.0, px)
                oy = rng.uniform(0.0, py)
                cand = []
                if nrow >= 2:
                    for j in range(nrow):
                        y = (oy + j * py) % H
                        off = 0.0 if (j & 1) == 0 else px * 0.5
                        for i in range(ncol):
                            x = (ox + off + i * px) % W
                            cand.append((x, y, d_big))
                else:
                    y = oy % H
                    for i in range(ncol):
                        x = (ox + i * px) % W
                        cand.append((x, y, d_big))

                kept = []
                for c in cand:
                    if _fits(c[0], c[1], c[2], kept, W, H, gap):
                        kept.append(c)

                # triangular interstitial sites (up and down triangles)
                if nrow >= 2:
                    for j in range(nrow):
                        y0 = (oy + j * py) % H
                        y1 = (oy + ((j + 1) % nrow) * py) % H
                        off0 = 0.0 if (j & 1) == 0 else px * 0.5
                        off1 = 0.0 if ((j + 1) & 1) == 0 else px * 0.5
                        for i in range(ncol):
                            x0 = (ox + off0 + i * px) % W
                            x0r = (ox + off0 + ((i + 1) % ncol) * px) % W
                            x1 = (ox + off1 + i * px) % W
                            x1r = (ox + off1 + ((i + 1) % ncol) * px) % W
                            # up triangle: (x0,y0),(x0r,y0),(x1,y1)
                            cx = (x0 + x0r + x1) / 3.0
                            cy = (y0 + y0 + y1) / 3.0
                            if _fits(cx % W, cy % H, d_small,
                                     kept, W, H, gap):
                                kept.append((cx % W, cy % H, d_small))
                            # down triangle: (x0,y0),(x1,y1),(x1r,y1)
                            cx2 = (x0 + x1 + x1r) / 3.0
                            cy2 = (y0 + y1 + y1) / 3.0
                            if _fits(cx2 % W, cy2 % H, d_small,
                                     kept, W, H, gap):
                                kept.append((cx2 % W, cy2 % H, d_small))

                a = _area(kept)
                if a > best_a:
                    best_a = a
                    best = kept
    return best


# ---------------------------------------------------------------------------
# MEC void-fill (two-stage: largest fit, then small-tube top-up)
# ---------------------------------------------------------------------------

def _void_fill(layout, W, H, gap, diameters_desc, rng, sweep_cap=1500):
    added = 0
    for _ in range(sweep_cap):
        sites = _mec_sites(layout, W, H, rng, n_seeds=10 + 2 * len(layout))
        placed = False
        for (sx, sy) in sites:
            d = _max_fit_diameter(sx, sy, layout, W, H, gap, diameters_desc)
            if d is not None:
                layout.append((sx, sy, d))
                added += 1
                placed = True
                break
        if not placed:
            break
    return added


def _void_fill_small(layout, W, H, gap, diameters_desc, rng,
                     sweep_cap=400):
    """Insert smallest-fit tubes into remaining narrow gaps (1.0/1.5 only)."""
    small = [d for d in sorted(diameters_desc) if d <= 1.5 + 1e-9]
    if not small:
        return 0
    added = 0
    for _ in range(sweep_cap):
        sites = _mec_sites(layout, W, H, rng, n_seeds=8 + len(layout),
                           extra_corner=False)
        placed = False
        for (sx, sy) in sites[:12]:
            for d in small:
                if _fits(sx, sy, d, layout, W, H, gap):
                    layout.append((sx, sy, d))
                    added += 1
                    placed = True
                    break
            if placed:
                break
        if not placed:
            break
    return added


# ---------------------------------------------------------------------------
# Local refinement (adaptive-step, slack + area non-decreasing)
# ---------------------------------------------------------------------------

def _refine(layout, W, H, gap, diameters_desc, rng, rounds=6,
            strength=0.30):
    best = list(layout)
    best_a = _area(best)
    cur = list(layout)
    for r in range(rounds):
        step_scale = strength * (0.72 ** r)
        n = len(cur)
        for i in range(n):
            xi, yi, di = cur[i]
            others = cur[:i] + cur[i + 1:]
            for _ in range(10):
                dx = rng.uniform(-step_scale, step_scale)
                dy = rng.uniform(-step_scale, step_scale)
                nx = (xi + dx) % W
                ny = (yi + dy) % H
                if _fits(nx, ny, di, others, W, H, gap):
                    old_slack = _slack(xi, yi, others, W, H)
                    new_slack = _slack(nx, ny, others, W, H)
                    if new_slack >= old_slack - 1e-6:
                        cur[i] = (nx, ny, di)
                        xi, yi = nx, ny
                        break
        # try to add a circle at a MEC site
        sites = _mec_sites(cur, W, H, rng, n_seeds=8, extra_corner=False)
        for (sx, sy) in sites[:6]:
            d = _max_fit_diameter(sx, sy, cur, W, H, gap, diameters_desc)
            if d is not None:
                cur.append((sx, sy, d))
                break
        a = _area(cur)
        if a > best_a + 1e-12:
            best_a = a
            best = list(cur)
        else:
            cur = list(best)
    return best


# ---------------------------------------------------------------------------
# Simulated annealing with diameter-exchange moves
# ---------------------------------------------------------------------------

def _anneal(layout, W, H, gap, diameters, rng,
            steps=7000, T0=1.2, T1=0.02):
    desc = sorted(diameters, reverse=True)
    cur = list(layout)
    cur_a = _area(cur)
    best, best_a = list(cur), cur_a
    for s in range(steps):
        T = T0 * (T1 / T0) ** (s / max(steps - 1, 1))
        op = rng.random()
        cand = list(cur)
        if op < 0.35 and cand:
            # reposition a tube
            i = rng.randrange(len(cand))
            x, y, d = cand.pop(i)
            nx_ = rng.uniform(0.0, W)
            ny_ = rng.uniform(0.0, H)
            if _fits(nx_, ny_, d, cand, W, H, gap):
                cand.append((nx_, ny_, d))
            else:
                cand.append((x, y, d))
        elif op < 0.60:
            # add a tube
            d = rng.choice(diameters)
            x = rng.uniform(0.0, W)
            y = rng.uniform(0.0, H)
            if _fits(x, y, d, cand, W, H, gap):
                cand.append((x, y, d))
        elif op < 0.85 and cand:
            # diameter exchange: swap a placed tube for a larger/smaller class
            i = rng.randrange(len(cand))
            x, y, d_old = cand[i]
            others = cand[:i] + cand[i + 1:]
            trial = rng.sample(desc, len(desc))
            swapped = False
            for d_new in trial:
                if d_new == d_old:
                    continue
                if _fits(x, y, d_new, others, W, H, gap):
                    cand = others + [(x, y, d_new)]
                    swapped = True
                    break
            if not swapped:
                cand = list(cur)
        else:
            # remove a tube, then try to add larger tubes
            if cand:
                cand.pop(rng.randrange(len(cand)))
            for _ in range(40):
                d = rng.choice(desc)
                x = rng.uniform(0.0, W)
                y = rng.uniform(0.0, H)
                if _fits(x, y, d, cand, W, H, gap):
                    cand.append((x, y, d))
        ca = _area(cand)
        if ca >= cur_a or rng.random() < math.exp((ca - cur_a) / max(T, 1e-9)):
            cur, cur_a = cand, ca
            if cur_a > best_a:
                best, best_a = list(cur), cur_a
    return best, best_a


# ---------------------------------------------------------------------------
# Top-level make_layout
# ---------------------------------------------------------------------------

def _try_add_best(cand_fn, best, best_a):
    try:
        lay = cand_fn()
    except Exception:
        return best, best_a
    if not lay:
        return best, best_a
    a = _area(lay)
    if a > best_a:
        return lay, a
    return best, best_a


def make_layout(W, H, gap, diameters, rng):
    diameters_desc = sorted(diameters, reverse=True)
    best = []
    best_a = -1.0

    # ---- Path A: triangular lattice backbones (all diameters) ----
    try:
        for d in diameters_desc:
            lat = _lattice_fill(d, W, H, gap, rng, tries=10)
            if not lat:
                continue
            lay = list(lat)
            _void_fill(lay, W, H, gap, diameters_desc, rng, sweep_cap=1400)
            _void_fill_small(lay, W, H, gap, diameters_desc, rng,
                             sweep_cap=300)
            lay = _refine(lay, W, H, gap, diameters_desc, rng,
                          rounds=5, strength=0.30)
            _void_fill(lay, W, H, gap, diameters_desc, rng, sweep_cap=500)
            a = _area(lay)
            if a > best_a:
                best_a = a
                best = lay
    except Exception:
        pass

    # ---- Path B: periodic-corrected hex backbones (all diameters) ----
    try:
        for dback in diameters_desc:
            lat = _hex_fill(dback, W, H, gap, rng, tries=14)
            if not lat:
                continue
            lay = list(lat)
            _void_fill(lay, W, H, gap, diameters_desc, rng, sweep_cap=1400)
            _void_fill_small(lay, W, H, gap, diameters_desc, rng,
                             sweep_cap=300)
            lay = _refine(lay, W, H, gap, diameters_desc, rng,
                          rounds=5, strength=0.30)
            _void_fill(lay, W, H, gap, diameters_desc, rng, sweep_cap=500)
            a = _area(lay)
            if a > best_a:
                best_a = a
                best = lay
    except Exception:
        pass

    # ---- Path C: explicit hex-pair interstitial backbones ----
    try:
        for dbig in diameters_desc:
            for dsmall in diameters_desc:
                if dsmall > dbig + 1e-9:
                    continue
                hp = _hex_pair_backbone(dbig, dsmall, W, H, gap, rng,
                                        tries=4)
                if not hp:
                    continue
                _void_fill(hp, W, H, gap, diameters_desc, rng,
                           sweep_cap=1200)
                _void_fill_small(hp, W, H, gap, diameters_desc, rng,
                                 sweep_cap=300)
                hp = _refine(hp, W, H, gap, diameters_desc, rng,
                             rounds=5, strength=0.28)
                _void_fill(hp, W, H, gap, diameters_desc, rng,
                           sweep_cap=500)
                a = _area(hp)
                if a > best_a:
                    best_a = a
                    best = list(hp)
    except Exception:
        pass

    # ---- Path D: MEC Apollonian growth from diameter-diverse seeds ----
    for r_try in range(6):
        try:
            lay = []
            d0 = diameters_desc[r_try % len(diameters_desc)]
            lay.append((rng.uniform(0.0, W), rng.uniform(0.0, H), d0))
            for _ in range(300):
                sites = _mec_sites(lay, W, H, rng, n_seeds=8 + 2 * len(lay))
                placed = False
                for (sx, sy) in sites:
                    d = _max_fit_diameter(sx, sy, lay, W, H, gap,
                                          diameters_desc)
                    if d is not None:
                        lay.append((sx, sy, d))
                        placed = True
                        break
                if not placed:
                    break
            _void_fill(lay, W, H, gap, diameters_desc, rng, sweep_cap=1200)
            _void_fill_small(lay, W, H, gap, diameters_desc, rng,
                             sweep_cap=300)
            lay = _refine(lay, W, H, gap, diameters_desc, rng,
                          rounds=5, strength=0.28)
            a = _area(lay)
            if a > best_a:
                best_a = a
                best = lay
        except Exception:
            continue

    # ---- Final polish: 3 annealing rounds then MEC refill ----
    try:
        for _round in range(3):
            lay, a = _anneal(best, W, H, gap, diameters, rng,
                             steps=7000, T0=1.2, T1=0.02)
            if a > best_a:
                best_a = a
                best = lay
            best = _refine(best, W, H, gap, diameters_desc, rng,
                           rounds=3, strength=0.20)
            _void_fill(best, W, H, gap, diameters_desc, rng, sweep_cap=800)
            _void_fill_small(best, W, H, gap, diameters_desc, rng,
                             sweep_cap=200)
            a2 = _area(best)
            if a2 > best_a:
                best_a = a2
    except Exception:
        pass

    # Final safety: drop any illegal circle.
    clean = []
    for (x, y, d) in best:
        x = x % W
        y = y % H
        if _fits(x, y, d, clean, W, H, gap):
            clean.append((x, y, d))
    return clean
# EVOLVE-BLOCK-END


# 供 EVOLVE-BLOCK 内部使用（定义在块外，避免被变异破坏）


def run_experiment(random_seed=None, **kwargs):
    """ShinkaEvolve 评估入口（签名与 evaluate.py 契约固定，勿改动）。"""
    W = float(kwargs.get("box_w", DEFAULT_BOX_W))
    H = float(kwargs.get("box_h", DEFAULT_BOX_H))
    gap = float(kwargs.get("safety_gap", DEFAULT_SAFETY_GAP))
    diameters = list(kwargs.get("diameters", DEFAULT_DIAMETERS))
    rng = random.Random(random_seed)
    return make_layout(W, H, gap, diameters, rng)
