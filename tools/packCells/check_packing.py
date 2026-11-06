#!/usr/bin/env python3
"""Check a packing written by packCells: no enclosing ellipsoids overlap each
other or a wall below a minimum scale.

usage: check_packing.py sX sY sZ [options] FILE.pos [FILE.pos ...]

The box and the walls are given as to packCells (--walls, --cylinder,
--offset). The cell type of a file is its name (RBC.pos: RBC, with the
ellipsoid packCells uses for it); --cell NAME DX DY DZ gives the ellipsoid of
a type of your own, or of a renamed file (NAME.pos).

For every pair of cells, also through the periodic images along the periodic
axes (a cell's own image included), the scale at which their nominal
ellipsoids touch is computed from the Perram-Wertheim contact function; for
the walls, the scale at which an ellipsoid touches a face of the box (exact)
or the cylinder (from 3600 points around its cross-section). The smallest of
these must be at least --min_scale (default: the largest mesh fit of the
types, as for packCells; 1 for --cell types).

Exit status: 0 if the packing passes, 1 if not, 2 for invalid arguments.
Needs numpy.
"""
import argparse
import itertools
import math
import os
import sys

import numpy as np

# Ellipsoid diameters [um] along the cell's x, y and z axes, and the scale at
# which the ellipsoid still contains the cell's mesh (the "mesh fit"), as in
# packCells.cpp (builtinTypes): keep them equal.
BUILTIN = {
    'RBC': (8.4, 4.4, 8.4, 0.931),
    'PLT': (2.7, 1.17, 2.7, 0.929),
    'WBC': (8.4, 8.4, 8.4, 0.952),
    'vRBC': (3.5, 6.0, 11.0, 1.089),
    'RBC_m': (6.23, 3.26, 6.23, 0.931),
    'PLT_m': (2.05, 1.17, 2.05, 0.927),
    'PLT_mko': (3.68, 3.68, 3.68, 0.929),
}


def rotation(rx, ry, rz):
    """Body-to-world rotation of a .pos line's angles [degrees]: about x,
    then y, then z, as packCells writes them and HemoCell applies them."""
    a, b, c = np.radians([rx, ry, rz])
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    Rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def touching_scales(MA, MB, r):
    """Perram-Wertheim: the scales at which pairs of ellipsoids
    {x : x^T M^-1 x <= 1}, B at A + r, touch: sqrt(max over lambda of
    lambda (1 - lambda) r^T ((1 - lambda) MA + lambda MB)^-1 r), the function
    being concave in lambda. MA, MB: (n, 3, 3); r: (n, 3)."""
    def F(lam):
        C = (1 - lam)[:, None, None] * MA + lam[:, None, None] * MB
        return lam * (1 - lam) * np.einsum('ni,ni->n', r, np.linalg.solve(C, r[:, :, None])[:, :, 0])
    n = len(r)
    grid = np.linspace(0, 1, 21)[1:-1]
    values = np.array([F(np.full(n, g)) for g in grid])
    best = values.argmax(axis=0)
    lo = np.where(best > 0, grid[np.maximum(best - 1, 0)], 0.0)
    hi = np.where(best < len(grid) - 1, grid[np.minimum(best + 1, len(grid) - 1)], 1.0)
    golden = (math.sqrt(5) - 1) / 2
    for _ in range(60):
        c = hi - golden * (hi - lo)
        d = lo + golden * (hi - lo)
        left = F(c) > F(d)
        hi = np.where(left, d, hi)
        lo = np.where(left, lo, c)
    return np.sqrt(np.maximum(F(0.5 * (lo + hi)), values.max(axis=0)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('size', nargs=3, type=float, metavar='S', help='box size sX sY sZ [um]')
    parser.add_argument('files', nargs='+', metavar='FILE.pos')
    parser.add_argument('--walls', default='', help='axes with walls, e.g. z or yz')
    parser.add_argument('--cylinder', nargs=2, metavar=('AXIS', 'R'))
    parser.add_argument('--offset', nargs=3, type=float, default=[0, 0, 0], metavar=('X', 'Y', 'Z'))
    parser.add_argument('--cell', nargs=4, action='append', default=[], metavar=('NAME', 'DX', 'DY', 'DZ'))
    parser.add_argument('--min_scale', type=float)
    args = parser.parse_intermixed_args()

    box = np.array(args.size)
    if not all(c in 'xyz' for c in args.walls):
        parser.error('--walls: expected axes x, y and/or z')
    cylinder = None
    if args.cylinder:
        if args.cylinder[0] not in ('x', 'y', 'z'):
            parser.error('--cylinder: expected the axis x, y or z')
        cylinder = ('xyz'.index(args.cylinder[0]), float(args.cylinder[1]))
    periodic = [('xyz'[a] not in args.walls) and (cylinder is None or a == cylinder[0]) for a in range(3)]
    types = dict(BUILTIN)
    for name, dx, dy, dz in args.cell:
        types[name] = (float(dx), float(dy), float(dz), 1.0)

    # The cells: type, centre (without the offset), nominal shape matrix M,
    # largest semi-axis.
    names, centres, shapes, reach = [], [], [], []
    counts = {}
    for path in args.files:
        name = os.path.basename(path)[:-4] if path.endswith('.pos') else os.path.basename(path)
        if name not in types:
            parser.error(f'{path}: unknown cell type {name} (give its ellipsoid with --cell {name} DX DY DZ)')
        semi = np.array(types[name][:3]) / 2
        with open(path) as f:
            values = f.read().split()
        n = int(values[0])
        rows = np.array(values[1:1 + 6 * n], dtype=float).reshape(n, 6)
        counts[name] = n
        for row in rows:
            R = rotation(*row[3:])
            names.append(name)
            centres.append(row[:3] - np.array(args.offset))
            shapes.append(R @ np.diag(semi ** 2) @ R.T)
            reach.append(semi.max())
    names = np.array(names)
    C = np.array(centres).reshape(-1, 3)
    M = np.array(shapes).reshape(-1, 3, 3)
    reach = np.array(reach)
    present = [t for t, n in counts.items() if n > 0]
    min_scale = args.min_scale if args.min_scale is not None else max(types[t][3] for t in present)
    cap = 1.05 * max(1.0, min_scale)  # pairs and walls beyond this scale are not computed

    print(f'{len(C)} cells ({", ".join(f"{t} {n}" for t, n in counts.items())}) in '
          f'{box[0]:g} x {box[1]:g} x {box[2]:g} um, periodic along '
          f'{"".join("xyz"[a] for a in range(3) if periodic[a]) or "none"}'
          + (f', cylinder along {"xyz"[cylinder[0]]} of radius {cylinder[1]:g} um' if cylinder else ''))

    # Centres the loader would not accept, or outside the walls
    bad = 0
    for a in range(3):
        if periodic[a]:
            bad += int(np.sum((C[:, a] < 0) | (C[:, a] >= box[a])))
        elif cylinder is None or a == cylinder[0]:
            bad += int(np.sum((C[:, a] < 0) | (C[:, a] > box[a])))
    if cylinder:
        b, c = [a for a in range(3) if a != cylinder[0]]
        bad += int(np.sum(np.hypot(C[:, b] - box[b] / 2, C[:, c] - box[c] / 2) > cylinder[1]))

    # Pairs within reach, through the periodic images (each contact once:
    # i < j, and a cell's own images from one side)
    ranges = [range(-k, k + 1) for k in
              (int(math.ceil(2 * cap * reach.max() / box[a])) if periodic[a] else 0 for a in range(3))]
    pi, pj, pr = [], [], []
    for s in itertools.product(*ranges):
        shift = np.array(s) * box
        positive = s > (0, 0, 0)
        for start in range(0, len(C), 512):
            i = np.arange(start, min(start + 512, len(C)))
            r = C[None, :, :] + shift - C[i, None, :]
            near = np.einsum('ijk,ijk->ij', r, r) < (cap * (reach[i, None] + reach[None, :])) ** 2
            keep = near & ((np.arange(len(C))[None, :] > i[:, None]) |
                           ((np.arange(len(C))[None, :] == i[:, None]) & positive))
            ii, jj = np.nonzero(keep)
            pi.append(i[ii])
            pj.append(jj)
            pr.append(r[ii, jj])
    pi = np.concatenate(pi) if pi else np.zeros(0, int)
    pj = np.concatenate(pj) if pj else np.zeros(0, int)
    pr = np.concatenate(pr) if pr else np.zeros((0, 3))
    pair_min = {}
    if len(pr):
        mu = np.concatenate([touching_scales(M[pi[k:k + 20000]], M[pj[k:k + 20000]], pr[k:k + 20000])
                             for k in range(0, len(pr), 20000)])
        for k in range(len(mu)):
            key = '-'.join(sorted((names[pi[k]], names[pj[k]])))
            pair_min[key] = min(pair_min.get(key, np.inf), mu[k])

    # Walls: the faces exactly (the ellipsoid reaches sqrt(M_aa) along axis
    # a), the cylinder from its cross-section (the projection of the
    # ellipsoid, an ellipse with the 2x2 block of M), bisected for the scale
    # at which its farthest point from the axis reaches R
    wall_min = np.inf
    for a in range(3):
        if 'xyz'[a] in args.walls:
            extent = np.sqrt(M[:, a, a])
            wall_min = min(wall_min, (np.maximum(C[:, a], 0) / extent).min(),
                           (np.maximum(box[a] - C[:, a], 0) / extent).min())
    if cylinder:
        b, c = [a for a in range(3) if a != cylinder[0]]
        q = np.stack([C[:, b] - box[b] / 2, C[:, c] - box[c] / 2], 1)
        near = np.hypot(q[:, 0], q[:, 1]) + cap * reach > cylinder[1]
        if near.any():
            P = M[near][:, [b, c]][:, :, [b, c]]
            L = np.linalg.cholesky(P)
            t = np.linspace(0, 2 * np.pi, 3600, endpoint=False)
            e = np.einsum('nij,jk->nik', L, np.stack([np.cos(t), np.sin(t)]))  # (n, 2, 3600)
            qn = q[near][:, :, None]

            def farthest(s):
                return np.hypot(qn[:, 0] + s[:, None] * e[:, 0], qn[:, 1] + s[:, None] * e[:, 1]).max(1)
            lo = np.zeros(near.sum())
            hi = np.full(near.sum(), cap)
            inside = farthest(np.zeros(near.sum())) < cylinder[1]
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                out = farthest(mid) > cylinder[1]
                hi = np.where(out, mid, hi)
                lo = np.where(out, lo, mid)
            wall_min = min(wall_min, np.where(inside, lo, 0).min())

    smallest = min([wall_min] + list(pair_min.values()))
    shown = ', '.join(f'{k} {v:.4f}' for k, v in sorted(pair_min.items())) or f'none below {cap:.3f}'
    print(f'smallest scale at which two ellipsoids touch: {shown}')
    if args.walls or cylinder:
        print('smallest scale at which an ellipsoid touches a wall: '
              + (f'{wall_min:.4f}' if wall_min < np.inf else f'none below {cap:.3f}'))
    if bad:
        print(f'FAILED: {bad} centres outside the box (or the cylinder), or on its upper face along a periodic axis')
        return 1
    if smallest < min_scale:
        print(f'FAILED: overlap at scale {smallest:.4f}, below the minimum {min_scale:.4f}')
        return 1
    print(f'OK: no overlap below scale {min(smallest, cap):.4f} (minimum {min_scale:.4f})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
