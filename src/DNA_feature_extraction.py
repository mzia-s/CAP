#!/usr/bin/env python
"""
DNA_feature_extraction.py
=========================

* Point-cloud generation  →  *.npz
* Facet-curve & rates extraction  →  *.csv

Only the DNA feature extraction is released in this code. Protein sequence embeddings can be generated with Transformer Protein language model ESM2. Run ``python DNA_feature_extraction.py --help``
for usage.
"""
from __future__ import annotations

import argparse
import bisect
import math
import os
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

try:
    import gudhi  # pip install gudhi
except ImportError as exc:  # pragma: no cover
    print("[!] Gudhi not installed — run `pip install gudhi`.", file=sys.stderr)
    raise exc

###############################################################################
# Legacy facet implementation (verbatim)                                      #
###############################################################################

class facet:  # noqa: N801 – keep original name
    def __init__(
        self,
        points: np.ndarray,
        max_dim: int = 0,
        min_edge_length: float = 0.0,
        max_edge_length: float = 24.5,
        num_samples: int = 50,
    ) -> None:
        self.points = points.astype(float)
        self.max_dim = max_dim
        self.min_edge_length = min_edge_length
        self.max_edge_length = max_edge_length
        self.num_samples = num_samples

    # ────────────────────────────────────────────────────────── Vietoris–Rips ──
    def vietoris_rips_simplices_birth_death(self):
        n = len(self.points)
        tree = gudhi.RipsComplex(
            points=self.points, max_edge_length=self.max_edge_length
        ).create_simplex_tree(max_dimension=self.max_dim + 1)

        subset_birth = {
            frozenset(s): f
            for s, f in tree.get_simplices()
            if (len(s) - 1) <= self.max_dim
        }

        subset_death = {}
        for subset in subset_birth:
            k = len(subset)
            if k == self.max_dim + 1:
                subset_death[subset] = math.inf
            else:
                cand = [
                    subset_birth[subset | {v}]
                    for v in range(n)
                    if v not in subset
                    and len(subset | {v}) == k + 1
                    and (subset | {v}) in subset_birth
                ]
                subset_death[subset] = min(cand) if cand else math.inf

        barcodes: Dict[int, List[tuple[float, float]]] = {
            d: [] for d in range(self.max_dim + 1)
        }
        for subset, birth in subset_birth.items():
            dim = len(subset) - 1
            death = subset_death[subset]
            if birth != death:
                barcodes[dim].append((birth, death))
        return barcodes

    # ───────────────────────────────────────────────────────────── helper I/O ──
    def prepare_intervals_from_barcodes(self, dimension: int | None = None):
        bc = self.vietoris_rips_simplices_birth_death()
        intervals = (
            bc.get(dimension, []) if dimension is not None else [iv for v in bc.values() for iv in v]
        )
        births = sorted(b for b, _ in intervals)
        deaths = sorted(d for _, d in intervals)
        return births, deaths

    def count_active_intervals_sorted(self, births, deaths, t):
        return bisect.bisect_right(births, t) - bisect.bisect_left(deaths, t)

    def compute_active_curve(self, dimension, t_values):
        births, deaths = self.prepare_intervals_from_barcodes(dimension)
        return [self.count_active_intervals_sorted(births, deaths, t) for t in t_values]

    def compute_active_rates(self, dimension, t_values):
        births, deaths = self.prepare_intervals_from_barcodes(dimension)
        return [
            self.count_active_intervals_sorted(births, deaths, t) / t if t else 0.0
            for t in t_values
        ]

    # ─────────────────────────────────────────────────────────── public curves ─
    def facet_curves(self):
        t = np.linspace(self.min_edge_length, self.max_edge_length, self.num_samples)
        return [self.compute_active_curve(d, t) for d in range(self.max_dim + 1)]

    def facet_rates(self):
        t = np.linspace(self.min_edge_length, self.max_edge_length, self.num_samples)
        return [self.compute_active_rates(d, t) for d in range(self.max_dim + 1)]


###############################################################################
# 1. K-mer position clouds                                                     #
###############################################################################

def extract_kmer_position_clouds(sequence: str, k: int = 1) -> Dict[str, np.ndarray]:
    seq = sequence.upper().replace("U", "T")
    allowed = {"A", "C", "G", "T"}
    kmers = ["".join(p) for p in product("ACGT", repeat=k)]
    clouds: Dict[str, List[int]] = {km: [] for km in kmers}

    for i in range(len(seq) - k + 1):
        window = seq[i : i + k]
        if set(window) <= allowed:
            clouds[window].append(i + 1)  # 1-based

    return {
        km: np.asarray(pos, dtype=np.float32).reshape(-1, 1) for km, pos in clouds.items()
    }


def write_point_clouds_from_fasta(fasta_dir: Path, out_dir: Path, k: int = 1) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for f in sorted(fasta_dir.glob("*.fasta")):
        seq = "".join(p.strip() for p in f.read_text().splitlines()[1:])
        clouds = extract_kmer_position_clouds(seq, k=k)
        np.savez_compressed(out_dir / f"{f.stem}_k{k}_clouds.npz", **clouds)
        print(f"[✓] clouds → {f.stem}_k{k}_clouds.npz")
    print(f"Generated {len(list(fasta_dir.glob('*.fasta')))} file(s) → {out_dir}\n")


###############################################################################
# 2. Modern Facet wrapper (cached)                                             #
###############################################################################

@dataclass
class FacetOpts:
    max_dim: int = 2
    min_edge: float = 0.0
    max_edge: float = 24.5
    num_samples: int = 50


class FacetCurves:
    def __init__(self, points: np.ndarray, opts: FacetOpts):
        self.points = points.astype(float)
        self.o = opts
        self.tree = (gudhi.RipsComplex(points=self.points, max_edge_length=self.o.max_edge)
            ).create_simplex_tree(max_dimension=self.o.max_dim + 1)
        self._bc: Dict[int, List[tuple[float, float]]] | None = None

    # cached barcodes
    def _barcodes(self):
        if self._bc is not None:
            return self._bc
        sb = {
            frozenset(s): f
            for s, f in self.tree.get_simplices()
            if (len(s) - 1) <= self.o.max_dim
        }
        sd, n = {}, len(self.points)
        for subset in sb:
            k = len(subset)
            if k == self.o.max_dim + 1:
                sd[subset] = math.inf
            else:
                cand = [
                    sb[subset | {v}]
                    for v in range(n)
                    if v not in subset and len(subset | {v}) == k + 1 and (subset | {v}) in sb
                ]
                sd[subset] = min(cand) if cand else math.inf
        bc: Dict[int, List[tuple[float, float]]] = {d: [] for d in range(self.o.max_dim + 1)}
        for subset, birth in sb.items():
            dim, death = len(subset) - 1, sd[subset]
            if birth != death:
                bc[dim].append((birth, death))
        self._bc = bc
        return bc

    # curves + rates
    def curves_and_rates(self):
        grid = np.linspace(self.o.min_edge, self.o.max_edge, self.o.num_samples)
        curves, rates = [], []
        for d in range(self.o.max_dim + 1):
            births, deaths = zip(*self._barcodes().get(d, [])) if self._barcodes().get(d) else ([], [])
            births, deaths = np.sort(births), np.sort(deaths)
            counts = [
                bisect.bisect_right(births, t) - bisect.bisect_left(deaths, t) for t in grid
            ]
            curves.append(counts)
            rates.append([c / t if t else 0.0 for c, t in zip(counts, grid)])
        return np.array(curves), np.array(rates)


###############################################################################
# 3. Feature matrix builder                                                    #
###############################################################################

def pad(arr: np.ndarray, length: int) -> np.ndarray:
    out = np.zeros(length, dtype=float)
    out[: min(len(arr), length)] = arr[:length]
    return out


def features_from_npz(
    npz_dir: Path,
    out_csv: Path,
    k: int = 1,
    opts: FacetOpts | None = None,
) -> None:
    opts = opts or FacetOpts()
    expected_block = (opts.max_dim + 1) * opts.num_samples * 2
    rows, ids = [], []

    for fn in sorted(npz_dir.glob(f"*_k{k}_clouds.npz")):
        clouds = np.load(fn)
        blocks = []
        for km in sorted(clouds.files):
            pts = clouds[km]
            if pts.shape[0] < 2:
                blocks.append(np.zeros(expected_block))
                continue
            curves, rates = FacetCurves(pts, opts).curves_and_rates()
            blocks.append(np.concatenate([curves.flatten(), rates.flatten()]))
        rows.append(np.concatenate(blocks))
        ids.append(fn.stem.split("_k", 1)[0])

    df = pd.DataFrame(rows, index=ids)
    df.index.name = "Sequence_ID"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv)
    print(f"[✓] feature matrix → {out_csv}")


###############################################################################
# 4. CLI wrapper                                                               #
###############################################################################

def main() -> None:  # pragma: no cover
    p = argparse.ArgumentParser(description="CAP DNA feature pipeline")
    sub = p.add_subparsers(dest="cmd", required=True)

    cl = sub.add_parser("clouds", help="Generate k-mer point clouds (*.npz) from FASTAs")
    cl.add_argument("fasta_dir", type=Path)
    cl.add_argument("out_dir", type=Path)
    cl.add_argument("--k", type=int, default=1)

    fe = sub.add_parser("features", help="Compute facet curves/rates → CSV")
    fe.add_argument("npz_dir", type=Path)
    fe.add_argument("out_csv", type=Path)
    fe.add_argument("--k", type=int, default=1)
    fe.add_argument("--max-dim", type=int, default=2)
    fe.add_argument("--num-samples", type=int, default=50)

    args = p.parse_args()

    if args.cmd == "clouds":
        write_point_clouds_from_fasta(args.fasta_dir, args.out_dir, k=args.k)
    elif args.cmd == "features":
        opts = FacetOpts(max_dim=args.max_dim, num_samples=args.num_samples)
        features_from_npz(args.npz_dir, args.out_csv, k=args.k, opts=opts)
    else:
        p.error("Unknown command")


if __name__ == "__main__":
    main()
