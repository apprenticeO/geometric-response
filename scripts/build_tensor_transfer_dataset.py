#!/usr/bin/env python3
"""
Build a tensor-channel dataset for the slow-band transfer identity check from the
fast delta-ladder TEBD script (6a_run_delta_ladder_tebd_fast.py / run_delta_ladder).

Outputs a CSV with:
  - t
  - R_info_<comp> columns built from triad components (√F_Q · S_A · I(A:Ā)) per A
  - optional sqrtF_<comp>, S_<comp>, I_<comp> diagnostic columns
  - optional R_<comp> columns if a response CSV is provided (either per-comp or scalar R)
Also saves a small JSON with metadata (L, ell, dt, steps, sample_k, projector).

Usage example:
  python3 scripts/build_tensor_transfer_dataset.py \
      --L 8 --ell 4 --dt 0.05 --steps 60 --chi 64 --cut 1e-9 \
      --sample-k 5 --out-csv TRIAD_Geometry/results/tensor_transfer.csv

To attach a scalar response series and replicate across components:
  python3 scripts/build_tensor_transfer_dataset.py \
      ... --response-csv TRIAD_Geometry/results/two_qubit_pipeline/two_qubit.csv \
      --response-col R --replicate-response
"""
from __future__ import annotations

import argparse
import json
import os
from typing import List, Dict, Any, Optional

import numpy as np
import pandas as pd


def load_run_delta_ladder(module_path: str):
    import importlib.util
    spec = importlib.util.spec_from_file_location('delta_ladder_mod', module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import run_delta_ladder from {module_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    if not hasattr(mod, 'run_delta_ladder'):
        raise AttributeError("run_delta_ladder not found in module")
    return getattr(mod, 'run_delta_ladder')


def main():
    ap = argparse.ArgumentParser(description="Build tensor-channel dataset from fast delta-ladder TEBD output")
    ap.add_argument("--L", type=int, default=8)
    ap.add_argument("--ell", type=int, default=4)
    ap.add_argument("--J", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.0)
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--steps", type=int, default=60)
    ap.add_argument("--chi", type=int, default=64)
    ap.add_argument("--cut", type=float, default=1e-9)
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-6], help="epsilon(s) for commutator ladder (unused here, passthrough)")
    ap.add_argument("--sample-k", type=int, default=5)
    ap.add_argument("--out-csv", type=str, required=True, help="Output CSV path")
    ap.add_argument("--out-json", type=str, default="", help="Optional metadata JSON path")
    # Response attachment (optional)
    ap.add_argument("--response-csv", type=str, default="", help="CSV with response series to attach")
    ap.add_argument("--response-col", type=str, default="R", help="Column name in response CSV for scalar response, or prefix 'R_' for per-comp")
    ap.add_argument("--replicate-response", action="store_true", help="If scalar response, replicate across components")
    # Diagnostics
    ap.add_argument("--include-legs", action="store_true", help="Include sqrtF_<comp>, S_<comp>, I_<comp> columns")
    args = ap.parse_args()

    # Locate sibling 6_run_delta_ladder_tebd.py
    here = os.path.dirname(os.path.abspath(__file__))
    ladder_module = os.path.join(os.path.dirname(here), "Scripts_Dyanmics_triad", "6_run_delta_ladder_tebd.py")
    run_delta_ladder = load_run_delta_ladder(ladder_module)

    # Run the fast ladder with provided params
    res = run_delta_ladder(
        L=args.L, ell=args.ell, J=args.J, h=args.h,
        dt=args.dt, steps=args.steps, chi=args.chi, cut=args.cut,
        epsilons=list(map(float, args.eps)),
        sample_k=args.sample_k,
        alpha=1.0, S_floor=0.0, I_floor=0.0, use_norm=False, report_gap=False
    )
    # res is a dataclass; convert to plain dict
    # Expect fields: triad_values [T x ell], sqrtF_values [T x ell], S_values [T x ell], I_values [T x ell]
    res_dict: Dict[str, Any] = {}
    for f in res.__dataclass_fields__:  # type: ignore[attr-defined]
        res_dict[f] = getattr(res, f)

    triad_vals: List[List[float]] = res_dict.get("triad_values", [])
    sqrtF_vals: List[List[float]] = res_dict.get("sqrtF_values", [])
    S_vals: List[List[float]] = res_dict.get("S_values", [])
    I_vals: List[List[float]] = res_dict.get("I_values", [])

    triad_arr = np.array(triad_vals, dtype=float) if triad_vals else np.zeros((0, args.ell), dtype=float)
    sqrtF_arr = np.array(sqrtF_vals, dtype=float) if sqrtF_vals else np.zeros_like(triad_arr)
    S_arr = np.array(S_vals, dtype=float) if S_vals else np.zeros_like(triad_arr)
    I_arr = np.array(I_vals, dtype=float) if I_vals else np.zeros_like(triad_arr)

    T = triad_arr.shape[0]
    if T == 0:
        raise RuntimeError("No samples returned by run_delta_ladder; increase steps or reduce sample_k")

    # Reconstruct sampled time grid (sampled every sample_k TEBD steps)
    t = np.arange(1, T + 1, dtype=float) * (args.dt * args.sample_k)

    # Build output dataframe
    cols: Dict[str, np.ndarray] = {"t": t}
    components = [f"A{i+1}" for i in range(args.ell)]
    for i, comp in enumerate(components):
        cols[f"R_info_{comp}"] = triad_arr[:, i]
        if args.include_legs:
            cols[f"sqrtF_{comp}"] = sqrtF_arr[:, i]
            cols[f"S_{comp}"] = S_arr[:, i]
            cols[f"I_{comp}"] = I_arr[:, i]
    # Provide a scalar Π as a simple systemic channel (sum of components)
    try:
        ri_stack = np.stack([cols[f"R_info_{c}"] for c in components], axis=0)
        cols["Pi"] = np.sum(ri_stack, axis=0)
    except Exception:
        pass

    # Optional response attachment
    if args.response_csv:
        df_resp = pd.read_csv(args.response_csv)
        # Align by length via min length (simple alignment)
        n_min = min(len(df_resp), T)
        if n_min <= 0:
            raise ValueError("Response CSV has no rows")
        # Per-component response present?
        per_comp_available = all((f"R_{c}" in df_resp.columns) for c in components)
        if per_comp_available:
            for c in components:
                cols[f"R_{c}"] = np.asarray(df_resp[f"R_{c}"].to_numpy(dtype=float)[:n_min])
            # Truncate time/others too
            for k in list(cols.keys()):
                cols[k] = cols[k][:n_min]
        else:
            # Scalar response column
            if args.response_col not in df_resp.columns:
                raise ValueError(f"Response column '{args.response_col}' not found in {args.response_csv}")
            R_scalar = np.asarray(df_resp[args.response_col].to_numpy(dtype=float)[:n_min])
            # Truncate time/others too
            for k in list(cols.keys()):
                cols[k] = cols[k][:n_min]
            if args.replicate_response:
                for c in components:
                    cols[f"R_{c}"] = R_scalar.copy()
            else:
                cols["R"] = R_scalar

    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True) if os.path.dirname(args.out_csv) else None
    pd.DataFrame(cols).to_csv(args.out_csv, index=False)
    print(f"Saved tensor dataset: {args.out_csv} (components={','.join(components)})")

    if args.out_json:
        meta = dict(
            L=args.L, ell=args.ell, J=args.J, h=args.h, dt=args.dt, steps=args.steps,
            chi=args.chi, cut=args.cut, sample_k=args.sample_k, eps=args.eps,
            components=components,
            response_attached=bool(args.response_csv),
            response_cols=[c for c in cols.keys() if c.startswith("R_")] if args.response_csv else [],
        )
        with open(args.out_json, "w") as f:
            json.dump(meta, f, indent=2)
        print(f"Saved metadata: {args.out_json}")


if __name__ == "__main__":
    main()


