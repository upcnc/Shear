#!/usr/bin/env python3
"""Plot history.csv from wang_ibm_detach."""
import csv
import sys
from pathlib import Path

def load(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append({k: float(v) if k != "step" else int(float(v)) for k, v in row.items()})
    return rows

def main():
    if len(sys.argv) < 2:
        print("usage: python plot_wang_detach.py result_dir [more dirs...]")
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; print table only")
        for d in sys.argv[1:]:
            rows = load(Path(d) / "history.csv")
            print(d, "final detached", rows[-1]["frac_detached"], "bonds", rows[-1]["n_inter_bonds"])
        return
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for d in sys.argv[1:]:
        rows = load(Path(d) / "history.csv")
        t = [r["time_s"] for r in rows]
        ax[0].plot(t, [r["frac_detached"] for r in rows], label=Path(d).name)
        ax[1].plot(t, [r["n_inter_bonds"] for r in rows], label=Path(d).name)
    ax[0].set_xlabel("t / s"); ax[0].set_ylabel("detached fraction")
    ax[1].set_xlabel("t / s"); ax[1].set_ylabel("living inter-unit bonds")
    ax[0].legend(); ax[1].legend()
    fig.tight_layout()
    fig.savefig("wang_detach_compare.png", dpi=140)
    print("wrote wang_detach_compare.png")

if __name__ == "__main__":
    main()
