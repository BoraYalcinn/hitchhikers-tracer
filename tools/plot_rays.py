"""Plot ray paths written by app/dump_rays.f90.

Reads output/rays.csv (columns: ray_id, b, r, phi), turns the polar
coordinates (r, phi) into x, y and draws every ray around the black hole.

Usage (from the repository root, with the virtual environment active):
    python tools/plot_rays.py                  # output/rays.csv -> output/rays.png
    python tools/plot_rays.py in.csv out.png   # custom file names
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # draw into a file, no window needed
import matplotlib.pyplot as plt
import numpy as np

R_HORIZON = 2.0                    # event horizon, 2M
R_PHOTON_SPHERE = 3.0              # photon sphere, 3M
B_CRITICAL = 3.0 * np.sqrt(3.0)    # critical impact parameter, about 5.196


def load(path):
    """Read the CSV into a table with named columns."""
    return np.genfromtxt(path, delimiter=",", names=True)


def plot(data, out_path):
    fig, ax = plt.subplots(figsize=(9, 7))

    # Black hole and photon sphere.
    ax.add_patch(plt.Circle((0, 0), R_HORIZON, color="black", zorder=3))
    ax.add_patch(plt.Circle((0, 0), R_PHOTON_SPHERE, fill=False,
                            linestyle="--", color="gray", zorder=2))

    # One line per ray, coloured by impact parameter b.
    ray_ids = np.unique(data["ray_id"])
    colors = plt.cm.viridis(np.linspace(0.0, 0.95, len(ray_ids)))

    for color, ray_id in zip(colors, ray_ids):
        rows = data[data["ray_id"] == ray_id]
        x = rows["r"] * np.cos(rows["phi"])
        y = rows["r"] * np.sin(rows["phi"])
        b = rows["b"][0]
        fate = "falls in" if b < B_CRITICAL else "escapes"
        ax.plot(x, y, color=color, linewidth=1.5,
                label=f"b = {b:g}  ({fate})")

    ax.set_aspect("equal")        # circles must look like circles
    ax.set_xlim(-14, 34)
    ax.set_ylim(-16, 16)
    ax.set_xlabel("x  (M)")
    ax.set_ylabel("y  (M)")
    ax.set_title("Light rays around a Schwarzschild black hole")
    ax.legend(loc="lower right", fontsize=8,
              title=f"b_c = 3*sqrt(3) = {B_CRITICAL:.3f}")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


def main():
    in_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("output/rays.csv")
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("output/rays.png")
    if not in_path.exists():
        print(f"{in_path} not found. Run dump_rays first.")
        sys.exit(1)
    plot(load(in_path), out_path)


if __name__ == "__main__":
    main()