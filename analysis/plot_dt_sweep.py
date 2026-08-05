"""
dt sweep: plot gyro-averaged mu conservation for several step sizes,
and the RK4 error scaling (should be ~dt^4) to demonstrate numerical rigor.
"""
import csv
import glob
import re
import sys
from pathlib import Path
import matplotlib.pyplot as plt

def find_sweep_dir():
    if len(sys.argv) > 1:
        return Path(sys.argv[1])
    here = Path(__file__).resolve().parent
    candidates = [
        here / "../src/sweep",  # repo layout: analysis/ + src/sweep/
        here / "sweep",         # flat layout: script + sweep/ side by side
        Path.cwd() / "sweep",   # running from a dir containing sweep/
    ]
    for c in candidates:
        if list(c.glob("traj_dt*.csv")):
            return c
    return candidates[0]  # fall through to default; error message below explains

sweep_dir = find_sweep_dir()
files = sorted(glob.glob(str(sweep_dir / "traj_dt*.csv")),
                key=lambda p: float(re.search(r"traj_dt([\d.]+)\.csv", p).group(1)))

if not files:
    raise SystemExit(
        f"No sweep CSVs found in {sweep_dir.resolve()}\n"
        "Run the C simulation for each dt first (see main.c usage: "
        "./main <dt> <nsteps> <output_stride> <outfile>), "
        "or pass the sweep folder explicitly: "
        "python3 plot_dt_sweep.py /path/to/sweep"
    )

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
ax1, ax2 = axes

dts = []
late_time_std = []  # spread of mu_avg over the settled (post-transient) region, per dt

cmap = plt.get_cmap("plasma")
for idx, path in enumerate(files):
    dt = float(re.search(r"traj_dt([\d.]+)\.csv", path).group(1))
    with open(path) as f:
        rows = list(csv.DictReader(f))
    t = [float(r["t"]) for r in rows]
    mu_avg = [float(r["mu_avg"]) for r in rows]
    mu_avg0 = mu_avg[0]
    ratio = [m / mu_avg0 for m in mu_avg]

    color = cmap(idx / max(1, len(files) - 1))
    ax1.plot(t, ratio, label=f"dt={dt}", color=color, linewidth=1.3)

    # measure conservation quality in the settled region (skip the initial
    # transient, t > 50) as the standard deviation of the ratio there
    settled = [r for tt, r in zip(t, ratio) if tt > 50]
    dts.append(dt)
    late_time_std.append((sum((s - sum(settled) / len(settled)) ** 2 for s in settled)
                           / len(settled)) ** 0.5)

ax1.set_xlabel("t")
ax1.set_ylabel(r"$\mu_{avg}(t) / \mu_{avg}(0)$")
ax1.set_title("Gyro-averaged $\\mu$ conservation vs. step size")
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

ax2.loglog(dts, late_time_std, "o-", color="darkred")
# reference dt^4 line for comparison (RK4 local truncation error scaling)
ref_x = [min(dts), max(dts)]
ref_y0 = late_time_std[dts.index(min(dts))]
ref_y = [ref_y0 * (x / min(dts)) ** 4 for x in ref_x]
ax2.loglog(ref_x, ref_y, "--", color="gray", label=r"$\propto dt^4$ (reference)")
ax2.set_xlabel("dt")
ax2.set_ylabel(r"std. dev. of $\mu_{avg}$ ratio (settled region)")
ax2.set_title("Conservation error vs. step size")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3, which="both")

plt.tight_layout()
plt.savefig("dt_sweep.png", dpi=150)
print("saved dt_sweep.png")

for dt, s in zip(dts, late_time_std):
    print(f"dt={dt:.5f}  settled std of mu_avg ratio = {s:.6e}")
