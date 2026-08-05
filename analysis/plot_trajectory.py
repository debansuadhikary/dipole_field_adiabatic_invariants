"""
Plot the 3D trajectory written by the C simulation (src/main -> trajectory.csv).

Usage: python3 plot_trajectory.py [path/to/trajectory.csv]
"""
import sys
import csv
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "../src/trajectory.csv"

t, x, y, z = [], [], [], []
with open(path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        t.append(float(row["t"]))
        x.append(float(row["x"]))
        y.append(float(row["y"]))
        z.append(float(row["z"]))

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Color the path by time so gyration/bounce/drift structure is visible
sc = ax.scatter(x, y, z, c=t, cmap="viridis", s=1)
ax.plot(x, y, z, color="gray", alpha=0.15, linewidth=0.5)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Charged particle trajectory in dipole field\n(gyration + bounce + drift)")
cbar = fig.colorbar(sc, ax=ax, shrink=0.6, label="time")

plt.tight_layout()
plt.savefig("trajectory_3d.png", dpi=150)
print("saved trajectory_3d.png")
