# Windows convenience script (equivalent of `make sweep`).
# Run from the repo root: .\scripts\sweep.ps1
# Requires main.exe already built, e.g.:
#   gcc -O2 -Wall -Wextra -o src\main.exe src\integrator.c src\fields.c src\gyro_avg.c src\main.c -lm

New-Item -ItemType Directory -Force -Path src\sweep | Out-Null

.\src\main.exe 0.02    15000  10  src\sweep\traj_dt0.02.csv
.\src\main.exe 0.01    30000  20  src\sweep\traj_dt0.01.csv
.\src\main.exe 0.005   60000  40  src\sweep\traj_dt0.005.csv
.\src\main.exe 0.002   150000 100 src\sweep\traj_dt0.002.csv
.\src\main.exe 0.001   300000 200 src\sweep\traj_dt0.001.csv
.\src\main.exe 0.0005  600000 400 src\sweep\traj_dt0.0005.csv

Write-Host "Sweep complete - 6 CSVs written to src\sweep\"
