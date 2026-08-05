CC = gcc
CFLAGS = -O2 -Wall -Wextra
LDFLAGS = -lm

SRC = src

.PHONY: all test run sweep clean

all: $(SRC)/main $(SRC)/test_sho $(SRC)/test_fields

$(SRC)/main: $(SRC)/integrator.c $(SRC)/fields.c $(SRC)/gyro_avg.c $(SRC)/main.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(SRC)/test_sho: $(SRC)/integrator.c $(SRC)/test_sho.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(SRC)/test_fields: $(SRC)/fields.c $(SRC)/test_fields.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Run both validation harnesses (RK4 vs. analytic SHO, dipole field vs. analytic formulas)
test: $(SRC)/test_sho $(SRC)/test_fields
	./$(SRC)/test_sho
	./$(SRC)/test_fields

# Single run at default dt, writes src/trajectory.csv
run: $(SRC)/main
	cd $(SRC) && ./main

# Reproduce the dt sweep (six step sizes, same total simulated time)
sweep: $(SRC)/main
	mkdir -p $(SRC)/sweep
	cd $(SRC) && ./main 0.02    15000  10  sweep/traj_dt0.02.csv
	cd $(SRC) && ./main 0.01    30000  20  sweep/traj_dt0.01.csv
	cd $(SRC) && ./main 0.005   60000  40  sweep/traj_dt0.005.csv
	cd $(SRC) && ./main 0.002   150000 100 sweep/traj_dt0.002.csv
	cd $(SRC) && ./main 0.001   300000 200 sweep/traj_dt0.001.csv
	cd $(SRC) && ./main 0.0005  600000 400 sweep/traj_dt0.0005.csv

clean:
	rm -f $(SRC)/main $(SRC)/test_sho $(SRC)/test_fields
	rm -f $(SRC)/*.csv $(SRC)/sweep/*.csv
	rm -f analysis/*.png
