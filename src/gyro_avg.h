#ifndef GYRO_AVG_H
#define GYRO_AVG_H

/*
 * Running gyro-average of a scalar (mu), using a window whose length
 * adapts to the *current* local gyro-period. This filters out the
 * gyro-frequency ripple in the raw per-step mu (a real finite-Larmor-
 * radius effect in non-uniform fields — see main.c comments) so we can
 * see the actual secular trend: how well mu is conserved over bounce
 * and drift timescales, which is the physically meaningful question.
 *
 * The window is a sliding average over the *last* N samples, where N is
 * recomputed every push() from the current gyro-period / dt. Because the
 * field varies slowly along the trajectory (that's the adiabatic
 * assumption we're testing!), N itself changes slowly, so trimming the
 * buffer by a few samples per step is enough to track it — no need to
 * rebuild the buffer from scratch.
 */

#define GYRO_AVG_MAX_WINDOW 200000

typedef struct {
    double buf[GYRO_AVG_MAX_WINDOW];
    int head;       /* index where the next sample will be written */
    int count;      /* number of valid samples currently in the window */
    double sum;      /* running sum of samples currently in the window */
} GyroAvg;

void gyro_avg_init(GyroAvg *g);

/*
 * Push one new raw sample, taken at gyro-period T_local (in the same time
 * units as dt) and step size dt. Returns the current windowed average.
 * target_window = round(T_local / dt), clamped to [1, GYRO_AVG_MAX_WINDOW].
 */
double gyro_avg_push(GyroAvg *g, double sample, double T_local, double dt);

#endif
