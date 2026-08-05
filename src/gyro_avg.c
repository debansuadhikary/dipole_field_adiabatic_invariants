#include <math.h>
#include "gyro_avg.h"

void gyro_avg_init(GyroAvg *g) {
    g->head = 0;
    g->count = 0;
    g->sum = 0.0;
}

double gyro_avg_push(GyroAvg *g, double sample, double T_local, double dt) {
    /* how many samples should the window hold right now */
    int target_window = (int)lround(T_local / dt);
    if (target_window < 1) target_window = 1;
    if (target_window > GYRO_AVG_MAX_WINDOW) target_window = GYRO_AVG_MAX_WINDOW;

    /* push the new sample */
    g->buf[g->head] = sample;
    g->sum += sample;
    g->head = (g->head + 1) % GYRO_AVG_MAX_WINDOW;
    g->count++;

    /* trim oldest samples until we're back down to target_window */
    while (g->count > target_window) {
        int tail = (g->head - g->count + GYRO_AVG_MAX_WINDOW) % GYRO_AVG_MAX_WINDOW;
        g->sum -= g->buf[tail];
        g->count--;
    }

    return g->sum / g->count;
}
