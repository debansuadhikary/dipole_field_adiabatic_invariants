#include <stdio.h>
#include <math.h>
#include "integrator.h"

/*
 * Simple harmonic oscillator: x'' = -omega^2 * x
 * state = {x, v}, dstate = {v, -omega^2 * x}
 * params = pointer to a single double: omega^2
 */
void sho_deriv(double t, const double *state, double *dstate, void *params) {
    (void)t;
    double omega2 = *(double *)params;
    dstate[0] = state[1];
    dstate[1] = -omega2 * state[0];
}

int main(void) {
    double omega = 2.0;
    double omega2 = omega * omega;
    double state[2] = {1.0, 0.0}; /* x(0)=1, v(0)=0 -> x(t) = cos(omega t) */
    double scratch[2 * 5];
    double dt = 0.001;
    int nsteps = 10000; /* t_final = 10 */
    int i;

    printf("t,x_numeric,x_analytic,error,energy\n");
    for (i = 0; i <= nsteps; i++) {
        double t = i * dt;
        double x_analytic = cos(omega * t);
        double energy = 0.5 * state[1] * state[1] + 0.5 * omega2 * state[0] * state[0];
        if (i % 1000 == 0) {
            printf("%f,%f,%f,%e,%f\n", t, state[0], x_analytic,
                   fabs(state[0] - x_analytic), energy);
        }
        rk4_step(sho_deriv, t, dt, state, 2, &omega2, scratch);
    }
    return 0;
}
