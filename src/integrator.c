#include "integrator.h"

/*
 * scratch layout (caller allocates 5*dim doubles):
 *   scratch[0*dim .. 1*dim)  = k1
 *   scratch[1*dim .. 2*dim)  = k2
 *   scratch[2*dim .. 3*dim)  = k3
 *   scratch[3*dim .. 4*dim)  = k4
 *   scratch[4*dim .. 5*dim)  = tmp (temporary state for midpoint/endpoint evals)
 */
void rk4_step(DerivFunc deriv, double t, double dt, double *state, int dim,
              void *params, double *scratch) {
    double *k1  = scratch + 0 * dim;
    double *k2  = scratch + 1 * dim;
    double *k3  = scratch + 2 * dim;
    double *k4  = scratch + 3 * dim;
    double *tmp = scratch + 4 * dim;

    int i;

    /* k1 = f(t, state) */
    deriv(t, state, k1, params);

    /* k2 = f(t + dt/2, state + dt/2 * k1) */
    for (i = 0; i < dim; i++) tmp[i] = state[i] + 0.5 * dt * k1[i];
    deriv(t + 0.5 * dt, tmp, k2, params);

    /* k3 = f(t + dt/2, state + dt/2 * k2) */
    for (i = 0; i < dim; i++) tmp[i] = state[i] + 0.5 * dt * k2[i];
    deriv(t + 0.5 * dt, tmp, k3, params);

    /* k4 = f(t + dt, state + dt * k3) */
    for (i = 0; i < dim; i++) tmp[i] = state[i] + dt * k3[i];
    deriv(t + dt, tmp, k4, params);

    /* state += dt/6 * (k1 + 2*k2 + 2*k3 + k4) */
    for (i = 0; i < dim; i++) {
        state[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}
