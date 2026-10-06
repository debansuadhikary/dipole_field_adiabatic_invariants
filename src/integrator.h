#ifndef INTEGRATOR_H
#define INTEGRATOR_H

/* Generic RK4 integrator.
   deriv signature:
    void deriv(double t, const double *state, double *dstate, void *params)
 
    t       - current time
    state   - current state vector, length `dim`
    dstate  - OUTPUT: derivative of state at (t, state), length `dim`
    params  - opaque pointer to whatever extra data deriv needs
              (e.g. charge, mass, field strength) — avoids global variables
 */

typedef void (*DerivFunc)(double t, const double *state, double *dstate, void *params);

/*  
Advance `state` in place by one RK4 step of size `dt`. `dim` is the length of the state vector. `scratch` must point to a caller-allocated buffer of size 4*dim doubles (k1,k2,k3,k4) plus dim doubles for the temporary state — see .c file for exact layout. Passing scratch in avoids malloc/free every step.
*/
void rk4_step(DerivFunc deriv, double t, double dt, double *state, int dim,
              void *params, double *scratch);

#endif
