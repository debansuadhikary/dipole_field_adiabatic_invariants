#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integrator.h"
#include "fields.h"
#include "gyro_avg.h"

/*
 * Charged particle state: {x, y, z, vx, vy, vz}
 * Equation of motion: dv/dt = (q/m) * (v x B)
 *
 * We carry q/m as a single ratio "qm" rather than separate q and m —
 * that's all the Lorentz force actually depends on.
 */
typedef struct {
    double qm;         /* charge-to-mass ratio */
    double m;          /* mass — kept explicit (=1 by default) so mu's
                           formula reads the same as the physics: only the
                           equation of motion needs q/m, mu needs m alone */
    DipoleParams dip;  /* dipole field parameters */
} ParticleParams;

/*
 * Diagnostic: first adiabatic invariant mu = m*v_perp^2 / (2*|B|), plus the
 * local field magnitude (useful on its own for sanity checks / later
 * bounce-point detection). v_perp is the component of v perpendicular to
 * the *local* B direction — NOT just vx,vy, since B's direction rotates
 * as the particle moves off the equatorial plane.
 *
 * This does an extra dipole_field() call at the *current* state, separate
 * from the ones rk4_step() does internally at intermediate RK4 stages —
 * we want mu evaluated exactly at the state we're about to record, not at
 * some k2/k3 midpoint.
 */
static void compute_mu(const double *state, const ParticleParams *pp,
                        double *mu_out, double *Bmag_out) {
    double x = state[0], y = state[1], z = state[2];
    double vx = state[3], vy = state[4], vz = state[5];
    double Bx, By, Bz;

    dipole_field(x, y, z, &pp->dip, &Bx, &By, &Bz);

    double Bmag = sqrt(Bx * Bx + By * By + Bz * Bz);
    double bx = Bx / Bmag, by = By / Bmag, bz = Bz / Bmag; /* unit vector b-hat */

    double vpar = vx * bx + vy * by + vz * bz;             /* v . b-hat */
    double vperp_x = vx - vpar * bx;
    double vperp_y = vy - vpar * by;
    double vperp_z = vz - vpar * bz;
    double vperp2 = vperp_x * vperp_x + vperp_y * vperp_y + vperp_z * vperp_z;

    *mu_out = pp->m * vperp2 / (2.0 * Bmag);
    *Bmag_out = Bmag;
}

static void particle_deriv(double t, const double *state, double *dstate, void *params_v) {
    (void)t;
    ParticleParams *pp = (ParticleParams *)params_v;
    double x = state[0], y = state[1], z = state[2];
    double vx = state[3], vy = state[4], vz = state[5];
    double Bx, By, Bz;

    dipole_field(x, y, z, &pp->dip, &Bx, &By, &Bz);

    /* dr/dt = v */
    dstate[0] = vx;
    dstate[1] = vy;
    dstate[2] = vz;

    /* dv/dt = qm * (v x B) */
    dstate[3] = pp->qm * (vy * Bz - vz * By);
    dstate[4] = pp->qm * (vz * Bx - vx * Bz);
    dstate[5] = pp->qm * (vx * By - vy * Bx);
}

int main(int argc, char **argv) {
    /* Defaults — all overridable via command line so we can sweep dt later
       without recompiling (see project plan: dt sweep is the key rigor result) */
    double dt = 0.001;
    long nsteps = 2000000;
    long output_stride = 200; /* write every Nth step, else the CSV is huge */
    const char *outpath = "trajectory.csv";

    if (argc > 1) dt = atof(argv[1]);
    if (argc > 2) nsteps = atol(argv[2]);
    if (argc > 3) output_stride = atol(argv[3]);
    if (argc > 4) outpath = argv[4];

    ParticleParams pp;
    pp.qm = 1.0;
    pp.m = 1.0;
    pp.dip.M = 1000.0;

    /* Initial conditions: start in the equatorial plane at L=5, where the
       field points purely along -z (Bz = -M/L^3 there, Bx=By=0), so
       "parallel" and "perpendicular" are simply z vs. the xy-plane at t=0.
       vy = perpendicular (gyration) speed, vz = parallel (bounce) speed. */
    double state[6] = {
        5.0, 0.0, 0.0,   /* x, y, z */
        0.0, 1.0, 0.5    /* vx, vy, vz */
    };

    double scratch[6 * 5];

    FILE *f = fopen(outpath, "w");
    if (!f) {
        fprintf(stderr, "error: could not open %s for writing\n", outpath);
        return 1;
    }
    fprintf(f, "t,x,y,z,vx,vy,vz,mu,mu_avg,Bmag\n");

    /* Record mu at t=0 so every later value in analysis can be normalized
       as mu(t)/mu(0) - the actual "how well is it conserved" metric. */
    double mu0, Bmag0;
    compute_mu(state, &pp, &mu0, &Bmag0);
    fprintf(stderr, "initial mu = %.8f (Bmag=%.6f)\n", mu0, Bmag0);

    GyroAvg gavg;
    gyro_avg_init(&gavg);

    long i;
    for (i = 0; i <= nsteps; i++) {
        double t = i * dt;

        /* mu must be pushed every step (not just on output rows) so the
           averaging window actually contains a full local gyro-period's
           worth of samples, regardless of output_stride. */
        double mu, Bmag;
        compute_mu(state, &pp, &mu, &Bmag);
        double T_local = 2.0 * M_PI / (pp.qm * Bmag); /* local gyro-period */
        double mu_avg = gyro_avg_push(&gavg, mu, T_local, dt);

        if (i % output_stride == 0) {
            fprintf(f, "%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.10f,%.10f,%.8f\n",
                    t, state[0], state[1], state[2],
                    state[3], state[4], state[5], mu, mu_avg, Bmag);
        }
        if (i < nsteps) {
            rk4_step(particle_deriv, t, dt, state, 6, &pp, scratch);
        }
    }

    fclose(f);
    fprintf(stderr, "wrote %ld points to %s (dt=%.6f, nsteps=%ld)\n",
            (nsteps / output_stride) + 1, outpath, dt, nsteps);
    return 0;
}
