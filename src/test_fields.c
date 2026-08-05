#include <stdio.h>
#include <math.h>
#include "fields.h"

int main(void) {
    DipoleParams p = {.M = 1.0};
    double Bx, By, Bz;

    /* On-axis (x=0,y=0,z): analytic B = 2M/z^3 along z (pure axial field) */
    double z = 2.0;
    dipole_field(0.0, 0.0, z, &p, &Bx, &By, &Bz);
    double B_onaxis_analytic = 2.0 * p.M / (z * z * z);
    printf("on-axis z=%.2f: numeric=(%.6f,%.6f,%.6f)  analytic Bz=%.6f\n",
           z, Bx, By, Bz, B_onaxis_analytic);

    /* Equatorial (z=0, r=x): analytic B = -M/r^3 along z (points opposite moment) */
    double x = 2.0;
    dipole_field(x, 0.0, 0.0, &p, &Bx, &By, &Bz);
    double B_eq_analytic = -p.M / (x * x * x);
    printf("equatorial x=%.2f: numeric=(%.6f,%.6f,%.6f)  analytic Bz=%.6f\n",
           x, Bx, By, Bz, B_eq_analytic);

    /* Sanity: |B| should fall off faster on-axis than equatorial is wrong way
       to check -- instead confirm ratio B_onaxis/B_eq = 2 at equal distance,
       a standard dipole property. */
    printf("ratio (on-axis Bz) / (equatorial |Bz|) at equal r = %.6f (expect 2.0)\n",
           fabs(B_onaxis_analytic) / fabs(B_eq_analytic));

    return 0;
}
