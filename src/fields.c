#include <math.h>
#include "fields.h"

void dipole_field(double x, double y, double z, const DipoleParams *p,
                   double *Bx, double *By, double *Bz) {
    double r2 = x * x + y * y + z * z;
    double r = sqrt(r2);
    double r5 = r2 * r2 * r; /* r^5 */
    double M = p->M;

    *Bx = M * 3.0 * x * z / r5;
    *By = M * 3.0 * y * z / r5;
    *Bz = M * (3.0 * z * z - r2) / r5;
}
