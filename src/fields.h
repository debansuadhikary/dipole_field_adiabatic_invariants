#ifndef FIELDS_H
#define FIELDS_H

/*
 * Axial magnetic dipole field, moment aligned with +z (like Earth's
 * approximate field, or the core field of a magnetic mirror/trap).
 *
 *   B(r) = (M / r^5) * (3xz, 3yz, 3z^2 - r^2)
 *
 * where r = sqrt(x^2+y^2+z^2), and M is the dipole strength constant
 * (rolled together with mu0/4pi so we don't carry SI constants around —
 * pick M so field magnitudes near your particle's orbit are O(1) in
 * your chosen units).
 *
 * NOTE: this diverges at r=0. Never place a particle at the origin.
 */
typedef struct {
    double M; /* dipole strength */
} DipoleParams;

/* Fills Bx,By,Bz with the field at (x,y,z) for the given dipole. */
void dipole_field(double x, double y, double z, const DipoleParams *p,
                   double *Bx, double *By, double *Bz);

#endif
