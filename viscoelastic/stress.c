#include <math.h>

#include "matrix.h"
#include "drying.h"
#include "material-data.h"
#include "visco.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */
#define NTERMS 100

double RelaxCummings(double t, double T, double X, int deriv)
{
    double G;
    maxwell *m;
    m = CreateMaxwell();
    if(deriv)
        G = DMaxwellRelax(m, t, T, X);
    else
        G = MaxwellRelax(m, t, T, X);
    DestroyMaxwell(m);
    return G;
}

double RelaxLaura(double t, double T, double X, int deriv)
{
    double J;
    if(deriv)
        J = DMaxwellRelaxLaura(t, T, X);
    else
        J = MaxwellRelaxLaura(t, T, X);
    return J;
}

double RelaxZhu(double t, double T, double X, int deriv)
{
    double K;
    maxwell *m;
    m = CreateMaxwellZhu();
    if(deriv)
        K = DMaxwellRelax(m, t, T, X);
    else
        K = MaxwellRelax(m, t, T, X);
    DestroyMaxwell(m);
    return K;
}

/**
 * Calculate the strain assuming that capillary pressure is the driving
 * force for shrinkage. Capillary pressure is determined from the Kelvin
 * equation, and the water activity is determined from the Oswin isotherm
 * model. The creep compliance function is from Gina's thesis, and has been
 * fitted to take into account the nonlinear variation with applied stress.
 *
 * @param t Simulation time [s]
 * @param x Length coordinate [m]
 * @param RH Relative humidity of the surrounding air. [-]
 * @param X0 Initial moisture content of the sample [kg/kg db]
 * @param Xe Equilibrium moisture content of the sample based on the supplied
 *      relative humidity
 * @param L Sample thickness [m]
 * @param T Drying temperature [K]
 *
 * @returns Infintesimal strain [-]
 */
double VEStress(double x, double t, drydat d,
                double (*E)(double, double, drydat, double),
                double (*G)(double, double, double, int))
{
    int i,
        nt = 100; /* Number of time steps to use */
    double Xdb, 
           s = 0,
           e,
           eta = .2,
           dt = t/nt; /* Size of each time step */
    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, d);
        e = E(i*dt, x, d, .2);

        s += G(t-i*dt, d.T, Xdb, 1) * e * dt;
    }

    Xdb = CrankEquationFx(x, t, d);
    e = E(t, x, d, eta);
    s += G(0, d.T, Xdb, 0)*e;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return s;
}

