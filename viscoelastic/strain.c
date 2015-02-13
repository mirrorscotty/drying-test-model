#include <math.h>

#include "matrix.h"
#include "drying.h"
#include "material-data.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */
#define NTERMS 100

double CreepGina(double t, double T, double X, double P, int deriv)
{
    double J;
    burgerse *b;
    b = CreateBurgersE();
    if(deriv)
        J = DBurgersECreep(b, t, T, X, P);
    else
        J = BurgersECreep(b, t, T, X, P);
    DestroyBurgersE(b);
    return J;
}

double CreepLaura(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DMaxwellCreepConverted(t, T, X);
    else
        J = MaxwellCreepConverted(t, T, X);
    return J;
}

double CreepLaura2(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DMaxwellCreepLaura(t, T, X);
    else
        J = MaxwellCreepLaura(t, T, X);
    return J;
}

double CreepZhu(double t, double T, double X, double P, int deriv)
{
    double J, Bc = 8.2e-14;
    maxwell *m;
    m = CreateMaxwellZhu();
    if(deriv)
        J = DMaxwellCreep(m, t, T, X);
    else
        J = MaxwellCreep(m, t, T, X);
    DestroyMaxwell(m);
    return J*1e0*Bc;
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
double strainpc(double x, double t, drydat d,
                double (*F)(double, double, drydat),
                double (*J)(double, double, double, double, int))
{
    int i,
        nt = 1000; /* Number of time steps to use */
    double Xdb, 
           P,
           e = 0,
           dt = t/nt; /* Size of each time step */

    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, d);
        P = F(x, t, d);
        /* Use a modified integral formula to calculate strain. This has been
         * integrated by parts to eliminate the numerical error associated with
         * approximating the pressure time derivative. */
        e += J(t-i*dt, d.T, Xdb, -1*P, 1) * P  * dt;
    }

    /* The other part of the integration formula */
    Xdb = CrankEquationFx(x, t, d);
    e += -1*J(t, d.T, Xdb, -1*P, t)*F(x, 0, d);
    Xdb = CrankEquationFx(x, 0, d);
    e += J(0, d.T, Xdb, -1*P, 0)*F(x, t, d);

    return e*6*(.5-.35);
}

/**
 * Calculate the equilibrium strain using Laura's creep compliance data.
 */
double EqStrainPc(double t, double x, drydat d)
{
    double Xdb, 
           Pc,
           Ea,
           e = 0;

    /* Calculate moisture content using the Crank equation */
    Xdb = CrankEquationFx(x, t, d);
    /* Pore pressure */
    Pc = pore_press(Xdb, d.T); 
    Ea = 68.18*(1/(1+exp((Xdb-250.92*exp(-0.0091*d.T))/2.19))+0.078) * 1e6;
    e = Pc/Ea;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return e;
}

/**
 * Integrate strain to find the displacement at any point throughout the slab.
 * @param xi Index of the x coordinate to find the displacement of.
 * @param strain Vector of strain values for the whole slab.
 * @param L Slab thickness [m]
 *
 * @returns Displacement [m]
 */
double displacement(int xi, vector* strain, double L)
{
    int i; /* Loop index */
    double dx = L/len(strain), /* Delta x between points */
           u = 0; /* Displacement */
    /* Integrate strain from x=0 to x=xi*dx with respect to x. */
    for(i=0; i<xi; i++) {
        u += valV(strain, i) * dx;
    }

    return u;
}

/**
 * Calculate a vector of displacement values at each location.
 * @param strain Vector of strain values across the thickness of the sample. [-]
 * @param L Sample thickness [m]
 * @returns Displacement [m]
 */
vector* displacementV(vector* strain, double L)
{
    int i;
    vector *u;
    u = CreateVector(len(strain));
    for(i=0; i<len(strain); i++)
        setvalV(u, i, displacement(i, strain, L));
    return u;
}
