#include <math.h>

#include "matrix.h"
#include "drying.h"
#include "material-data.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */
#define NTERMS 100

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
double strainpc(double t, double x, drydat d)
{
    int i,
        nt = 1000; /* Number of time steps to use */
    double Xdb, 
           Pc,
           e = 0,
           dt = t/nt; /* Size of each time step */
    burgerse *b;

    b = CreateBurgersE();
    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, d);
        //Pc = pore_press(Xdb, T) - pore_press(.3, T);
        /* Pore pressure */
        Pc = pore_press(Xdb, d.T);
        /* Use a modified integral formula to calculate strain. This has been
         * integrated by parts to eliminate the numerical error associated with
         * approximating the pressure time derivative. */
        e += DBurgersECreep(b, t-i*dt, d.T, Xdb, -1*Pc) * Pc  * dt;
    }

    Xdb = CrankEquationFx(x, t, d);
    Pc = pore_press(Xdb, d.T);
    /* The other part of the integration formula */
    e += BurgersECreep(b, 0, d.T, Xdb, -1*Pc)*Pc;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return e;
}

double MaxwellStrainPc(double t, double x, drydat d)
{
    int i,
        nt = 1000; /* Number of time steps to use */
    double Xdb, 
           Pc,
           e = 0,
           dt = t/nt; /* Size of each time step */

    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, d);
        /* Pore pressure */
        Pc = pore_press(Xdb, d.T);
        /* Use a modified integral formula to calculate strain. This has been
         * integrated by parts to eliminate the numerical error associated with
         * approximating the pressure time derivative. */
        e += DMaxwellCreepConverted(t-i*dt, d.T, Xdb) * Pc  * dt;
    }

    Xdb = CrankEquationFx(x, t, d);
    Pc = pore_press(Xdb, d.T);
    /* The other part of the integration formula */
    e += MaxwellCreepConverted(0, d.T, Xdb)*Pc;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return .06*e;
}

double ZhuMaxwellStrain(double t, double x, drydat d)
{
    int i,
        nt = 1000; /* Number of time steps to use */
    double Bc = 8.2e-14,
           Xdb, 
           Pc,
           e = 0,
           dt = t/nt; /* Size of each time step */
    maxwell *m;

    m = CreateMaxwellZhu();
    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, d);
        //Pc = pore_press(Xdb, T) - pore_press(.3, T);
        /* Pore pressure */
        Pc = pore_press(Xdb, d.T);
        /* Use a modified integral formula to calculate strain. This has been
         * integrated by parts to eliminate the numerical error associated with
         * approximating the pressure time derivative. */
        e += DMaxwellCreep(m, t-i*dt, d.T, Xdb) * Pc  * dt;
    }

    Xdb = CrankEquationFx(x, t, d);
    Pc = pore_press(Xdb, d.T);
    /* The other part of the integration formula */
    e += MaxwellCreep(m, 0, d.T, Xdb)*Pc;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return 10000*e*Bc;
}

double EqStrainPc(double t, double x, drydat d)
{
    double Xdb, 
           Pc,
           Ea,
           e = 0;

    /* Calculate moisture content using the Crank equation */
    Xdb = CrankEquationFx(x, t, d);
    //if(x>=9.4e-4)
    //    printf("x = %g, t = %g, Xdb = %g\n", x, t, Xdb);
    /* Pore pressure */
    Pc = pore_press(Xdb, d.T); 
    Ea = 68.18*(1/(1+exp((Xdb-250.92*exp(-0.0091*d.T))/2.19))+0.078) * 1e6;
    e = Pc/Ea;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return e*.06;
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

vector* displacementV(vector* strain, double L)
{
    int i;
    vector *u;
    u = CreateVector(len(strain));
    for(i=0; i<len(strain); i++)
        setvalV(u, i, displacement(i, strain, L));
    return u;
}
