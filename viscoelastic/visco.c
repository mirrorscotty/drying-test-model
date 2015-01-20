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
 * Calculate the strain assuming that binding energy gradient is the driving
 * force for shrinkage. Binding energy is calculated using the Clausius-
 * Cleyperon equation from supplied moisture profile. The time derivative of
 * binding energy gradient is then numerically integrated to calculate
 * viscoelastic strain.
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
double strain(double t, double x, double RH, double D, double X0, double Xe, double L, double T)
{
    int i, /* Loop index */
        nt = 100; /* Number of time steps to use */
    double Xdb1, Xdb2,
           Xdb3, Xdb4,
           Eb1, Eb2,
           Eb3, Eb4,
           GradEb, GradEbph,
           h = 1e-5, /* Used for numeric differentiation */
           dt = t/nt, /* Time step size */
           e = 0;

    oswin *d;
    maxwell *m;

    /* Point 1: x-h, t-h
     * Point 2: x+h, t-h
     * Point 3: x-h, t+h
     * Point 4: x+h, t+h
     */

    d = CreateOswinData();
    m = CreateMaxwell();
    
    for(i=0; i<nt; i++) {
        /* Calculate Xdb at each of four points */
        Xdb1 = CrankEquationFx(x-h, t-h, L, D, Xe, X0, NTERMS);
        Xdb2 = CrankEquationFx(x+h, t-h, L, D, Xe, X0, NTERMS);
        Xdb3 = CrankEquationFx(x-h, t+h, L, D, Xe, X0, NTERMS);
        Xdb4 = CrankEquationFx(x+h, t+h, L, D, Xe, X0, NTERMS);

        /* Do the same for binding energy */
        Eb1 = BindingEnergyOswin(d, Xdb1, T);
        Eb2 = BindingEnergyOswin(d, Xdb2, T);
        Eb3 = BindingEnergyOswin(d, Xdb3, T);
        Eb4 = BindingEnergyOswin(d, Xdb4, T);

        /* Calculate binding energy gradient at t-h and t+h */
        GradEb = (Eb2 - Eb1)/(2*h);
        GradEbph = (Eb4 - Eb3)/(2*h);

        /* Integrate to find creep */
        e += MaxwellCreep(m, t-i*dt, T, Xdb1) * (GradEbph-GradEb)/(2*h) * dt;
    }

    /* Debugging stuff */
    puts("----------------------------------------");
    printf("Xdb(x-h, t-h) = %g, Xdb(x+h, t-h) = %g\n", Xdb1, Xdb2);
    printf("Xdb(x-h, t+h) = %g, Xdb(x+h, t+h) = %g\n", Xdb3, Xdb4);
    printf("Eb(x-h, t-h) = %g, Eb(x+h, t-h) = %g\n", Eb1, Eb2);
    printf("Eb(x-h, t+h) = %g, Eb(x+h, t+h) = %g\n", Eb3, Eb4);
    printf("GradEb(t+h) = %g, GradEb(t-h) = %g\n",
            GradEbph, GradEb);
    printf("DGradEbDt = %g\n", (GradEbph-GradEb)/(2*h));

    return e;
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
double strainpc(double t, double x, double RH, double D, double X0, double Xe, double L, double T)
{
    int i,
        nt = 100; /* Number of time steps to use */
    double Xdb, 
           Pc,
           e = 0,
           dt = t/nt; /* Size of each time step */
    burgerse *b;

    b = CreateBurgersE();
    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, L, D, Xe, X0, NTERMS);
        //Pc = pore_press(Xdb, T) - pore_press(.3, T);
        /* Pore pressure */
        Pc = pore_press(Xdb, T);
        /* Use a modified integral formula to calculate strain. This has been
         * integrated by parts to eliminate the numerical error associated with
         * approximating the pressure time derivative. */
        e += DBurgersECreep(b, t-i*dt, T, Xdb, -1*Pc) * Pc  * dt;
    }

    Xdb = CrankEquationFx(x, t, L, D, Xe, X0, NTERMS);
    Pc = pore_press(Xdb, T);
    /* The other part of the integration formula */
    e += BurgersECreep(b, 0, T, Xdb, -1*Pc)*Pc;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return 0.06*e;
}

double MaxwellStrainPc(double t, double x, double RH, double D, double X0, double Xe, double L, double T)
{
    int i,
        nt = 1000; /* Number of time steps to use */
    double Xdb, 
           Pc,
           e = 0,
           dt = t/nt; /* Size of each time step */

    for(i=0; i<nt; i++) {
        /* Calculate moisture content using the Crank equation */
        Xdb = CrankEquationFx(x, i*dt, L, D, Xe, X0, NTERMS);
        /* Pore pressure */
        Pc = pore_press(Xdb, T);
        /* Use a modified integral formula to calculate strain. This has been
         * integrated by parts to eliminate the numerical error associated with
         * approximating the pressure time derivative. */
        e += DMaxwellCreepConverted(t-i*dt, T, Xdb) * Pc  * dt;
    }

    Xdb = CrankEquationFx(x, t, L, D, Xe, X0, NTERMS);
    Pc = pore_press(Xdb, T);
    /* The other part of the integration formula */
    e += MaxwellCreepConverted(0, T, Xdb)*Pc;

    /* Multiply strain (or, more accurately, stress) by porosity to get
     * effective stress (hopefully) */
    return e*.06;
}

double EqStrainPc(double t, double x, double RH, double D, double X0, double Xe, double L, double T)
{
    double Xdb, 
           Pc,
           Ea,
           e = 0;

    /* Calculate moisture content using the Crank equation */
    Xdb = CrankEquationFx(x, t, L, D, Xe, X0, NTERMS);
    /* Pore pressure */
    Pc = pore_press(Xdb, T); 
    Ea = 68.18*(1/(1+exp((Xdb-250.92*exp(-0.0091*T))/2.19))+0.078) * 1e6;
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

