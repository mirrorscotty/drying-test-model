#include "pasta.h"
#include "matrix.h"
#include "drying.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */
#define NTERMS 1000

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

    //Pc = pore_press(Xdb, T) - pore_press(.3, T);
    Pc = pore_press(Xdb, T);
    /* The other part of the integration formula */
    e += BurgersECreep(b, 0, T, Xdb, -1*Pc)*Pc;

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

int main(int argc, char *argv[])
{
    double L = 1e-3, /* Length [m] */
           D, /* Diffusivity [m^2/s] */
           X0 = .3, /* Initial moisture content */
           Xe, /* Equilibrium moisture content */
           T = 60+273.15, /* Drying temperature */
           h = 1e-7,
           Xdbph,
           t,
           RH;

    oswin *d; /* Isotherm data */

    int npts = 51, /* Number of points to subdivide the thickness of the slab
                      into */
        i; /* Loop index */
    
    vector *z, *Xdb, *Eb, *GradEb, *str, *u, *Pc;
    matrix *output;

    /* Print out a usage statement if needed */
    if(argc != 3) {
        puts("Usage:");
        puts("dry <aw> <t>");
        puts("aw: Final water activity to dry to.");
        puts("t: Time to output a profile for.");
        exit(0);
    }

    /* Store command line arguments */
    RH = atof(argv[1]);
    t = atof(argv[2]);

    /* Calculate equilibrium moisture content and average diffusivity */
    d = CreateOswinData();
    Xe = OswinIsotherm(d, RH, T);
    D = DiffCh10((Xe+X0)/2, T);

    /* Make a bunch of vectors */
    z = CreateVector(npts);
    Xdb = CreateVector(npts);
    Pc = CreateVector(npts);
    Eb = CreateVector(npts);
    GradEb = CreateVector(npts);
    str = CreateVector(npts);
    u = CreateVector(npts);

    /* Print equilibrium moisture content */
    printf("Xe = %g\n", Xe);

    /* Set up the domain by assigning the z value of each point to a slot in
     * the vector */
    for(i=0; i<npts; i++)
        setvalV(z, i, L/npts * i);

    for(i=0; i<npts; i++) {
        /* Moisture content */
        setvalV(Xdb, i, CrankEquationFx(valV(z, i), t, L, D, Xe, X0, NTERMS));
        /* Binding energy */
        setvalV(Eb, i, BindingEnergyOswin(d, valV(Xdb, i), T));
        /* Calculate binding energy at x+h so that it can be used to determine
         * the binding energy gradient */
        Xdbph = CrankEquationFx(valV(z,i)+h, t, L, D, Xe, X0, NTERMS);
        /* Binding energy gradient */
        setvalV(GradEb, i, (BindingEnergyOswin(d, Xdbph, T) - BindingEnergyOswin(d, valV(Xdb, i), T))/h);
        /* Capillary Pressure */
        setvalV(Pc, i, pore_press(valV(Xdb, i), T));
        /* Strain */
        setvalV(str, i, strainpc(t, valV(z, i), RH, D, X0, Xe, L, T));
    }
    /* Find the displacement at all points */
    for(i=0;i<npts;i++)
        setvalV(u, i, displacement(i, str, L));

    /* Smash together the desired vectors and print the resulting matrix */
    output = CatColVector(5, z, Xdb, Pc, str, u);
    mtxprntfilehdr(output, "out.csv","x [m],Xdb [kg/kg db],Pc [Pa],strain,u [m]\n");

    return 0;
}

