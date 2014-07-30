#include "pasta.h"
#include "matrix.h"
#include "drying.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */
#define NTERMS 40

/**
 * Calculate the strain assuming that binding energy gradient is the driving
 * force for shrinkage. Binding energy is calculated using the Clausius-
 * Cleyperon equation from supplied moisture profile.
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
    double Xdb1, Xdb2,
           Eb1, Eb2,
           GradEb,
           h = 1e-5, /* Used for numeric differentiation */
           e = 0,
           E = 5e9; /* From Guinea et al. 2004 */

    oswin *d;

    /* Point 1: x-h, t
     * Point 2: x+h, t
     */

    d = CreateOswinData();
    
    /* Calculate the moisture content at two points */
    Xdb1 = CrankEquationFx(x-h, t, L, D, Xe, X0, NTERMS);
    Xdb2 = CrankEquationFx(x+h, t, L, D, Xe, X0, NTERMS);

    /* Determine the binding energies at those two points */
    Eb1 = BindingEnergyOswin(d, Xdb1, T);
    Eb2 = BindingEnergyOswin(d, Xdb2, T);

    /* Take the derivative using the central difference method */
    GradEb = (Eb2 - Eb1)/(2*h);
    /* Calculate elastic strain */
    e = -GradEb / E;

    return e;
}

/**
 * Calculate the strain due to capillary pressure.
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
    double Xdb, Pc,
           e = 0,
           E = 5e9; /* From Guinea et al. 2004 */

    /* Determine moisture content */
    Xdb = CrankEquationFx(x, t, L, D, Xe, X0, NTERMS);

    /* Calculate pore pressure */
    Pc = pore_press(Xdb, T);
    /* Find elastic strain */
    e = Pc / E;

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
                    * into. */
        i; /* Loop index */
    
    vector *z, *Xdb, *Eb, *GradEb, *str, *strpc, *u, *upc;
    matrix *output;

    /* Print out a usage statement if an incorrect number of arguments is
     * supplied. */
    if(argc != 3) {
        puts("Usage:");
        puts("dry <aw> <t>");
        puts("aw: Final water activity to dry to.");
        puts("t: Time to output a profile for.");
        exit(0);
    }

    /* Store the command line arguments in the appropriate variables. */
    RH = atof(argv[1]);
    t = atof(argv[2]);

    /* Calculate equilibrium moisture content and average diffusivity */
    d = CreateOswinData();
    Xe = OswinIsotherm(d, RH, T);
    D = DiffCh10((Xe+X0)/2, T);

    /* Make a bunch of vectors */
    z = CreateVector(npts);
    Xdb = CreateVector(npts);
    Eb = CreateVector(npts);
    GradEb = CreateVector(npts);
    str = CreateVector(npts);
    strpc = CreateVector(npts);
    u = CreateVector(npts);
    upc = CreateVector(npts);

    /* Print out the equilibrium moisture content for fun */
    printf("Xe = %g\n", Xe);
    
    /* Set up the domain by assigning the z value of each point to a slot in
     * the vector */
    for(i=0; i<npts; i++)
        setvalV(z, i, L/(npts-1) * i);

    for(i=0; i<npts; i++) {
        /* Moisture content */
        setvalV(Xdb, i, CrankEquationFx(valV(z, i), t, L, D, Xe, X0, NTERMS));
        /* Binding energy */
        setvalV(Eb, i, BindingEnergyOswin(d, valV(Xdb, i), T));
        /* Calculate the binding energy at x+h so that it can be used to
         * determine binding energy gradient */
        Xdbph = CrankEquationFx(valV(z,i)+h, t, L, D, Xe, X0, NTERMS);
        /* Binding energy gradient */
        setvalV(GradEb, i, (BindingEnergyOswin(d, Xdbph, T) - BindingEnergyOswin(d, valV(Xdb, i), T))/h);
        /* Strain due to binding energy gradient */
        setvalV(str, i, strain(t, valV(z, i), RH, D, X0, Xe, L, T));
        /* Strain due to capillary pressure */
        setvalV(strpc, i, strainpc(t, valV(z, i), RH, D, X0, Xe, L, T));
    }
    /* Find the displacement at all points. */
    for(i=0;i<npts;i++) {
        setvalV(u, i, displacement(i, str, L));
        setvalV(upc, i, displacement(i, strpc, L));
    }

    /* Smash the desired vectors together and print the resulting matrix */
    output = CatColVector(5, z, Xdb, GradEb, strpc, upc);
    mtxprntfilehdr(output, "out.csv", "z [m],Xdb [kg/kg db],GradEb,e,u [m]\n");

    return 0;
}

