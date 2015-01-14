#include "material-data.h"
#include "matrix.h"
#include "drying.h"
#include "visco.h"

#include <stdlib.h>
#include <stdio.h>

#define NTERMS 100

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
    
    vector *z, *Xdb, *Eb, *GradEb, *str, *u, *Pc, *Mstr, *Mu;
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
    Mstr = CreateVector(npts);
    u = CreateVector(npts);
    Mu = CreateVector(npts);

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
        setvalV(Mstr, i, MaxwellStrainPc(t, valV(z, i), RH, D, X0, Xe, L, T));
    }
    /* Find the displacement at all points */
    for(i=0;i<npts;i++)
        setvalV(u, i, displacement(i, str, L));
        setvalV(Mu, i, displacement(i, Mstr, L));

    /* Smash together the desired vectors and print the resulting matrix */
    output = CatColVector(7, z, Xdb, Pc, str, u, Mstr, Mu);
    mtxprntfilehdr(output, "out.csv","x [m],Xdb [kg/kg db],Pc [Pa],strain,u [m],Maxwell Strain,Maxwell Displacement\n");

    return 0;
}

