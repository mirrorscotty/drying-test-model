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
           t,
           RH;
    drydat cond;

    int npts = 51, /* Number of points to subdivide the thickness of the slab
                      into */
        i; /* Loop index */
    
    vector *z, *Xdb, *Eb, *GradEb, *str, *u, *Pc;
    matrix *output;

    /* Print out a usage statement if needed */
    if(argc != 5) {
        puts("Usage:");
        puts("dry <T> <X0> <aw> <t>");
        puts("aw: Final water activity to dry to.");
        puts("t: Time to output a profile for.");
        exit(0);
    }

    /* Store command line arguments */
    cond.T = atof(argv[1]) + 237.15;
    cond.X0 = atof(argv[2]);
    RH = atof(argv[3]);
    t = atof(argv[4]);

    /* Calculate equilibrium moisture content and average diffusivity */
    //Xe = OswinIsotherm(d, RH, T);
    cond.Xe = RH;
    cond.D = DiffCh10((cond.Xe+cond.X0)/2, cond.T);
    cond.L = L;
    cond.nterms = NTERMS;

    /* Make a bunch of vectors */
    z = CreateVector(npts);
    Xdb = CreateVector(npts);
    Pc = CreateVector(npts);
    Eb = CreateVector(npts);
    GradEb = CreateVector(npts);
    str = CreateVector(npts);
    u = CreateVector(npts);

    /* Set up the domain by assigning the z value of each point to a slot in
     * the vector */
    for(i=0; i<npts; i++)
        setvalV(z, i, L/npts * i);

    for(i=0; i<npts; i++) {
        /* Moisture content */
        setvalV(Xdb, i, CrankEquationFx(valV(z, i), t, cond));
        /* Binding energy */
        setvalV(Eb, i, Esurf(valV(Xdb, i), cond.T));
        /* Binding energy gradient */
        setvalV(GradEb, i, GradEsurf(valV(z, i), t, cond));
        /* Capillary Pressure */
        setvalV(Pc, i, PoreP(valV(z, i), t, cond));
        /* Strain */
        setvalV(str, i, strainpc(valV(z, i), t, cond, &GradEsurf, &CreepZhu));
    }
    /* Find the displacement at all points */
    u = displacementV(str, L);
    PrintVector(Pc);

    /* Smash together the desired vectors and print the resulting matrix */
    output = CatColVector(7, z, Xdb, Pc, Eb, GradEb, str, u);
    mtxprntfilehdr(output, "out.csv","x [m],Xdb [kg/kg db],Pc [Pa],Eb [J/m^3],GradEb [Pa/m^3],strain,u [m]\n");

    return 0;
}

