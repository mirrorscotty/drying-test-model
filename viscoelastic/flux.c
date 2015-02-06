#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "material-data.h"
#include "visco.h"
#include "drying.h"

#define NX 20
#define NTERMS 100

double SurfDisplace(double t, drydat cond)
{
    int nx = NX;
    double dx = cond.L/nx,
           d;
    int i;

    vector *x, *Xdb, *str;
    x = CreateVector(nx);
    Xdb = CreateVector(nx);
    str = CreateVector(nx);

    for(i=0; i<nx; i++)
        setvalV(x, i, dx*i);
    for(i=0; i<nx; i++)
        setvalV(str, i, strainpc(valV(x,i),
                                 t,
                                 cond,
                                 &GradEsurf,
                                 &CreepGina));

    d = displacement(nx-1, str, cond.L);
    DestroyVector(x);
    DestroyVector(Xdb);
    DestroyVector(str);

    return d;
}

double SurfDisplaceMax(double t, drydat cond)
{
    int nx = NX;
    double dx = cond.L/nx, d;
    int i;

    vector *x, *Xdb, *str;
    x = CreateVector(nx);
    Xdb = CreateVector(nx);
    str = CreateVector(nx);

    for(i=0; i<nx; i++)
        setvalV(x, i, dx*i);
    for(i=0; i<nx; i++)
        setvalV(str, i, maxstrain(t, valV(x,i), cond));

    d = displacement(nx-1, str, cond.L);
    DestroyVector(x);
    DestroyVector(Xdb);
    DestroyVector(str);

    return d;
}

double SurfDisplaceEta(double t, drydat cond)
{
    int nx = NX;
    double dx = cond.L/nx, d;
    int i;

    vector *x, *Xdb, *str;
    x = CreateVector(nx);
    Xdb = CreateVector(nx);
    str = CreateVector(nx);

    for(i=0; i<nx; i++)
        setvalV(x, i, dx*i);
    for(i=0; i<nx; i++)
        setvalV(str, i, etastrain(t, valV(x,i), cond, .3));

    d = displacement(nx-1, str, cond.L);
    DestroyVector(x);
    DestroyVector(Xdb);
    DestroyVector(str);

    return d;
}

double SurfDisplaceEq(double t, drydat cond)
{
    int nx = NX;
    double dx = cond.L/nx,
           d;
    int i;

    vector *x, *Xdb, *str;
    x = CreateVector(nx);
    Xdb = CreateVector(nx);
    str = CreateVector(nx);

    for(i=0; i<nx; i++)
        setvalV(x, i, dx*i);
    for(i=0; i<nx; i++)
        setvalV(str, i, EqStrainPc(t, valV(x,i), cond));

    //PrintVector(str);
    d = displacement(nx-1, str, cond.L);
    DestroyVector(x);
    DestroyVector(Xdb);
    DestroyVector(str);

    return d;
}

double SurfMoistureFlux(double t, drydat cond)
{
    double J, h = 1e-7, Xdbs, Xdbsmh;

    Xdbs = CrankEquationFx(cond.L, t, cond);
    Xdbsmh = CrankEquationFx(cond.L-h, t, cond);

    J = cond.D*(Xdbs-Xdbsmh) / h;

    return J;
}

double StressCummings(double t, drydat cond)
{
    int nslice = 10, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += VEStress(i*dx, t, cond, &etastrain, &RelaxCummings)/nslice;
    return avg;
}
double StressLaura(double t, drydat cond)
{
    int nslice = 10, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += VEStress(i*dx, t, cond, &etastrain, &RelaxLaura)/nslice;
    return avg;
}
double StressZhu(double t, drydat cond)
{
    int nslice = 10, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += VEStress(i*dx, t, cond, &etastrain, &RelaxZhu)/nslice;
    return avg;
}
double StressEsurf(double t, drydat cond)
{
    int nslice = 10, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += GradEsurf(i*dx, t, cond)/nslice;
    return avg;
}

int main(int argc, char *argv[])
{
    double L = 1e-3, /* Length [m] */
           t,
           nt,
           RH,
           dt;

    drydat cond;
    //oswin *d; /* Isotherm data */

    int i; /* Loop index */
    
    vector *Xdb, /* Slab moisture content */
           *Pc,
           *PL, *PC, *PZ, *PE,
           *Js, /* Moisture flux at the surface */
           *tv, /* Time vector */
           *VV0,
           *VV0Max,
           *VV0Eq,
           *MaxDisp,
           *EtaDisp,
           *EqDisp,
           *Disp;

    matrix *out;

    /* Print out a usage statement if needed */
    if(argc != 6) {
        puts("Usage:");
        puts("dry <T> <X0> <aw> <t> <nt>");
        puts("T: Drying temperature [C]");
        puts("X0: Initial moisture content of the pasta. [kg/kg db]");
        puts("aw: Final water activity to dry to. [-]");
        puts("t: Length of time to simulate. [s]");
        puts("nt: Number of time steps to use.");
        exit(0);
    }

    /* Store command line arguments */
    cond.T = atof(argv[1]) + 273.15;
    cond.X0 = atof(argv[2]);
    RH = atof(argv[3]);
    t = atof(argv[4]);
    nt = atof(argv[5]);

    dt = t/nt;

    /* Calculate equilibrium moisture content and average diffusivity */
    //d = CreateOswinData();
    cond.Xe = RH; //OswinIsotherm(d, RH, cond.T);
    cond.D = DiffCh10((cond.Xe+cond.X0)/2, cond.T);
    cond.L = L;
    cond.nterms = NTERMS;

    tv = CreateVector(nt);
    Xdb = CreateVector(nt);
    Pc = CreateVector(nt);
    PL = CreateVector(nt);
    PC = CreateVector(nt);
    PZ = CreateVector(nt);
    PE = CreateVector(nt);
    Js = CreateVector(nt);
    Disp = CreateVector(nt);
    MaxDisp = CreateVector(nt);
    EqDisp = CreateVector(nt);
    EtaDisp = CreateVector(nt);
    VV0 = CreateVector(nt);
    VV0Max = CreateVector(nt);
    VV0Eq = CreateVector(nt);

    for(i=0; i<nt; i++) {
        setvalV(tv, i, dt*i);

        //setvalV(Xdb, i, CrankEquationFx(.9, i*dt, L, D, Xe, X0, NTERMS));
        setvalV(Xdb, i, CrankEquation(i*dt, cond));
        //setvalV(Xdb, i, EqStrainPc(i*dt, .9, 0, D, X0, Xe, L, T));

        setvalV(Pc, i, pore_press(valV(Xdb, i), cond.T));
        setvalV(PL, i, StressLaura(i*dt, cond));
        setvalV(PC, i, StressCummings(i*dt, cond));
        setvalV(PZ, i, StressZhu(i*dt, cond));
        setvalV(PE, i, StressEsurf(i*dt, cond));

        setvalV(Js, i, SurfMoistureFlux(i*dt, cond));

        setvalV(Disp, i, SurfDisplace(dt*i, cond));
        setvalV(MaxDisp, i, SurfDisplaceMax(dt*i, cond));
        setvalV(EqDisp, i, SurfDisplaceEq(dt*i, cond));
        setvalV(EtaDisp, i, SurfDisplaceEta(dt*i, cond));

        setvalV(VV0, i, (1e-3+valV(Disp, i))/1e-3);
        setvalV(VV0Max, i, (1e-3+valV(MaxDisp, i))/1e-3);
        setvalV(VV0Eq, i, (1e-3+valV(EqDisp, i))/1e-3);
    }

    out = CatColVector(15, tv, Xdb, Disp, EtaDisp, MaxDisp, EqDisp, Pc, PL, PC, PZ, PE, Js, VV0, VV0Max, VV0Eq);

    mtxprntfilehdr(out, "out.csv", "Time [s],Moisture Content [kg/kg db],Surface Displacement [m],EtaDisp [m],Max Surf Disp [m],Eq Disp [m],Pressure [Pa],PL,PC,PZ,PE,Surface Water Flux [kg/m^2],V/V0,V/V0 Max,V/V0 Eq\n");

    return 0;
}

