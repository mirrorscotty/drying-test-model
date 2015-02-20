#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "material-data.h"
#include "visco.h"
#include "drying.h"

#define NX 20
#define NTERMS 100

double SurfDisplace(double t, drydat cond,
                    double (*F)(double, double, drydat),
                    double (*J)(double, double, double, double, int))
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
                                 F,
                                 J));

    d = displacement(nx-1, str, cond.L);
    DestroyVector(x);
    DestroyVector(Xdb);
    DestroyVector(str);

    return -1*d;
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
        setvalV(str, i, etastrain(t, valV(x,i), cond, ETA));

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

double AvgVEStress(double t, drydat cond,
                   double (*E)(double, double, drydat, double),
                   double (*G)(double, double, double, int))
{
    int nslice = 11, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += VEStress(i*dx, t, cond, E, G)/nslice;
    return avg;
}

double StressEsurf(double t, drydat cond)
{
    int nslice = 11, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += GradEsurf(i*dx, t, cond)/nslice;
    return avg;
}

double StressPoreP(double t, drydat cond)
{
    int nslice = 11, i;
    double avg = 0,
           dx = cond.L/(nslice-1);

    for(i=0; i<nslice; i++)
        avg += -1*EffPoreP(i*dx, t, cond)/nslice;
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
    oswin *d; /* Isotherm data */

    int i; /* Loop index */
    
    vector *Xdb, /* Slab moisture content */
           *Pc,
           *PL, *PC, *PZ, *PE, *PG,
           *Js, /* Moisture flux at the surface */
           *tv, /* Time vector */
           *VV0,
           *VV0Max,
           *VV0Eq,
           *MaxDisp,
           *EtaDisp,
           *EqDisp,
           *DPPG, *DEPPG, *DGEG,
           *DPPL, *DEPPL, *DGEL;

    matrix *out;
    gordontaylor *gt;

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

    /* Print out the moisture content where the glass transition occurs. */
    gt = GTSemolina();
    printf("Tg = %g at Xdb = %g\n", cond.T, GordonTaylorInv(gt, cond.T));
    DestroyGT(gt);

    dt = t/nt;

    /* Calculate equilibrium moisture content and average diffusivity */
    d = CreateOswinData();
    cond.Xe = OswinIsotherm(d, RH, cond.T);
    cond.D = DiffCh10((cond.Xe+cond.X0)/2, cond.T);
    cond.L = L;
    cond.nterms = NTERMS;

    tv = CreateVector(nt);
    Xdb = CreateVector(nt);
    Pc = CreateVector(nt);
    PL = CreateVector(nt);
    PC = CreateVector(nt);
    PZ = CreateVector(nt);
    PG = CreateVector(nt);
    PE = CreateVector(nt);
    Js = CreateVector(nt);
    DPPG = CreateVector(nt);
    DEPPG = CreateVector(nt);
    DGEG = CreateVector(nt);
    DPPL = CreateVector(nt);
    DEPPL = CreateVector(nt);
    DGEL = CreateVector(nt);
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

        setvalV(Pc, i, StressPoreP(i*dt, cond));
        setvalV(PL, i, AvgVEStress(i*dt, cond, &etastrain, &RelaxLaura));
        setvalV(PC, i, AvgVEStress(i*dt, cond, &etastrain, &RelaxCummings)); 
        setvalV(PG, i, AvgVEStress(i*dt, cond, &etastrain, &RelaxGina));
        setvalV(PZ, i, AvgVEStress(i*dt, cond, &etastrain, &RelaxZhu));
        setvalV(PE, i, StressEsurf(i*dt, cond)); /*

        setvalV(Js, i, SurfMoistureFlux(i*dt, cond));

        setvalV(DPPG, i, SurfDisplace(dt*i, cond, &PoreP, &CreepGina));
        setvalV(DEPPG, i, SurfDisplace(dt*i, cond, &EffPoreP, &CreepGina));
        setvalV(DGEG, i, SurfDisplace(dt*i, cond, &GradEsurf, &CreepGina));
        setvalV(DPPL, i, SurfDisplace(dt*i, cond, &PoreP, &CreepLaura2));
        setvalV(DEPPL, i, SurfDisplace(dt*i, cond, &EffPoreP, &CreepLaura2));
        setvalV(DGEL, i, SurfDisplace(dt*i, cond, &GradEsurf, &CreepLaura2));
        setvalV(MaxDisp, i, SurfDisplaceMax(dt*i, cond));
        setvalV(EqDisp, i, SurfDisplaceEq(dt*i, cond)); */
        setvalV(EtaDisp, i, SurfDisplaceEta(dt*i, cond));

        //setvalV(VV0, i, (1e-3+valV(Disp, i))/1e-3);
        //setvalV(VV0Max, i, (1e-3+valV(MaxDisp, i))/1e-3);
        //setvalV(VV0Eq, i, (1e-3+valV(EqDisp, i))/1e-3);
    }

    out = CatColVector(20, tv, Xdb, DPPG, DEPPG, DGEG, DPPL, DEPPL, DGEL, EtaDisp, MaxDisp, EqDisp, Pc, PL, PC, PZ, PE, Js, VV0, VV0Max, VV0Eq);

    mtxprntfilehdr(out, "out.csv", "Time [s],Moisture Content [kg/kg db],Surface Disp Pc G [m],Disp EPc G,Disp Eb G,Disp Pc L, Disp EPc L,Disp Eb L,EtaDisp [m],Max Surf Disp [m],Eq Disp [m],Pressure [Pa],PL,PC,PZ,PE,Surface Water Flux [kg/m^2],V/V0,V/V0 Max,V/V0 Eq\n");

    return 0;
}

