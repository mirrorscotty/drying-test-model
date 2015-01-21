#include <math.h>

#include "matrix.h"
#include "drying.h"
#include "material-data.h"

#include <stdlib.h>
#include <stdio.h>

/* Number of terms to use for the Crank equation when solving for moisture
 * profile. */
#define NTERMS 1000

double maxstrain(double t, double x, drydat d)
{
    double Xdb = CrankEquationFx(x, t, d),
           rhow, rhos,
           vs, vw;
    choi_okos *co;

        co = CreateChoiOkos(WATERCOMP);
        rhow = rho(co, d.T);
        DestroyChoiOkos(co);

        co = CreateChoiOkos(PASTACOMP);
        rhos = rho(co, d.T);
        DestroyChoiOkos(co);

        vs = 1/(1+rhos*d.X0/rhow);
        vw = Xdb*rhos/rhow * vs;

        return (vs+vw)-1;
}

