#include "pasta.h"
#include "choi-okos.h"

double Esurf(double Xdb, double T)
{
    double Mw = 18,
           rhos,
           Xmin = 1e-9,
           dX,
           X,
           E = 0;
    int nslice = 100,
        i;
    oswin *d;
    choi_okos *co;

    d = CreateOswinData();
    co = CreateChoiOkos(PASTACOMP);
    
    dX = (Xdb-Xmin)/nslice;
    for(i=0; i<nslice; i++) {
        X = Xmin+i*dX;
        E += BindingEnergyOswin(d, X, T) * dX;
    }

    rhos = rho(co, T);
    return Mw*Xdb*rhos*E;
}

