#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/solver.h"

void solve(double *Told, double *Tval, int Ntot, double *aE, double *aW,
            double *aP, double *bP, double **T, int n){
    for(int i=0;i<Ntot;i++)
        Told[i]= T[n+1][i];
    Tval[0] = ( aE[0] * Told[1] + bP[0]) / aP[0];
    for(int i=1;i<Ntot-1;i++)
        Tval[ i ] = (aE[i]*Told[i+1]+aW[i]*Tval[i-1]+bP[i])/aP[i];
    Tval[Ntot-1] = (aW[Ntot-1] * Tval[Ntot-2]+bP[Ntot-1])/aP[Ntot-1];
}

void update(double **T, double *Tval, int n, int Ntot){
    for(int i=0; i<Ntot;i++)
        T[n+1][i] = Tval[i];
}