#ifndef SOLVER_H
#define SOLVER_H

void solve(double *Told, double *Tval, int Ntot, double *aE, double *aW,
            double *aP, double *bP, double **T,int n);

void update(double **T, double *Tval, int n, int Ntot);

#endif