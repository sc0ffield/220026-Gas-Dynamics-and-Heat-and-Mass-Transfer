#ifndef COEFF_H
#define COEFF_H

void clear(double *x);

void coeff_boundary(double H, double W, 
                double beta, double lambda_A, double Text, double t, double lambda_D,
                int Ntot, double Pext, double **T, double *x_nodes, int n,
                double *aE, double *aW, double *aP, double *bP, double told);

void coeff_int(int Ntot, int Na, int Nb, int Nc, int Nd, int Nf,double told, double t,
                double *x, double *x_nodes,double lambda_A, double lambda_B,double lambda_C,
                double lambda_D, double lambda_F, double H, double W,double **T,
                double rhoA, double rhoB, double rhoC, double rhoD, double rhoF, int n,
                double ef,int numfin, double Text, double Tsky,double Pext,double vext,
                double epsilon, double sigma, double deltat,double alpha_ext,double *Qrad,
                double *QconvNat,double beta, double *aE, double *aW, double *aP, 
                double *bP);
                
#endif 