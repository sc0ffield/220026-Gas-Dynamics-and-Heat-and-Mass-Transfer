#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "headers/presolver.h"
#include "headers/coeff.h"
#include "headers/solver.h"

int main(){
    //1.ENTRADA DE DADES
        //fisiques
    const int Na = 30, Nb = 10, Nc = 10, Nd = 30, Nf = 50;
    const double ea = 0.01, eb = 0.005, ec = 0.005, ed = 0.01, ef = 0.0;
    const double Lf = 0.4, H = 0.5, W = 0.5, h=0.05, rhoF = 9100.0;
    const double vext = 20.0, alpha_ext = 70.0, Text = 293.0, Pext = 1e5;
    const double Tsky = 268.0, epsilon = 0.9, sigma = 5.67e-8, lambda_F = 220.0; 
    const double rhoA = 6200.0, rhoB = 7800.0, rhoC = 8800, rhoD = 6700;
    const double lambda_A = 50.0, lambda_B = 220.0, lambda_C = 220.0, lambda_D = 50.0;
        //numeriques
    const int nmax = 1200;
    const double Tstart = 273.0;
    const double delta = 5e-5, deltat = 0.5, beta = 0.5;
    double maxdiff = 1e20, critFinal = 1e9, Mdiff=0, cFi=0;

    //2.DEFINICIO DE VECTORS I CALCULS PREVIS
    int  crit = 0, n = 0, numfin = 100;
    double t = 0, told,diff,cFinal;
    int n_total = Na + Nb + Nc + Nd + Nf, Ntot = n_total + 2;
    double *Qrad = (double*)malloc((Ntot)*sizeof(double));
    double *QconvNat = (double*)malloc((Ntot)*sizeof(double));
    double *aE = (double*)malloc((Ntot)*sizeof(double));
    double *aW = (double*)malloc((Ntot)*sizeof(double));
    double *aP = (double*)malloc((Ntot)*sizeof(double));
    double *bP = (double*)malloc((Ntot)*sizeof(double));
    double *Told = (double*)malloc((Ntot)*sizeof(double));       
    double *Tval = (double*)malloc((Ntot)*sizeof(double));
        //malla
    double delta_xA = ea/Na, delta_xB = eb/Nb, delta_xC = ec/Nc;
    double delta_xD = ed/Nd, delta_xF = ef/Nf;
    double *x = (double*)malloc((n_total+1)*sizeof(double));
    double *x_nodes = (double*)malloc((Ntot)*sizeof(double));
    double **T = (double**)malloc((nmax)*sizeof(double));
    for(int i=0;i<nmax;i++)
        T[i]=malloc(Ntot*sizeof(double));
    double *Tcomp = (double*)malloc((Ntot)*sizeof(double));
    mesh(x,n_total,Na,Nb,Nc,Nd,Nf,delta_xA,delta_xB,delta_xC,delta_xD,delta_xF,x_nodes);
        //factors de radiacio 
    double F23, F24, F2sky, F32, F34, F3sky, F42, F43, F4sky;
    FV(W,h,Lf,&F23,&F24,&F2sky,&F32,&F34,&F3sky,&F42,&F43,&F4sky);
    
    //3.
        //3.1.MAPA INICIAL DE TEMPERATURES
    for(int i=0;i<Ntot;i++){
        Tcomp[i]=Tstart;
        T[0][i]=Tcomp[i];
    }
        //3.2.COEFICIENTS
    while(crit==0)
    {
        told=t;
        t=t+deltat;
        while(maxdiff>delta){
            clear(aE); clear(aW); clear(aP); clear(bP);
            //coefficients dels nodes externs
            coeff_boundary(H, W, beta, lambda_A, Text, 
                            t, lambda_D, Ntot, Pext, T, x_nodes, n,aE,aW,aP,bP,told);
            //coefficients dels nodes interns
            coeff_int( Ntot, Na, Nb, Nc, Nd, Nf, told, t,x, x_nodes, lambda_A, lambda_B,
                         lambda_C, lambda_D, lambda_F, H, W,T, rhoA, rhoB, rhoC, rhoD, rhoF,
                        n,ef,numfin,Text,Tsky, Pext,vext,epsilon,sigma,deltat,alpha_ext,
                        Qrad,QconvNat,beta,aE,aW,aP,bP);
            //3.3.SOLVER
            solve(Told, Tval, Ntot, aE, aW,aP, bP, T, n);
            for (int i=0; i<Ntot;i++) {
                diff = fabs(Tval[i] - Told[i]);
                if (diff>Mdiff)
                    Mdiff=diff;
            }
            if(Mdiff<maxdiff)
                maxdiff=Mdiff;
            update(T, Tval, n, Ntot);
        }
        n++;
        for(int i=1;i<nmax;i++){
            cFinal = fabs(T[i]-T[i-1]);
            if(cFinal>cFi)
                cFi = cFinal;
        }
        if(cFi<critFinal)
            critFinal=cFi;
        //4.CRITERI DE CONVERGENCIA
        if(n>=nmax||critFinal<delta)  
            crit = 1;
    }

    //5.IMPRESSIO DE RESULTATS I CALCULS FINALS
    //el codi d'impressio no esta inclos aqui ja que no s'ha fet en C, sino en Python

    free(x); free(x_nodes); free(T); free(Tcomp); free(aE); free(aW); free(bP);
    free(aP); free(Tval); free(Told); free(Qrad); free(QconvNat);

    //6.FI
    return 0;
}