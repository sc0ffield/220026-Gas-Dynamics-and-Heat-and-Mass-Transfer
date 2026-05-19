#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/solver.h"

double lambda_arm(double lambdai, double lambdaj){ 
    return ((2*lambdai*lambdaj)/(lambdai+lambdaj));
}

void coeff(int N, double *aW, double *aP, double *aE, double *lambdas){
    for(int i=0;i<N+2;i++){
        aW[i]=(i==0)?0:lambda_arm(lambdas[i],lambdas[i-1]);
        aE[i]=(i==N+1)?0:lambda_arm(lambdas[i],lambdas[i+1]);
        aP[i]=aW[i]+aE[i];
    }
}

void GaussSeidel(int N, double *aW, double *aP, double *aE,
                    double *T_c, double *T_s, double *max_error){
    double error, Merror=0;
    T_c[0]=T_s[0]; //mantenir condicions de contorn
    T_c[N+1]=T_s[N+1];
    for(int i=1;i<N+1;i++){ //nodes interns
        T_c[i]=(aE[i]*T_s[i+1]+aW[i]*T_c[i-1])/aP[i];
        error=fabs(T_s[i]-T_c[i]);
        if(error>Merror)
            Merror = error;
    }
    *max_error=Merror;
}
