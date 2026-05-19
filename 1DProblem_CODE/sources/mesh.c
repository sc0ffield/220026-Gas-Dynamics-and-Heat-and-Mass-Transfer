#include <stdlib.h>
#include <stdio.h>
#include "../headers/mesh.h"

void generate_mesh(double *x, double dx, double xchange, 
                    int N, double L, double *lambdas,
                    double lambda1, double lambda2,
                    double *x_vc){
    x[0]=0; x[N+1]=L;
    x_vc[0]=0; lambdas[0]=lambda1; lambdas[N+1]=lambda2;
    for(int i=1;i<N+1;i++){  
        x_vc[i]=dx*i;
        x[i]=x_vc[i]-dx/2;
        if(x_vc[i]<xchange) 
            lambdas[i]=lambda1;
        else
            lambdas[i]=lambda2;
    }
}