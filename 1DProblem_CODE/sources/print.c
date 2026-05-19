#include <stdlib.h>
#include <stdio.h>
#include "../headers/print.h"

void print_state(double L, double xchange, double lambda1, 
                double lambda2, double T0, double TL){
    printf("1D HEAT CONDUCTION IN A COMPOSITE WALL\n\n"
            "PROBLEM STATEMENT\n\n");
    printf("A wall of total length %fm is composed of "
            "two materials joined at x_c = %fm. The first segment "
            "(0<x<x_c) is made of a material A with thermal conductivity "
            "lambda_1 = %fW/mK, and the second segment (x_c<x<L) is made "
            "of material B with thermal conductivity lambda_1 = %fW/mK.\n\n"
            "The system is in steady-state, and heat conduction is "
            "one-dimensional along the x-axis. The temperature at "
            "the left end of the wall is maintained at T_0 = %fK, and the "
            "temperature at the right end is maintained at T_L = %fK.\n\n",
            L,xchange,lambda1,lambda2,T0,TL);
    printf("Assuming:\n-Steady-state heat conduction\n-No internal heat " 
            "generation\n-Perfect thermal contact between the two materials "
            "\n-Constant thermal conductivity in each material\n-One-dimensional "
            "conduction only.\n\nDetermine the temperature distribution T(x) along "
            "the lenght of the wall.\n");
}

void print_mesh(int N, double dx, double* x, double *lambdas, double *x_vc){
    printf("\n\n1. DISCRETIZATION OF THE MESH\n\n");
    printf("Contemplating %d control volumes:\ndx = %f \nNode positions [m] = {",N,dx);
    for(int i=0;i<N+2;i++){
        printf("%f",x[i]);
        if(i<(N+1))
            printf(", ");
    }
    printf("}\nFace positions [m] = {");
    for(int i=0;i<N+1;i++){
        printf("%f",x_vc[i]);
        if(i<N)
            printf(", ");
    }
    printf("}\nThermal conductivities per node [W/mK] = {");
    for(int i=0;i<N+1;i++){
        printf("%f",lambdas[i]);
        if(i<N)
            printf(", ");
    }
    printf("}\n");
}

void print_initialf(double Tstart, double T0, double TL){
    printf("\n\n2. INITIAL MAP DEFINITION (internal nodes)\n\n");
    printf("T_supposed[i]=T_initial=%fK\n",Tstart);
    printf("\n\n3. BOUNDARY CONDITIONS (external nodes)\n\n");
    printf("T[0]=T_0=%fK\n",T0);
    printf("T[N+1]=T_L=%fK\n",TL);
}

void print_coeff(){
    printf("\n\n4. COEFFICIENT CALCULATION\n\n");
    printf("For this case:\n  aW = lambda_w\n  aE = lambda_e\n  aP=aW+aE "
            "\n  bP = 0\nbeing lambda_w and lambda_e the harmonic lambda "
            "between the nodes W-P and P-E.\n  lambda_harm = (2*lambdai*lambdaj)"
            "/(lambdai+lambdaj) for equal distances dx and surfaces S.\n");
}

void print_solver(){
    printf("\n\n5. GAUSS-SEIDEL SOLVER\n\n");
    printf("T_calculated[i]=(aE[i]*T[i+1]+aW[i]*T[i-1]+bP[i])/aP[i]\n"
            "\n\n6. Is max|T_calculated[i]-T_supposed[i]|<delta?\n\n  NO: "
            " T_s[i]=T_c[i] and back to point 5.\n  YES: point 7.\n");
}

void print_results(double delta, int iter, int N, double *T){
    printf("\n\n7. RESULTS\n\n");
    printf("Contemplating delta = %e, and with %d iterations:\n",delta, iter);
    printf("T[i] [K] = {");
    for(int i=0;i<N+2;i++){
        printf("%f",T[i]);
        if(i<(N+1))
            printf(", ");
    }
    printf("}\n\n\n8. END\n\n");
}

