#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "headers/print.h"
#include "headers/mesh.h"
#include "headers/solver.h"

int main(){
    //1.ENTRADA DE DADES
        //físiques
    const double L = 5.0, xchange = 2.8; 
    const double lambda1 = 120, lambda2 = 180;
    const double T0 = 100, TL = 200;  
        //numériques
    const int N = 100; const double delta = 1e-9, Tstart = (T0+TL)/2;

    //2.DEFINICIÓ DE VECTORS I CÀLCULS PREVIS
    double *x = (double*)malloc((N+2)*sizeof(double)); //posició dels nodes
    double *x_vc = (double*)malloc((N+1)*sizeof(double)); //posició de les cares
    double *lambdas = (double*)malloc((N+2)*sizeof(double)); //lambdas dels nodes
    const double dx = L/N; //distància entre cada node
    generate_mesh(x, dx, xchange, N, L, lambdas, lambda1, lambda2, x_vc);

    //3.MAPA INICIAL DE TEMPERATURA
    double *T_s = (double*)malloc((N+2)*sizeof(double)); //temperatura suposada
    double *T_c = (double*)malloc((N+2)*sizeof(double)); //temperatura calculada

    for(int i=0;i<N+2;i++)
        T_s[i] = Tstart;
    T_s[0]=T0; T_s[N+1]=TL;

    //4.CÀLCUL DE COEFICIENTS DE DISCRETITZACIÓ
    double *aE = (double*)malloc((N+2)*sizeof(double));
    double *aP = (double*)malloc((N+2)*sizeof(double));
    double *aW = (double*)malloc((N+2)*sizeof(double));
    //bP = 0 (no hi ha fonts internes ni convecció)

    coeff(N, aW, aP, aE, lambdas);

    //5.SOLVER (GAUSS-SEIDEL)
    int iter=0; double max_error;
    while(1){
        iter++;
        GaussSeidel( N, aW, aP, aE, T_c, T_s, &max_error);
        if(max_error<delta) //criteri de convergència
            break;
        for(int i=0;i<N+2;i++)
            T_s[i]=T_c[i]; //actualitzar
    }
    
    free(T_s); free(aE); free(aW); free(aP);
    
    //7.CÀLCULS FINALS I IMPRESSIÓ DE RESULTATS
    print_state(L,xchange,lambda1,lambda2,T0,TL);
    print_mesh(N, dx, x, lambdas, x_vc);
    print_initialf(Tstart, T0, TL);
    print_coeff();
    print_solver();
    print_results(delta, iter, N, T_c);
    
    //PLOTING
    //FILE *file = fopen("results.txt", "w");
    //for (int i = 0; i < N+2; i++) //escriure .txt
    //    fprintf(file, "%f %f\n", x[i], T_c[i]);
    //fclose(file);

    //FILE *pipeForGNUPlot=NULL;
    //pipeForGNUPlot = popen("gnuplot -p","w"); //plot amb gnuplot
    //fprintf(pipeForGNUPlot,"set title \"T(x) results\"\n");
    //fprintf(pipeForGNUPlot,"set xlabel \"x [m]\"\n");
    //fprintf(pipeForGNUPlot,"set ylabel \"T(x) [K]\"\n");
    //fprintf(pipeForGNUPlot,"plot 'results.txt' w l linewidth 3 \n");
    //pclose(pipeForGNUPlot);

    free(T_c); free(x); free(x_vc); free(lambdas);
    
    //8.FI
    return 0;
}
