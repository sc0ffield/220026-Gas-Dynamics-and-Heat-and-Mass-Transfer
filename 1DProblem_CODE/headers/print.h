#ifndef PRINT_H
#define PRINT_H

void print_state(double L, double xchange, double lambda1, 
                double lambda2, double T0, double TL);

void print_mesh(int N, double dx, double* x, double *lambdas, double *x_vc);

void print_initialf(double Tstart, double T0, double TL);

void print_coeff();

void print_solver();

void print_results(double delta, int iter, int N, double *T);

#endif