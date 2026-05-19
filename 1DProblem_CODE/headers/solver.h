#ifndef SOLVER_H
#define SOLVER_H

double lambda_arm(double lambdai, double lambdaj);

void coeff(int N, double *aW, double *aP, double *aE, 
            double *lambdas);

void GaussSeidel(int N, double *aW, double *aP, double *aE,
                    double *T_c, double *T_s, double *max_error);

#endif