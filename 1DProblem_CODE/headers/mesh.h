#ifndef MESH_H
#define MESH_H

void generate_mesh(double *x, double dx, double xchange, 
                    int N, double L, double *lambdas,
                    double lambda1, double lambda2,
                    double *x_vc);

#endif