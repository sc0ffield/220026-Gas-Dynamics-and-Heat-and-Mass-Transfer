#ifndef PRESOLVER_H
#define PRESOLVER_H

void FV(double W, double h, double Lf, double *F23,
        double *F24, double *F2sky, double *F32,
        double *F34, double *F3sky, double *F42,
        double *F43, double *F4sky);

void mesh(double *x, int n_total, int Na, int Nb, int Nc,
            int Nd, int Nf, double delta_xA, double delta_xB,
            double delta_xC, double delta_xD, double delta_xF,
            double *x_nodes);

double ConvecF(double Tw,double Text,double Pext,double x,
                double vext);

double ConvecNatu(double Tw,double Text,double Pext,double x);

double q_A(double t); double q_D(double t);
double q_VA(double t); double q_VB(double t);
double q_VC(double t); double q_VD(double t);

double Cpa(double T); double Cpb(double T);
double Cpc(double T); double Cpd(double T);
double Cpf(double T);

#endif