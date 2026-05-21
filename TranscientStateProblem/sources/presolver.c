#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/presolver.h"

void FV(double W, double h, double Lf, double *F23,
        double *F24, double *F2sky, double *F32,
        double *F34, double *F3sky, double *F42,
        double *F43, double *F4sky){
    double A2 = W*h, A3 = W*Lf, A4 = W*h;
    double theta = atan(h/Lf);
    
    *F23 = ((Lf+h)/(Lf/cos(theta)))/(2*h);
    *F24 = ((Lf/cos(theta)*2)-(2*Lf))/(2*h);
    *F2sky = 1.0-*F23-*F24;
    *F32 = (A2/A3)*(*F23);
    *F34 = (((Lf+h)-(Lf/cos(theta)))/(2*Lf));
    *F3sky = 1.0-*F32-*F34;
    *F42 = (A2/A4)*(*F24);
    *F43 = (A3/A4)*(*F34);
    *F4sky = 1.0-*F42-*F43;
}


void mesh(double *x, int n_total, int Na, int Nb, int Nc,
            int Nd, int Nf, double delta_xA, double delta_xB,
            double delta_xC, double delta_xD, double delta_xF,
            double *x_nodes){
    double dx;
    x[0] = 0;
    x_nodes[0] = 0;
    for(int i=1;i<n_total+1;i++){ 
        if(i<=Na)
            dx = delta_xA;
        else if(i<=Na+Nb)
            dx = delta_xB;
        else if(i<=Na+Nb+Nf)
            dx = delta_xF;
        else if(i<=Na+Nb+Nf+Nc)
            dx = delta_xC;
        else
            dx = delta_xD;
        x[i]=x[i-1]+dx;
        x_nodes[i] = (x[i]+x[i-1])/2;
    }
    x_nodes[n_total+1]=x[n_total];
}

double ConvecF(double Tw,double Text,double Pext,double x,
                double vext){
    double mu, Pr, alpha_F, Re, Nu;
    double Tm = (Tw+Text)/2;
    double cp=1034.09-2.849e-1*Tm+7.817e-4*pow(Tm,2)-4.971e-7
                *pow(Tm,3)+1.077e-10*pow(Tm,4);
    double rho = Pext/(287*Tm);
    double lamb = (2.648e-3*pow(Tm,0.5))/(1+((245.4/Tm))*(pow(10,-12/Tm)));
    
    if(Tm<1100){
        double C1 = 1.458e-6, S = 110.4;
        mu = (C1*pow(Tm,1.5))/(Tm+S);
        Pr = mu*cp/lamb;
    }

    else{
        mu =(2.5393e-5*(sqrt(Tm/273.15)))/(1+(122/Tm));
        Pr =0.71;
    }

    Re =(rho*vext*x)/mu;
    if(Re<5e5)
        Nu =0.332*pow(Re,1/2)*pow(Pr,1/3);
    else
        Nu =0.029*pow(Re,4/5)*pow(Pr,1/3);

    alpha_F =(Nu*lamb)/x;
    return (alpha_F);
}

double ConvecNatu(double Tw,double Text,double Pext,double x){
    double g=9.81, mu,Pr,Gr,Ra,C,n,k, alpha_nat,Nu;
    double Tm =(Tw+Text)/2;
    double cp=1034.09-2.849e-1*Tm+7.817e-4*pow(Tm,2)-4.971e-7
                *pow(Tm,3)+1.077e-10*pow(Tm,4);
    double rho = Pext/(287*Tm);
    double lamb = (2.648e-3*pow(Tm,0.5))/(1+((245.4/Tm))*(pow(10,-12/Tm)));

    if(Tm<1100){
        double C1 =1.458e-6, S = 110.4;
        mu =(C1*pow(Tm,1.5))/(Tm+S);
        Pr =(mu)*cp/lamb;
    }
    else{
        mu =(2.5393e-5*(sqrt(Tm/273.15)))/(1+(122/Tm));
        Pr =0.71;
    }
    Gr = (g*(1/Tm)*pow(rho,2)*abs(Tw-Text)*pow(x,3))/pow(mu,2);
    Ra = Pr*Gr;
    if(Ra<pow(10,9)){
        C =0.8;
        n =1/4;
        k =pow((1+pow((1+1/(sqrt(Pr))),2)),-1/4);
    }
    else{
        C =0.0246;
        n =2/5;
        k = pow((pow(Pr,1/6)/(1+0.494*pow(Pr,2/3))),2/5);
    }
    Nu =(C*pow(Ra,n))*k;
    alpha_nat =(Nu*lamb)/x;
    return(alpha_nat);
}

double q_A(double t){
    return 1000*sin(M_PI*t/60);
    //return 0;
}

double q_VA(double t){
    return 1e5;
    //return 0;
}

double q_VB(double t){
    return 5e4;
    //return 0;
}

double q_VC(double t){
    return 2e4;
    //return 0;
}

double q_VD(double t){
    return 1e4;
    //return 0;
}

double q_D(double ){
    return 5e4;
    //return 0;
}

//Cp -> 0.117*T+0.0004*T^2
double Cpa(double T){
    return 500;
}

double Cpb(double T){
    return 700;
}

double Cpc(double T){
    return 250;
}

double Cpd(double T){
    return 200;
}

double Cpf(double T){
    return 600;
}
//cp = 417.1+1.55*T-0.001183*T^2+3.86e-7*T^3