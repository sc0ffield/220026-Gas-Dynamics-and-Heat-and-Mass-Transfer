#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/coeff.h"
#include "../headers/presolver.h"

void clear(double *x){
    for(int i=0;i<sizeof(x);i++)
        x[i]=0;
}

void coeff_boundary(double H, double W,
                double beta, double lambda_A, double Text, double t, double lambda_D,
                int Ntot, double Pext, double **T, double *x_nodes, int n,
                double *aE, double *aW, double *aP, double *bP, double told){
    //node inicial paret
    double dPEini = x_nodes[1] - x_nodes[0];
    double Seini = H*W;
    double alpha_natin = ConvecNatu(T[n+1][0], Text, Pext, H);
    aE[0] = ((beta*lambda_A)/dPEini)*Seini;
    aW[0] = 0;
    aP[0] = ((beta*lambda_A )/dPEini)*Seini+beta*alpha_natin*Seini;
    bP[0] = beta*q_A(t)*Seini+(1-beta)*(q_A(told)+lambda_A*((T[n][1]-T[n][0])/dPEini)
            +alpha_natin*(Text-T[n][0]))*Seini+alpha_natin*beta*Text*Seini;
    
    //node final paret
    double dPW = x_nodes[Ntot-1]-x_nodes[Ntot-2];
    double alpha_natout = ConvecNatu(T[n+1][Ntot-1],Text,Pext,H);
    aE[Ntot-1] = 0;
    aW[Ntot-1] = beta*(lambda_D/dPW)*H*W;
    aP[Ntot-1] = beta*(lambda_D/dPW)*H*W+alpha_natout*beta*H*W;
    bP[Ntot-1] = (1-beta) *((-(lambda_D/dPW)*(T[n][Ntot-1]-T[n][Ntot-2]))+q_D(told)
                    -alpha_natout*(T[n][Ntot -1]-Text))*H*W+q_D(t)*H*W*beta
                    +alpha_natout*beta*Text*H*W;
}

void coeff_int(int Ntot, int Na, int Nb, int Nc, int Nd, int Nf,double told, double t,
                double *x, double *x_nodes,double lambda_A, double lambda_B,double lambda_C,
                double lambda_D, double lambda_F, double H, double W,double **T,
                double rhoA, double rhoB, double rhoC, double rhoD, double rhoF, int n,
                double ef,int numfin, double Text, double Tsky, double Pext, double vext,
                double epsilon, double sigma,double deltat,double alpha_ext,double *Qrad,
                double *QconvNat, double beta, double *aE, double *aW, double *aP, 
                double *bP){
    //calcul de distancies, superficies, propietats,... corresponents al node
    double dPE, dPW, dPe, deE, dPw, dwW, Se, Sw, Vp, qv, qv_old, Sconv, Cpp, rho, nat;
    double lambdaE, lambdaW, Sconv2, alpha_convNat;
    for(int i=1;i=Ntot-1;i++){
        dPE = x_nodes[i+1]-x_nodes[i];
        dPW = x_nodes[i]-x_nodes[i-1];
        if(i==Na||i==Na+1||i==Na+Nb||i==Na+Nb+1||i==Na+Nb+Nf+1||i==Na+Nb+Nf||
            i==Na+Nb+Nc+Nf+1||i==Na+Nb+Nc+Nf){
            dPe=x[i]-x_nodes[i];
            deE=x_nodes[i+1]-x[i];
            dPw =x_nodes[i]-x[i-1];
            dwW =x[i-1]-x_nodes[i-1];
        }
        if(i<Na){
            lambdaE = lambda_A;
            lambdaW = lambda_A;
            Se = H*W;
            Sw = H*W;
            Vp = (x[i]-x[i-1])*W*H;
            qv = q_VA(t);
            qv_old = q_VA(told);
            Sconv = 0;
            Cpp = Cpa((T[n+1][i]+T[n][i])/2);
            rho = rhoA;
            nat = 1;
        }
        else if(i==Na){
            lambdaE = dPE/((dPe/lambda_A)+(deE/lambda_B));
            lambdaW = lambda_A;
            Se = H*W;
            Sw = H*W;
            Vp = (x[i]-x[i-1])*W*H;
            qv = q_VA(t);
            qv_old = q_VA(told);
            Sconv = 0;
            Cpp = Cpa((T[n+1][i]+T[n][i])/2);
            rho = rhoA;
            nat =1;
        }
        else if(i==Na+1){
            lambdaE = lambda_B;
            lambdaW = dPW/((dPw/lambda_B)+(dwW/lambda_A));
            Se = H*W;
            Sw = H*W;
            Vp = (x[i]-x[i-1])*W*H;
            qv = q_VB(t);
            qv_old = q_VB(told);
            Sconv = 0;
            Cpp = Cpb((T[n+1][i]+T[n][i])/2);
            rho = rhoB;
            nat = 1;
        }
        else if(Na+1<i<Na+Nb){
            lambdaE = lambda_B;
            lambdaW = lambda_B;
            Se = H*W;
            Sw = H*W;
            Vp = (x[i]-x[i-1])*W*H;
            qv = q_VB(t);
            qv_old = q_VB(told);
            Sconv = 0;
            Cpp = Cpb((T[n+1][i]+T[n][i])/2);
            rho = rhoB;
            nat = 1;
        }
        else if(i==Na+Nb){
            lambdaE = dPE/((dPe/lambda_B)+(deE/lambda_F));
            lambdaW = lambda_B;
            Se = ef*W*numfin;
            Sw = H*W;
            Sconv = H*W-numfin*ef*W;
            Vp = (x[i]-x[i-1])*W*H;
            qv = q_VB(t);
            qv_old = q_VB(told);
            Cpp = Cpb((T[n+1][i]+T[n][i])/2);
            rho = rhoB;
            nat = 1;
        }
        else if(i==Na+Nb+1){
            lambdaE = lambda_F;
            lambdaW = dPW/((dPw/lambda_F)+(dwW/lambda_B));
            Se = ef * W * numfin;
            Sw = ef * W * numfin;
            Sconv = 2* numfin *(x[i]-x[i-1])*W;
            Vp =(x[i] - x[i-1]) * W * ef * numfin;
            qv = 0;
            qv_old = 0;
            Cpp = Cpf((T[n+1][i]+T[n][i])/2);
            rho = rhoF;
            nat =0;
        }
        else if(Na + Nb +1 <i < Na + Nb + Nf){
            lambdaE = lambda_F;
            lambdaW = lambda_F;
            Se = ef * W * numfin;
            Sw = ef * W * numfin;
            Sconv = 2 * numfin * ( x[ i ] - x[ i - 1]) * W;
            Vp = ( x [ i ] - x [ i - 1]) * W * ef * numfin;
            qv = 0;
            qv_old = 0;
            Cpp = Cpf((T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoF;
            nat =0;
        }
        else if(i == Na + Nb + Nf){
            lambdaE = dPE /(( dPe / lambda_F ) +( deE / lambda_C ) );
            lambdaW = lambda_F;
            Se = ef * W * numfin;
            Sw = ef * W * numfin;
            Sconv = 2 * numfin * ( x[ i ] - x[ i - 1]) * W;
            Vp = ( x[ i ] - x[ i - 1]) * W * ef * numfin;
            qv = 0;
            qv_old = 0;
            Cpp = Cpf(( T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoF;
            nat =0;
        }
        else if (i == Na + Nb + Nf +1){
            lambdaE = lambda_C;
            lambdaW = dPW /(( dPw / lambda_C ) +( dwW / lambda_F ) );
            Se = H * W;
            Sw = ef * W * numfin;
            Sconv = H *W - numfin * W * ef;
            Vp = ( x[ i ] - x[ i - 1]) * W * H;
            qv = q_VC ( t );
            qv_old = q_VC ( told );
            Cpp = Cpc (( T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoC;
            nat =1;
        }
        else if(Na + Nb + Nf +1 <i < Na + Nb + Nf + Nc){
            lambdaE = lambda_C;
            lambdaW = lambda_C;
            Se = H * W;
            Sw = H * W;
            Vp =( x[ i ] - x[ i - 1]) * W * H;
            qv = q_VC ( t );
            qv_old = q_VC ( told );
            Sconv =0;
            Cpp = Cpc (( T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoC;
            nat =1;
        }
        else if (i == Na + Nb + Nf + Nc){
            lambdaE = dPE /(( dPe / lambda_C ) +( deE / lambda_D ) );
            lambdaW = lambda_C;
            Se = H * W;
            Sw = H * W;
            Vp =( x[ i ] - x[ i - 1]) * W * H;
            qv = q_VC ( t );
            qv_old = q_VC ( told );
            Sconv =0;
            Cpp = Cpc (( T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoC;
            nat =1;
        }
        else if (i == Na + Nb + Nf + Nc +1){
            lambdaE = lambda_D;
            lambdaW = dPW /(( dPw / lambda_D ) +( dwW / lambda_C ) );
            Se = H * W;
            Sw = H * W;
            Vp = ( x[ i ] - x[ i - 1]) * W * H;
            qv = q_VD ( t );
            qv_old = q_VD ( told );
            Sconv =0;
            Cpp = Cpd (( T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoD;
            nat =1;
        }
        else{
            lambdaE = lambda_D;
            lambdaW = lambda_D;
            Se = H * W;
            Sw = H * W;
            Vp = ( x[ i ] - x[ i - 1]) * W * H;
            qv = q_VD ( t );
            qv_old = q_VD ( told );
            Sconv =0;
            Cpp = Cpd (( T[ n +1][ i ]+ T[ n ][ i ]) /2);
            rho = rhoD;
            nat =1;
        }
        //radiacio i conveccio natural
        Sconv2 =( x[ i ] - x[i -1])*W*nat;
        alpha_convNat = ConvecF ( T[ n + 1][ i ] ,Text ,Pext ,x_nodes[i] , vext );
        Qrad[ i ] = Sconv2 * epsilon * sigma * ( pow(T[ n + 1][ i ],4)-pow(Tsky,4));
        QconvNat[ i ] = Sconv2 * alpha_convNat * ( T[ n + 1][ i ]- Text );

        //calcul de coeficients dels nodes interns
        aE[i]= beta *( Se / dPE ) * lambdaE;
        aW[ i ]= beta *( Sw / dPW ) * lambdaW;
        aP[ i ]= aE[ i ]+ aW[ i ]+( rho * Vp * Cpp ) / deltat + alpha_ext *Sconv * beta + 
                alpha_convNat * beta * Sconv2 + beta *epsilon * sigma * pow(T[ n +1][ i ],3)* Sconv2;
        bP[ i ]= beta *( qv * Vp ) +((( rho * Vp * Cpp ) / deltat ) * T[ n ][ i ])+
                (1 - beta ) *( qv_old *Vp +(( lambdaE *(T[ n ][ i +1] - T[ n ][ i])*Se)/ 
                dPE)-(( lambdaW *( T[ n ][ i ] - T[ n ][ i -1]) * Sw ) /dPW ) -
                ( alpha_ext *(T[ n ][ i ] - Text ) * Sconv ) -(alpha_convNat *(T[ n ][ i ] - 
                Text ) * Sconv2 ) -( epsilon *sigma *( pow(T[ n ][ i ],4) - pow(Tsky,4)) * Sconv2 ) ) +
                alpha_ext *Text * Sconv * beta + alpha_convNat * beta * Text * Sconv2 +
                beta * epsilon * sigma * pow(Tsky,4)* Sconv2;
    }
}