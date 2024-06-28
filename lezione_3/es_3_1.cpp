#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

#include "random.h"
#include "motions.h"
#include "FunzioneBase.h"
#include "Integral.h"

using namespace std;

// erf function definition
class erf_func: public FunzioneBase{

  public:
    
   double Eval(double x) const override{

    Random rndi;
    rndi.Initialize();

    gaussiana g;    
    MonteCarloMedia I(0,x,rndi);

    return 2./sqrt(M_PI)*I.Integra(10000,g);

   }

};


int main(){

    int N=100; // number of blocks
    int M=50000; // elements per block

    // setting parameters
    double r=0.1; //risk-free interest rate
    double sigma=0.25; // volatitlity
    double K=100; //strike price
    double t=1; // delivery time
    double S_0=100; // initial asset price

    // class for sorting geometric brownian motion
    GBM gbm;

    double C=0,P=0,C_fin=0,P_fin=0,C2_fin=0,P2_fin=0; // call and put values direct
    double C_d=0,P_d=0,C_d_fin=0,P_d_fin=0,C2_d_fin=0,P2_d_fin=0; // call and put values discrete

    //output files
    ofstream fileout("results_3_1_direct.csv");
    fileout << "call,err_call,put,err_put" << endl;

    ofstream fileout2("results_3_1_discrete.csv");
    fileout2 << "call,err_call,put,err_put" << endl;

    //help variable
    int count =1;

    for(int i=0; i<N*M; i++){

        //direct step
        double S_T= gbm.step(S_0,t,0,r,sigma);

        //discrete step
        double S_i=S_0;

        for(int i=0;i<100;i++){
            S_i=gbm.step(S_i,(i+2)/100,(i+1)/100,r,sigma);
        }

        //direct calculation of C(T) and P(T)

        double max=0;
        double min=0;

        if(S_T-K>0){
            max=S_T-K;
        }

        if(S_T-K<0){
            min=K-S_T;
        }

        C+=exp(-r)*max/M;
        P+=exp(-r)*min/M;

        //discrete calculation of C(T) and P(T)

        double max_d=0;
        double min_d=0;

        if(S_i-K>0){
            max_d=S_i-K;
        }

        if(S_i-K<0){
            min_d=K-S_i;
        }

        C_d+=exp(-r)*max_d/M;
        P_d+=exp(-r)*min_d/M;

        if(i/M==count){

            //direct
            C_fin+=C;
            P_fin+=P;
            C2_fin+=pow(C,2);
            P2_fin+=pow(P,2);
            fileout << C_fin/count << "," << sqrt((C2_fin/count-pow(C_fin/count,2))/(count)) << "," << P_fin/count << "," << sqrt((P2_fin/count-pow(P_fin/count,2))/(count)) << endl;
            C=0;
            P=0;

            //discrete
            C_d_fin+=C_d;
            P_d_fin+=P_d;
            C2_d_fin+=pow(C_d,2);
            P2_d_fin+=pow(P_d,2);
            fileout2 << C_d_fin/count << "," << sqrt((C2_d_fin/count-pow(C_d_fin/count,2))/(count)) << "," << P_d_fin/count << "," << sqrt((P2_d_fin/count-pow(P_d_fin/count,2))/(count)) << endl;
            C_d=0;
            P_d=0;

            count++;
        }


    }

    fileout.close();
    fileout2.close();

    return 0;
}