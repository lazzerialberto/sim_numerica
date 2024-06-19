#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "random.h"

using namespace std;

int main(){

    vector<int> N={1,2,10,100};
    int M=10000;
    vector<double> appo2;
    Random rnd;
    rnd.Initialize();

    //Loading Random Numbers
    for(int i=0;i<10000000;i++){
        appo2.push_back(rnd.Rannyu());
    }


    //Loading data on files for each N for each distribution
    ofstream fileout("risultati_2_1.csv");
    fileout<<"1-unif,1-expo,1-lor,2-unif,2-expo,2-lor,10-unif,10-expo,10-lor,100-unif,100-expo,100-lor"<<endl;

    for(int j=0; j<M; j++){

        double sum1_unif=0;
        double sum1_expo=0;
        double sum1_lor=0;
        double sum2_unif=0;
        double sum2_expo=0;
        double sum2_lor=0;
        double sum3_unif=0;
        double sum3_expo=0;
        double sum3_lor=0;
        double sum4_unif=0;
        double sum4_expo=0;
        double sum4_lor=0;

        for(int i=0; i<N.size(); i++){
            for(int k=0; k<N[i]; k++){
                if(N[i]==1){
                    sum1_unif+=appo2[j+k*N[i]]/N[i];
                    sum1_expo+=-log(1-appo2[j+k*N[i]])/N[i];
                    sum1_lor+=tan(M_PI * (appo2[j+k*N[i]] - 0.5))/N[i];
                }
                else if(N[i]==2){
                    sum2_unif+=appo2[j+k*N[i]]/N[i];
                    sum2_expo+=-log(1-appo2[j+k*N[i]])/N[i];
                    sum2_lor+=tan(M_PI * (appo2[j+k*N[i]] - 0.5))/N[i];
                }
                else if(N[i]==10){
                    sum3_unif+=appo2[j+k*N[i]]/N[i];
                    sum3_expo+=-log(1-appo2[j+k*N[i]])/N[i];
                    sum3_lor+=tan(M_PI * (appo2[j+k*N[i]] - 0.5))/N[i];
                }
                else if(N[i]==100){
                    sum4_unif+=appo2[j+k*N[i]]/N[i];
                    sum4_expo+=-log(1-appo2[j+k*N[i]])/N[i];
                    sum4_lor+=tan(M_PI * (appo2[j+k*N[i]] - 0.5))/N[i];
                }
            }
        }
        fileout<<sum1_unif<<","<<sum1_expo<<","<<sum1_lor<<","<<sum2_unif<<","<<sum2_expo<<","<<sum2_lor<<","<<sum3_unif<<","<<sum3_expo<<","<<sum3_lor<<","<<sum4_unif<<","<<sum4_expo<<","<<sum4_lor<<endl;
    }

    fileout.close();

    rnd.SaveSeed();

    return 0;
}