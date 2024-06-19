#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

int main(){

    int M=20000;
    double L = 1.3; //needle lenght
    double d = 2.0; //distance between lines

    int N_hit=0;
    int N_tot=0;
    double pi,pi2,err;

    Random rnd;
    rnd.Initialize();

    ofstream fileout("results_3_1.csv");

    fileout << "pi,err" << endl;

    for(int j=0; j<100; j++){

        int i=0;

        while(N_tot<M){

            double appo,appo1,appo2,theta,pos;
            appo=rnd.Rannyu()*d; //random number between 0 and d
            appo1=rnd.Rannyu(-1,1); //random number between -1 and 1
            appo2=rnd.Rannyu(); //Random number between 0 and 1

            //sorting angles from 0 to 180° in order to the simmetry of the problem
            if(sqrt(pow(appo1,2)+pow(appo2,2))<=1){
                theta=appo1/sqrt(pow(appo1,2)+pow(appo2,2));

                //Genero un angolo e la posizione di un estremo tra 0 e d, poi per capire se ha sorpassato utilizzo la proiezione sull'asse x.

                pos=appo+L*theta;

                if(appo==0 or pos<=0 or pos>=d){
                    N_hit++;
                }
                N_tot++;
            }
            i++;

        }

        if(j==0 or j==1){
            cout << N_hit << "   " << N_tot << "  " << (2*L*N_tot/(N_hit*d)) << endl;
        }

        pi+=(2*L*N_tot/(N_hit*d));
        pi2+=pow((2*L*N_tot/(N_hit*d)),2);
        err=sqrt((pi2/(j+1)-pow(pi/(j+1),2))/(j+1));

        fileout << pi/(j+1) << "," << err << endl;
        N_hit=0;
        N_tot=0;


    }
    fileout.close();

    rnd.SaveSeed();

    return 0;
}