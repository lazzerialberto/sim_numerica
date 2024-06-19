#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "random.h"

using namespace std;

int main(){

    int M=10000; //elements each block
    int N=100; //number of blocks

    double func_value=0;
    double func_value2=0;
    int count=1;
    double I1=0;
    double I2=0;
    double error=0;

    double I1_new=0;
    double I2_new=0;
    double error_new=0;
    double func_value_new=0;
    double func_value2_new=0;

    //initializing random numbers generator
    Random rnd;
    rnd.Initialize();

    ofstream fileout("risultati_2_1.csv");
    fileout << "I-value,I-unif,error,I-value-new,I-im-sampl,error-new" << endl;

    for(int i=1; i<=M;i++){

        double appo;
        double appo2;

        //uniform sampling
        appo=rnd.Rannyu();
        func_value+=M_PI*0.5*cos(M_PI*0.5*appo);
        func_value2+=pow(M_PI*0.5*cos(M_PI*0.5*appo),2);


        //importance sampling
        appo2=2/M_PI*acos(1-appo);
        func_value_new+=1./tan(M_PI*0.5*appo2);
        func_value2_new+=pow(1./tan(M_PI*0.5*appo2),2);


        if(i==count*N){
            I1+=func_value/(N);
            I2+=func_value2/N;
            error=sqrt((I2/(count+1)-pow(I1/(count+1),2))/(count+1));
            I1_new+=func_value_new/N;
            I2_new+=func_value2_new/N;
            error_new=sqrt((I2_new/(count+1)-pow(I1_new/(count+1),2))/(count+1));
            cout << error << endl;
            fileout <<func_value/(N)<<"," << I1/(count+1) << "," << error << "," << func_value_new/N << "," << I1_new/(count+1) << "," << error_new << endl;
            count++;
            func_value=0;
            func_value2=0;
            func_value_new=0;
            func_value2_new=0;
        }
    }

    fileout.close();

    rnd.SaveSeed();

    return 0;
}