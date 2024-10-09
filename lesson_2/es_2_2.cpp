#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "posizione.h"

using namespace std;

int main(){

//DISCRETE 3D RW

    ifstream filein1("trajs.dat");

    vector<posizione> pos1;
    double x,y,z;
    
    while(!filein1.eof()){
        filein1 >> x;
        filein1 >> y;
        filein1 >> z;
        pos1.push_back(posizione(x,y,z));
    }

    filein1.close();

    ofstream fileout1("risultati_2_2.csv");
    fileout1 << "rw,err" << endl;

    //blocking average and error

    for(int i=0; i<100; i++){

        double sum=0;
        double sum2=0;
        double error=0; 
        double appo=0;

        for(int j=0;j<10000;j++){
            appo+=pow(pos1[i+j*100].Distanza(posizione(0,0,0)),2);
        
            if((j+1)%100==0){
                sum+=sqrt(appo/100.);
                sum2+=appo/100.;
                appo=0;
            }
        }

        error=sqrt((sum2/(100)-pow(sum/(100),2))/(100));
        fileout1 << sum/100 << "," << error << endl;
    }

    fileout1.close();

//CONTINUUM 3D RW

    ifstream filein("trajs_con.dat");

    vector<posizione> pos;
    double x,y,z;
    
    while(!filein.eof()){
        filein >> x;
        filein >> y;
        filein >> z;
        pos.push_back(posizione(x,y,z));
    }

    filein.close();

    ofstream fileout("risultati_2_2_con.csv");
    fileout << "rw,err" << endl;

    //blocking average and error

    for(int i=0; i<100; i++){

        double sum=0;
        double sum2=0;
        double error=0; 
        double appo=0;

        for(int j=0;j<10000;j++){
            appo+=pow(pos[i+j*100].Distanza(posizione(0,0,0)),2);
        
            if((j+1)%100==0){
                sum+=sqrt(appo/100.);
                sum2+=appo/100.;
                appo=0;
            }
        }

        error=sqrt((sum2/(100)-pow(sum/(100),2))/(100));
        fileout << sum/100 << "," << error << endl;
    }

    fileout.close();


    return 0;
}