#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "posizione.h"
#include "random.h"

using namespace std;

int main(){

    int n_steps=100; //number of steps of each trajectory
    int M=10000; //number of different trajectories
    int count=0;

    Random rnd;
    rnd.Initialize();

    //discrete trajectories files
    ofstream fileout("trajectories.csv");
    ofstream fileout2("trajs.dat");
    fileout << "x,y,z" << endl;

    for(int i=0; i<M; i++){

        vector<posizione> pos;
        pos.push_back(posizione(0,0,0));

        for(int j=1; j<=n_steps; j++){

            fileout << pos[j-1].GetX() << "," << pos[j-1].GetY() << "," << pos[j-1].GetZ() << endl;
            fileout2 << pos[j-1].GetX() << endl << pos[j-1].GetY() << endl << pos[j-1].GetZ() << endl;

            double appo, appo1;
            appo=rnd.Rannyu();
            appo1=rnd.Rannyu();

            if(appo1>=0.5){
                appo1=1; //step forward
            }
            if(appo1<0.5){
                appo1=-1; //step backward
            }

            //sorting the step direction

            if(static_cast<int>(appo*3)==0){
                pos.push_back(posizione((pos[j-1].GetX()+appo1),pos[j-1].GetY(),pos[j-1].GetZ()));
            }
            
            if(static_cast<int>(appo*3)==1){
                pos.push_back(posizione(pos[j-1].GetX(),(pos[j-1].GetY()+appo1),pos[j-1].GetZ()));
            }

            if(static_cast<int>(appo*3)==2){
                pos.push_back(posizione(pos[j-1].GetX(),pos[j-1].GetY(),(pos[j-1].GetZ()+appo1)));
            }

            else{
                cout << "not doing anything, any problem?" << endl;
                cout << appo1 << " , "<< static_cast<int>(appo*3) << endl;
                count ++;
            }

        }

    }

    cout << count << endl;

    fileout.close();
    fileout2.close();

    rnd.SaveSeed();

    return 0;
}