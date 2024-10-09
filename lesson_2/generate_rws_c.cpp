#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "posizione.h"
#include "random.h"

using namespace std;

int main(){

    int n_steps=100; //number of steps
    int M=10000; //number of different trajectories

    Random rnd;
    rnd.Initialize();

    //output files
    ofstream fileout("trajectories_cont.csv");
    ofstream fileout2("trajs_con.dat");

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

            double theta=0;
            double phi=0;

            //sorting angle between 0 e 2pi
            theta=2*M_PI*appo;
            phi=acos(1-2*appo1);

            pos.push_back(posizione((pos[j-1].GetX()+sin(phi)*cos(theta)),(pos[j-1].GetY()+sin(phi)*sin(theta)),(pos[j-1].GetZ()+cos(phi))));
        }

    }

    fileout.close();
    fileout2.close();

    rnd.SaveSeed();

    return 0;
}