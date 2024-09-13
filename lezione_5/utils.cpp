#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "utils.h"

using namespace std;

Brownian :: Brownian(Random rnd){m_rand=rnd;}
Brownian :: ~Brownian(){}

//return x_i+1 for a brownian motion with constant mu and sigma
double Brownian :: step(double x_i,double t_ii,double t_i,double mu,double sigma){

    double x_ii=x_i+mu*sqrt(t_ii-t_i)+sigma*m_rand.Gauss(mu,sigma)*sqrt(t_ii-t_i);

    return x_ii;

}

//return x_i+1 for a brownian motion with time-dependent mu and sigma_squared
/*double Brownian :: step(double x_i,double t_ii,double t_i,FunzioneBase & mu,FunzioneBase & sigma_squared){

    MonteCarloMedia I(t_i,t_ii,m_rand);

    double x_ii=x_i+I.Integra(10000,mu)+m_rand.Gauss(mu.Eval(t_ii),sqrt(sigma_squared.Eval(t_ii)))*sqrt(I.Integra(10000,sigma_squared));

    return x_ii;

}*/


//constructor and destructor for GBM
GBM :: GBM(){m_rand.Initialize();};
GBM :: ~GBM(){m_rand.SaveSeed();};

//return x_i+1 for a geometric brownian motion with constant mu and sigma
double GBM :: step(double x_i,double t_ii,double t_i,double mu,double sigma){

    double x_ii=x_i*exp((mu-0.5*pow(sigma,2))*sqrt(t_ii-t_i)+sigma*m_rand.Gauss(0,1)*sqrt(t_ii-t_i));

    return x_ii;

}



//constructor and destructor for Metropolis
Metropolis::Metropolis(double delta,bool gauss,Random &rnd): n_rand(rnd) {n_accept=0,n_total=0,n_delta=delta, n_gauss=gauss;}
Metropolis::~Metropolis(){}

double Metropolis::AccRate(){

    cout << "Transition rate = " << double(double(n_accept)/n_total)*100<< "%"<< endl;

    return double(double(n_accept)/n_total)*100;
}

void Metropolis::Step(posizione &r,FunzioneBase & prob){

    posizione r_c=r;

    if(n_gauss){
        r_c.AddPos(n_rand.Gauss(0,n_delta),n_rand.Gauss(0,n_delta),n_rand.Gauss(0,n_delta));
    }
    else{
        r_c.AddPos(n_rand.Rannyu(-n_delta,n_delta),n_rand.Rannyu(-n_delta,n_delta),n_rand.Rannyu(-n_delta,n_delta));
    }

    double ratio=(pow(prob.Eval(r_c.GetX(),r_c.GetY(),r_c.GetZ()),2))/(pow(prob.Eval(r.GetX(),r.GetY(),r.GetZ()),2));
    double A=std::min(1.,ratio);

    if(n_rand.Rannyu()<=A){
        r=r_c;
        n_accept++;
        //cout << r.GetR() << endl;
    }
    n_total++;
}

void Metropolis::SetStepLenght(posizione &r,FunzioneBase & prob){

    for(int i=0;i<1000;i++){
        this->Step(r,prob);
    }

    cout << this->AccRate()<< endl;;

    while(this->AccRate()>50.01 or this->AccRate()<49.99){

        if(this->AccRate()>50.1){
            n_delta=n_delta+0.005;
            cout << "increasing step lenght by 0.005" << endl;
        }

        else{
            n_delta=n_delta-0.01;
            cout << "decreasing step lenght by 0.005" << endl;
        }

        n_accept=0;
        n_total=0;

        for(int j=0;j<10000;j++){
            this->Step(r,prob);
        }

    }

    cout << "Final step lenght: " << n_delta << endl;
    n_accept=0;
    n_total=0;

}


blockingaverage::blockingaverage(char* argv[],int n_steps){

    m_av=0,m_av2=0,m_err=0;
    std::string path_out="OUTPUT/";
    
    ofstream fileout;

    if(std::string(argv[3])=="g"){
        fileout.open(path_out+"out_ground.dat");
        m_fileout=path_out+"out_ground.dat";
        fileout << "***************************************************"<< endl << "            HYDROGEN ATOM SIMULATION\n***************************************************" << endl << endl;
        fileout << "Estimation of <r> in hydrogen atom ground state."<< endl;
    }
    else if(std::string(argv[3])=="e"){
        fileout.open(path_out+"out_excited.dat");
        m_fileout=path_out+"out_excited.dat";
        fileout << "***************************************************"<< endl << "            HYDROGEN ATOM SIMULATION\n***************************************************" << endl << endl;
        fileout << "Estimation of <r> in hydrogen atom first excited state."<< endl;
    }

    fileout << "Number of steps for each block: " << n_steps << "." << endl;
    fileout << "Transition probability: " << std::string(argv[1]) << "." << endl;;
    fileout << "Starting position: (" << atof(argv[4]) << ", " << atof(argv[5]) << ", " << atof(argv[6]) << ")a_0." << endl;

    fileout << endl << "BLOCK:      ACTUAL VALUE:       AVERAGE VALUE:      ERROR:" << endl;

}

blockingaverage::~blockingaverage(){}

double blockingaverage::measure(int n_steps, Metropolis metro,FunzioneBase & prob,posizione &r){

    double sum=0;

    for(int i=0; i<n_steps; i++){
        metro.Step(r,prob);
        sum+=r.GetR();
    }

    return sum/n_steps;

}

void blockingaverage::averages(int i_block ,int n_steps,Metropolis metro, FunzioneBase & prob,posizione & r){

    if(i_block==0){
        m_av=0;
        m_av2=0;
        m_err=0;
    }

    double average=this->measure(n_steps,metro, prob, r);
    m_av+=average;
    m_av2+=pow(average,2);
    m_err=sqrt( fabs(m_av/double(i_block+1) - pow( m_av/double(i_block+1) ,2) )/double(i_block+1) );

    ofstream outfile(m_fileout,ios::app);

    outfile << i_block+1 << setw(20) << average << setw(20) << m_av/double(i_block+1) << setw(20) << m_err << endl;

}

double blockingaverage::GetAv(){return m_av;}

double blockingaverage::GetErr(){return m_err;}

