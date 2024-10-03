#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>

#include "utils.h"

using namespace std;

Brownian :: Brownian(Random rnd){m_rand=rnd;}
Brownian :: ~Brownian(){}

//return x_i+1 for a brownian motion with constant mu and sigma
double Brownian :: step(double x_i,double t_ii,double t_i,double mu,double sigma){

    double x_ii=x_i+mu*sqrt(t_ii-t_i)+sigma*m_rand.Gauss(mu,sigma)*sqrt(t_ii-t_i);

    return x_ii;

}


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

    //cout << "Transition rate = " << double(double(n_accept)/n_total)*100<< "%"<< endl; if you want output remove comment

    return double(double(n_accept)/n_total)*100;
}


void Metropolis::Step(double &x,FunzioneBase & prob, double mu, double sigma){

    double x_new=x;

    if(n_gauss){
        x_new=x+n_rand.Gauss(0,n_delta);
    }
    else{
        x_new=x+n_rand.Rannyu(-n_delta,n_delta);
    }

    double ratio=(pow(prob.Eval(x_new,mu,sigma),2))/(pow(prob.Eval(x,mu,sigma),2));
    double A=std::min(1.,ratio);

    if(n_rand.Rannyu()<=A){
        x=x_new;
        n_accept++;
    }
    n_total++;
}

void Metropolis::SetStepLenght(double &x,FunzioneBase & prob, double mu, double sigma){

    for(int i=0;i<1000;i++){
        this->Step(x,prob,mu,sigma);
    }

    //cout << this->AccRate()<< endl;;

    while(this->AccRate()>50.01 or this->AccRate()<49.99){

        if(this->AccRate()>50.01){
            n_delta=n_delta+0.005;
            //cout << "increasing step lenght by 0.005" << endl; remove comment if you want terminal output
        }

        else{
            n_delta=n_delta-0.005;
            //cout << "decreasing step lenght by 0.005" << endl;
        }

        n_accept=0;
        n_total=0;

        for(int j=0;j<1000;j++){
            this->Step(x,prob,mu,sigma);
        }

    }

    //cout << "Final step lenght: " << n_delta << endl;
    n_accept=0;
    n_total=0;

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
    }
    n_total++;
}

void Metropolis::SetStepLenght(posizione &r,FunzioneBase & prob){

    for(int i=0;i<1000;i++){
        this->Step(r,prob);
    }

    cout << this->AccRate()<< endl;;

    while(this->AccRate()>50.01 or this->AccRate()<49.99){

        if(this->AccRate()>50.01){
            n_delta=n_delta+0.005;
            cout << "increasing step lenght by 0.005" << endl;
        }

        else{
            n_delta=n_delta-0.005;
            cout << "decreasing step lenght by 0.005" << endl;
        }

        n_accept=0;
        n_total=0;

        for(int j=0;j<1000;j++){
            this->Step(r,prob);
        }

    }

    cout << "Final step lenght: " << n_delta << endl;
    n_accept=0;
    n_total=0;

}


blockingaverage::blockingaverage(){m_av=0,m_av2=0,m_err=0;}

blockingaverage::blockingaverage(char* argv[],int n_steps){

    m_av=0,m_av2=0,m_err=0;
    std::string path_out="OUTPUT/";

    ofstream fileout(path_out+"output.dat");
    
    ofstream fileout2(path_out+"energy_mean.dat");
    m_fileout=path_out+"energy_mean.dat";

    fileout << "***************ENERGY***************" << endl;
    fileout << "#  Steps for each block: " << n_steps << "." << endl;
    fileout << "#  Transition probability: " << std::string(argv[1]) << "." << endl;;
    fileout << "#  Starting position: " << atof(argv[3]) << endl;
    fileout << "#  mu: " << atof(argv[4]) << endl;
    fileout << "#  sigma: " << atof(argv[5]) << endl;
    fileout.close();

    fileout2 << endl << "#   BLOCK:      ACTUAL_E:       AVERAGE_E:      ERROR:" << endl;

}

blockingaverage::~blockingaverage(){}

double blockingaverage::measure_E(int n_steps, Metropolis metro,FunzioneBase & prob, double &x,double mu,double sigma){

    double sum=0;

    for(int i=0; i<n_steps; i++){
        metro.Step(x,prob,mu,sigma);
        sum+=-(0.5*((pow(x-mu,2)/(pow(sigma,4)))-1./(sigma*sigma))*exp(-(pow(x-mu,2)/(2.*sigma*sigma)))+((pow(x+mu,2)/(pow(sigma,4)))-1./(sigma*sigma))*exp(-(pow(x+mu,2)/(2.*sigma*sigma))))/(prob.Eval(x,mu,sigma))+pow(x,4)-2.5*x*x;
    }

    return sum/n_steps;

}

double blockingaverage::measure(int n_steps, Metropolis metro,FunzioneBase & prob,posizione &r){

    double sum=0;

    for(int i=0; i<n_steps; i++){
        metro.Step(r,prob);
        sum+=r.GetR();
    }

    return sum/n_steps;

}

void blockingaverage::measure_savepos(int n_steps, Metropolis metro,FunzioneBase & prob,double &x,double mu,double sigma){

    ofstream fileout("./OUTPUT/positions.dat");

    fileout << "#  x" << endl;

    for(int i=0; i<n_steps; i++){
        metro.Step(x,prob,mu,sigma);
        fileout << x << endl;
    }

}

void blockingaverage::measure_savepos_2(int n_steps, Metropolis metro,FunzioneBase & prob,double &x,double mu,double sigma){

    ofstream fileout("./OUTPUT_2/positions.dat");

    fileout << "#  x" << endl;

    for(int i=0; i<n_steps; i++){
        metro.Step(x,prob,mu,sigma);
        fileout << x << endl;
    }

}

void blockingaverage::measure_savepos(int n_steps, Metropolis metro,FunzioneBase & prob,posizione &r){

    ofstream fileout("./OUTPUT/positions.dat");

    fileout << "x,y,z,r" << endl;

    for(int i=0; i<n_steps; i++){
        metro.Step(r,prob);
        fileout << r.GetX() << "," << r.GetY() << "," << r.GetZ() << "," << r.GetR() << endl;
    }

}

void blockingaverage::averages(bool final,int i_block ,int n_steps,Metropolis metro, FunzioneBase & prob,double & x,double mu,double sigma){

    if(i_block==0){
        m_av=0;
        m_av2=0;
        m_err=0;
    }

    double average=this->measure_E(n_steps,metro, prob, x,mu,sigma);
    m_av+=average;
    m_av2+=pow(average,2);
    m_err=sqrt( fabs(m_av/double(i_block+1) - pow( m_av/double(i_block+1) ,2) )/double(i_block+1) );

    if(final){
        ofstream outfile(m_fileout,ios::app);
        outfile << i_block+1 << setw(20) << average << setw(20) << m_av/double(i_block+1) << setw(20) << m_err << endl;
    }

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

double blockingaverage::GetAv(int block){return m_av/block;}

double blockingaverage::GetErr(){return m_err;}

void blockingaverage::change_fileout(std::string fileout){
    m_fileout=fileout;
    ofstream out6(fileout);
    out6 << "#   BLOCK       ACTUAL      AVE     ERROR" << endl;
}


SAnnealing::SAnnealing(int N_blocks,int M_steps,double mu,double sigma,double T,double delta_mu,double delta_sigma,Random &rnd): _rnd(rnd){
    _mu=mu;
    _sigma=sigma;
    _T=T;
    _T_start=T;
    _delta_mu=delta_mu;
    _delta_sigma=delta_sigma;
    _step=0.9;
    _N_block=N_blocks;
    _M_steps=M_steps;
    _fileout_out="./OUTPUT_2/output.dat";
    _fileout_E="./OUTPUT_2/E_min.dat";
    _fileout_pos="./OUTPUT_2/positions.dat";
    _fileout_mu="./OUTPUT_2/mu.dat";
    _fileout_sigma="./OUTPUT_2/sigma.dat";

    ofstream out1(_fileout_out);
    ofstream out2(_fileout_E);
    ofstream out3(_fileout_mu);
    ofstream out4(_fileout_sigma);
    ofstream out5(_fileout_pos);

    out1 << "***************ENERGY MINIMIZATION***************"<< endl;
    out1 << "Simulation type: Simulated Annealing + Variational MonteCarlo" << endl;
    out1 << "Number of block to estimate <H>: " << _N_block << endl;
    out1 << "Number of steps for each block: " << _M_steps << endl;
    out1 << "Starting mu: " << mu << endl;
    out1 << "Starting sigma: " << sigma << endl;
    out1 << "Starting Temperature: " << T << endl;
    out1 << "Temperature moltiplication factor in each step: " << _step << endl;

    out2 << "#      TEMP         BETA       E_AVE       ERROR" << endl;
    out3 << "#      TEMP         BETA       MU" << endl;
    out4 << "#      TEMP         BETA       SIGMA" << endl;
    out5 << "#   x" << endl;

    out1.close();
    out2.close();
    out3.close();
    out4.close();
    out5.close();
}

SAnnealing::~SAnnealing(){;}

void SAnnealing::move(FunzioneBase &boltz){

    ofstream out2(_fileout_E,ios::app);
    ofstream out3(_fileout_mu,ios::app);
    ofstream out4(_fileout_sigma,ios::app);

    double mu_purp=_mu+_rnd.Rannyu(-_delta_mu,_delta_mu);
    double sigma_purp=_sigma+_rnd.Rannyu(-_delta_sigma,_delta_sigma);

    _E=this->calculate_E("unif",0.1,mu_purp,sigma_purp,2.);

    _E_min=_E;

    if(_E[0]<=_E_old[0]){
        _mu=mu_purp;
        _sigma=sigma_purp;
    }
    else{
        if(boltz.Eval((_E[0]-_E_old[0]),_T,0)>_rnd.Rannyu()){
            _mu=mu_purp;
            _sigma=sigma_purp;
        }
        else{
            _E=_E_old;
        }
    }

    if(_E[0]<_E_min[0]){
        _E_min=_E;
        _mu_min=_mu;
        _sigma_min=_sigma;
    }

    _E_old=_E;
}

void SAnnealing::Step(int i_step,int n_steps,FunzioneBase &boltz){

    ofstream out2(_fileout_E,ios::app);
    ofstream out3(_fileout_mu,ios::app);
    ofstream out4(_fileout_sigma,ios::app);

    if(i_step==0){

        _E_old=this->calculate_E("unif",0.1,_mu,_sigma,2);

        out2 << setw(12) << _T << setw(12) << 1./_T << setw(12) << _E_old[0] << setw(12) << _E_old[1] << endl;
        out3 << setw(12) << _T << setw(12) << 1./_T << setw(12) << _mu << endl;
        out4 << setw(12) << _T << setw(12) << 1./_T << setw(12) << _sigma << endl;

        _T=_T*_step;
    }

    cout << "Step: " << i_step << endl;
    
    for(int m=0; m<n_steps; m++){
        this->move(boltz);
        Progress_Bar(m,n_steps);
    }

    cout << endl;

    _E=_E_min;
    _E_old=_E_min;
    _mu=_mu_min;
    _sigma=_sigma_min;

    out2 << setw(12) << _T << setw(12) << 1./_T << setw(12) << _E[0] << setw(12) << _E[1] << endl;
    out3 << setw(12) << _T << setw(12) << 1./_T << setw(12) << _mu << endl;
    out4 << setw(12) << _T << setw(12) << 1./_T << setw(12) << _sigma << endl;

    _T=_T*_step;

}

void SAnnealing::finalize(int count){

    ofstream out1(_fileout_out,ios::app);

    out1 << "Final Temperature: " << this->get_T() << endl;
    out1 << "Min Energy: " << this->get_E() << endl;
    out1 << "Final mu: " << this->get_mu() << endl;
    out1 << "Final sigma: " << this->get_sigma() << endl;
    out1 << "Number of iterations needed: " << count << endl;
}

double SAnnealing::get_mu()const {
    return _mu;
}

double SAnnealing::get_sigma() const {
    return _sigma;
}

double SAnnealing::get_T() const {
    return _T;
}

double SAnnealing::get_E() const {
    return _E[0];
}

double SAnnealing::get_E_old() const {
    return _E_old[0];
}

//EXCERCISE 8.2 - SPECIFIC FUNCTION TO CALCULATE ENERGY
vector<double> SAnnealing::calculate_E(string unif_gauss, double delta,double mu,double sigma,double x_start){
    gausswf prob;
    double x=x_start;
    bool gauss;
    vector<double> v;

    if(unif_gauss=="gauss"){
        gauss=true;
    }
    else if(unif_gauss=="unif"){
        gauss=false;
    }

    Metropolis metro(delta,gauss,_rnd);

    //setting delta in order to have ~50% acceptance
    metro.SetStepLenght(x,prob,mu,sigma);

    blockingaverage block;
    block.change_fileout("./OUTPUT_2/fin_energy.dat");

    for(int i=0;i<_N_block; i++){
        block.averages(true,i,_M_steps,metro,prob,x,mu,sigma);
        if(i==_N_block-1){
            v.push_back(block.GetAv(i));
            v.push_back(block.GetErr());
        }
    }

    block.measure_savepos_2(_M_steps,metro,prob,x,_mu,_sigma);

    cout << endl;
    return v;
}


void SAnnealing::Progress_Bar(int progress, int total) {

    float percentage = static_cast<float>(progress) / total;
    int pos = static_cast<int>(50 * percentage);

    std::string bar;
    for (int i = 0; i < 50; ++i) {
        if (i < pos) bar += "=";
        else if(i==pos) bar += ">";
        else bar += " ";
    }

    std::cout<< "|" << bar << "| " << setprecision(3) << int(percentage * 100.0)<< "%\r";
    std::fflush(stdout);
}

