#ifndef __utils__
#define __utils__

#include <string>
#include <vector>
#include <algorithm>
#include "random.h"
#include "FunzioneBase.h"
#include "posizione.h"

//this class contains functions to sort a Brownian motion
class Brownian{

private:
    Random m_rand;

public:

    //Default constructor
    Brownian(Random rnd);

    //Default destructor
    ~Brownian();

    //Brownian step with constant mean and standard deviation
    double step(double x_i,double t_ii,double t_i,double mu,double sigma);

    //Brownian step with time-dependent mean and variance
    //double step(double x_i,double t_ii,double t_i,FunzioneBase & mu,FunzioneBase & sigma);
};


//this class contains functions to sort a Geometric Brownian Motion (GBM)
class GBM{

    private:

        Random m_rand;

    public:

        //Default constructor
        GBM();

        //Default destructor
        ~GBM();

        //Brownian step with constant mean and standard deviation
        double step(double x_i,double t_ii,double t_i,double mu,double sigma);

};


//this class is for the Metropolis algorithm
class Metropolis{

    private:
        int n_accept,n_total;
        double n_delta;
        bool n_gauss;
        Random &n_rand;
    
    public:

        //default constructor
        Metropolis(double,bool,Random &);

        //default destructor
        ~Metropolis();

        //return the acceptance rate
        double AccRate();

        //metrpolis step (1D)
        void Step(double &x,FunzioneBase & prob, double mu, double sigma);

        //find delta to have 50% acceptance rate (1D)
        void SetStepLenght(double &x,FunzioneBase &prob,double mu,double sigma);

        //metrpolis step (3D)
        void Step(posizione &,FunzioneBase &);

        //find delta to have 50% acceptance rate (3D)
        void SetStepLenght(posizione &,FunzioneBase &);

};


//class to compute and print averages
class blockingaverage{

    private:
        double m_av,m_av2,m_err;
        std::string m_fileout;

    public:

        blockingaverage();

        //constructor
        blockingaverage(char*argv[],int);

        //destructor
        ~blockingaverage();

        //measuring
        double measure(int,Metropolis,FunzioneBase&,posizione &);

        //measuring energy
        double measure_E(int n_steps,Metropolis metro,FunzioneBase&prob,double &x,double mu, double sigma);

        //measuring and saving positions
        void measure_savepos(int,Metropolis,FunzioneBase&,posizione &);

        //measuring and saving positions
        void measure_savepos(int n_steps,Metropolis metro,FunzioneBase& prob,double & x,double mu,double sigma);

        //measuring and saving positions ex_8_2
        void measure_savepos_2(int n_steps,Metropolis metro,FunzioneBase& prob,double & x,double mu,double sigma);

        //computing averages
        void averages(int,int,Metropolis, FunzioneBase & ,posizione &);
        void averages(bool final,int i_block ,int n_steps,Metropolis metro, FunzioneBase & prob,double & x,double mu,double sigma);

        //results
        double GetAv(int block);
        double GetErr();

        void change_fileout(std::string fileout);
};

class SAnnealing{

    public:

        SAnnealing(int N_blocks,int M_steps,double mu, double sigma,double T,double delta_mu,double delta_sigma,Random &rnd);

        ~SAnnealing();

        // move of Metropolis algorithm
        void move(FunzioneBase &boltz);

        void Step(int i_step,int n_steps,FunzioneBase &boltz);

        std::vector<double> calculate_E(std::string unif_gauss, double delta,double mu, double sigma,double x_start);

        void finalize(int);

        double get_T()const;
        double get_mu() const;
        double get_sigma() const;
        double get_E()const;
        double get_E_old()const;

        void Progress_Bar(int progress, int total);

    private:
        Random _rnd;
        double _mu,_sigma,_T,_T_start,_step,_mu_min,_sigma_min;
        double _delta_mu,_delta_sigma;
        std::vector<double> _E_old;
        std::vector<double> _E;
        std::vector<double> _E_min;
        int _N_block,_M_steps;
        std::string _fileout_out;
        std::string _fileout_E;
        std::string _fileout_pos;
        std::string _fileout_mu;
        std::string _fileout_sigma;

};



#endif // __utils__