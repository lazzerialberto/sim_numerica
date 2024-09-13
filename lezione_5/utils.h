#ifndef __utils__
#define __utils__

#include <string>
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

        //metrpolis step
        void Step(posizione &,FunzioneBase &);

        //find delta to have 50% acceptance rate
        void SetStepLenght(posizione &,FunzioneBase &);

};


//class to compute and print averages
class blockingaverage{

    private:
        double m_av,m_av2,m_err;
        std::string m_fileout;

    public:

        //constructor
        blockingaverage(char*argv[],int);

        //destructor
        ~blockingaverage();

        //measuring
        double measure(int,Metropolis,FunzioneBase&,posizione &);

        //computing averages
        void averages(int,int,Metropolis, FunzioneBase & ,posizione &);

        //results
        double GetAv();
        double GetErr();

};



#endif // __utils__