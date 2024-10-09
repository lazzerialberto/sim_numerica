#ifndef __motions__
#define __motions__

#include "random.h"
#include "FunzioneBase.h"
#include "Integral.h"

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
    double step(double x_i,double t_ii,double t_i,FunzioneBase & mu,FunzioneBase & sigma);
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




#endif // __motions__