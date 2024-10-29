#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "motions.h"
#include "random.h"
#include "FunzioneBase.h"
#include "Integral.h"

using namespace std;

Brownian :: Brownian(Random rnd){m_rand=rnd;}
Brownian :: ~Brownian(){}

//return x_i+1 for a brownian motion with constant mu and sigma
double Brownian :: step(double x_i,double t_ii,double t_i,double mu,double sigma){

    double x_ii=x_i+mu*sqrt(t_ii-t_i)+sigma*m_rand.Gauss(mu,sigma)*sqrt(t_ii-t_i);

    return x_ii;

}

//return x_i+1 for a brownian motion with time-dependent mu and sigma_squared
double Brownian :: step(double x_i,double t_ii,double t_i,FunzioneBase & mu,FunzioneBase & sigma_squared){

    MonteCarloMedia I(t_i,t_ii,m_rand);

    double x_ii=x_i+I.Integra(10000,mu)+m_rand.Gauss(mu.Eval(t_ii),sqrt(sigma_squared.Eval(t_ii)))*sqrt(I.Integra(10000,sigma_squared));

    return x_ii;

}


//constructor and destructor for GBM
GBM :: GBM(){m_rand.Initialize();};
GBM :: ~GBM(){m_rand.SaveSeed();};

//return x_i+1 for a geometric brownian motion with constant mu and sigma
double GBM :: step(double x_i,double t_ii,double t_i,double mu,double sigma){

    double x_ii=x_i*exp((mu-0.5*pow(sigma,2))*(t_ii-t_i)+sigma*m_rand.Gauss(0,1)*sqrt(t_ii-t_i));

    return x_ii;

}