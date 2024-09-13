#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__

#include <cmath>

//#include "Integral.h"
#include "random.h"

class FunzioneBase {

  public:

  virtual double Eval (double,double,double) const = 0;
  virtual ~FunzioneBase() {;};

};

class gaussiana: public FunzioneBase{

  public:

    double Eval(double x,double y,double z) const override{
      double r=sqrt(x*x+y*y+z*z);
      return exp(-pow(r,2));
    }

};

//ground state wavefunction (a_0 esxpressed in nm=10^-9 m)
class groundwf: public FunzioneBase{

  public:

    double Eval(double x,double y,double z) const override{
      double r=sqrt(x*x+y*y+z*z);
      return (1./sqrt(M_PI*pow(0.0529,3)))*exp(-r);
    }

};

//excited state wavefunction (a_0 esxpressed in nm=10^-9 m)
class excitedwf: public FunzioneBase{

  public:

    double Eval(double x,double y,double z) const override{
      double r=sqrt(x*x+y*y+z*z);
      return (1./sqrt(32.*M_PI*pow(0.0529,3)))*exp(-r/2.)*r*(z/r);
    }

};

#endif // __FunzioneBase__