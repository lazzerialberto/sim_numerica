#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__

#include <cmath>

//#include "Integral.h"
#include "random.h"

class FunzioneBase {

  public:

  virtual double Eval (double) const = 0;
  virtual ~FunzioneBase() {;};

};

class gaussiana: public FunzioneBase{

  public:

    double Eval(double x) const override{
      return exp(-pow(x,2));
    }

};

//Implement specific functions for mean and standard deviation

//... class mean: public FunzioneBase ... virtual methods requieres const override

//... class std_dev: public FunzioneBase ... virtual methods requieres const override


#endif // __FunzioneBase__