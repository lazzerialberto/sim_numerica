#ifndef __INTEGRAL_H__
#define __INTEGRAL_H__

#include <iostream>
#include <cmath>
#include <vector>
#include "random.h"
#include "FunzioneBase.h"

using namespace std;

class Integral {

 public:
  
  Integral (double a, double b){
    checkInterval (a,b);
    m_nstep = 0;
    m_h = 0; 
    m_sum = 0;
    m_integral =0;
  }

  virtual double Integra(unsigned int N, FunzioneBase& f, double fmax) = 0;

 protected:

  void checkInterval (double a, double b){
		if(a<b){
			m_sign=1;
			m_a=a;
			m_b=b;
		}
		if(a>b){
			m_sign=-1;
			m_a=b;
			m_b=a;
		}
		if(a==b){
			cerr << "Errore: intervallo degenere. ESCO" << endl;
			exit(-3);
		}
  }

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
	double m_prec, m_err;
};

class MonteCarloMedia: public Integral{

public:

	MonteCarloMedia( double a, double b, Random rnd): Integral(a,b) {m_rand=rnd;};
	~MonteCarloMedia(){};

	virtual double Integra(unsigned int N, FunzioneBase& f, double fmax=0){

			double sum=0;
			vector<double> v;
			for(int i=0; i<N; i++){
				double x=m_rand.Rannyu(m_a,m_b);
				v.push_back(f.Eval(x));
				sum+=v[i];
			}
			return ((m_b-m_a)*sum/N);
		};


private:

	Random m_rand;

};

class MonteCarloHoM: public Integral{

public:

	MonteCarloHoM(double a, double b,Random rnd): Integral(a,b) {m_rand=rnd;};

	virtual double Integra(unsigned int N, FunzioneBase& f, double fmax){

			int N_tot=0;
			int N_hit=0;
			for(int i=0; i<N; i++){
				double x=m_rand.Rannyu(m_a,m_b);
				double y=m_rand.Rannyu(0,fmax);
				if(y<=f.Eval(x)){
					N_hit++;
				}
				N_tot++;
			}
			return ((m_b-m_a)*fmax*(static_cast<double>(N_hit)/static_cast<double>(N_tot)));

		};


private:

	Random m_rand;

};


#endif //__INTEGRAL_H__