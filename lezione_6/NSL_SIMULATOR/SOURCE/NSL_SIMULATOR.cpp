/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <cmath>
#include "system.h"

using namespace std;

void Progress_Bar(int progress, int total, int bar_width=50) {

    float percentage = static_cast<float>(progress) / total;
    int pos = static_cast<int>(bar_width * percentage);

    std::string bar;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) bar += "=";
        else if(i==pos) bar += ">";
        else bar += " ";
    }

    std::cout<< "|" << bar << "| " << setprecision(3) << int(percentage * 100.0)<< "%\r";
    std::fflush(stdout);
}

int main (int argc, char * argv[]){

  if(argc!=2){
    cout << "Error: temperature step is required" << endl;
    exit(-1);
  }

  double step_lenght= atof(argv[1]);

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  for(int k=0; k<1000; k++){
    SYS.step();
    SYS.measure();
  }

  for(int l=0;l<= static_cast<int>((2.0-0.5)/step_lenght);l++){

    SYS.read_configuration();

    if(l!=0){
      double temp=SYS.get_temp()+step_lenght;
      SYS.set_temp(temp);
    }

    SYS.block_reset(0);

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);


    }
    SYS.finalize();
    //adding progress bar
    Progress_Bar(l,static_cast<int>((2.0-0.5)/step_lenght));
  }
  cout << endl;

  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
