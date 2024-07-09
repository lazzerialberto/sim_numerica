#include <iostream>
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

int main (int argc, char *argv[]){

  ofstream fileout;
  fileout.open("../OUTPUT/parameters.dat",ios::app);

  //fileout << "CUTOFF RAD:     DENSITY:     FINAL TEMP:     INITIAL TEMP:" << endl;

  vector<double> temp={0.6,0.9,1.1,1.6,2.0,2.2}; //riprendo da temp=2.0
  vector<double> r_cut={2.0,2.2,2.5,2.7,4.5,5.0,5.5};
  vector<double> rho={0.02,0.05,0.1,0.7,0.8,0.9,1.0,1.1,1.2};

  for(int l=4;l<temp.size(); l++){

    for(int m=0;m<r_cut.size(); m++){

      for(int k=0;k<rho.size();k++){

        int nconf = 1;
        System SYS;
        SYS.initialize_for_equilibration(temp[l],r_cut[m],rho[k]);
        SYS.initialize_properties();
        SYS.block_reset(0);

        for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
          for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
            SYS.step();
            SYS.measure();
          }
          SYS.averages(i+1);
          SYS.block_reset(i+1);
        }

        SYS.print_parameters(temp[l]);
        SYS.finalize();
      
      }
    }
    //adding progress bar
    Progress_Bar(l,temp.size());
  }
  cout << endl;

  return 0;
}

