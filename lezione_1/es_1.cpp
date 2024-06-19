#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "random.h"

using namespace std;

double error(double AV, double AV2, int n) {
    if (n == 0)
        return 0;
    else
        return sqrt((AV2 - pow(AV, 2)) / n);
}

int main() {
    const int N = 100;    // Number of  blocks
    const int M = 10000;   // Numer of iterations per block
    double av=0;
    double av_2=0;
    double sum2=0;
    double sum2_2=0;
    double sum_copy=0;
    double sum_copy_2=0;
    vector<double> appo2;
    int count=1;
    double appo;
    Random rnd;
    rnd.Initialize();

    //output files
    ofstream fileout("risultati_1.csv");
    ofstream fileout2("risultati_1_2.csv");

    fileout<<"sum_prog,sum2_prog,err_prog"<<endl;
    fileout2<<"sum_prog,sum2_prog,err_prog"<<endl;

    for(int i=0; i<1000000;i++){

        appo=rnd.Rannyu();
        av+=appo;
        av_2+=pow((appo-0.5),2);

        if(i/1000==count){

            //<r>
            av/=N;
            sum_copy+=av;
            sum2+=pow(av,2);

            //<r-1/2>^2
            av_2/=N;
            sum_copy_2+=av_2;
            sum2_2+=pow(av_2,2);

            //Loading up files
            fileout<<sum_copy/i*100<<","<<sum2/i*10<<","<<error(sum_copy/i*100,sum2/i*10,count)<<endl;
            fileout2<<sum_copy_2/i*100<<","<<sum2_2/i*10<<","<<error(sum_copy_2/i*100,sum2_2/i*10,count)<<endl;

            av=0;
            av_2=0;
            count++;
        }

    }

    fileout.close();
    fileout2.close();

    //Loading chi squared

    for(int i=0;i<10000000; i++){
        appo2.push_back(rnd.Rannyu());
    }
    std::cout << appo2[0] << " " << appo2.size() << endl;


    //chi2 for uniform distribution
    ofstream fileout3("risultati_1_4.csv");
    fileout3<<"chi2"<<endl;

    for(int j=0; j<1000; j++){

        int n_j=0;
        double chi2=0;

        for(int l=1; l<=N; l++){

            for(int k=0; k<M; k++){
                if(appo2[k+j*M]>=(static_cast<double>(l-1)/N) && appo2[k+j*M]<(static_cast<double>(l)/N)){
                    n_j++;
                }
            }

            std::cout << n_j << " " << static_cast<double>(l)/N << endl;
            chi2+=pow((n_j-M/N),2)/(M/N);
            n_j=0;
        }

        fileout3<<chi2<<endl;
    }

    fileout3.close();

    rnd.SaveSeed();

    return 0;
}


