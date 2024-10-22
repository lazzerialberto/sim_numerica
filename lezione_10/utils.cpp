#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <mpi.h>

#include "posizione.h"
#include "random.h"
#include "utils.h"

using namespace std;



element::element(int cities) {
    _ncities=cities;
    _r.resize(_ncities);
    _labels.resize(_ncities);
}

posizione element::get_pos(int i)const {
    return _r[i];
}

int element::get_label(int i)const{
    return _labels[i];
}

std::vector<int> element::get_label() const {
    return _labels;
}

void element::calculate_dist() {

    double d=0;

    for(int i=0; i<_ncities; i++){
        if(i+1!=_ncities){
            d+=_r[_labels[i]-1].Distanza(_r[_labels[i+1]-1]);
        }
        else{
            d+=_r[_labels[i]-1].Distanza(_r[_labels[0]-1]);
        }
    }

    _dist=d;
}

double element::get_dist() const {
    return _dist;
}

int element::get_n_cities() const {
    return _ncities;
}

void element::set_label(std::vector<int> new_labels){
    for(int i=0; i<_ncities; i++){
        _labels[i]=new_labels[i];
    }
}

void element::generate_initial_config(Random &_rnd,std::string circle){

    //generate initial configuration choosing random point on unitary circonference
    if(circle=="circle"){

        _labels[0]=1;
        _r[0].SetPos(1.,0.,0.);

        for(int i=1; i<_ncities; i++){
            double appo=_rnd.Rannyu();
            _labels[i]=i+1;
            _r[i].SetPos(cos(appo*2.*M_PI),sin(appo*2.*M_PI),0.);
        }
    }
    else if(circle=="square"){ //generate initial configuration inside a 1x1 square

        _labels[0]=1;
        _r[0].SetPos(0.8,0.,0.);

        for(int i=1; i<_ncities; i++){
            _labels[i]=i+1;
            _r[i].SetPos(_rnd.Rannyu(),_rnd.Rannyu(),0.);
        }

    }
    else if(circle=="italy"){ //read initial configuration from file

        ifstream infile("./INPUT/cap_prov_ita.dat");

        double x,y;

        for(int i=0; i<_ncities; i++){
            _labels[i]=i+1;
            infile >> x >> y;
            _r[i].SetPos(x,y,0);
        }
        infile.close();

    }
}

bool element::check_elem(){

    bool appo=true;

    if(_labels[0]!=1){
        std::cout << "Error: different starter point\nReject this configuration" << endl;
        appo=false;
    }
    else{
        for(int i=1; i<_ncities-1; i++){
            for(int j=i+1;j<_ncities;j++){
                if(_labels[j]==_labels[i]){
                    std::cout << "Error: city: " << _labels[j] << " in position " << i << " and " << j << "\nReject this configuration" << endl;
                    appo=false;
                }
                if(!appo){
                    break;
                }
            }
            if(!appo){
                break;
            }
        }
    }
    return appo;
}

void element::print_config(std::string filename){

    ofstream out(filename);
    out << "#         City         x            y" << endl;

    for(int i=0; i<_ncities; i++){
        out << setw(12) << _labels[i] << setw(15) << _r[_labels[i]-1].GetX() << setw(15) << _r[_labels[i]-1].GetY() << endl;
    }
}

element element::crossover(Random &_rnd, element mother){

    std::vector<int> son_labels; //final offspring labels vector
    std::vector<int> daughter_labels; //final offspring labels vector

    //help vectors
    std::vector<int> missing_son;
    std::vector<int> missing_daughter;
    std::vector<int> son_sorted;
    std::vector<int> daughter_sorted;

    int cut=static_cast<int>(_rnd.Rannyu(2,_ncities));

    for(int i=0; i<_ncities;i++){
        if(i<cut){
            son_labels.push_back(_labels[i]);
            daughter_labels.push_back(mother.get_label(i));
        }
        else{
            missing_son.push_back(_labels[i]);
            missing_daughter.push_back(mother.get_label(i));
        }
    }

    for(int i=0; i<_ncities; i++){
        for(int j=0; j<_ncities-cut; j++){
            if(missing_son[j]==mother.get_label(i)){
                son_sorted.push_back(missing_son[j]);
            }
            if(missing_daughter[j]==_labels[i]){
                daughter_sorted.push_back(_labels[i]);
            }
        }
    }

    son_labels.insert(son_labels.end(),son_sorted.begin(),son_sorted.end());
    daughter_labels.insert(daughter_labels.end(),daughter_sorted.begin(),daughter_sorted.end());

    element offspring2(mother);
    offspring2.set_label(daughter_labels);

    this->set_label(son_labels);
    bool ok=this->check_elem();
    ok=offspring2.check_elem();

    return offspring2;

}

void element::pair_permutation(Random &_rnd){

    //choose random elements to switch
    int elem1= static_cast<int>(_rnd.Rannyu(1,_ncities));
    int elem2= static_cast<int>(_rnd.Rannyu(1,_ncities));

    while(elem2==elem1){ //if they're the same element do it again
        elem2=static_cast<int>(_rnd.Rannyu(1,_ncities));
    }

    int appo=_labels[elem1];
    _labels[elem1]=_labels[elem2];
    _labels[elem2]=appo;

}

void element::shift(Random &_rnd){

    int n=static_cast<int>(_rnd.Rannyu(1,_ncities-1)); //step lenght
    int m=static_cast<int>(_rnd.Rannyu(1,_ncities)); //number of elements
    int l=static_cast<int>(_rnd.Rannyu(1,_ncities-1)); //starting point

    std::vector<int> subVector(_labels.begin() + 1, _labels.end()); //not touching first city

    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){

            //claculate index for periodicity
            int pbc1= (l+m+i-j-1)%(_ncities-1);
            int pbc2= (pbc1+1)%(_ncities-1);

            //switch element
            int appo=subVector[pbc1];
            subVector[pbc1]=subVector[pbc2];
            subVector[pbc2]=appo;
        }
    }

    subVector.insert(subVector.begin(), 1); //first city in the first position
    _labels=subVector;
}

void element::permutation(Random &_rnd) {
    int m = static_cast<int>(_rnd.Rannyu(1,_ncities/2)); // number of elements to permute
    int l = static_cast<int>(_rnd.Rannyu(1,_ncities-m));       // starting point 1 (esclude la prima città)
    int s = static_cast<int>(_rnd.Rannyu(1,_ncities-m));       // starting point 2 (esclude la prima città)

    while (l == s || (l + m > s && l < s + m)){ //checking not going over vector size
        l = static_cast<int>(_rnd.Rannyu(1,_ncities-m));
        s = static_cast<int>(_rnd.Rannyu(1,_ncities-m));
    }

    std::vector<int> temp(_labels.begin() + l, _labels.begin() + l + m);
    std::copy(_labels.begin() + s, _labels.begin() + s + m, _labels.begin() + l);
    std::copy(temp.begin(), temp.end(), _labels.begin() + s);
}

void element::inversion(Random &_rnd){

    int m=static_cast<int>(_rnd.Rannyu(1,_ncities));//number of element to invert
    int l=static_cast<int>(_rnd.Rannyu(1,_ncities));//starting point

    while(l+m>_ncities-1){
        m=static_cast<int>(_rnd.Rannyu(1,_ncities));
        l=static_cast<int>(_rnd.Rannyu(1,_ncities));
    }

    std::reverse(_labels.begin() + l, _labels.begin() + l + m);
}



Population::Population(int individuals,int cities) {

    element e(cities);
    _pop.resize(individuals,e);
    _nindividuals=individuals;
    _ncities=cities;
}

void Population::initial_pop(Random &_rnd,string circle){

    _circle=circle;
    element e(_ncities);
    e.generate_initial_config(_rnd,circle); //generate the starting configuration

    _pop[0]=e;

    std::string path_name;

    if(_circle=="circle"){
        path_name="./OUTPUT/CIRCLE/PATHS/";
    }
    else if(_circle=="square"){
        path_name="./OUTPUT/SQUARE/PATHS/";
    }
    else if(_circle=="italy"){
        path_name="./OUTPUT/ITALY/PATHS/";
    }
    _pop[0].print_config(path_name+"init_conf.dat");

    for(int i=1; i<_nindividuals; i++){

        //generating the first generation by copying the initial configuration with some mutation each with equal probability
        for(int j=0; j<5; j++){

            double a=_rnd.Rannyu();
            bool appo=true;

            //mutation of starting configuration with equal probability
            if(a<=0.25){
                e.pair_permutation(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.pair_permutation(_rnd);
                    appo=e.check_elem();
                }
            }
            else if(a>0.25 and a<=0.5){
                e.shift(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.shift(_rnd);
                    appo=e.check_elem();
                }
            }
            else if(a>0.5 and a<=0.75){
                e.inversion(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.inversion(_rnd);
                    appo=e.check_elem();
                }
            }
            else if(a>0.75){
                e.permutation(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.permutation(_rnd);
                    appo=e.check_elem();
                }
            }
        }
        
        _pop[i]=e;
        //_pop[i].print_config("./OUTPUT/config_"+std::to_string(i+1));
    }
}

void Population::initial_pop_mpi(Random &_rnd,string circle,element e){

    _circle=circle;

    _pop[0]=e;

    for(int i=1; i<_nindividuals; i++){

        //generating the first generation by copying the initial configuration with some mutation each with equal probability
        for(int j=0; j<5; j++){

            double a=_rnd.Rannyu();
            bool appo=true;

            //mutation of starting configuration with equal probability
            if(a<=0.25){
                e.pair_permutation(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.pair_permutation(_rnd);
                    appo=e.check_elem();
                }
            }
            else if(a>0.25 and a<=0.5){
                e.shift(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.shift(_rnd);
                    appo=e.check_elem();
                }
            }
            else if(a>0.5 and a<=0.75){
                e.inversion(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.inversion(_rnd);
                    appo=e.check_elem();
                }
            }
            else if(a>0.75){
                e.permutation(_rnd);
                appo=e.check_elem();
                while(!appo){
                    e.permutation(_rnd);
                    appo=e.check_elem();
                }
            }
        }
        
        _pop[i]=e;
        //_pop[i].print_config("./OUTPUT/config_"+std::to_string(i+1));
    }
}

void Population::Order_pop(){

    for(int k=0; k<_nindividuals; k++){
        _pop[k].calculate_dist();
    }

    for(int i=0; i<_nindividuals-1; i++){
        for(int j=i+1; j<_nindividuals; j++){
            if(_pop[j].get_dist()<_pop[i].get_dist()){
                element appo(_pop[i]);
                _pop[i]=_pop[j];
                _pop[j]=appo;
            }
        }
    }

    _ordered=true;
}

int Population::get_n_individuals() const {
    return _nindividuals;
}

int Population::get_n_cities() const {
    return _ncities;
}

element Population::get_element(int i)const{
    return _pop[i];
}

void Population::Set_element(int i, element e){
    _pop[i]=e;
}

void Population::print_Losses(int generation){

    ofstream out_1;
    ofstream out_2;

    if(_circle=="circle"){
        out_1.open("./OUTPUT/CIRCLE/best_loss.dat",ios::app);
        out_2.open("./OUTPUT/CIRCLE/half_loss.dat",ios::app);
    }
    else if(_circle=="square"){
        out_1.open("./OUTPUT/SQUARE/best_loss.dat",ios::app);
        out_2.open("./OUTPUT/SQUARE/half_loss.dat",ios::app);
    }
    else if(_circle=="italy"){
        out_1.open("./OUTPUT/ITALY/best_loss.dat",ios::app);
        out_2.open("./OUTPUT/ITALY/half_loss.dat",ios::app);
    }

    if(_ordered){

        double loss=0;

        for(int i=0; i<static_cast<int>(double(_nindividuals)/2); i++){
            _pop[i].calculate_dist();
            loss+=_pop[i].get_dist();
        }

        loss/=static_cast<int>(double(_nindividuals)/2);


        out_1 << setw(12) << generation << setw(20) << _pop[0].get_dist() << endl;
        out_2 << setw(12) << generation << setw(20) << loss << endl;
    }
    else{
        cout << "Error: population of the " << generation << "-th generation is not ordered" << endl;
    }

    out_1.close();
    out_2.close();
}

void Population::print_best_path(int generation){

    std::string path_name;

    if(_circle=="circle"){
        path_name="./OUTPUT/CIRCLE/PATHS";
    }
    else if(_circle=="square"){
        path_name="./OUTPUT/SQUARE/PATHS";
    }
    else if(_circle=="italy"){
        path_name="./OUTPUT/ITALY/PATHS";
    }

    string filename=path_name+"/path_"+to_string(generation)+".dat";

    if(_ordered){
        _pop[0].print_config(filename);
    }
    else{
        cout << "Error: population of the " << generation << "-th generation is not ordered" << endl;
    }
}

void Population::Mutation(Random &_rnd,double _p_pp,double _p_s,double _p_p, double _p_i){

    for(int i=0; i<_nindividuals; i++){
        if(_rnd.Rannyu()<=_p_pp){
            _pop[i].pair_permutation(_rnd);
        }
        if(_rnd.Rannyu()<=_p_s){
            _pop[i].shift(_rnd);
        }
        if(_rnd.Rannyu()<=_p_p){
            _pop[i].permutation(_rnd);
        }
        if(_rnd.Rannyu()<=_p_i){
            _pop[i].inversion(_rnd);
        }
        bool ok=_pop[i].check_elem();
    }
}

void Population::selection_and_crossover(Random &rnd, double _p_sel,double _p_cross){

    Population new_pop(_nindividuals,_ncities);

    this->Order_pop();

    for(int i=0; i<_nindividuals; i+=2){

        long double r_p=pow(rnd.Rannyu(),_p_sel);
        int sel=static_cast<int>(_nindividuals*r_p);
        
        if(rnd.Rannyu()<=_p_cross){

            element offspring1(_pop[i]);
            element parent(_pop[sel]);

            element offspring2(_pop[sel]);
            offspring2=offspring1.crossover(rnd,parent);
            new_pop.Set_element(i,offspring1);
            new_pop.Set_element(i+1,offspring2);
        }
        else{
            new_pop.Set_element(i,_pop[i]);
            new_pop.Set_element(i+1,_pop[sel]);
        }
    }

    new_pop.Order_pop();
    for(int k=0; k<_nindividuals; k++){
        _pop[k]=new_pop.get_element(k);
        bool ok=_pop[k].check_elem();
    }
}



TSP::TSP(){
    
    ifstream fin("./INPUT/input.dat");
    
    std::string property;

    std::string disp;

    while(!fin.eof()){

        fin>>property;

        if(property=="Configuration"){
            fin >> disp;
            if(disp=="circle"){
                _circle="circle";
                _path_name="./OUTPUT/CIRCLE/";
            }
            else if(disp=="square"){
                _circle="square";
                _path_name="./OUTPUT/SQUARE/";
            }
            else if(disp=="italy"){
                _circle="italy";
                _path_name="./OUTPUT/ITALY/";
            }
        }
        else if(property=="Individuals"){
            fin >> _nindivid;
        }
        else if(property=="Cities"){
            if(_circle=="circle" or _circle=="square"){
                fin >> _ncities;
            }
            else if(_circle=="italy"){
                ifstream fin("./INPUT/prov_ita.txt");

                std::string appo;

                while(std::getline(fin, appo)){
                    _ncities++;
                }
            }
        }
        else if(property=="Generations"){
            fin >> _ngen;
        }
        else if(property=="Loss"){
            fin >> _norm;
        }
        else if(property=="Selection_factor"){
            fin >> _p_sel;
        }
        else if(property=="Crossover_prob"){
            fin >> _p_cross;
        }
        else if(property=="Pair_perm_prob"){
            fin >> _p_pp;
        }
        else if(property=="Perm_prob"){
            fin >> _p_p;
        }
        else if(property=="Shift_prob"){
            fin >> _p_s;
        }
        else if(property=="Inversion_prob"){
            fin >> _p_i;
        }
    }
    
    fin.close();
}

void TSP::print_output(){

    _out1.open(_path_name+"output.dat");

    _out1 << "**************** TRAVELING SALESMAN PROBLEM *********************" << endl;
    _out1 << "Disposition: " << _circle << endl;
    _out1 << "Number of individuals: " << _nindivid << endl;
    _out1 << "Number of cities: " << _ncities << endl;
    _out1 << "Number of generations: " << _ngen << endl;
    _out1 << "Loss used for minimization: " << _norm << endl;
    _out1 << "Selection factor: " << _p_sel << endl;
    _out1 << "Crossover probability: " << _p_cross << endl;
    _out1 << "Pair permutation probability: " << _p_pp << endl;
    _out1 << "Permutation probability: " << _p_p << endl;
    _out1 << "Shift probability: " << _p_s << endl;
    _out1 << "Inversion probability: " << _p_i << endl;

    _out1.close();


    ofstream out2(_path_name+"best_loss.dat");
    ofstream out3(_path_name+"half_loss.dat");

    out2 << "#        GENERATION          LOSS" << endl;
    out3 << "#        GENERATION          LOSS" << endl;

    out2.close();
    out3.close();

}

void TSP::print_output(int rank,int processes){

    if(rank==0){

        _out1.open(_path_name+"output.dat");

        _out1 << "**************** TRAVELING SALESMAN PROBLEM *********************" << endl;
        _out1 << "Disposition: " << _circle << endl;
        _out1 << "Number of individuals: " << _nindivid << endl;
        _out1 << "Number of cities: " << _ncities << endl;
        _out1 << "Number of generations: " <<  _ngen << endl;
        _out1 << "Loss used for minimization: " <<  _norm << endl;
        _out1 << "Selection factor: " << _p_sel << endl;
        _out1 << "Crossover probability: " <<  _p_cross << endl;
        _out1 << "Pair permutation probability: " <<  _p_pp << endl;
        _out1 << "Permutation probability: " << _p_p << endl;
        _out1 << "Shift probability: " << _p_s << endl;
        _out1 << "Inversion probability: " << _p_i << endl;
        _out1 << "Number of processes: " << processes << endl;
        _out1 << "Printing only rank: " << rank+1 << endl;

        _out1.close();


        ofstream out2(_path_name+"best_loss.dat");
        ofstream out3(_path_name+"half_loss.dat");

        out2 << "#        GENERATION          LOSS" << endl;
        out3 << "#        GENERATION          LOSS" << endl;

        out2.close();
        out3.close();
    }

}

void TSP::finalize(int rank, double dt){

    _out1.open(_path_name+"output.dat",ios::app);
    _out1 << "Time for process " << rank+1 << ": " << dt << " s" << endl;
    _out1.close();
}

std::string TSP::get_circle() const {
    return _circle;
}

int TSP::get_individuals() const {
    return _nindivid;
}

int TSP::get_cities() const {
    return _ncities;
}

int TSP::get_generations() const {
    return _ngen;
}

double TSP::get_p_sel() const {
    return _p_sel;
}

double TSP::get_p_cross() const {
    return _p_cross;
}

double TSP::get_p_pp() const {
    return _p_pp;
}

double TSP::get_p_p() const {
    return _p_p;
}

double TSP::get_p_s() const {
    return _p_s;
}

double TSP::get_p_i() const {
    return _p_i;
}



