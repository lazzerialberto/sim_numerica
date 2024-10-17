#include "posizione.h"
#include "random.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>


class element{

    public:

        element(int cities); //constructor

        posizione get_pos(int i) const;
        int get_label(int i)const;
        double get_dist() const;
        void calculate_dist(); //compute distance in this configuration
        int get_n_cities() const;

        void set_label(std::vector<int> new_labels); //setting manually the order of the cities

        void generate_initial_config(Random &_rnd, bool circle);

        bool check_elem(); //check function
        void print_config(std::string filename); //printing configuration

        element crossover(Random &_rnd,element mother);

        //mutations
        void pair_permutation(Random &_rnd);
        void shift(Random &_rnd);
        void permutation(Random &_rnd);
        void inversion(Random &_rnd);

    private:

        std::vector<posizione> _r; //contains the position of cities ordered for the entire simulation by the initial configuration
        //to access an element use _r[_labels[i]-1]
        std::vector<int> _labels; //contains the order of the cities, this is the vector which is being mutated during the simulation
        double _dist; //total distance of this specific configuration
        int _ncities;
};


class Population{

    public:

        Population(int individuals, int cities); //constructor

        void initial_pop(Random &_rnd,bool circle); //create the first generation

        void Order_pop(); //ordering _pop depending on L(2)

        int get_n_individuals() const;
        int get_n_cities() const;

        element get_element(int i)const; //getting i-th element from _pop
        void Set_element(int i,element e); //modifiyng i-th element from _pop

        void print_Losses(int generation);
        void print_best_path(int generation);

        void Mutation(Random &_rnd,double _p_pp,double _p_s,double _p_p, double _p_i);
        void selection_and_crossover(Random &rnd, double _p_sel,double _p_cross);

    private:
        std::vector<element> _pop; //vector with type <element> describing a population
        int _nindividuals,_ncities;
        bool _ordered; // tells if the population has been already ordered
        bool _circle;
};


class TSP {

    public:

        TSP();

        //methods to get variables
        bool get_circle() const;
        int get_individuals() const;
        int get_cities() const;
        int get_generations() const;
        double get_p_sel() const;
        double get_p_cross() const;
        double get_p_pp() const;
        double get_p_p() const;
        double get_p_s() const;
        double get_p_i() const;
    
    private:

        std::ofstream _out1;
        int _nindivid,_ncities,_ngen;
        bool _circle;
        double _p_sel,_p_cross,_p_s, _p_pp, _p_p,_p_i; // probabilities: selection, crossover, shift, pair permutation, permutation, inversion

};