#ifndef SSE_H
#define SSE_H

#include <iostream>
#include <vector>
#include <string>
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"

using namespace std;

class mc {
private:
    int M;
    int L;
    double T;
    double filling;
    double a;
    double U;
    double t;
    double epsilon;
    int therm;
    vector<int> state;
    vector<int> sm;
    int sweep;
    int init_n_max;
    int loop_term;
    int n;
    vector<double> weight;
    vector<int> vtx_type;
    vector<double> prob;

public:    
    parser param;
    void param_init(string dir) {param.read_file(dir);}
    randomnumbergenerator * rng;
    void random_init() {
        if (param.defined("SEED"))
            rng = new randomnumbergenerator(param.value_of<luint>("SEED"));
        else
            rng = new randomnumbergenerator();
    }
    void random_write(odump& d) {rng->write(d);}
    void seed_write(string fn) {
        ofstream s;
        s.open(fn.c_str());
        s << rng->seed()<<endl;
        s.close();
    }
    void random_read(idump& d) {
        rng = new randomnumbergenerator();
        rng->read(d);
    }
    void random_clear() {delete rng;}
    double random01() {return rng->d();}
    int random0N(int N) {return rng->i(N);}

    void init();
    void do_update();
    void do_measurement();
    void write(string);
    bool read(string);
    void write_output(string);
    
    bool is_thermalized();
    measurements measure;    
    
    mc(string);
    ~mc();
};

#endif
