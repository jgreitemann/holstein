#ifndef SSE_H
#define SSE_H

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"

using namespace std;

#ifndef NDEBUG
template class vector<int>;
template class vector<double>;
#endif

class mc {
private:
    uint M;
    uint L;
    double T;
    uint N_el_up, N_el_down;
    double a;
    double U;
    double t;
    double mu;
    double epsilon;
    uint therm;
    vector<int> state;
    vector<int> sm;
    uint sweep;
    uint init_n_max;
    uint loop_term;
    uint n;
    double N_loop;
    bool dublon_rejected;
    vector<double> weight;
    vector<int> vtx_type;
    vector<double> prob;
    vector<int> ns;

public:    
    parser param;
    void param_init(string dir) {param.read_file(dir);}
    randomnumbergenerator *rng;
    void random_init() {
        if (param.defined("SEED")) {
            rng = new randomnumbergenerator(
                      param.value_of<luint>("SEED")
                  );
        } else {
            rng = new randomnumbergenerator();
        }
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
    int random0N(int N) {
        assert(N > 0);
        return (int)(N*rng->d());
    }

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
