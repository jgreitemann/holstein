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

typedef unsigned char byte;

enum el_state {
    empty,
    up,
    down,
    dublon
};

enum leg {
    bottom_left,
    bottom_right,
    top_right,
    top_left
};

enum worm_type {
    up_worm,
    down_worm,
    dublon_worm
};

enum operator_type {
    electron_diag,
    up_hopping,
    down_hopping,
    creator_left,
    creator_right,
    annihilator_left,
    annihilator_right,
    phonon_diag
};

struct bond_operator {
    operator_type type  : 3;
    unsigned short bond : 13;
};

const bond_operator identity = {electron_diag, 0};

union vertex {
    struct {
        el_state bottom_left  : 2;
        el_state bottom_right : 2;
        el_state top_right    : 2;
        el_state top_left     : 2;
    };
    byte int_repr;
    el_state get_state(leg l) { return static_cast<el_state>(int_repr>>(2*l) & 3); }
};

union assignment {
    struct {
        vertex vtx;
        leg ent_leg    : 2;
        leg exit_leg   : 2;
        worm_type worm : 4;
    };
    unsigned short int_repr;
    vertex flipped_vtx() {
        vertex res;
        res.int_repr = vtx.int_repr ^ ((worm+1) << (2*ent_leg))
                                    ^ ((worm+1) << (2*exit_leg));
        return res;
    };
};

struct subseq_node {
    int i;
    int Nd;
    int m;
    double r;
};

enum lock_flag {
    unlocked,
    left_lock,
    right_lock,
    total_lock
};

union list_position {
    struct {
        leg vtx_leg : 2;
        int index   : 30;
    };
    int int_repr;
    list_position() {}
    list_position(leg vtx_leg, int index) : vtx_leg(vtx_leg), index(index) {}
    bool operator==(list_position& lp) { return int_repr==lp.int_repr; }
};

const list_position invalid_pos = {bottom_left, -1};

#ifndef NDEBUG
template class vector<int>;
template class vector<double>;
template class vector<el_state>;
template class vector<bond_operator>;
template class vector<lock_flag>;
template class vector<list_position>;
#endif

class mc {
private:
    uint M;
    uint Np;
    uint L;
    double T;
    uint N_el_up, N_el_down;
    double a;
    double U;
    double mu;
    double omega;
    double g;
    double epsilon;
    uint therm;
    vector<el_state> state;
    vector<int> occ;
    vector<bond_operator> sm;
    uint sweep;
    uint init_n_max;
    uint loop_term;
    uint n;
    double vtx_visited;
    uint N_loop;
    bool dublon_rejected;
    vector<double> weight;
    vector<int> vtx_type;
    vector<double> prob;
    vector<int> ns;
    double avg_worm_len;
    uint worm_len_sample_size;

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
