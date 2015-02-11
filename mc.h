#ifndef SSE_H
#define SSE_H

#define N_BOND 8                // number of bond operator flavors
#define NB (L)                  // number of bonds
#define COORD 2                 // coordination number
#define LEFT_SITE(b) (b-1)      // lattice -
#define RIGHT_SITE(b) ((b)%L)   //           geometry
#define LEFT_BOND(s) ((s+L-1)%L+1)
#define RIGHT_BOND(s) (s+1)
#define N_WORM 3                // number of worm types

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
    inline bool operator== (const bond_operator& other) {
        return bond == other.bond && type == other.type;
    }
    inline bool operator!= (const bond_operator& other) {
        return !(*this == other);
    }
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
    inline bool operator== (const list_position& lp) {
        return int_repr == lp.int_repr;
    }
    inline bool operator!= (const list_position& lp) {
        return int_repr != lp.int_repr;
    }
};

const list_position invalid_pos = list_position(bottom_left, -1);

#ifndef NDEBUG
template class vector<int>;
template class vector<double>;
template class vector<el_state>;
template class vector<bond_operator>;
template class vector<lock_flag>;
template class vector<list_position>;
template class vector<operator_type>;
#endif

class mc {
private:
    uint M;
    uint Np;
    uint L;
    double T;
    uint N_el_up, N_el_down;
    double enlargement_factor;
    double U;
    double mu;
    double omega;
    double g;
#ifdef MCL_PT
    int label;
    int myrep;
    int pt_spacing;
    vector<double> gvec;
    vector<double> muvec;
#endif
    double epsilon;
    double q_S, q_chi;
    uint therm, total_therm;
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
#ifdef MCL_PT
    vector<double> other_weight;
#endif
    vector<operator_type> vtx_type;
    vector<double> prob;
    double avg_worm_len;
    uint worm_len_sample_size;
    bool mu_adjust;
    double mu_adjust_range;
    int mu_adjust_N;
    int mu_adjust_therm;
    int mu_adjust_sweep;
    int mu_index;

    // workspace variables
    double a[5][6];
    vector<vector<subseq_node> > subseq;
    vector<int> initial_Nd;
    vector<el_state> current_state;
    vector<int> current_occ;
    vector<vertex> vtx;
    vector<lock_flag> lock;
    vector<list_position> link;
    vector<list_position> first;
    vector<list_position> last;
    vector<int> sum_n, sum_s, sum_nn, sum_ss, sum_m;
    vector<double> S_rho_r, S_sigma_r;
    vector<double> chi_rho_r, chi_sigma_r;
    vector<double> mean_m;
    vector<double> cos_q_S;
    vector<double> sin_q_S;
    vector<double> cos_q_chi;
    vector<double> sin_q_chi;
    vector<double> mus;
    vector<double> N_mus;

    inline vertex diag_vertex_at_bond (vector<el_state>& state, unsigned short b) {
        vertex v;
        v.bottom_left = state[LEFT_SITE(b)];
        v.bottom_right = state[RIGHT_SITE(b)];
        v.top_right = v.bottom_right;
        v.top_left = v.bottom_left;
        return v;
    }

    void recalc_weights(vector<double> &weight, double mu);
    void recalc_directed_loop_probs();

public:    
    parser param;
#ifdef MCL_PT
    vector<measurements> measure;
#else
    measurements measure;
#endif

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
#ifdef MCL_PT
    void write_output(string,int);
#else
    void write_output(string);
#endif

    bool is_thermalized();
#ifdef MCL_PT
    bool request_global_update();
    void change_parameter(int);
    void change_to(int);
    double get_weight(int);
    int get_label();
#endif
    mc(string);
    ~mc();
};

#ifdef MCL_PT
typedef mc mc_pt;
#endif

#endif
