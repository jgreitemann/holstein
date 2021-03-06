// holstein
// Copyright (C) 2015  Jonas Greitemann <j.greitemann@lmu.de>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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
#define N_GROUP 10              // number of assignment groups

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"
#include <endian.hpp>

using namespace std;

typedef unsigned char byte;

enum el_state {
    empty,
    up,
    down,
    dublon
};

inline byte number_of_electrons (el_state s) {
    return (s & 1) + (s >> 1);
}

inline signed char local_magnetization (el_state s) {
    return (s & 1) - (s >> 1);
}

enum worm_type {
    up_worm,
    down_worm,
    dublon_worm
};

inline el_state flipped_state (el_state to_flip, worm_type what) {
    return static_cast<el_state>(to_flip ^ (what+1));
}

enum leg {
    bottom_left,
    bottom_right,
    top_right,
    top_left
};

inline leg straight(leg l) {
    return static_cast<leg>(l ^ 3);
}

enum operator_type {
    electron_diag,
    up_hopping,
    down_hopping,
    creator_left,
    creator_right,
    annihilator_left,
    annihilator_right,
    phonon_diag,
    disallowed
};

enum vtx_type {
    W_1m,
    W_10,
    W_1p,
    W_2m,
    W_2p,
    W_3,
    W_4,
    W_invalid
};

enum assignment_role {
    role_b1,
    role_b2,
    role_a,
    role_b,
    role_c,
    no_role
};

enum thermalization_stage {
    initial_stage,
    lower_stage,
    upper_stage,
    convergence_stage,
    tempering_stage,
    bdoub_therm_stage,
    final_stage,
    thermalized
};

struct thermalization_state {
    thermalization_stage stage;
    int sweeps;
    void set_stage (thermalization_stage s) {
        stage = s;
        sweeps = 0;
    }
    bool in_logging_stage () {
        switch (stage) {
            case lower_stage:
            case upper_stage:
            case convergence_stage:
                return true;
            default:
                return false;
        }
    } 
};

struct bond_operator {
    operator_type type  : 4;
    unsigned short bond : 12;
    inline bool operator== (const bond_operator& other) {
        return bond == other.bond && type == other.type;
    }
    inline bool operator!= (const bond_operator& other) {
        return !(*this == other);
    }
};

const bond_operator identity = {electron_diag, 0};

union vertex {
#ifdef BOOST_BIG_ENDIAN
    struct {
        el_state top_left     : 2;
        el_state top_right    : 2;
        el_state bottom_right : 2;
        el_state bottom_left  : 2;
    };
#else
#ifndef BOOST_LITTLE_ENDIAN
#warning Endianness unknown, defaulting for little endian!
#endif
    struct {
        el_state bottom_left  : 2;
        el_state bottom_right : 2;
        el_state top_right    : 2;
        el_state top_left     : 2;
    };
#endif
    byte int_repr;
    bool operator== (const vertex& other) {
        return int_repr == other.int_repr;
    }
    el_state get_state(leg l) {
        return static_cast<el_state>(int_repr>>(2*l) & 3);
    }
    vtx_type type() {
        int n_left = number_of_electrons(bottom_left);
        int n_right = number_of_electrons(bottom_right);
        switch (op_type()) {
            case electron_diag:
                switch (n_left + n_right) {
                    case 0: return W_1p;
                    case 1: return W_2p;
                    case 2:
                        if (n_left == 1 && n_right == 1)
                            return W_3;
                        else
                            return W_10;
                    case 3: return W_2m;
                    case 4: return W_1m;
                }
                break;
            case up_hopping:
            case down_hopping:
                return W_4;
        }
        return W_invalid;
    }
    operator_type op_type() {
        bool down_diag =    ((bottom_left  & down) == (top_left  & down))
                         && ((bottom_right & down) == (top_right & down));
        bool up_diag =    ((bottom_left  & up) == (top_left  & up))
                       && ((bottom_right & up) == (top_right & up));
        bool down_hop =    ((bottom_left  & down) == (top_right & down))
                        && ((bottom_right & down) == (top_left  & down));
        bool up_hop =    ((bottom_left  & up) == (top_right & up))
                      && ((bottom_right & up) == (top_left  & up));
        if (down_diag && up_diag) return electron_diag;
        if (down_diag && up_hop)  return up_hopping;
        if (up_diag && down_hop)  return down_hopping;
        return disallowed;
    }
};

union assignment {
#ifdef BOOST_BIG_ENDIAN
    struct {
        worm_type worm : 4;
        leg exit_leg   : 2;
        leg ent_leg    : 2;
        vertex vtx;
    };
#else
#ifndef BOOST_LITTLE_ENDIAN
#warning Endianness unknown, defaulting for little endian!
#endif
    struct {
        vertex vtx;
        leg ent_leg    : 2;
        leg exit_leg   : 2;
        worm_type worm : 4;
    };
#endif
    unsigned short int_repr;
    bool operator== (const assignment& other) {
        return int_repr == other.int_repr;
    }
    vertex flipped_vtx() {
        vertex res;
        res.int_repr = vtx.int_repr ^ ((worm+1) << (2*ent_leg))
                                    ^ ((worm+1) << (2*exit_leg));
        return res;
    }
    assignment flipped_assign() {
        assignment flipped;
        flipped.vtx = flipped_vtx();
        flipped.ent_leg = exit_leg;
        flipped.exit_leg = ent_leg;
        flipped.worm = worm;
        return flipped;
    }
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
#ifdef BOOST_BIG_ENDIAN
    struct {
        int index   : 30;
        leg vtx_leg : 2;
    };
#else
#ifndef BOOST_LITTLE_ENDIAN
#warning Endianness unknown, defaulting for little endian!
#endif
    struct {
        leg vtx_leg : 2;
        int index   : 30;
    };
#endif
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

struct fourier_mode {
    double q;
    double n_q_re;
    double n_q_im;
    double s_q_re;
    double s_q_im;
    double sum_n_q_re;
    double sum_n_q_im;
    double sum_s_q_re;
    double sum_s_q_im;
    vector<double> cos_q;
    vector<double> sin_q;
};

#ifndef NDEBUG
template class vector<int>;
template class vector<double>;
template class vector<el_state>;
template class vector<bond_operator>;
template class vector<lock_flag>;
template class vector<list_position>;
template class vector<operator_type>;
template class vector<fourier_mode>;
#endif

class mc {
private:
    uint M;
    uint Np;
    uint N_tau;
    uint L;
    double beta, init_beta, final_beta;
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
    double delta;
    double epsilon;
    double q_S;
    int matsubara;
    uint therm;
    uint tempering_therm;
    int bdoub_level;
    uint bdoub_therm;
    double tempering_exp;
    vector<el_state> state;
    vector<int> occ;
    vector<bond_operator> sm;
    uint sweep;
    uint init_n_max;
    uint loop_term;
    uint n;
    uint n_hop;
    double vtx_visited;
    uint N_loop;
    bool dublon_rejected;
    vector<double> weight;
#ifdef MCL_PT
    vector<double> other_weight;
#endif
    vector<double> prob;
    double avg_worm_len;
    uint worm_len_sample_size;
    bool mu_adjust;
    int mu_adjust_therm;
    int mu_adjust_sweep;
    double mu_adjust_range;
    double mu_adjust_tol;
    bool mus_file;
    thermalization_state therm_state;
    bool calc_dyn;
    uint bin_length;
    uint thermlog_interval;

    // vertex/assignment classification
    vector<assignment_role> role;
    vector<int> assign_group;
    vector<vtx_type> v_type;
    vector<operator_type> op_type;

    // workspace variables
    double a[no_role][N_GROUP];
    vector<vector<subseq_node> > subseq;
    vector<int> initial_Nd;
    vector<el_state> current_state;
    vector<int> current_occ;
    vector<vertex> vtx;
    vector<lock_flag> lock;
    vector<list_position> link;
    vector<list_position> first;
    vector<list_position> last;
    vector<double> S_rho_r, S_sigma_r;
    vector<double> C_rho_q, C_sigma_q;
    vector<double> tau;
    vector<int> n_p, s_p;
    vector<double> cos_q_S;
    vector<fourier_mode> ns_q;
    double lower_mu, upper_mu;
    double lower_N, upper_N;
    int N_mu;
    stringstream mu_data, bisection_protocol, thermlog;

    inline vertex diag_vertex_at_bond (vector<el_state>& state,
                                       unsigned short b) {
        vertex v;
        v.bottom_left = state[LEFT_SITE(b)];
        v.bottom_right = state[RIGHT_SITE(b)];
        v.top_right = v.bottom_right;
        v.top_left = v.bottom_left;
        return v;
    }

    void recalc_weights(vector<double> &weight, double mu, double &delta);
    void recalc_directed_loop_probs();
    void init_assignments();
    void init_vertices();

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
