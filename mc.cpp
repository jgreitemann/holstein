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

#include "mc.h"
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>

void bubble_sort_perm (int *a, int *p, uint k) {
    // initialize permutation
    for (uint i = 0; i < k; ++i) {
        p[i] = i;
    }

    // bubble sort
    for (uint i = 0; i < k-1; ++i) {
        for (uint j = 0; j < k-i-1; ++j) {
            if (a[j] > a[j+1]) {
                int tmp1 = a[j];
                int tmp2 = p[j];
                a[j] = a[j+1];
                p[j] = p[j+1];
                a[j+1] = tmp1;
                p[j+1] = tmp2;
            }
        }
    }
}

double default_mu(double g, double omega) {
    return g*g/omega;
}

mc :: mc (string dir) {
    vector<int> qvec;
    // initialize job parameters
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    bdoub_level = param.value_or_default<int>("BETA_DOUBLING_LEVEL", 0);
    final_beta = param.value_or_default<double>("BETA", L);
    init_beta = param.value_or_default<double>("INIT_BETA",
                                        final_beta * pow(0.5, bdoub_level));
    beta = init_beta;
    epsilon = param.value_or_default<double>("EPSILON", 0.01);
    N_el_up = param.value_or_default<int>("N_el_up", L/2);
    N_el_down = param.value_or_default<int>("N_el_down", N_el_up);
    enlargement_factor = param.value_or_default<double>("ENLARGEMENT_FACTOR",
                                                        1.3);
    U = param.value_or_default<double>("U", 0.);
    omega = param.value_or_default<double>("OMEGA", 1.);
#ifdef MCL_PT
    gvec = param.return_vector<double>("@G");
    pt_spacing = param.value_or_default<int>("PT_SPACING", 100);
#else
    g = param.value_or_default<double>("G", 0.);
#endif
    mu = param.value_or_default<double>("MU", default_mu(g, omega));
    q_S = param.value_or_default<double>("Q_S", 2*M_PI/L);
    calc_dyn = param.value_or_default<int>("DYNAMICAL_CORRELATIONS", 1);
    bin_length = param.value_or_default<int>("BINLENGTH", 1);
    if (calc_dyn)
        qvec = param.return_vector<int>("@Q");
    matsubara = param.value_or_default<int>("MATSUBARA", 0);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 50000);
    bdoub_therm = param.value_or_default<int>("BETA_DOUBLING_THERM", therm);
    tempering_therm = param.value_or_default<int>("TEMPERING_THERM", 0);
    tempering_exp = param.value_or_default<int>("TEMPERING_EXP", 1.);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    vtx_visited = param.value_or_default<double>("VTX_VISITED", 2.0);
    Np = param.value_or_default<int>("N_P", 20);
    mu_adjust = param.value_or_default<bool>("MU_ADJUST", 0);
    mu_adjust_range = param.value_or_default<double>("MU_ADJUST_RANGE", 0.5);
    mu_adjust_therm = param.value_or_default<int>("MU_ADJUST_THERM", 5000);
    mu_adjust_sweep = param.value_or_default<int>("MU_ADJUST_SWEEP", 10000);
    mu_adjust_tol = param.value_or_default<double>("MU_ADJUST_TOLERANCE", 0.01);
    bool pbc = param.value_or_default<int>("PERIODIC_BOUNDARY", 1);
    bool mus_file = param.value_or_default<int>("MUS_FILE", 1);
    thermlog_interval = param.value_or_default<int>("THERMLOG_INTERVAL", 0);
    assert(N_el_up <= L && N_el_down <= L);
    assert(N_el_up % 2 == pbc && N_el_down % 2 == pbc);
#ifdef MCL_PT
    // mu adjustment must not be used in conjunction with PT
    assert(!mu_adjust);
#endif

    // initialize vectors
    init_vertices();
    init_assignments();
    weight.resize(256);
    prob.resize(N_WORM<<12);
    subseq.resize(L, vector<subseq_node>());
    initial_Nd.resize(L);
    first.resize(L);
    last.resize(L);
    S_rho_r.resize(L);
    S_sigma_r.resize(L);
    n_p.resize(L);
    s_p.resize(L);
    C_rho_q.resize(qvec.size());
    C_sigma_q.resize(qvec.size());
    cos_q_S.resize(L);
#ifdef MCL_PT
    other_weight.resize(256);
    measure.resize(gvec.size());
#endif

    // calculate trigonometric factors for Fourier transforms
    for (uint s = 0; s < L; ++s)
        cos_q_S[s] = cos(q_S*s);

    if (calc_dyn) {
        vector<int>::iterator qit;
        fourier_mode mode;
        for (qit = qvec.begin(); qit != qvec.end(); ++qit) {
            mode.q = *qit * (2*M_PI/L);
            mode.cos_q.resize(L);
            mode.sin_q.resize(L);
            for (uint j = 0; j < L; ++j) {
                mode.cos_q[j] = cos(mode.q*j);
                mode.sin_q[j] = sin(mode.q*j);
            }
            ns_q.push_back(mode);
        }
    }
}

void mc :: recalc_weights(vector<double> &weight, double mu, double &delta) {
    fill(weight.begin(), weight.end(), 0.);

    // define weights
    double C = (U/4 > -abs(mu)) ? (U/4 + 2*abs(mu)) : (-U/4);
    double W[] = {
                     C - U/4 - 2*mu,
                     C - U/4,
                     C - U/4 + 2*mu,
                     C - mu,
                     C + mu,
                     C + U/4,
                     1.
                 };
    for (uint i = 0; i < 7; ++i) {
        assert(W[i] >= -1e-14
               || ((cerr << "W[" << i << "]=" << W[i] << endl) && false));
        if (abs(W[i]) < 1e-14)
            W[i] = 0.;
    }

    // calculate loop segment weights
    vtx_type groups[N_GROUP][3] = {
                                    {W_1p, W_2p, W_4},      // 0
                                    {W_1m, W_2m, W_4},      // 1
                                    {W_10, W_2p, W_4},      // 2
                                    {W_10, W_2m, W_4},      // 3
                                    {W_2p, W_3,  W_4},      // 4
                                    {W_2m, W_3,  W_4},      // 5
                                    {W_4,  W_4,  W_invalid},// 6
                                    {W_10, W_1p, W_invalid},// 7
                                    {W_1m, W_10, W_invalid},// 8
                                    {W_2m, W_2p, W_invalid} // 9
                                  };
    for (uint gr = 0; gr < 10; ++gr) {
        if (groups[gr][2] == W_invalid) {   // 2x2 group
            bool bigger = W[groups[gr][0]] < W[groups[gr][1]];
            a[role_a][gr] = W[groups[gr][!bigger]];
            a[role_b1][gr] = (!bigger) ? (W[groups[gr][bigger]]
                                          - W[groups[gr][!bigger]]) : 0.;
            a[role_b2][gr] = (bigger)  ? (W[groups[gr][bigger]]
                                          - W[groups[gr][!bigger]]) : 0.;
        } else {                            // 3x3 group
            a[role_a][gr] = 0.5 * (  W[groups[gr][0]]
                                   + W[groups[gr][1]]
                                   - W[groups[gr][2]]);
            a[role_b][gr] = 0.5 * (  W[groups[gr][0]]
                                   - W[groups[gr][1]]
                                   + W[groups[gr][2]]);
            a[role_c][gr] = 0.5 * ( -W[groups[gr][0]]
                                   + W[groups[gr][1]]
                                   + W[groups[gr][2]]);
            a[role_b1][gr] = (a[role_c][gr] < 0) ? -2.*a[role_c][gr] : 0.0;
            a[role_b2][gr] = (a[role_b][gr] < 0) ? -2.*a[role_b][gr] : 0.0;
            a[role_a][gr] += 0.5 * (-a[role_b1][gr] - a[role_b2][gr]);
            a[role_b][gr] += 0.5 * (-a[role_b1][gr] + a[role_b2][gr]);
            a[role_c][gr] += 0.5 * ( a[role_b1][gr] - a[role_b2][gr]);
        }
    }

    // determine delta
    delta = 0.0;
    for (uint gr = 0; gr < N_GROUP; ++gr) {
        if (W[groups[gr][2]] != W_invalid && a[role_a][gr] < -delta) {
            delta = -a[role_a][gr];
        }
    }
    for (uint gr = 0; gr < N_GROUP; ++gr) {
        if (W[groups[gr][2]] != W_invalid) {
            a[role_a][gr] += delta + epsilon;
        }
    }
    for (uint i = 0; i < 6; ++i) {
        W[i] += delta + epsilon;
    }

    // assign vertex weights
    for (int vtx_i = 0; vtx_i < 256; ++vtx_i)
        if (v_type[vtx_i] != W_invalid)
            weight[vtx_i] = W[v_type[vtx_i]];
}

void mc :: recalc_directed_loop_probs() {
    recalc_weights(weight, mu, delta);
    fill(prob.begin(), prob.end(), 0.);

    // calculate transition probabilities
    assignment assign;
    for (int assign_i = 0; assign_i < (N_WORM << 12); ++assign_i) {
        if (role[assign_i] == no_role)
            continue;
        assign.int_repr = assign_i;
        prob[assign_i] = a[role[assign_i]][assign_group[assign_i]]
                         / weight[assign.vtx.int_repr];
    }

    // cumulate transition probabilities
    assignment assign2;
    for (unsigned short vertex_i = 0; vertex_i < 256; ++vertex_i) {
        assign.vtx.int_repr = vertex_i;
        for (byte ent_leg_i = bottom_left; ent_leg_i <= top_left; ++ent_leg_i) {
            assign.ent_leg = static_cast<leg>(ent_leg_i);
            for (byte worm = 0; worm < N_WORM; ++worm) {
                assign.worm = static_cast<worm_type>(worm);
                assign.exit_leg = bottom_left;
                assign2 = assign;
                assign2.exit_leg = bottom_right;
                prob[assign2.int_repr] += prob[assign.int_repr];
                assign.exit_leg = top_right;
                prob[assign.int_repr] += prob[assign2.int_repr];
                assign2.exit_leg = top_left;
                prob[assign2.int_repr] += prob[assign.int_repr];
            }
        }
    }
}

mc :: ~mc() {
    random_clear();
    state.clear();
    occ.clear();
    sm.clear();
    weight.clear();
    v_type.clear();
    op_type.clear();
    role.clear();
    assign_group.clear();
    prob.clear();
    subseq.clear();
    initial_Nd.clear();
    current_state.clear();
    current_occ.clear();
    vtx.clear();
    lock.clear();
    link.clear();
    first.clear();
    last.clear();
    S_rho_r.clear();
    S_sigma_r.clear();
    C_rho_q.clear();
    C_sigma_q.clear();
    n_p.clear();
    s_p.clear();
    ns_q.clear();
    cos_q_S.clear();
}

void mc :: do_update() {
    // log particle numbers during thermalization if so desired
    if (thermlog_interval && therm_state.stage != thermalized
            && therm_state.sweeps % thermlog_interval == 0) {
        int N = 0;
        for (uint s = 0; s < L; ++s) {
            N += number_of_electrons(state[s]);
        }
        thermlog << therm_state.stage << " "
                 << therm_state.sweeps << " "
                 << N << " "
                 << beta << endl;
    }

    // change thermalization stage as necessary
    switch (therm_state.stage) {
        case initial_stage:
            beta = init_beta;
            if (therm_state.sweeps >= therm) {
                N_loop = (uint)(vtx_visited / avg_worm_len * M);
                therm_state.set_stage(mu_adjust ? lower_stage
                                                : bdoub_therm_stage);
                N_mu = 0;
            }
            break;
        case lower_stage:
            beta = init_beta;
            if (therm_state.sweeps >= mu_adjust_therm+mu_adjust_sweep) {
                lower_N = 1. * N_mu / mu_adjust_sweep;
                mu_data << mu << " " << lower_N << endl;
                if (lower_N < N_el_up+N_el_down) {
                    lower_mu = lower_mu - 0.5 * (upper_mu - lower_mu);
                    therm_state.set_stage(lower_stage);
                    mu = lower_mu;
                } else {
                    therm_state.set_stage(upper_stage);
                    mu = upper_mu;
                }
                recalc_directed_loop_probs();
                N_mu = 0;
            }
            break;
        case upper_stage:
            beta = init_beta;
            if (therm_state.sweeps >= mu_adjust_therm+mu_adjust_sweep) {
                upper_N = 1. * N_mu / mu_adjust_sweep;
                mu_data << mu << " " << upper_N << endl;
                if (upper_N > N_el_up+N_el_down) {
                    upper_mu = upper_mu + 0.5 * (upper_mu - lower_mu);
                    therm_state.set_stage(upper_stage);
                    mu = upper_mu;
                } else {
                    therm_state.set_stage(convergence_stage);
                    mu = 0.5 * (lower_mu + upper_mu);
                }
                recalc_directed_loop_probs();
                N_mu = 0;
            }
            break;
        case convergence_stage:
            beta = init_beta;
            if (therm_state.sweeps >= mu_adjust_therm+mu_adjust_sweep) {
                double center_N = 1. * N_mu / mu_adjust_sweep;
                mu_data << mu << " " << center_N << endl;
                if (center_N < N_el_up+N_el_down) {
                    upper_mu = 0.5 * (lower_mu + upper_mu);
                    upper_N = center_N;
                } else {
                    lower_mu = 0.5 * (lower_mu + upper_mu);
                    lower_N = center_N;
                }
                bisection_protocol << lower_mu << " " << upper_mu << endl;
                if (upper_mu - lower_mu > mu_adjust_tol) {
                    therm_state.set_stage(convergence_stage);
                } else {
                    therm_state.set_stage(tempering_stage);
                }
                mu = 0.5 * (lower_mu + upper_mu);
                recalc_directed_loop_probs();
                N_mu = 0;
            }
            break;
        case tempering_stage:
            beta = init_beta + (final_beta*pow(0.5, bdoub_level) - init_beta)
                   * pow(1.*therm_state.sweeps/tempering_therm, tempering_exp);
            if (therm_state.sweeps >= tempering_therm) {
                therm_state.set_stage(final_stage);
            }
            break;
        case final_stage:
            beta = final_beta * pow(0.5, bdoub_level);
            if (therm_state.sweeps >= therm) {
                therm_state.set_stage(bdoub_therm_stage);
                if (!mus_file)
                    break;
                stringstream fname;
                fname << "../mus/" << setprecision(4) << U << "_" << g << "_"
                      << omega << ".mu";
                ofstream fstr(fname.str().c_str());
                if (fstr.is_open()) {
                    fstr << mu << " " << U << " " << g << " " << omega << endl
                         << "# (above) mu_best U g omega" << endl
                         << endl << endl
                         << "# (below) mu N_mu" << endl
                         << mu_data.str()
                         << endl << endl
                         << "# (below) lower_mu upper_mu" << endl
                         << bisection_protocol.str();
                    fstr.close();
                }
            }
            break;
        case bdoub_therm_stage:
            beta = final_beta * pow(0.5, bdoub_level);
            if (therm_state.sweeps == 1) {
                if (bdoub_therm <= 1 || bdoub_level <= 0) {
                    therm_state.set_stage(thermalized);
                    break;
                }
                // append operator string to itself
                sm.resize(2 * M);
                for (uint i = 0; i < M; ++i) {
                    sm[i+M] = sm[i];
                }
                M *= 2;
                n *= 2;
                n_hop *= 2;
                --bdoub_level;
                beta = final_beta * pow(0.5, bdoub_level);
            }
            if (therm_state.sweeps == bdoub_therm) {
                therm_state.set_stage(bdoub_level > 0 ? bdoub_therm_stage
                                                      : thermalized);
            }
            break;
        case thermalized:
            beta = final_beta;
            break;
    }

    // diagonal update & subsequence construction
    for_each(subseq.begin(), subseq.end(),
             mem_fun_ref(&vector<subseq_node>::clear));
    fill(initial_Nd.begin(), initial_Nd.end(), 0);
    current_state = state;
    current_occ = occ;
    for (uint i = 0; i < M; ++i) {
        // electronic diagonal update
        if (sm[i] == identity) {
            bond_operator b;
            b.type = electron_diag;
            b.bond = random0N(NB)+1;
            vertex vtx = diag_vertex_at_bond(current_state, b.bond);
            if (random01() < NB*beta*weight[vtx.int_repr]/(M-n)) {
                sm[i] = b;
                n++;
            }
        } else if (sm[i].type == electron_diag) {
            vertex vtx = diag_vertex_at_bond(current_state, sm[i].bond);
            if (random01() < (M-n+1)/(NB*beta*weight[vtx.int_repr])) {
                sm[i] = identity;
                n--;
            }
        }

        // phononic diagonal update
        if (Np && sm[i] == identity) {
            bond_operator b;
            b.type = phonon_diag;
            b.bond = random0N(NB)+1;
            int cocc = current_occ[LEFT_SITE(b.bond)];
            if (random01() < NB*beta*omega*(Np-cocc)/(M-n)) {
                sm[i] = b;
                n++;
            }
        } else if (sm[i].type == phonon_diag) {
            int cocc = current_occ[LEFT_SITE(sm[i].bond)];
            if (random01() < (M-n+1)/(NB*beta*omega*(Np-cocc))) {
                sm[i] = identity;
                n--;
            }
        }

        bond_operator b = sm[i];
        if (b != identity) {
            // append to the phonon subsequences
            switch (b.type) {
                case phonon_diag:
                    if (subseq[LEFT_SITE(b.bond)].empty()) {
                        ++initial_Nd[LEFT_SITE(b.bond)];
                    } else {
                        ++(subseq[LEFT_SITE(b.bond)].back().Nd);
                    }
                    break;
                case up_hopping:
                case down_hopping:
                    break;
                default:
                    if (!Np)
                        break;
                    unsigned short site;
                    switch (b.type) {
                        case electron_diag:
                            site = random0N(2) ? LEFT_SITE(b.bond)
                                               : RIGHT_SITE(b.bond);
                            break;
                        case creator_left:
                        case annihilator_left:
                            site = LEFT_SITE(b.bond);
                            break;
                        case creator_right:
                        case annihilator_right:
                            site = RIGHT_SITE(b.bond);
                            break;
                    }
                    vertex vtx = diag_vertex_at_bond(current_state, b.bond);
                    byte n_el = number_of_electrons(current_state[site]);
                    subseq_node newNode;
                    newNode.i = i;
                    newNode.Nd = 0;
                    newNode.m = current_occ[site];
                    newNode.r = weight[vtx.int_repr]/n_el;
                    subseq[site].push_back(newNode);
                    break;
            }

            // propagation of state
            switch (b.type) {
                case up_hopping:
                    current_state[LEFT_SITE(b.bond)] =
                        flipped_state(current_state[LEFT_SITE(b.bond)],
                                      up_worm);
                    current_state[RIGHT_SITE(b.bond)] =
                        flipped_state(current_state[RIGHT_SITE(b.bond)],
                                      up_worm);
                    break;
                case down_hopping:
                    current_state[LEFT_SITE(b.bond)] =
                        flipped_state(current_state[LEFT_SITE(b.bond)],
                                      down_worm);
                    current_state[RIGHT_SITE(b.bond)] =
                        flipped_state(current_state[RIGHT_SITE(b.bond)],
                                      down_worm);
                    break;
                case creator_left:
                    current_occ[LEFT_SITE(b.bond)]++;
                    break;
                case creator_right:
                    current_occ[RIGHT_SITE(b.bond)]++;
                    break;
                case annihilator_left:
                    current_occ[LEFT_SITE(b.bond)]--;
                    break;
                case annihilator_right:
                    current_occ[RIGHT_SITE(b.bond)]--;
                    break;
            }
        }
    }

    // transfer Nd to the back of the list & push front to the back
    for (uint s = 0; s < L; ++s) {
        if (subseq[s].size() >= 2) {
            subseq[s].back().Nd += initial_Nd[s];
        } else {
            subseq[s].clear();
        }
    }

    // subsequence phonon update
    for (uint s = 0; Np && s < L; ++s) {
        vector<subseq_node>::iterator i1, i2;
        if (subseq[s].size() < 2)
            continue;
        for (uint i = 0; i < M/L; ++i) {
            i2 = subseq[s].begin() + random0N(subseq[s].size());
            i1 = i2++;
            if (i2 == subseq[s].end()) {
                i2 = subseq[s].begin();
            }
            if (   sm[i1->i].type == electron_diag
                && sm[i2->i].type == electron_diag) {
                if (random0N(2)) { // (H_1, H_1) -> (H_5, H_4)
                    double prob = g*g * i1->m / (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m+1)/(Np-i1->m), i1->Nd);
                    if (random01() < prob) {
                        sm[i1->i].type = sm[i1->i].bond == LEFT_BOND(s)
                                         ? annihilator_right
                                         : annihilator_left;
                        sm[i2->i].type = sm[i2->i].bond == LEFT_BOND(s)
                                         ? creator_right
                                         : creator_left;
                        --(i2->m);
                    }
                } else { // (H_1, H_1) -> (H_4, H_5)
                    double prob = g*g * (i1->m+1) / (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m-1)/(Np-i1->m), i1->Nd);
                    if (i1->m < Np && random01() < prob) {
                        sm[i1->i].type = sm[i1->i].bond == LEFT_BOND(s)
                                         ? creator_right
                                         : creator_left;
                        sm[i2->i].type = sm[i2->i].bond == LEFT_BOND(s)
                                         ? annihilator_right
                                         : annihilator_left;
                        ++(i2->m);
                    }
                }
            } else if (   (   sm[i1->i].type == creator_left
                           || sm[i1->i].type == creator_right)
                       && (   sm[i2->i].type == annihilator_left
                           || sm[i2->i].type == annihilator_right)) {
                if (random0N(2)) { // (H_4, H_5) -> (H_5, H_4)
                    double prob = 1. * (i1->m) / (i1->m+1)
                                  * pow(1.*(Np-i1->m+1)/(Np-i1->m-1), i1->Nd);
                    if (random01() < prob) {
                        sm[i1->i].type = sm[i1->i].type == creator_right
                                         ? annihilator_right
                                         : annihilator_left;
                        sm[i2->i].type = sm[i2->i].type == annihilator_right
                                         ? creator_right
                                         : creator_left;
                        i2->m -= 2;
                    }
                } else { // (H_4, H_5) -> (H_1, H_1)
                    double prob = 1./g/g / (i1->m+1) * (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m)/(Np-i1->m-1), i1->Nd);
                    if (random01() < prob) {
                        sm[i1->i].type = electron_diag;
                        sm[i2->i].type = electron_diag;
                        --(i2->m);
                    }
                }
            } else if (   (   sm[i1->i].type == annihilator_left
                           || sm[i1->i].type == annihilator_right)
                       && (   sm[i2->i].type == creator_left
                           || sm[i2->i].type == creator_right)) {
                if (random0N(2)) { // (H_5, H_4) -> (H_4, H_5)
                    double prob = 1. / (i1->m) * (i1->m+1)
                                  * pow(1.*(Np-i1->m-1)/(Np-i1->m+1), i1->Nd);
                    if (i1->m < Np && random01() < prob) {
                        sm[i1->i].type = sm[i1->i].type == annihilator_right
                                         ? creator_right
                                         : creator_left;
                        sm[i2->i].type = sm[i2->i].type == creator_right
                                         ? annihilator_right
                                         : annihilator_left;
                        i2->m += 2;
                    }
                } else { // (H_5, H_4) -> (H_1, H_1)
                    double prob = 1./g/g / i1->m * (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m)/(Np-i1->m+1), i1->Nd);
                    if (random01() < prob) {
                        sm[i1->i].type = electron_diag;
                        sm[i2->i].type = electron_diag;
                        ++(i2->m);
                    }
                }
            }
        }
    }

    // mapping back to phonon state
    for (uint s = 0; s < L; ++s) {
        if (!subseq[s].empty()) {
            occ[s] = subseq[s].front().m;
        }
    }

    // adjust M during thermalization
    if (!is_thermalized() && enlargement_factor*n > M) {
        M = enlargement_factor*n;
        sm.resize(M, identity);
    }
    assert(n <= M); // You might need to increase "a" or the
                    // thermalization time if this assertion fails.

    // directed loops electron update
    if (n > 0) {
        // linked list construction
        vtx.resize(n);
        lock.resize(n);
        link.resize(4*n);
        fill(first.begin(), first.end(), invalid_pos);
        fill(last.begin(), last.end(), invalid_pos);
        current_state = state;
        uint p = 0;
        for (uint i = 0; p < n; ++i) {
            assert(i < M);
            if (sm[i] == identity)
                continue;
            // establish links
            if (first[LEFT_SITE(sm[i].bond)] == invalid_pos) {
                first[LEFT_SITE(sm[i].bond)] = list_position(bottom_left, p);
                last[LEFT_SITE(sm[i].bond)] = list_position(top_left, p);
            } else {
                link[list_position(bottom_left, p).int_repr] =
                    last[LEFT_SITE(sm[i].bond)];
                link[last[LEFT_SITE(sm[i].bond)].int_repr] =
                    list_position(bottom_left, p);
                last[LEFT_SITE(sm[i].bond)] = list_position(top_left, p);
            }
            if (first[RIGHT_SITE(sm[i].bond)] == invalid_pos) {
                first[RIGHT_SITE(sm[i].bond)] = list_position(bottom_right, p);
                last[RIGHT_SITE(sm[i].bond)] = list_position(top_right, p);
            } else {
                link[list_position(bottom_right, p).int_repr] =
                    last[RIGHT_SITE(sm[i].bond)];
                link[last[RIGHT_SITE(sm[i].bond)].int_repr] =
                    list_position(bottom_right, p);
                last[RIGHT_SITE(sm[i].bond)] = list_position(top_right, p);
            }

            // determine vertex type
            vtx[p].bottom_left = current_state[LEFT_SITE(sm[i].bond)];
            vtx[p].bottom_right = current_state[RIGHT_SITE(sm[i].bond)];
            switch (sm[i].type) {
                case up_hopping:
                    current_state[LEFT_SITE(sm[i].bond)] =
                        flipped_state(current_state[LEFT_SITE(sm[i].bond)],
                                      up_worm);
                    current_state[RIGHT_SITE(sm[i].bond)] =
                        flipped_state(current_state[RIGHT_SITE(sm[i].bond)],
                                      up_worm);
                    lock[p] = unlocked;
                    break;
                case down_hopping:
                    current_state[LEFT_SITE(sm[i].bond)] =
                        flipped_state(current_state[LEFT_SITE(sm[i].bond)],
                                      down_worm);
                    current_state[RIGHT_SITE(sm[i].bond)] =
                        flipped_state(current_state[RIGHT_SITE(sm[i].bond)],
                                      down_worm);
                    lock[p] = unlocked;
                    break;
                case creator_left:
                case annihilator_left:
                    lock[p] = left_lock; // require at least one electron on
                                         // left site
                    break;
                case creator_right:
                case annihilator_right:
                    lock[p] = right_lock; // require at least one electron on
                                          // right site
                    break;
                case phonon_diag:
                    lock[p] = total_lock; // require at least one electron on
                                          // both sites
                    break;
                default:
                    lock[p] = unlocked;
                    break;
            }
            vtx[p].top_right = current_state[RIGHT_SITE(sm[i].bond)];
            vtx[p].top_left = current_state[LEFT_SITE(sm[i].bond)];

            ++p;
        }
        for (uint s = 0; s < L; ++s) {
            if (last[s] != invalid_pos) {
                link[first[s].int_repr] = last[s];
                link[last[s].int_repr] = first[s];
            }
        }

        // directed loop construction
        int R;
        list_position j, j0;
        assignment assign;
        worm_type worm;
        double r;
        for (uint i = 0; i < N_loop; ++i) {
            R = random0N(N_WORM*4*n);
            worm = static_cast<worm_type>(R / (4*n));
            j0.int_repr = R % (4*n);
            j = j0;

            // check if dublon worm can start from here
            el_state ent_state = vtx[j0.index].get_state(j0.vtx_leg);
            if (worm == dublon_worm) {
                dublon_rejected = ent_state != empty && ent_state != dublon;
                if (dublon_rejected) {
                    // IMPORTANT: this has to be counted as a loop
                    continue; // not an empty or fully occupied site
                }
            }

            int k;
            for (k = 0; ; ++k) {
                if (k == loop_term*M) {
                    do_update();
                    return;
                }
                assign.vtx = vtx[j.index];
                assign.ent_leg = j.vtx_leg;
                assign.worm = worm;
                switch (lock[j.index]) {
                    case total_lock: // continue straight
                        assign.exit_leg = straight(j.vtx_leg);
                        break;
                    case left_lock:
                    case right_lock:
                        if (   (   lock[j.index] == left_lock
                                && (   j.vtx_leg == bottom_left
                                    || j.vtx_leg == top_left))
                            || (   lock[j.index] == right_lock
                                && (   j.vtx_leg == bottom_right
                                    || j.vtx_leg == top_right))) {
                            if (flipped_state(vtx[j.index].get_state(j.vtx_leg),
                                              worm) == empty) {
                                assign.exit_leg = j.vtx_leg; // bounce
                            } else {
                                if (vtx[j.index].get_state(j.vtx_leg)==dublon) {
                                    if (random0N(2)) {
                                         // continue straight
                                         assign.exit_leg = straight(j.vtx_leg);
                                    } else {
                                        assign.exit_leg = j.vtx_leg; // bounce
                                    }
                                } else {
                                    // continue straight
                                    assign.exit_leg = straight(j.vtx_leg);
                                }
                            }
                        } else {
                            // continue straight
                            assign.exit_leg = straight(j.vtx_leg);
                        }
                        break;
                    case unlocked:
                        r = random01();
                        byte exit_leg_i;
                        for (exit_leg_i = bottom_left;
                                exit_leg_i <= top_left;
                                ++exit_leg_i) {
                            assign.exit_leg = static_cast<leg>(exit_leg_i);
                            if (r < prob[assign.int_repr])
                                break;
                        }
                        assert(exit_leg_i < 4); // assert that break was called
                        break;
                }
                // do not count bounces into the worm length
                if (assign.ent_leg == assign.exit_leg) {
                    --k;
                }
                // flip the vertex:
                vtx[j.index] = assign.flipped_vtx();
                j.vtx_leg = assign.exit_leg;// exit leg position in linked list
                if (j == j0)    // loop closed (c.f. [SS02], Fig. 4b)
                    break;
                j = link[j.int_repr];
                if (j == j0)    // loop closed (c.f. [SS02], Fig. 4a)
                    break;
            }
            // logging worm length
            if (therm_state.stage == initial_stage
                    && therm_state.sweeps < therm
                    && therm_state.sweeps >= therm/2) {
                avg_worm_len *= 1.*worm_len_sample_size
                                / (worm_len_sample_size+1);
                avg_worm_len += 1.*k / (++worm_len_sample_size);
            }
        }

        // mapping back to operator sequence
        p = 0;
        operator_type new_type;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == identity)
                continue;
            switch (sm[i].type) {
                case electron_diag:
                case up_hopping:
                case down_hopping:
                    new_type = op_type[vtx[p].int_repr];
#ifdef MEASURE_KIN_ENERGY
                    if (sm[i].type == electron_diag && new_type != electron_diag) {
                        n_hop += 1;
                    } else if (sm[i].type != electron_diag && new_type == electron_diag) {
                        n_hop -= 1;
                    }
#endif
                    sm[i].type = new_type;
                    break;
            }
            ++p;
        }

        // updating the state, flipping states randomly on sites not
        // affected by the bond operators
        for (uint s = 0; s < L; ++s) {
            if (first[s] == invalid_pos) {
                state[s] = static_cast<el_state>(random0N(4));
            } else {
                state[s] = vtx[first[s].index].get_state(first[s].vtx_leg);
            }
        }
    } else { // empty operator sequence
        for (uint s = 0; s < L; ++s) {
            state[s] = static_cast<el_state>(random0N(4));
        }
    }

    // log the number of electrons if necessary
    if (mu_adjust && therm_state.in_logging_stage()
                  && therm_state.sweeps >= mu_adjust_therm) {
        for (uint s = 0; s < L; ++s) {
            N_mu += number_of_electrons(state[s]);
        }
    }
    ++(therm_state.sweeps);
}

void mc :: do_measurement() {
    uint N_up = 0, N_down = 0;
    for (uint s = 0; s < L; ++s) {
        N_up += state[s] == up || state[s] == dublon;
        N_down += state[s] == down || state[s] == dublon;
    }
#ifdef MCL_PT
    measure[myrep].add("N_up", N_up);
    measure[myrep].add("N_down", N_down);
    measure[myrep].add("dublon_rejection_rate", dublon_rejected);
#else
    measure.add("N_up", N_up);
    measure.add("N_down", N_down);
    measure.add("dublon_rejection_rate", dublon_rejected);
#endif

    // skip measurement if particle numbers are not right
    if (N_up != N_el_up || N_down != N_el_down)
        return;

    double C = (U/4 > -abs(mu)) ? (U/4 + 2*abs(mu)) : (-U/4);
    double energy = -1./beta*n + NB*(C+delta+epsilon) + L*omega*Np
                    + 2*mu*(L-N_el_up-N_el_down);

    // add data to measurement
#ifdef MCL_PT
    measure[myrep].add("Energy", energy);
    #ifdef MEASURE_KIN_ENERGY
    measure[myrep].add("kinetic_Energy", -1./beta*n_hop);
    #endif
#else
    measure.add("Energy", energy);
    #ifdef MEASURE_KIN_ENERGY
    measure.add("kinetic_Energy", -1./beta*n_hop);
    #endif
#endif

    // calculate equal-time real space correlation functions (at p = 0)
    int n_staggered = 0, s_staggered = 0;
    int m = 0;
    for (uint j = 0; j < L; ++j) {
        n_p[j] = number_of_electrons(state[j]);
        s_p[j] = local_magnetization(state[j]);
        if (j & 1) {
            n_staggered -= n_p[j];
            s_staggered -= s_p[j];
        } else {
            n_staggered += n_p[j];
            s_staggered += s_p[j];
        }
        m += occ[j];
    }
#ifdef MCL_PT
    measure[myrep].add("ph_density", 1.*m/L);
#else
    measure.add("ph_density", 1.*m/L);
#endif
    for (uint r = 0; r < L; ++r) {
        double sum_nn = 0, sum_ss = 0;
        for (uint j = 0; j < L; ++j) {
            sum_nn += n_p[(j+r)%L] * n_p[j];
            sum_ss += s_p[(j+r)%L] * s_p[j];
        }
        S_rho_r[r] = 1./L*sum_nn;
        S_sigma_r[r] = 1./L*sum_ss;
    }

    // Fourier transform to obtain structure factors
    double S_rho_q = 0.0;
    double S_sigma_q = 0.0;
    for (uint s = 0; s < L; ++s) {
        S_rho_q += cos_q_S[s] * S_rho_r[s];
        S_sigma_q += cos_q_S[s] * S_sigma_r[s];
    }
#ifdef MCL_PT
    measure[myrep].add("S_rho_q", S_rho_q);
    measure[myrep].add("S_sigma_q", S_sigma_q);
#else
    measure.add("S_rho_q", S_rho_q);
    measure.add("S_sigma_q", S_sigma_q);
#endif

    double omega_mats;
    if (calc_dyn) {
        // generate imaginary times and sort
        omega_mats = 2*M_PI / beta * matsubara;
        tau.resize(n+2);
        tau[0] = 0;
        for (uint p = 1; p <= n; ++p) {
            tau[p] = random01() * beta;
        }
        tau[n+1] = beta;
        sort(tau.begin(), tau.end());

        // initialize dynamical correlation functions
        vector<fourier_mode>::iterator qit;
        for (qit = ns_q.begin(); qit != ns_q.end(); ++qit) {
            qit->n_q_re = 0.;
            qit->n_q_im = 0.;
            qit->s_q_re = 0.;
            qit->s_q_im = 0.;
            qit->sum_n_q_re = 0.;
            qit->sum_n_q_im = 0.;
            qit->sum_s_q_re = 0.;
            qit->sum_s_q_im = 0.;
            for (uint j = 0; j < L; ++j) {
                double n_j = number_of_electrons(state[j]);
                double s_j = local_magnetization(state[j]);
                qit->n_q_re += n_j * qit->cos_q[j];
                qit->n_q_im += n_j * qit->sin_q[j];
                qit->s_q_re += s_j * qit->cos_q[j];
                qit->s_q_im += s_j * qit->sin_q[j];
            }
        }
    }

    // initialize staggered susceptibilities 
    current_state = state;
    uint p = 0;
    int W_up = 0, W_down = 0;
    int sum_n_staggered = 0, sum_n_staggered_sq = n_staggered * n_staggered;
    int sum_s_staggered = 0, sum_s_staggered_sq = s_staggered * s_staggered;
    
    // traverse the operator string
    for (uint i = 0; p < n; ++i) {
        // fast-forward to next proper operator
        if (sm[i] == identity)
            continue;
        
        // sample staggered density & spin & their squares
        sum_n_staggered += n_staggered;
        sum_n_staggered_sq += n_staggered * n_staggered;
        sum_s_staggered += s_staggered;
        sum_s_staggered_sq += s_staggered * s_staggered;
        
        if (calc_dyn) {
            // sample dynamical correlations
            double mats_cos = matsubara
                ? (cos(omega_mats * tau[p+1]) - cos(omega_mats * tau[p]))
                : (tau[p+1] - tau[p]);
            double mats_sin = matsubara
                ? (-sin(omega_mats * tau[p+1]) + sin(omega_mats * tau[p]))
                : 0;
            vector<fourier_mode>::iterator qit;
            for (qit = ns_q.begin(); qit != ns_q.end(); ++qit) {
                qit->sum_n_q_re += qit->n_q_re*mats_cos - qit->n_q_im*mats_sin;
                qit->sum_n_q_im += qit->n_q_re*mats_sin + qit->n_q_im*mats_cos;
                qit->sum_s_q_re += qit->s_q_re*mats_cos - qit->s_q_im*mats_sin;
                qit->sum_s_q_im += qit->s_q_re*mats_sin + qit->s_q_im*mats_cos;
            }
        }

        // imaginary time propagation of state
        bond_operator b = sm[i];
        int from, to;
        switch (b.type) {
            case up_hopping:
                current_state[LEFT_SITE(b.bond)] =
                    flipped_state(current_state[LEFT_SITE(b.bond)], up_worm);
                current_state[RIGHT_SITE(b.bond)] =
                    flipped_state(current_state[RIGHT_SITE(b.bond)], up_worm);
                if (current_state[LEFT_SITE(b.bond)] & up) {
                    if (b.bond & 1) {
                        n_staggered += 2;
                        s_staggered += 2;
                    } else {
                        n_staggered -= 2;
                        s_staggered -= 2;
                    }
                    if (b.bond == L) {
                        ++W_up;
                    }
                } else {
                    if (b.bond & 1) {
                        n_staggered -= 2;
                        s_staggered -= 2;
                    } else {
                        n_staggered += 2;
                        s_staggered += 2;
                    }
                    if (b.bond == L) {
                        --W_up;
                    }
                }
                if (calc_dyn) {
                    to = (current_state[LEFT_SITE(b.bond)] & up)
                                ? LEFT_SITE(b.bond) : RIGHT_SITE(b.bond);
                    from = (current_state[LEFT_SITE(b.bond)] & up)
                                ? RIGHT_SITE(b.bond) : LEFT_SITE(b.bond);
                    vector<fourier_mode>::iterator qit;
                    for (qit = ns_q.begin(); qit != ns_q.end(); ++qit) {
                        qit->n_q_re += qit->cos_q[to] - qit->cos_q[from];
                        qit->n_q_im += qit->sin_q[to] - qit->sin_q[from];
                        qit->s_q_re += qit->cos_q[to] - qit->cos_q[from];
                        qit->s_q_im += qit->sin_q[to] - qit->sin_q[from];
                    }
                }
                break;
            case down_hopping:
                current_state[LEFT_SITE(b.bond)] =
                    flipped_state(current_state[LEFT_SITE(b.bond)], down_worm);
                current_state[RIGHT_SITE(b.bond)] =
                    flipped_state(current_state[RIGHT_SITE(b.bond)], down_worm);
                if (current_state[LEFT_SITE(b.bond)] & down) {
                    if (b.bond & 1) {
                        n_staggered += 2;
                        s_staggered -= 2;
                    } else {
                        n_staggered -= 2;
                        s_staggered += 2;
                    }
                    if (b.bond == L) {
                        ++W_down;
                    }
                } else {
                    if (b.bond & 1) {
                        n_staggered -= 2;
                        s_staggered += 2;
                    } else {
                        n_staggered += 2;
                        s_staggered -= 2;
                    }
                    if (b.bond == L) {
                        --W_down;
                    }
                }
                if (calc_dyn) {
                    to = (current_state[LEFT_SITE(b.bond)] & down)
                                ? LEFT_SITE(b.bond) : RIGHT_SITE(b.bond);
                    from = (current_state[LEFT_SITE(b.bond)] & down)
                                ? RIGHT_SITE(b.bond) : LEFT_SITE(b.bond);
                    vector<fourier_mode>::iterator qit;
                    for (qit = ns_q.begin(); qit != ns_q.end(); ++qit) {
                        qit->n_q_re += qit->cos_q[to] - qit->cos_q[from];
                        qit->n_q_im += qit->sin_q[to] - qit->sin_q[from];
                        qit->s_q_re -= qit->cos_q[to] - qit->cos_q[from];
                        qit->s_q_im -= qit->sin_q[to] - qit->sin_q[from];
                    }
                }
                break;
        }
        ++p;
    }

    // measure winding number fluctuations
#ifdef MCL_PT
    measure[myrep].add("W", W_up + W_down);
    measure[myrep].add("W_sq", W_up*W_up + W_down*W_down);
#else
    measure.add("W", W_up + W_down);
    measure.add("W_sq", W_up*W_up + W_down*W_down);
#endif
    
    // cf. [DT01]
    double chi_rho_pi = beta/L * (1./n/(n+1) * sum_n_staggered * sum_n_staggered
                                  + 1./(n+1)/(n+1) * sum_n_staggered_sq);
    double chi_sigma_pi = beta/L * (1./n/(n+1) * sum_s_staggered * sum_s_staggered
                                    + 1./(n+1)/(n+1) * sum_s_staggered_sq);
#ifdef MCL_PT
    measure[myrep].add("chi_rho_pi", chi_rho_pi);
    measure[myrep].add("chi_sigma_pi", chi_sigma_pi);

    // measure equal-time real space correlations
    measure[myrep].add("S_rho_r", S_rho_r);
    measure[myrep].add("S_sigma_r", S_sigma_r);
#else
    measure.add("chi_rho_pi", chi_rho_pi);
    measure.add("chi_sigma_pi", chi_sigma_pi);

    // measure equal-time real space correlations
    measure.add("S_rho_r", S_rho_r);
    measure.add("S_sigma_r", S_sigma_r);
#endif

    // measure dynamical correlations
    if (calc_dyn) {
        // final p = n term
        double mats_cos = matsubara
            ? (cos(omega_mats * tau[p+1]) - cos(omega_mats * tau[p]))
            : (tau[p+1] - tau[p]);
        double mats_sin = matsubara
            ? (-sin(omega_mats * tau[p+1]) + sin(omega_mats * tau[p]))
            : 0;
        vector<fourier_mode>::iterator qit;
        for (qit = ns_q.begin(); qit != ns_q.end(); ++qit) {
            qit->sum_n_q_re += qit->n_q_re * mats_cos - qit->n_q_im * mats_sin;
            qit->sum_n_q_im += qit->n_q_re * mats_sin + qit->n_q_im * mats_cos;
            qit->sum_s_q_re += qit->s_q_re * mats_cos - qit->s_q_im * mats_sin;
            qit->sum_s_q_im += qit->s_q_re * mats_sin + qit->s_q_im * mats_cos;
        }

        // collect data
        vector<double>::iterator C_rho_it, C_sigma_it;
        double prefactor = 1./beta/(matsubara ? omega_mats*omega_mats : 1);
        for (qit = ns_q.begin(), C_rho_it = C_rho_q.begin(),
                C_sigma_it = C_sigma_q.begin(); qit != ns_q.end();
                ++qit, ++C_rho_it, ++C_sigma_it) {
            *C_rho_it = prefactor
                        * (qit->sum_n_q_re * qit->sum_n_q_re
                           + qit->sum_n_q_im * qit->sum_n_q_im);
            *C_sigma_it = prefactor
                          * (qit->sum_s_q_re * qit->sum_s_q_re
                             + qit->sum_s_q_im * qit->sum_s_q_im);
        }
#ifdef MCL_PT
        measure[myrep].add("C_rho_q", C_rho_q);
        measure[myrep].add("C_sigma_q", C_sigma_q);
#else
        measure.add("C_rho_q", C_rho_q);
        measure.add("C_sigma_q", C_sigma_q);
#endif
    }
}


bool mc :: is_thermalized() {
    return therm_state.stage == thermalized;
}

void mc :: init() {
    random_init();

    // initialize states randomly
    state.resize(L);
    occ.resize(L, 0);
    bool place_holes_up = N_el_up > L/2;
    bool place_holes_down = N_el_down > L/2;
    el_state initial_state =
        static_cast<el_state>((place_holes_up ? up : empty)
                              | (place_holes_down ? down : empty));
    for (uint i = 0; i < L; i++) {
        state[i] = initial_state;
    }
    // place down spins
    uint things_to_place = place_holes_down ? (L-N_el_down) : N_el_down;
    while (things_to_place > 0) {
        int site = random0N(L);
        if ((static_cast<el_state>(state[site] & down) == down)
                == place_holes_down) {
            state[site] = static_cast<el_state>(state[site] ^ down);
            things_to_place--;
        }
    }
    // place up spins
    things_to_place = place_holes_up ? (L-N_el_up) : N_el_up;
    while (things_to_place > 0) {
        int site = random0N(L);
        if ((static_cast<el_state>(state[site] & up) == up)
                == place_holes_up) {
            state[site] = static_cast<el_state>(state[site] ^ up);
            things_to_place--;
        }
    }

    n = 0;
    n_hop = 0;
    M = (uint)(enlargement_factor * init_n_max);
    dublon_rejected = true;
    avg_worm_len = 0;
    worm_len_sample_size = 0;
    N_loop = vtx_visited * M;
    sm.resize(M, identity);

    // read mu value from database if available & desired
#ifdef MCL_PT
    muvec.resize(gvec.size());
    if (mus_file) {
        for (uint i = 0; i < gvec.size(); ++i) {
            stringstream fname;
            fname << "../mus/" << setprecision(4) << U << "_" << gvec[i] << "_"
                  << omega << ".mu";
            ifstream fstr(fname.str().c_str());
            if (fstr.is_open()) {
                fstr >> muvec[i];
            } else {
                muvec[i] = default_mu(gvec[i], omega);
            }
        }
    } else {
        for (uint i = 0; i < gvec.size(); ++i) {
            muvec[i] = default_mu(gvec[i], omega);
        }
    }
#else
    if (mus_file) {
        stringstream fname;
        fname << "../mus/" << setprecision(4) << U << "_" << g << "_" << omega
              << ".mu";
        ifstream fstr(fname.str().c_str());
        if (fstr.is_open()) {
            fstr >> mu;
        }
    }
#endif

    // set up adjustment of mu if desired
    if (mu_adjust) {
        lower_mu = mu - 0.5*mu_adjust_range;
        upper_mu = mu + 0.5*mu_adjust_range;
        mu = lower_mu;
    }

    therm_state.set_stage(initial_stage);
#ifndef MCL_PT
    recalc_directed_loop_probs();
#endif

    // add observables
#ifdef MCL_PT
    for (uint i = 0; i < gvec.size(); ++i) {
        measure[i].add_observable("N_up");
        measure[i].add_observable("N_down");
        measure[i].add_observable("dublon_rejection_rate");
        measure[i].add_observable("Energy");
    #ifdef MEASURE_KIN_ENERGY
        measure[i].add_observable("kinetic_Energy");
    #endif
        measure[i].add_observable("ph_density");
        measure[i].add_observable("W");
        measure[i].add_observable("W_sq");
        measure[i].add_observable("S_rho_q");
        measure[i].add_observable("S_sigma_q");
        measure[i].add_observable("chi_rho_pi");
        measure[i].add_observable("chi_sigma_pi");
        measure[i].add_vectorobservable("S_rho_r", L, bin_length);
        measure[i].add_vectorobservable("S_sigma_r", L, bin_length);
        if (calc_dyn) {
            measure[i].add_vectorobservable("C_rho_q", ns_q.size(), bin_length);
            measure[i].add_vectorobservable("C_sigma_q", ns_q.size(), bin_length);
        }
    }
#else
    measure.add_observable("N_up");
    measure.add_observable("N_down");
    measure.add_observable("dublon_rejection_rate");
    measure.add_observable("Energy");
    #ifdef MEASURE_KIN_ENERGY
    measure.add_observable("kinetic_Energy");
    #endif
    measure.add_observable("ph_density");
    measure.add_observable("W");
    measure.add_observable("W_sq");
    measure.add_observable("S_rho_q");
    measure.add_observable("S_sigma_q");
    measure.add_observable("chi_rho_pi");
    measure.add_observable("chi_sigma_pi");
    measure.add_vectorobservable("S_rho_r", L, bin_length);
    measure.add_vectorobservable("S_sigma_r", L, bin_length);
    if (calc_dyn) {
        measure.add_vectorobservable("C_rho_q", ns_q.size(), bin_length);
        measure.add_vectorobservable("C_sigma_q", ns_q.size(), bin_length);
    }
#endif
}

void mc :: write(string dir) {
    odump d(dir + "dump");
    random_write(d);
    d.write(therm_state);
    d.write(state);
    d.write(occ);
    d.write(sm);
    d.write(n);
    d.write(n_hop);
    d.write(dublon_rejected);
    d.write(avg_worm_len);
    d.write(worm_len_sample_size);
    d.write(N_loop);
    d.write(mu);
    if (mu_adjust) {
        d.write(lower_mu);
        d.write(upper_mu);
        d.write(N_mu);
        string mu_data_str(mu_data.str());
        string bisection_protocol_str(bisection_protocol.str());
        d.write(mu_data_str);
        d.write(bisection_protocol_str);
    }
    if (thermlog_interval > 0) {
        string thermlog_str(thermlog.str());
        d.write(thermlog_str);
    }
#ifdef MCL_PT
    d.write(muvec);
#endif
    d.write(bdoub_level);
    d.close();
    seed_write(dir + "seed");
    dir += "bins";
    ofstream f;
    f.open(dir.c_str());
    f << ( (is_thermalized()) ? therm_state.sweeps : 0 ) << endl;
    f.close();
}

bool mc :: read(string dir) {
    idump d(dir+"dump");
    if (!d) {
        return false;
    } else {
        random_read(d);
        d.read(therm_state);
        d.read(state);
        d.read(occ);
        d.read(sm);
        M = sm.size();
        d.read(n);
        d.read(n_hop);
        d.read(dublon_rejected);
        d.read(avg_worm_len);
        d.read(worm_len_sample_size);
        d.read(N_loop);
        d.read(mu);
        if (mu_adjust) {
            d.read(lower_mu);
            d.read(upper_mu);
            d.read(N_mu);
            string mu_data_str, bisection_protocol_str;
            d.read(mu_data_str);
            mu_data << mu_data_str;
            d.read(bisection_protocol_str);
            bisection_protocol << bisection_protocol_str;
        }
        if (thermlog_interval > 0) {
            string thermlog_str;
            d.read(thermlog_str);
            thermlog << thermlog_str;
        }
#ifdef MCL_PT
        d.read(muvec);
#endif
        d.read(bdoub_level);
        d.close();
#ifndef MCL_PT
        recalc_directed_loop_probs();
#endif
        return true;
    }
}

#ifdef MCL_PT
void mc :: write_output(string dir, int para) {
#else
void mc :: write_output(string dir) {
#endif
    // add evalables
    ofstream f;
    f.open(dir.c_str());
    f << "PARAMETERS" << endl;
#ifdef MCL_PT
    param.get_all_with_one_from_specified_array("@G",para,f);
    measure[para].get_statistics(f);
#else
    param.get_all(f);
    measure.get_statistics(f);
#endif
    double C = (U/4 > -abs(mu)) ? (U/4 + 2*abs(mu)) : (-U/4);
    f << "SIMULATION PROPERTIES" << endl
      << "C+delta+epsilon = " << C << " + " << delta << " + "
      << epsilon << " = " << C+delta+epsilon << endl
      << "operator string max. length: " << M << endl
      << "average worm length: " << avg_worm_len << endl
      << "number of loops per MCS: " << N_loop << endl
      << "BISECTION PROTOCOL" << endl
      << "mu_best U g omega" << endl
      << mu << " " << U << " " << g << " " << omega << endl
      << endl << endl
      << "mu N_mu" << endl
      << mu_data.str()
      << endl << endl
      << "lower_mu upper_mu" << endl
      << bisection_protocol.str() << endl
      << "THERMLOG" << endl
      << thermlog.str();
}

void mc :: init_assignments() {
    fill(assign_group.begin(), assign_group.end(), 0);
    assign_group.resize(N_WORM << 12, 0);
    fill(role.begin(), role.end(), no_role);
    role.resize(N_WORM << 12, no_role);
    assignment mat[3][3];
    for (int worm_i = 0; worm_i < N_WORM; ++worm_i) {
        mat[0][0].worm = static_cast<worm_type>(worm_i);
        for (int vtx_i = 0; vtx_i <= 256; ++vtx_i) {
            if (v_type[vtx_i] == W_invalid)
                continue;
            mat[0][0].vtx.int_repr = vtx_i;
            for (int ent_leg_i = bottom_left;
                    ent_leg_i <= top_left;
                    ++ent_leg_i) {
                // find suitable top left element as anchor for new group
                mat[0][0].ent_leg = static_cast<leg>(ent_leg_i);
                if (mat[0][0].worm == dublon_worm) {
                    el_state ent_state =
                        mat[0][0].vtx.get_state(mat[0][0].ent_leg);
                    if (ent_state != empty && ent_state != dublon) {
                        continue;
                    }
                }
                mat[0][0].exit_leg = mat[0][0].ent_leg;
                if (assign_group[mat[0][0].int_repr] != 0)
                    continue;

                // spawn first row and column of this group and fill in the
                // main diagonal
                int k = 1;
                int W[3];
                W[0] = v_type[mat[0][0].vtx.int_repr];
                for (int exit_leg_i = (ent_leg_i+1)%4;
                        exit_leg_i != ent_leg_i && k < 3;
                        exit_leg_i=(exit_leg_i+1)%4) {
                    mat[0][k] = mat[0][0];
                    mat[0][k].exit_leg = static_cast<leg>(exit_leg_i);
                    mat[k][0] = mat[0][k].flipped_assign();
                    W[k] = v_type[mat[k][0].vtx.int_repr];
                    if (W[k] != W_invalid) {
                        mat[k][k] = mat[k][0];
                        mat[k][k].exit_leg = mat[k][k].ent_leg;
                        ++k;
                    }
                }

                // find remaining off-diagonals in case of a 3x3 group
                if (k == 3) {
                    mat[1][2] = mat[1][1];
                    for (int exit_leg_i = bottom_left;
                            exit_leg_i <= top_left;
                            ++exit_leg_i) {
                        mat[1][2].exit_leg = static_cast<leg>(exit_leg_i);
                        if (mat[1][2].worm == dublon_worm) {
                            el_state ent_state =
                                mat[1][2].vtx.get_state(mat[1][2].ent_leg);
                            if (ent_state != empty && ent_state != dublon) {
                                continue;
                            }
                        }
                        if (mat[1][2] == mat[1][0] || mat[1][2] == mat[1][1])
                            continue;
                        mat[2][1] = mat[1][2].flipped_assign();
                        if (v_type[mat[2][1].vtx.int_repr] == W_invalid)
                            continue;
                        assert(mat[2][1].vtx == mat[2][2].vtx);
                        break;
                    }
                } else if (k == 2) {
                    W[2] = W_invalid;
                } else {
                    cerr << "Found an assignment group of size " << k << "!"
                         << endl;
                }

                int p[3];
                bubble_sort_perm(W, p, k);

                // identify group
                int group = -1;
                if (k == 3) {
                    switch (W[0]) {
                        case W_1p: group = 0; break;
                        case W_1m: group = 1; break;
                        case W_10:
                            switch (W[1]) {
                                case W_2p: group = 2; break;
                                case W_2m: group = 3; break;
                            }
                            break;
                        case W_2p: group = 4; break;
                        case W_2m: group = 5; break;
                    }
                    assert(W[2] == W_4);
                } else if (k == 2) {
                    switch (W[0]) {
                        case W_4:  group = 6; break;
                        case W_10: group = 7; break;
                        case W_1m: group = 8; break;
                        case W_2m: group = 9; break;
                    }
                }
                assert(group >= 0);

                // save results
                role[mat[p[0]][p[1]].int_repr] = role_a;
                assign_group[mat[p[0]][p[1]].int_repr] = group;
                role[mat[p[1]][p[0]].int_repr] = role_a;
                assign_group[mat[p[1]][p[0]].int_repr] = group;
                role[mat[p[0]][p[0]].int_repr] = role_b1;
                assign_group[mat[p[0]][p[0]].int_repr] = group;
                role[mat[p[1]][p[1]].int_repr] = role_b2;
                assign_group[mat[p[1]][p[1]].int_repr] = group;
                if (k == 3) {
                    role[mat[p[0]][p[2]].int_repr] = role_b;
                    assign_group[mat[p[0]][p[2]].int_repr] = group;
                    role[mat[p[2]][p[0]].int_repr] = role_b;
                    assign_group[mat[p[2]][p[0]].int_repr] = group;
                    role[mat[p[1]][p[2]].int_repr] = role_c;
                    assign_group[mat[p[1]][p[2]].int_repr] = group;
                    role[mat[p[2]][p[1]].int_repr] = role_c;
                    assign_group[mat[p[2]][p[1]].int_repr] = group;
                }
            }
        }
    }
}

void mc :: init_vertices() {
    v_type.resize(256);
    op_type.resize(256);
    vertex vtx;
    for (int vtx_i = 0; vtx_i < 256; ++vtx_i) {
        vtx.int_repr = vtx_i;
        v_type[vtx_i] = vtx.type();
        if (v_type[vtx_i] != W_invalid) {
            op_type[vtx_i] = vtx.op_type();
        }
    }
}

#ifdef MCL_PT
void mc :: change_to(int i) {
    change_parameter(i);
}

void mc :: change_parameter(int i) {
    myrep = i;
    g = gvec[myrep];
    mu = muvec[myrep];
    recalc_directed_loop_probs();
    if (myrep == 0)
        label = 1;
    else if (myrep == (int) (gvec.size() - 1))
        label = -1;
}

bool mc :: request_global_update() {
    return    beta == final_beta
           && therm_state.sweeps % pt_spacing == 0;
}

double mc :: get_weight(int f) {
    double other_delta;
    recalc_weights(other_weight, muvec[f], other_delta);
    double log_weight = 0.;
    int n_phonon = 0;
    bond_operator b;
    vertex vtx;
    current_state = state;
    for (uint i = 0; i < M; ++i) {
        if (sm[i] == identity)
            continue;
        b = sm[i];
        switch (b.type) {
            case electron_diag:
                vtx = diag_vertex_at_bond(current_state, b.bond);
                log_weight += log(other_weight[vtx.int_repr]
                                  / weight[vtx.int_repr]);
                break;
            case up_hopping:
                current_state[LEFT_SITE(b.bond)] =
                    flipped_state(current_state[LEFT_SITE(b.bond)], up_worm);
                current_state[RIGHT_SITE(b.bond)] =
                    flipped_state(current_state[RIGHT_SITE(b.bond)], up_worm);
                break;
            case down_hopping:
                current_state[LEFT_SITE(b.bond)] =
                    flipped_state(current_state[LEFT_SITE(b.bond)], down_worm);
                current_state[RIGHT_SITE(b.bond)] =
                    flipped_state(current_state[RIGHT_SITE(b.bond)], down_worm);
                break;
            case creator_left:
            case creator_right:
            case annihilator_left:
            case annihilator_right:
                n_phonon++;
                break;
        }
    }
    return n_phonon * log(gvec[f] / g) + log_weight;
}

int mc :: get_label() {
    return label;
}
#endif
