#include "mc.h"
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cmath>
#include <gsl/gsl_fit.h>

inline byte number_of_electrons (el_state s) {
    return (s & 1) + (s >> 1);
}

inline signed char local_magnetization (el_state s) {
    return (s & 1) - (s >> 1);
}

inline el_state flipped_state (el_state to_flip, worm_type what) {
    return static_cast<el_state>(to_flip ^ (what+1));
}

inline leg straight(leg l) {
    return static_cast<leg>(l ^ 3);
}

mc :: mc (string dir) {
    // initialize job parameters
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    T = param.value_or_default<double>("T", 1.);
    N_el_up = param.value_or_default<int>("N_el_up", L/2);
    N_el_down = param.value_or_default<int>("N_el_down", N_el_up);
    a = param.value_or_default<double>("A", 1.3);
    U = param.value_or_default<double>("U", 1.);
    omega = param.value_or_default<double>("OMEGA", 1.);
    g = param.value_or_default<double>("G", 0.);
    mu = param.value_or_default<double>("MU", g*g/omega);
    q_chi = param.value_or_default<double>("Q_CHI", M_PI);
    q_S = param.value_or_default<double>("Q_S", 2*M_PI/L);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 50000);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    vtx_visited = param.value_or_default<double>("VTX_VISITED", 2.0);
    Np = param.value_or_default<int>("N_P", 10);
    mu_adjust = param.value_or_default<bool>("MU_ADJUST", 1);
    mu_adjust_range = param.value_or_default<double>("MU_ADJUST_RANGE", 0.1*abs(mu));
    mu_adjust_N = param.value_or_default<int>("MU_ADJUST_N", 20);
    mu_adjust_therm = param.value_or_default<int>("MU_ADJUST_THERM", 1000);
    mu_adjust_sweep = param.value_or_default<int>("MU_ADJUST_SWEEP", 10000);
    assert(N_el_up <= L && N_el_down <= L);
    assert(N_el_up % 2 == 1 && N_el_down % 2 == 1);

    if (mu_adjust) {
        total_therm = 2*therm+2*mu_adjust_N*mu_adjust_sweep + (2*mu_adjust_N-2)*mu_adjust_therm;
    } else {
        total_therm = therm;
    }

    // initialize vectors
    weight.resize(256);
    vtx_type.resize(256);
    prob.resize(N_WORM<<12);
    subseq.resize(L, vector<subseq_node>());
    initial_Nd.resize(L);
    first.resize(L);
    last.resize(L);
    sum_n.resize(L);
    sum_s.resize(L);
    sum_nn.resize(L);
    sum_ss.resize(L);
    sum_m.resize(L);
    S_rho_r.resize(L);
    S_sigma_r.resize(L);
    chi_rho_r.resize(L);
    chi_sigma_r.resize(L);
    mean_m.resize(L);
    cos_q_chi.resize(L);
    sin_q_chi.resize(L);
    cos_q_S.resize(L);
    sin_q_S.resize(L);

    // calculate trigonometric factors for Fourier transforms
    for (uint s = 0; s < L; ++s) {
        cos_q_chi[s] = cos(q_chi*s);
        sin_q_chi[s] = sin(q_chi*s);
        cos_q_S[s] = cos(q_S*s);
        sin_q_S[s] = sin(q_S*s);
    }
}

void mc :: recalc_directed_loop_probs() {
    fill(weight.begin(), weight.end(), 0.);
    fill(vtx_type.begin(), vtx_type.end(), electron_diag);
    fill(prob.begin(), prob.end(), 0.);

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
        if (W[i] < 0)
            W[i] = 0.;
    }

    // calculate loop segment weights
    double a[][6] = {
                        {0, 0, 0, 0, 0, 0}, // b_1
                        {0, 0, 0, 0, 0, 0}, // b_2
                        { // a
                            0.5 * (W[2] + W[4] - W[6]),
                            0.5 * (W[0] + W[3] - W[6]),
                            0.5 * (W[1] + W[4] - W[6]),
                            0.5 * (W[1] + W[3] - W[6]),
                            0.5 * (W[4] + W[5] - W[6]),
                            0.5 * (W[3] + W[5] - W[6])
                        },
                        { // b
                            0.5 * (W[2] - W[4] + W[6]),
                            0.5 * (W[0] - W[3] + W[6]),
                            0.5 * (W[1] - W[4] + W[6]),
                            0.5 * (W[1] - W[3] + W[6]),
                            0.5 * (W[4] - W[5] + W[6]),
                            0.5 * (W[3] - W[5] + W[6])
                        },
                        { // c
                            0.5 * (-W[2] + W[4] + W[6]),
                            0.5 * (-W[0] + W[3] + W[6]),
                            0.5 * (-W[1] + W[4] + W[6]),
                            0.5 * (-W[1] + W[3] + W[6]),
                            0.5 * (-W[4] + W[5] + W[6]),
                            0.5 * (-W[3] + W[5] + W[6])
                        }
                    };
    for (uint i = 0; i < 6; ++i) {
        a[1][i] = (a[3][i] < 0) ? -2*a[3][i] : 0;
        a[0][i] = (a[4][i] < 0) ? -2*a[4][i] : 0;
        a[3][i] += -a[0][i]/2 + a[1][i]/2;
        a[4][i] += a[0][i]/2 - a[1][i]/2;
        a[2][i] -= a[0][i]/2 + a[1][i]/2;
    }

    // determine epsilon
    epsilon = param.value_or_default<double>("EPSILON", -1.);
    double epsilon_min = 0.0;
    for (uint i = 0; i < 6; ++i) {
        if (a[2][i] < -epsilon_min) {
            epsilon_min = -a[2][i];
        }
    }
    if (epsilon < 0) {
        epsilon = epsilon_min;
    } else {
        assert(epsilon >= epsilon_min);
    }
    for (uint i = 0; i < 6; ++i) {
        a[2][i] += epsilon;
    }

    // parse vertex weights
    int vtx, j, i;
    for (uint i = 0; i < 6; ++i) {
        W[i] += epsilon;
    }
    ifstream file1("../vertex_types.txt");
    if (!file1.is_open()) {
        cerr << "Could not open file vertex_types.txt" << endl;
        exit(1);
    }
    while (file1 >> vtx >> j >> i) {
        weight[vtx] = W[i-1];
        vtx_type[vtx] = static_cast<operator_type>(j-1);
    }
    file1.close();

    // calculate transition probabilities
    ifstream file2("../assignments.txt");
    if (!file2.is_open()) {
        cerr << "Could not open file assignments.txt" << endl;
        exit(1);
    }
    assignment assign;
    while (file2 >> assign.int_repr >> i >> j) {
        if (assign.worm >= N_WORM) {
            continue;
        }
        if (i >= 6 || weight[assign.vtx.int_repr] == 0.) {    // 2x2 group
            prob[assign.int_repr] = 1.;
        } else {
            prob[assign.int_repr] = a[j][i] / weight[assign.vtx.int_repr];
        }
    }
    file2.close();

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
    vtx_type.clear();
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
    sum_n.clear();
    sum_s.clear();
    sum_nn.clear();
    sum_ss.clear();
    sum_m.clear();
    S_rho_r.clear();
    S_sigma_r.clear();
    chi_rho_r.clear();
    chi_sigma_r.clear();
    mean_m.clear();
}

void mc :: do_update() {
    // switch mu value as necessary
    if (sweep < total_therm && mu_adjust) {
        if (sweep == therm + ((mu_index<0) ? 1 : 2)*mu_adjust_sweep + (mu_index+mu_adjust_N-1)*(mu_adjust_sweep+mu_adjust_therm)) {
            if (sweep == total_therm-therm) {
                for (uint k = 0; k < mu_adjust_N; ++k) {
                    N_mus[k] = N_mus[k]/2/mu_adjust_sweep - (N_el_up+N_el_down);
                }
                double m, b, c00, c01, c11, sumsq;
                gsl_fit_linear(&mus[0], 1, &N_mus[0], 1, mu_adjust_N, &b, &m, &c00, &c01, &c11, &sumsq);
                mu = -b / m;
                recalc_directed_loop_probs();
            } else {
                mu = mus[abs(++mu_index)];
                recalc_directed_loop_probs();
            }
        }
    }

    // exiting thermalization stage
    if (sweep == total_therm) {
        N_loop = (uint)(vtx_visited / avg_worm_len * M);
    }

    // diagonal update & subsequence construction
    for_each(subseq.begin(), subseq.end(), mem_fun_ref(&vector<subseq_node>::clear));
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
            if (random01() < NB/T*weight[vtx.int_repr]/(M-n)) {
                sm[i] = b;
                n++;
            }
        } else if (sm[i].type == electron_diag) {
            vertex vtx = diag_vertex_at_bond(current_state, sm[i].bond);
            if (random01() < (M-n+1)/(NB/T*weight[vtx.int_repr])) {
                sm[i] = identity;
                n--;
            }
        }

        // phononic diagonal update
        if (sm[i] == identity) {
            bond_operator b;
            b.type = phonon_diag;
            b.bond = random0N(NB)+1;
            int cocc = current_occ[LEFT_SITE(b.bond)];
            if (random01() < NB/T*omega*(Np-cocc)/(M-n)) {
                sm[i] = b;
                n++;
            }
        } else if (sm[i].type == phonon_diag) {
            int cocc = current_occ[LEFT_SITE(sm[i].bond)];
            if (random01() < (M-n+1)/(NB/T*omega*(Np-cocc))) {
                sm[i] = identity;
                n--;
            }
        }

        bond_operator b = sm[i];
        if (b != identity) {
            // append to the phonon subsequences
            if (b.type == phonon_diag) {
                if (subseq[LEFT_SITE(b.bond)].empty()) {
                    ++initial_Nd[LEFT_SITE(b.bond)];
                } else {
                    ++(subseq[LEFT_SITE(b.bond)].back().Nd);
                }
            } else if (b.type != up_hopping && b.type != down_hopping) {
                unsigned short site;
                if (b.type == electron_diag) {
                    site = random0N(2) ? LEFT_SITE(b.bond) : RIGHT_SITE(b.bond);
                } else if (b.type == creator_left || b.type == creator_right) {
                    site = (b.type == creator_left) ? LEFT_SITE(b.bond)
                                                    : RIGHT_SITE(b.bond);
                } else if (b.type == annihilator_left || b.type == annihilator_right) {
                    site = (b.type == annihilator_left) ? LEFT_SITE(b.bond)
                                                        : RIGHT_SITE(b.bond);
                }
                vertex vtx = diag_vertex_at_bond(current_state, b.bond);
                byte n_el = number_of_electrons(current_state[site]);
                subseq_node newNode;
                newNode.i = i;
                newNode.Nd = 0;
                newNode.m = current_occ[site];
                newNode.r = weight[vtx.int_repr]/n_el;
                subseq[site].push_back(newNode);
            }

            // propagation of state
            if (b.type == up_hopping) {
                current_state[LEFT_SITE(b.bond)] = flipped_state(current_state[LEFT_SITE(b.bond)], up_worm);
                current_state[RIGHT_SITE(b.bond)] = flipped_state(current_state[RIGHT_SITE(b.bond)], up_worm);
            } else if (b.type == down_hopping) {
                current_state[LEFT_SITE(b.bond)] = flipped_state(current_state[LEFT_SITE(b.bond)], down_worm);
                current_state[RIGHT_SITE(b.bond)] = flipped_state(current_state[RIGHT_SITE(b.bond)], down_worm);
            } else if (b.type == creator_left) {
                current_occ[LEFT_SITE(b.bond)]++;
            } else if (b.type == creator_right) {
                current_occ[RIGHT_SITE(b.bond)]++;
            } else if (b.type == annihilator_left) {
                current_occ[LEFT_SITE(b.bond)]--;
            } else if (b.type == annihilator_right) {
                current_occ[RIGHT_SITE(b.bond)]--;
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
    for (uint s = 0; s < L; ++s) {
        vector<subseq_node>::iterator i1, i2;
        if (subseq[s].size() < 2)
            continue;
        for (uint i = 0; i < M/L; ++i) {
            i2 = subseq[s].begin() + random0N(subseq[s].size());
            i1 = i2++;
            if (i2 == subseq[s].end()) {
                i2 = subseq[s].begin();
            }
            if (sm[i1->i].type == electron_diag && sm[i2->i].type == electron_diag) {
                if (random0N(2)) { // (H_1, H_1) -> (H_5, H_4)
                    double prob = g*g * i1->m / (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m+1)/(Np-i1->m), i1->Nd);
                    if (random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_5^R
                            sm[i1->i].type = annihilator_right;
                            sm[i1->i].bond = LEFT_BOND(s);
                        } else { // H_5^L
                            sm[i1->i].type = annihilator_left;
                            sm[i1->i].bond = RIGHT_BOND(s);
                        }
                        if (type % 2) { // H_4^R
                            sm[i2->i].type = creator_right;
                            sm[i2->i].bond = LEFT_BOND(s);
                        } else { // H_4^L
                            sm[i2->i].type = creator_left;
                            sm[i2->i].bond = RIGHT_BOND(s);
                        }
                        --(i2->m);
                    }
                } else { // (H_1, H_1) -> (H_4, H_5)
                    double prob = g*g * (i1->m+1) / (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m-1)/(Np-i1->m), i1->Nd);
                    if (i1->m < Np && random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_4^R
                            sm[i1->i].type = creator_right;
                            sm[i1->i].bond = LEFT_BOND(s);
                        } else { // H_4^L
                            sm[i1->i].type = creator_left;
                            sm[i1->i].bond = RIGHT_BOND(s);
                        }
                        if (type % 2) { // H_5^R
                            sm[i2->i].type = annihilator_right;
                            sm[i2->i].bond = LEFT_BOND(s);
                        } else { // H_5^L
                            sm[i2->i].type = annihilator_left;
                            sm[i2->i].bond = RIGHT_BOND(s);
                        }
                        ++(i2->m);
                    }
                }
            } else if ((sm[i1->i].type == creator_left || sm[i1->i].type == creator_right)
                    && (sm[i2->i].type == annihilator_left || sm[i2->i].type == annihilator_right)) {
                if (random0N(2)) { // (H_4, H_5) -> (H_5, H_4)
                    double prob = 1. * (i1->m) / (i1->m+1)
                                  * pow(1.*(Np-i1->m+1)/(Np-i1->m-1), i1->Nd);
                    if (random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_5^R
                            sm[i1->i].type = annihilator_right;
                            sm[i1->i].bond = LEFT_BOND(s);
                        } else { // H_5^L
                            sm[i1->i].type = annihilator_left;
                            sm[i1->i].bond = RIGHT_BOND(s);
                        }
                        if (type % 2) { // H_4^R
                            sm[i2->i].type = creator_right;
                            sm[i2->i].bond = LEFT_BOND(s);
                        } else { // H_4^L
                            sm[i2->i].type = creator_left;
                            sm[i2->i].bond = RIGHT_BOND(s);
                        }
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
            } else if ((sm[i1->i].type == annihilator_left || sm[i1->i].type == annihilator_right)
                    && (sm[i2->i].type == creator_left || sm[i2->i].type == creator_right)) {
                if (random0N(2)) { // (H_5, H_4) -> (H_4, H_5)
                    double prob = 1. / (i1->m) * (i1->m+1)
                                  * pow(1.*(Np-i1->m-1)/(Np-i1->m+1), i1->Nd);
                    if (i1->m < Np && random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_4^R
                            sm[i1->i].type = creator_right;
                            sm[i1->i].bond = LEFT_BOND(s);
                        } else { // H_4^L
                            sm[i1->i].type = creator_left;
                            sm[i1->i].bond = RIGHT_BOND(s);
                        }
                        if (type % 2) { // H_5^R
                            sm[i2->i].type = annihilator_right;
                            sm[i2->i].bond = LEFT_BOND(s);
                        } else { // H_5^L
                            sm[i2->i].type = annihilator_left;
                            sm[i2->i].bond = RIGHT_BOND(s);
                        }
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
    if (!is_thermalized() && a*n > M) {
        M = a*n;
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
                link[list_position(bottom_left, p).int_repr] = last[LEFT_SITE(sm[i].bond)];
                link[last[LEFT_SITE(sm[i].bond)].int_repr] = list_position(bottom_left, p);
                last[LEFT_SITE(sm[i].bond)] = list_position(top_left, p);
            }
            if (first[RIGHT_SITE(sm[i].bond)] == invalid_pos) {
                first[RIGHT_SITE(sm[i].bond)] = list_position(bottom_right, p);
                last[RIGHT_SITE(sm[i].bond)] = list_position(top_right, p);
            } else {
                link[list_position(bottom_right, p).int_repr] = last[RIGHT_SITE(sm[i].bond)];
                link[last[RIGHT_SITE(sm[i].bond)].int_repr] = list_position(bottom_right, p);
                last[RIGHT_SITE(sm[i].bond)] = list_position(top_right, p);
            }

            // determine vertex type
            vtx[p].bottom_left = current_state[LEFT_SITE(sm[i].bond)];
            vtx[p].bottom_right = current_state[RIGHT_SITE(sm[i].bond)];
            if (sm[i].type == up_hopping) {
                current_state[LEFT_SITE(sm[i].bond)] = flipped_state(current_state[LEFT_SITE(sm[i].bond)], up_worm);
                current_state[RIGHT_SITE(sm[i].bond)] = flipped_state(current_state[RIGHT_SITE(sm[i].bond)], up_worm);
            } else if (sm[i].type == down_hopping) {
                current_state[LEFT_SITE(sm[i].bond)] = flipped_state(current_state[LEFT_SITE(sm[i].bond)], down_worm);
                current_state[RIGHT_SITE(sm[i].bond)] = flipped_state(current_state[RIGHT_SITE(sm[i].bond)], down_worm);
            }
            vtx[p].top_right = current_state[RIGHT_SITE(sm[i].bond)];
            vtx[p].top_left = current_state[LEFT_SITE(sm[i].bond)];

            if (sm[i].type == creator_left || sm[i].type == annihilator_left) {
                lock[p] = left_lock; // require at least one electron on left site
            } else if (sm[i].type == creator_right || sm[i].type == annihilator_right) {
                lock[p] = right_lock; // require at least one electron on right site
            } else if (sm[i].type == phonon_diag) {
                lock[p] = total_lock;
            } else {
                lock[p] = unlocked;
            }

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
                if (lock[j.index] == total_lock) {
                    assign.exit_leg = straight(j.vtx_leg); // continue straight
                } else if (lock[j.index] != unlocked) {
                    if ((lock[j.index] == left_lock && (j.vtx_leg == bottom_left || j.vtx_leg == top_left))
                        || (lock[j.index] == right_lock && (j.vtx_leg == bottom_right || j.vtx_leg == top_right))) {
                        if (flipped_state(vtx[j.index].get_state(j.vtx_leg), worm) == empty) {
                            assign.exit_leg = j.vtx_leg; // bounce
                        } else {
                            if (vtx[j.index].get_state(j.vtx_leg) == dublon) {
                                if (random0N(2)) {
                                    assign.exit_leg = straight(j.vtx_leg); // continue straight
                                } else {
                                    assign.exit_leg = j.vtx_leg; // bounce
                                }
                            } else {
                                assign.exit_leg = straight(j.vtx_leg); // continue straight
                            }
                        }
                    } else {
                        assign.exit_leg = straight(j.vtx_leg); // continue straight
                    }
                } else {
                    r = random01();
                    byte exit_leg_i;
                    for (exit_leg_i = bottom_left; exit_leg_i <= top_left; ++exit_leg_i) {
                        assign.exit_leg = static_cast<leg>(exit_leg_i);
                        if (r < prob[assign.int_repr])
                            break;
                    }
                    assert(exit_leg_i < 4); // assert that break was called
                }
                // do not count bounces into the worm length
                if (assign.ent_leg == assign.exit_leg) {
                    --k;
                }
                // flip the vertex:
                vtx[j.index] = assign.flipped_vtx();
                j.vtx_leg = assign.exit_leg; // exit leg position in linked list
                if (j == j0)    // loop closed (SS02, Fig. 4b)
                    break;
                j = link[j.int_repr];
                if (j == j0)    // loop closed (SS02, Fig. 4a)
                    break;
            }
            // logging worm length
            if (sweep < therm && sweep >= therm/2) {
                avg_worm_len *= 1.*worm_len_sample_size / (worm_len_sample_size+1);
                avg_worm_len += 1.*k / (++worm_len_sample_size);
            }
        }

        // mapping back to operator sequence
        p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == identity)
                continue;
            if (sm[i].type == electron_diag || sm[i].type == up_hopping || sm[i].type == down_hopping)
                sm[i].type = vtx_type[vtx[p].int_repr];
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
    if (sweep < total_therm && mu_adjust && sweep >= therm + ((mu_index<=0) ? 0 : 1)*mu_adjust_sweep + (mu_index+mu_adjust_N-1)*(mu_adjust_sweep+mu_adjust_therm) && sweep < therm + ((mu_index<0) ? 1 : 2)*mu_adjust_sweep + (mu_index+mu_adjust_N-1)*(mu_adjust_sweep+mu_adjust_therm)) {
        for (uint s = 0; s < L; ++s) {
            N_mus[abs(mu_index)] += number_of_electrons(state[s]);
        }
    }
    ++sweep;
}

void mc :: do_measurement() {
    uint N_up = 0, N_down = 0;
    for (uint s = 0; s < L; ++s) {
        N_up += state[s] == up || state[s] == dublon;
        N_down += state[s] == down || state[s] == dublon;
    }
    measure.add("N_up", N_up);
    measure.add("N_down", N_down);
    measure.add("dublon_rejection_rate", dublon_rejected);
    
    // skip measurement if particle numbers are not right
    if (N_up != N_el_up || N_down != N_el_down)
        return;

    double C = (U/4 > -abs(mu)) ? (U/4 + 2*abs(mu)) : (-U/4);
    double energy = -T * n + NB*(C+epsilon) + L*omega*Np;

    // add data to measurement
    measure.add("Energy", energy);

    // calculate correlation functions and susceptibilities
    fill(sum_n.begin(), sum_n.end(), 0);
    fill(sum_s.begin(), sum_s.end(), 0);
    fill(sum_nn.begin(), sum_nn.end(), 0);
    fill(sum_ss.begin(), sum_ss.end(), 0);
    fill(sum_m.begin(), sum_m.end(), 0);
    current_state = state;
    current_occ = occ;
    uint p = 0;
    int n_s, n_0, s_s, s_0;
    for (uint i = 0; p < n; ++i) {
        if (sm[i] == identity)
            continue;

        n_0 = number_of_electrons(current_state[0]);
        s_0 = local_magnetization(current_state[0]);
        for (uint s = 0; s < L; ++s) {
            n_s = number_of_electrons(current_state[s]);
            s_s = local_magnetization(current_state[s]);
            sum_n[s] += n_s;
            sum_nn[s] += n_s * n_0;
            sum_s[s] += s_s;
            sum_ss[s] += s_s * s_0;
            sum_m[s] += current_occ[s];
        }

        bond_operator b = sm[i];
        // propagation of state
        if (b.type == up_hopping) {
            current_state[LEFT_SITE(b.bond)] = flipped_state(current_state[LEFT_SITE(b.bond)], up_worm);
            current_state[RIGHT_SITE(b.bond)] = flipped_state(current_state[RIGHT_SITE(b.bond)], up_worm);
        } else if (b.type == down_hopping) {
            current_state[LEFT_SITE(b.bond)] = flipped_state(current_state[LEFT_SITE(b.bond)], down_worm);
            current_state[RIGHT_SITE(b.bond)] = flipped_state(current_state[RIGHT_SITE(b.bond)], down_worm);
        } else if (b.type == creator_left) {
            current_occ[LEFT_SITE(b.bond)]++;
        } else if (b.type == creator_right) {
            current_occ[RIGHT_SITE(b.bond)]++;
        } else if (b.type == annihilator_left) {
            current_occ[LEFT_SITE(b.bond)]--;
        } else if (b.type == annihilator_right) {
            current_occ[RIGHT_SITE(b.bond)]--;
        }
        ++p;
    }

    // calculate real space correlation functions and susceptibilities
    for (uint s = 0; s < L; ++s) {
        n_s = number_of_electrons(current_state[s]);
        s_s = local_magnetization(current_state[s]);
        sum_nn[s] += n_s * n_0;
        sum_ss[s] += s_s * s_0;
        S_rho_r[s] = 1./(n+1)*sum_nn[s];
        S_sigma_r[s] = 1./(n+1)*sum_ss[s];
        chi_rho_r[s] = 1./T/n/(n+1)*sum_n[s]*sum_n[0] + 1./T/(n+1)/(n+1)*sum_nn[s];
        chi_sigma_r[s] = 1./T/n/(n+1)*sum_s[s]*sum_s[0] + 1./T/(n+1)/(n+1)*sum_ss[s];
        mean_m[s] = 1./n*sum_m[s];
        // accumulate ms
        if (s > 0) {
            mean_m[s] += mean_m[s-1];
        }
    }
    mean_m[L-1] /= L;

    // Fourier transform
    double S_rho_q_re = 0.0;
    double S_rho_q_im = 0.0;
    double S_sigma_q_re = 0.0;
    double S_sigma_q_im = 0.0;
    double chi_rho_q_re = 0.0;
    double chi_rho_q_im = 0.0;
    double chi_sigma_q_re = 0.0;
    double chi_sigma_q_im = 0.0;
    for (uint s = 0; s < L; ++s) {
        S_rho_q_re += cos_q_S[s] * S_rho_r[s];
        S_rho_q_im += sin_q_S[s] * S_rho_r[s];
        S_sigma_q_re += cos_q_S[s] * S_sigma_r[s];
        S_sigma_q_im += sin_q_S[s] * S_sigma_r[s];
        chi_rho_q_re += cos_q_chi[s] * chi_rho_r[s];
        chi_rho_q_im += sin_q_chi[s] * chi_rho_r[s];
        chi_sigma_q_re += cos_q_chi[s] * chi_sigma_r[s];
        chi_sigma_q_im += sin_q_chi[s] * chi_sigma_r[s];
    }

    measure.add("S_rho_q_re", S_rho_q_re);
    measure.add("S_rho_q_im", S_rho_q_im);
    measure.add("S_sigma_q_re", S_sigma_q_re);
    measure.add("S_sigma_q_im", S_sigma_q_im);
    measure.add("chi_rho_q_re", chi_rho_q_re);
    measure.add("chi_rho_q_im", chi_rho_q_im);
    measure.add("chi_sigma_q_re", chi_sigma_q_re);
    measure.add("chi_sigma_q_im", chi_sigma_q_im);
    measure.add("ph_density", mean_m[L-1]);
}


bool mc :: is_thermalized() {
    return (sweep > total_therm);
}

void mc :: init() {
    random_init();

    // initialize states randomly
    state.resize(L);
    occ.resize(L, 0);
    bool place_holes_up = N_el_up > L/2;
    bool place_holes_down = N_el_down > L/2;
    el_state initial_state = static_cast<el_state>((place_holes_up ? up : empty)
                                                   | (place_holes_down ? down : empty));
    for (uint i = 0; i < L; i++) {
        state[i] = initial_state;
    }
    // place down spins
    uint things_to_place = place_holes_down ? (L-N_el_down) : N_el_down;
    while (things_to_place > 0) {
        int site = random0N(L);
        if ((static_cast<el_state>(state[site] & down) == down) == place_holes_down) {
            state[site] = static_cast<el_state>(state[site] ^ down);
            things_to_place--;
        }
    }
    // place up spins
    things_to_place = place_holes_up ? (L-N_el_up) : N_el_up;
    while (things_to_place > 0) {
        int site = random0N(L);
        if ((static_cast<el_state>(state[site] & up) == up) == place_holes_up) {
            state[site] = static_cast<el_state>(state[site] ^ up);
            things_to_place--;
        }
    }

    n = 0;
    sweep=0;
    M = (uint)(a * init_n_max);
    dublon_rejected = true;
    avg_worm_len = 0;
    worm_len_sample_size = 0;
    N_loop = vtx_visited * M;
    sm.resize(M, identity);

    // set up adjustment of mu if desired
    if (mu_adjust) {
        mus.resize(mu_adjust_N);
        N_mus.resize(mu_adjust_N, 0);

        // calculate mu values for adjustment
        for (uint i = 0; i < mu_adjust_N; i++) {
            mus[i] = (mu-mu_adjust_range) + 2*mu_adjust_range/(mu_adjust_N-1)*i;
        }
        mu_index = -(mu_adjust_N-1);
        mu = mus[abs(mu_index)];
    }
    recalc_directed_loop_probs();

    // add observables
    measure.add_observable("N_up");
    measure.add_observable("N_down");
    measure.add_observable("dublon_rejection_rate");
    measure.add_observable("Energy");
    measure.add_observable("S_rho_q_re");
    measure.add_observable("S_rho_q_im");
    measure.add_observable("S_sigma_q_re");
    measure.add_observable("S_sigma_q_im");
    measure.add_observable("chi_rho_q_re");
    measure.add_observable("chi_rho_q_im");
    measure.add_observable("chi_sigma_q_re");
    measure.add_observable("chi_sigma_q_im");
    measure.add_observable("ph_density");
}

void mc :: write(string dir) {
    odump d(dir + "dump");
    random_write(d);
    d.write(sweep);
    d.write(state);
    d.write(occ);
    d.write(sm);
    d.write(n);
    d.write(dublon_rejected);
    d.write(avg_worm_len);
    d.write(worm_len_sample_size);
    d.write(N_loop);
    d.write(mu);
    if (mu_adjust) {
        d.write(mu_index);
        d.write(mus);
        d.write(N_mus);
    }
    d.close();
    seed_write(dir + "seed");
    dir += "bins";
    ofstream f;
    f.open(dir.c_str());
    f << ( (is_thermalized()) ? sweep-total_therm : 0 ) << endl;
    f.close();
}

bool mc :: read(string dir) {
    idump d(dir+"dump");
    if (!d) {
        return false;
    } else {
        random_read(d);
        d.read(sweep);
        d.read(state);
        d.read(occ);
        d.read(sm);
        M = sm.size();
        d.read(n);
        d.read(dublon_rejected);
        d.read(avg_worm_len);
        d.read(worm_len_sample_size);
        d.read(N_loop);
        d.read(mu);
        if (mu_adjust) {
            d.read(mu_index);
            d.read(mus);
            d.read(N_mus);
        }
        d.close();
        recalc_directed_loop_probs();
        return true;
    }
}

void mc :: write_output(string dir) {
    // add evalables
    ofstream f;
    f.open(dir.c_str());
    f << "PARAMETERS" << endl;
    param.get_all(f);
    measure.get_statistics(f);
    f << "SIMULATION PROPERTIES" << endl;
    double C = (U/4 > -abs(mu)) ? (U/4 + 2*abs(mu)) : (-U/4);
    f << "C+epsilon = " << C+epsilon << endl;
    f << "operator string max. length: " << M << endl;
    f << "average worm length: " << avg_worm_len << endl;
    f << "number of loops per MCS: " << N_loop << endl;
    f << "mu = " << mu << endl;
    if (mu_adjust) {
        f << "table of mu values and mean particle number deviations from desired values:" << endl;
        for (uint i = 0; i < mus.size(); ++i) {
            f << mus[i] << " " << N_mus[i] << endl;
        }
    }
}
