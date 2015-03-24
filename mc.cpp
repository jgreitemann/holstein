#include "mc.h"
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>
#include <gsl/gsl_fit.h>

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

mc :: mc (string dir) {
    // initialize job parameters
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    T = param.value_or_default<double>("T", .05);
    epsilon = param.value_or_default<double>("EPSILON", 0.01);
    N_el_up = param.value_or_default<int>("N_el_up", L/2);
    N_el_down = param.value_or_default<int>("N_el_down", N_el_up);
    a = param.value_or_default<double>("A", 1.3);
    U = param.value_or_default<double>("U", 0.);
    omega = param.value_or_default<double>("OMEGA", 1.);
    g = param.value_or_default<double>("G", 0.);
    mu = param.value_or_default<double>("MU", g*g/omega);
    q_chi = param.value_or_default<double>("Q_CHI", M_PI);
    q_S = param.value_or_default<double>("Q_S", 2*M_PI/L);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 50000);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    vtx_visited = param.value_or_default<double>("VTX_VISITED", 2.0);
    Np = param.value_or_default<int>("N_P", 20);
    mu_adjust = param.value_or_default<bool>("MU_ADJUST", 0);
    mu_adjust_range = param.value_or_default<double>("MU_ADJUST_RANGE", 0.1);
    mu_adjust_N = param.value_or_default<int>("MU_ADJUST_N", 10);
    mu_adjust_therm = param.value_or_default<int>("MU_ADJUST_THERM", 5000);
    mu_adjust_sweep = param.value_or_default<int>("MU_ADJUST_SWEEP", 10000);
    assert(N_el_up <= L && N_el_down <= L);
    assert(N_el_up % 2 == 1 && N_el_down % 2 == 1);

    if (mu_adjust) {
        total_therm = 2*therm + 2*mu_adjust_N*mu_adjust_sweep
                      + (2*mu_adjust_N-2) * mu_adjust_therm;
    } else {
        total_therm = therm;
    }

    // initialize vectors
    init_vertices();
    init_assignments();
    weight.resize(256);
    prob.resize(N_WORM<<12);
    subseq.resize(L, vector<subseq_node>());
    initial_Nd.resize(L);
    first.resize(L);
    last.resize(L);
    sum_n.resize(L);
    sum_s.resize(L);
    sum_m.resize(L);
    sum_nn.resize(L);
    sum_ss.resize(L);
    n_p.resize(L);
    s_p.resize(L);
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
    double a[no_role][N_GROUP];
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
        int nextstop = therm                           // initial thermalization
            + ((mu_index<0) ? 1 : 2) * mu_adjust_sweep // reversal-point sweep
            + (mu_index+mu_adjust_N-1) * (mu_adjust_sweep+mu_adjust_therm);
        if (sweep == nextstop) {
            if (sweep == total_therm-therm) { // final stop, determine mu
                for (uint k = 0; k < mu_adjust_N; ++k) {
                    N_mus[k] = N_mus[k]/2/mu_adjust_sweep - (N_el_up+N_el_down);
                }
                double m, b, c00, c01, c11, sumsq;
                gsl_fit_linear(&mus[0], 1, &N_mus[0], 1, mu_adjust_N,
                               &b, &m, &c00, &c01, &c11, &sumsq);
                mu = -b / m;

                // save results to database
                stringstream fname;
                fname << "../mus/" << setprecision(4) << U << "_" << g << "_"
                      << omega << ".mu";
                ofstream fstr(fname.str().c_str());
                if (fstr.is_open()) {
                    fstr << mu << " " << U << " " << g << " " << omega << endl
                         << "# (above) mu_best U g omega" << endl
                         << endl << endl
                         << "# (below) mu DeltaN" << endl;
                    for (uint i = 0; i < mu_adjust_N; ++i) {
                        fstr << mus[i] << " " << N_mus[i] << endl;
                    }
                    fstr << "# L = " << L << endl
                         << "# T = " << T << endl
                         << "# mu_adjust_therm = " << mu_adjust_therm << endl
                         << "# mu_adjust_sweep = " << mu_adjust_sweep << endl
                         << "# linear model: f(x) = " << m << " * x + "
                         << b << endl;
                } else {
                    cerr << "Warning: could not write results from mu"
                            "adjustment to file " << fname.str() << endl;
                }

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
            if (   sm[i1->i].type == electron_diag
                && sm[i2->i].type == electron_diag) {
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
            } else if (   (   sm[i1->i].type == creator_left
                           || sm[i1->i].type == creator_right)
                       && (   sm[i2->i].type == annihilator_left
                           || sm[i2->i].type == annihilator_right)) {
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
            } else if (   (   sm[i1->i].type == annihilator_left
                           || sm[i1->i].type == annihilator_right)
                       && (   sm[i2->i].type == creator_left
                           || sm[i2->i].type == creator_right)) {
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
            if (sweep < therm && sweep >= therm/2) {
                avg_worm_len *= 1.*worm_len_sample_size
                                / (worm_len_sample_size+1);
                avg_worm_len += 1.*k / (++worm_len_sample_size);
            }
        }

        // mapping back to operator sequence
        p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == identity)
                continue;
            switch (sm[i].type) {
                case electron_diag:
                case up_hopping:
                case down_hopping:
                    sm[i].type = op_type[vtx[p].int_repr];
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
    if (   sweep < total_therm && mu_adjust
        && sweep >= therm + ((mu_index<=0) ? 0 : 1) * mu_adjust_sweep
                    + (mu_index+mu_adjust_N-1)*(mu_adjust_sweep+mu_adjust_therm)
        && sweep < therm + ((mu_index<0) ? 1 : 2) * mu_adjust_sweep
                   + (mu_index+mu_adjust_N-1)
                        * (mu_adjust_sweep+mu_adjust_therm)) {
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
    double energy = -T * n + NB*(C+delta+epsilon) + L*omega*Np;

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
    for (uint i = 0; p < n; ++i) {
        if (sm[i] == identity)
            continue;

        for (uint j = 0; j < L; ++j) {
            n_p[j] = number_of_electrons(current_state[j]);
            s_p[j] = local_magnetization(current_state[j]);
            sum_n[j] += n_p[j];
            sum_s[j] += s_p[j];
            sum_m[j] += current_occ[j];
        }

        for (uint r = 0; r < L; ++r) {
            for (uint j = 0; j < L; ++j) {
                sum_nn[r] += n_p[j+r] * n_p[j];
                sum_ss[r] += s_p[j+r] * s_p[j];
            }
        }

        bond_operator b = sm[i];
        // propagation of state
        switch (b.type) {
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
        ++p;
    }

    // calculate real space correlation functions and susceptibilities
    for (uint j = 0; j < L; ++j) {
        n_p[j] = number_of_electrons(current_state[j]);
        s_p[j] = local_magnetization(current_state[j]);
    }
    for (uint r = 0; r < L; ++r) {
        chi_rho_r[r] = 0.;
        chi_sigma_r[r] = 0.;
        for (uint j = 0; j < L; ++j) {
            sum_nn[r] += n_p[j+r] * n_p[j];
            sum_ss[r] += s_p[j+r] * s_p[j];
            chi_rho_r[r] += 1./T/L/n/(n+1) * sum_n[j+r] * sum_n[j];
            chi_sigma_r[r] += 1./T/L/n/(n+1) * sum_s[j+r] * sum_s[j];
        }
        S_rho_r[r] = 1./L/(n+1)*sum_nn[r];
        S_sigma_r[r] = 1./L/(n+1)*sum_ss[r];
        // cf. [DT01]
        chi_rho_r[r] += 1./T/L/(n+1)/(n+1)*sum_nn[r];
        chi_sigma_r[r] += 1./T/L/(n+1)/(n+1)*sum_ss[r];
        mean_m[r] = 1./n*sum_m[r];
        // accumulate ms
        if (r > 0) {
            mean_m[r] += mean_m[r-1];
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
    sweep=0;
    M = (uint)(a * init_n_max);
    dublon_rejected = true;
    avg_worm_len = 0;
    worm_len_sample_size = 0;
    N_loop = vtx_visited * M;
    sm.resize(M, identity);

    // read mu value from database if available
    stringstream fname;
    fname << "../mus/" << setprecision(4) << U << "_" << g << "_" << omega
          << ".mu";
    ifstream fstr(fname.str().c_str());
    if (fstr.is_open()) {
        fstr >> mu;
    }

    // set up adjustment of mu if desired
    if (mu_adjust) {
        mus.resize(mu_adjust_N);
        N_mus.resize(mu_adjust_N, 0);

        // calculate mu values for adjustment
        for (uint i = 0; i < mu_adjust_N; i++) {
            mus[i] = (mu-mu_adjust_range)
                     + 2*mu_adjust_range/(mu_adjust_N-1)*i;
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
    double C = (U/4 > -abs(mu)) ? (U/4 + 2*abs(mu)) : (-U/4);
    f << "SIMULATION PROPERTIES" << endl
      << "C+delta+epsilon = " << C << " + " << delta << " + "
      << epsilon << " = " << C+delta+epsilon << endl
      << "operator string max. length: " << M << endl
      << "average worm length: " << avg_worm_len << endl
      << "number of loops per MCS: " << N_loop << endl
      << "mu = " << mu << endl;
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
