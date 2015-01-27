#include "mc.h"

inline byte number_of_electrons (el_state s) {
    return (s & 1) + (s >> 1);
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
    epsilon = param.value_or_default<double>("EPSILON", -1.);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 10000);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    vtx_visited = param.value_or_default<double>("VTX_VISITED", 2.0);
    Np = param.value_or_default<int>("N_P", 10);
    assert(N_el_up <= L && N_el_down <= L);
    assert(N_el_up % 2 == 1 && N_el_down % 2 == 1);

    // resize vectors
    state.resize(L);
    occ.resize(L, 0);
    weight.resize(256, 0);
    vtx_type.resize(256, electron_diag);
    prob.resize(N_WORM<<12, 0);
    ns.resize(L);


    // define weights
    double C = (U > -abs(mu)/4) ? (U/4 + 2*abs(mu)) : (-U/4);
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
        assert(W[i] >= -1e-14);
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
    ns.clear();
}

void mc :: do_update() {
    // exiting thermalization stage
    if (sweep == therm) {
        N_loop = (uint)(vtx_visited / avg_worm_len * M);
    }

    // diagonal update & subsequence construction
    vector<vector<subseq_node> > subseq(L, vector<subseq_node>());
    vector<int> initial_Nd(L, 0);
    vector<el_state> current_state(state);
    vector<int> current_occ(occ);
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
    current_state.clear();
    current_occ.clear();

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
    subseq.clear();

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
        vector<vertex> vtx(n, vertex());
        vector<lock_flag> lock(n, unlocked);
        vector<list_position> link(4*n, list_position());
        vector<list_position> first(L, invalid_pos);
        vector<list_position> last(L, invalid_pos);
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
            }

            ++p;
        }
        for (uint s = 0; s < L; ++s) {
            if (last[s] != invalid_pos) {
                link[first[s].int_repr] = last[s];
                link[last[s].int_repr] = first[s];
            }
        }
        current_state.clear();

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
        link.clear();
        last.clear();
        lock.clear();

        // mapping back to operator sequence
        p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == identity)
                continue;
            if (sm[i].type == up_hopping || sm[i].type == down_hopping)
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
        vtx.clear();
        first.clear();
    } else { // empty operator sequence
        for (uint s = 0; s < L; ++s) {
            state[s] = static_cast<el_state>(random0N(4));
        }
    }

    ++sweep;
}

void mc :: do_measurement() {
    uint N_up = 0, N_down = 0;
    for (uint s = 0; s < L; ++s) {
        ns[s] = number_of_electrons(state[s]);
        N_up += state[s] == up || state[s] == dublon;
        N_down += state[s] == down || state[s] == dublon;
    }
    measure.add("N_up", N_up);
    measure.add("N_down", N_down);
    measure.add("dublon_rejection_rate", dublon_rejected);
    
    // skip measurement if particle numbers are not right
    if (N_up != N_el_up || N_down != N_el_down)
        return;

    double C = (U > -abs(mu)/4) ? (U/4 + 2*abs(mu)) : (-U/4);
    double energy = -T * n + NB*(C+epsilon) + L*omega*Np;

    // add data to measurement
    measure.add("Energy", energy);
    measure.add("n_i", ns);

    // calculate and measure density-density correlation
    for (uint s = L; s > 0; --s) {
        ns[s-1] *= ns[0];
    }
    measure.add("n_1n_i", ns);

    measure.add("m_i", occ);
}


bool mc :: is_thermalized() {
    return (sweep>therm);
}

void mc :: init() {
    random_init();

    // initialize states randomly
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

    // add observables
    measure.add_observable("N_up");
    measure.add_observable("N_down");
    measure.add_observable("dublon_rejection_rate");
    measure.add_observable("Energy");
    measure.add_vectorobservable("n_i", L);
    measure.add_vectorobservable("n_1n_i", L);
    measure.add_vectorobservable("m_i", L);
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
    d.close();
    seed_write(dir + "seed");
    dir += "bins";
    ofstream f;
    f.open(dir.c_str());
    f << ( (is_thermalized()) ? sweep-therm : 0 ) << endl;
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
        d.close();
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
    double C = (U > -abs(mu)/4) ? (U/4 + 2*abs(mu)) : (-U/4);
    f << "C+epsilon = " << C+epsilon << endl;
    f << "operator string max. length: " << M << endl;
    f << "average worm length: " << avg_worm_len << endl;
    f << "number of loops per MCS: " << N_loop << endl;
}
