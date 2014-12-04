#include "mc.h"

#define N_BOND 8                // number of bond operator flavors
#define NB (L)                  // number of bonds
#define COORD 2                 // coordination number
#define LEFT_SITE(b) (b-1)      // lattice -
#define RIGHT_SITE(b) ((b)%L)   //           geometry
#define N_WORM 3                // number of worm types

#ifdef HEAT_BATH
inline int flipped_vtx(int vtx) {
    int worm = vtx >> 12;
    int ent_leg = (vtx >> 8) & 3;
    int exit_leg = (vtx >> 10) & 3;
    return (vtx & 255) ^ ((worm+1) << (2*ent_leg))
                       ^ ((worm+1) << (2*exit_leg));
}
#endif

mc :: mc (string dir) {
    // initialize job parameters
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    T = param.value_or_default<double>("T", 1.);
    N_el_up = param.value_or_default<int>("N_el_up", L/2);
    N_el_down = param.value_or_default<int>("N_el_down", N_el_up);
    a = param.value_or_default<double>("A", 1.3);
    U = param.value_or_default<double>("U", 1.);
    t = param.value_or_default<double>("HOPPING", 1.);
    mu = param.value_or_default<double>("MU", 0.);
    omega = param.value_or_default<double>("OMEGA", 1.);
    g = param.value_or_default<double>("ELECTRON_PHONON_COUPLING", 0.);
    epsilon = param.value_or_default<double>("EPSILON", -1.);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 10000);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    vtx_visited = param.value_or_default<double>("VTX_VISITED", 2.0);
    Np = param.value_or_default<int>("N_P", 10);
    assert(N_el_up <= L && N_el_down <= L);
    assert(N_el_up % 2 == 1 && N_el_down % 2 == 1);
    assert(t > 0);

    // resize vectors
    state.resize(L);
    occ.resize(L, 0);
    weight.resize(256, 0);
    vtx_type.resize(256, 0);
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
                     t
                 };
    for (uint i = 0; i < 7; ++i) {
        assert(W[i] >= 0);
    }

#ifndef HEAT_BATH
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
#endif

    // determine epsilon
    double epsilon_min = 0.0;
#ifndef HEAT_BATH
    for (uint i = 0; i < 6; ++i) {
        if (a[2][i] < -epsilon_min) {
            epsilon_min = -a[2][i];
        }
    }
#endif
    if (epsilon < 0) {
        epsilon = epsilon_min;
    } else {
        assert(epsilon >= epsilon_min);
    }
#ifndef HEAT_BATH
    for (uint i = 0; i < 6; ++i) {
        a[2][i] += epsilon;
    }
#endif

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
        vtx_type[vtx] = j;
    }
    file1.close();

    // calculate transition probabilities
#ifndef HEAT_BATH
    ifstream file2("../assignments.txt");
    if (!file2.is_open()) {
        cerr << "Could not open file assignments.txt" << endl;
        exit(1);
    }
    while (file2 >> vtx >> i >> j) {
        if (vtx >= (N_WORM << 12)) {
            continue;
        }
        if (i >= 6) {    // 2x2 group
            prob[vtx] = 1.;
        } else {
            prob[vtx] = a[j][i] / weight[vtx & 255];
        }
    }
    file2.close();
#else
    for (i = 0; i < 1024; ++i) {
        for (uint worm = 0; worm < N_WORM; ++worm) {
            double total = 0.;
            for (int vtx = i; vtx < 4096; vtx += 1024) {
                prob[(worm<<12)+vtx]=weight[flipped_vtx((worm<<12)+vtx)];
                total += prob[(worm<<12)+vtx];
            }
            for (int vtx = i; vtx < 4096; vtx += 1024) {
                prob[(worm<<12)+vtx] /= total;
            }
        }
    }
#endif

    // cumulate transition probabilities
    for (i = 0; i < 1024; ++i) {
        for (uint worm = 0; worm < N_WORM; ++worm) {
            prob[(worm<<12)+(1<<10)+i] += prob[(worm<<12)+i];
            prob[(worm<<12)+(2<<10)+i] += prob[(worm<<12)+(1<<10)+i];
            prob[(worm<<12)+(3<<10)+i] += prob[(worm<<12)+(2<<10)+i];
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

    vector<vector<int> > i1(L, vector<int>());
    vector<vector<int> > i2(L, vector<int>());
    vector<vector<int> > Nd(L, vector<int>());
    vector<vector<int> > m(L, vector<int>());
    vector<vector<double> > r(L, vector<double>());

    // diagonal update & subsequence construction
    vector<int> current_state(state);
    vector<int> current_occ(occ);
    for (uint i = 0; i < M; ++i) {
        if (sm[i] == 0) {   // identity operator
            int b = random0N(NB)+1;
            if (random0N(2)) {  // try inserting an H_1
                int vtx = current_state[LEFT_SITE(b)]
                          + (current_state[RIGHT_SITE(b)]<<2)
                          + (current_state[LEFT_SITE(b)]<<6)
                          + (current_state[RIGHT_SITE(b)]<<4);
                if (random01() < NB/T*weight[vtx]/(M-n)) {
                    sm[i] = N_BOND*b;
                    n++;
                    n_Hubb++;
                }
            } else {            // try inserting an H_6
                int cocc = current_occ[LEFT_SITE(b)];
                if (random01() < NB/T*omega*(Np-cocc)/(M-n)) {
                    sm[i] = N_BOND*b + 7;
                    n++;
                }
            }
        } else {
            int b = sm[i] / N_BOND;
            int vtx = current_state[LEFT_SITE(b)]
                      + (current_state[RIGHT_SITE(b)]<<2)
                      + (current_state[LEFT_SITE(b)]<<6)
                      + (current_state[RIGHT_SITE(b)]<<4);
            if (sm[i] % N_BOND == 0) {   // diagonal Hubbard U
                if (random01() < (M-n+1)/(NB/T*weight[vtx])) {
                    sm[i] = 0;
                    n--;
                    n_Hubb--;
                }
            } else if (sm[i] % N_BOND == 7) { // diagonal phonon operator
                int cocc = current_occ[LEFT_SITE(b)];
                if (random01() < (M-n+1)/(NB/T*omega*(Np-cocc))) {
                    sm[i] = 0;
                    n--;
                }
            } else if (sm[i] % N_BOND == 1) {   // spin up hopping
                int left_up = current_state[LEFT_SITE(b)] & 1;
                current_state[LEFT_SITE(b)] =
                    (current_state[LEFT_SITE(b)] & 2)
                    + (current_state[RIGHT_SITE(b)] & 1);
                current_state[RIGHT_SITE(b)] =
                    (current_state[RIGHT_SITE(b)] & 2) + left_up;
            } else if (sm[i] % N_BOND == 2) {   // spin down hopping
                int left_down = (current_state[LEFT_SITE(b)] & 2);
                current_state[LEFT_SITE(b)] =
                    (current_state[LEFT_SITE(b)] & 1)
                    + (current_state[RIGHT_SITE(b)] & 2);
                current_state[RIGHT_SITE(b)] =
                    (current_state[RIGHT_SITE(b)] & 1) + left_down;
            } else if (sm[i] % N_BOND == 3) {
                current_occ[LEFT_SITE(b)]++;
            } else if (sm[i] % N_BOND == 4) {
                current_occ[RIGHT_SITE(b)]++;
            } else if (sm[i] % N_BOND == 5) {
                current_occ[LEFT_SITE(b)]--;
            } else if (sm[i] % N_BOND == 6) {
                current_occ[RIGHT_SITE(b)]--;
            }
            // append to the phonon subsequences
            if (sm[i] % N_BOND == 7) {
                Nd[LEFT_SITE(b)][Nd[LEFT_SITE(b)].size()-1]++;
            } else {
                int site;
                bool close;
                if (sm[i] % N_BOND == 0) {
                    site = random0N(2) ? LEFT_SITE(b) : RIGHT_SITE(b);
                    close = !i1[site].empty() && sm[i1[site].back()] % N_BOND == 0;
                } else if (sm[i] % N_BOND == 3 || sm[i] % N_BOND == 4) {
                    site = (sm[i]%N_BOND==3) ? LEFT_SITE(b)
                                             : RIGHT_SITE(b);
                    close = !i1[site].empty() && (sm[i1[site].back()] % N_BOND == 5
                            || sm[i1[site].back()] % N_BOND == 6);
                } else if (sm[i] % N_BOND == 5 || sm[i] % N_BOND == 6) {
                    site = (sm[i]%N_BOND==5) ? LEFT_SITE(b)
                                             : RIGHT_SITE(b);
                    close = !i1[site].empty() && (sm[i1[site].back()] % N_BOND == 3
                            || sm[i1[site].back()] % N_BOND == 4);
                } else {
                    continue;
                }
                int n_el = (current_state[site] & 1)
                           + (current_state[site] >> 1);
                if (close) {
                    i2[site].push_back(i);
                    r[site][r[site].size()-1] *= weight[vtx]/n_el;
                } else if (!i1[site].empty()){
                    i1[site].pop_back();
                    Nd[site].pop_back();
                    m[site].pop_back();
                    r[site].pop_back();
                }
                i1[site].push_back(i);
                Nd[site].push_back(0);
                m[site].push_back(current_occ[site]);
                r[site].push_back(weight[vtx]/n_el);
            }
        }
    }
    current_state.clear();
    current_occ.clear();

    // adjust M during thermalization
    if (!is_thermalized() && a*n > M) {
        M = a*n;
        sm.resize(M, 0);
    }
    assert(n <= M); // You might need to increase "a" or the
                    // thermalization time if this assertion fails.

    // subsequence phonon update
    i1.clear();
    i2.clear();
    Nd.clear();
    m.clear();
    r.clear();

    // directed loops electron update
    vector<bool> ph_el_op(L, false);
    if (n_Hubb > 0) {
        // linked list construction
        vector<int> vtx(n_Hubb, -1);
        vector<int> link(4*n_Hubb, 0);
        vector<int> first(L, -1);
        vector<int> last(L, -1);
        current_state = state;
        uint p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == 0)
                continue;
            if (sm[i] % N_BOND == 3 || sm[i] % N_BOND == 5) {
                ph_el_op[LEFT_SITE(sm[i]/N_BOND)] = true;
                p++;
                continue;
            }
            if (sm[i] % N_BOND == 4 || sm[i] % N_BOND == 6) {
                ph_el_op[RIGHT_SITE(sm[i]/N_BOND)] = true;
                p++;
                continue;
            }
            if (sm[i] % N_BOND == 7) {
                p++;
                continue;
            }
            // establish links
            if (first[LEFT_SITE(sm[i]/N_BOND)] == -1) {
                first[LEFT_SITE(sm[i]/N_BOND)] = 4 * p;
                last[LEFT_SITE(sm[i]/N_BOND)] = 4 * p + 3;
            } else {
                link[4*p] = last[LEFT_SITE(sm[i]/N_BOND)];
                link[last[LEFT_SITE(sm[i]/N_BOND)]] = 4 * p;
                last[LEFT_SITE(sm[i]/N_BOND)] = 4 * p + 3;
            }
            if (first[RIGHT_SITE(sm[i]/N_BOND)] == -1) {
                first[RIGHT_SITE(sm[i]/N_BOND)] = 4 * p + 1;
                last[RIGHT_SITE(sm[i]/N_BOND)] = 4 * p + 2;
            } else {
                link[4*p+1] = last[RIGHT_SITE(sm[i]/N_BOND)];
                link[last[RIGHT_SITE(sm[i]/N_BOND)]] = 4 * p + 1;
                last[RIGHT_SITE(sm[i]/N_BOND)] = 4 * p + 2;
            }

            // determine vertex type
            vtx[p] = current_state[LEFT_SITE(sm[i]/N_BOND)]
                     + (current_state[RIGHT_SITE(sm[i]/N_BOND)] << 2);
            if (sm[i] % N_BOND == 1) {   // spin up hopping
                int left_up = current_state[LEFT_SITE(sm[i]/N_BOND)] & 1;
                current_state[LEFT_SITE(sm[i]/N_BOND)] =
                    (current_state[LEFT_SITE(sm[i]/N_BOND)] & 2)
                    + (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 1);
                current_state[RIGHT_SITE(sm[i]/N_BOND)] =
                    (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 2)
                    + left_up;
            } else if (sm[i] % N_BOND == 2) {   // spin down hopping
                int left_down =
                    (current_state[LEFT_SITE(sm[i]/N_BOND)] & 2);
                current_state[LEFT_SITE(sm[i]/N_BOND)] =
                    (current_state[LEFT_SITE(sm[i]/N_BOND)] & 1)
                    + (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 2);
                current_state[RIGHT_SITE(sm[i]/N_BOND)] =
                    (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 1)
                    + left_down;
            }
            vtx[p] += (current_state[RIGHT_SITE(sm[i]/N_BOND)] << 4)
                      + (current_state[LEFT_SITE(sm[i]/N_BOND)] << 6);

            ++p;
        }
        for (uint s = 0; s < L; ++s) {
            if (last[s] != -1) {
                link[first[s]] = last[s];
                link[last[s]] = first[s];
            }
        }
        current_state.clear();

        // directed loop construction
        int j, j0, ent_vtx, exit_leg, worm;
        double r;
        for (uint i = 0; i < N_loop; ++i) {
            j0 = random0N(N_WORM*4*n_Hubb);
            worm = j0 / (4*n_Hubb);
            j0 %= 4*n_Hubb;
            j = j0;

            // check if dublon worm can start from here
            if (worm == 2) {
                int ent_state = (vtx[j0/4] >> (2*(j0%4))) & 3;
                dublon_rejected = (ent_state & 1) ^ (ent_state >> 1);
                if (dublon_rejected) {
                    // IMPORTANT: this has to be counted as a loop
                    continue; // not an empty or fully occupied site
                }
            }

            uint k;
            for (k = 0; ; ++k) {
                if (k == loop_term*M) {
                    do_update();
                    return;
                }
                assert(j/4 < (int)n_Hubb);
                ent_vtx = (worm << 12) | ((j%4) << 8) | vtx[j/4];
                r = random01();
                for (exit_leg = 0; exit_leg < 4; ++exit_leg)
                    if (r < prob[(exit_leg << 10) | ent_vtx])
                        break;
                assert(exit_leg < 4); // assert that break was called
                // flip the vertex:
                vtx[j/4] ^= ((worm+1) << 2*(j%4))
                            ^ ((worm+1) << 2*exit_leg);
                j += exit_leg-(j%4); // exit leg position in linked list
                if (j == j0)    // loop closed (SS02, Fig. 4b)
                    break;
                j = link[j];
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

        // mapping back to operator sequence
        p = 0;
        for (uint i = 0; p < n_Hubb; ++i) {
            if (sm[i] == 0 || sm[i] > 2)
                continue;
            sm[i] = N_BOND*(sm[i]/N_BOND) + vtx_type[vtx[p]] - 1;
            ++p;
        }

        // updating the state, flipping states randomly on sites not
        // affected by the bond operators
        for (uint s = 0; s < L; ++s) {
            if (first[s] == -1 && ph_el_op[s] == false) {
                state[s] = random0N(4);
            } else {
                uint leg = first[s] % 4;
                state[s] = (vtx[first[s]/4] & (3<<(2*leg))) >> (2*leg);
            }
        }
        vtx.clear();
        first.clear();
    } else { // no electron operators in the sequence
        uint p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == 0)
                continue;
            if (sm[i] % N_BOND == 3 || sm[i] % N_BOND == 5) {
                ph_el_op[LEFT_SITE(sm[i]/N_BOND)] = true;
                p++;
                continue;
            }
            if (sm[i] % N_BOND == 4 || sm[i] % N_BOND == 6) {
                ph_el_op[RIGHT_SITE(sm[i]/N_BOND)] = true;
                p++;
                continue;
            }
            if (sm[i] % N_BOND == 7) {
                p++;
                continue;
            }
        }
        for (uint s = 0; s < L; ++s) {
            if (ph_el_op[s] == false) {
                state[s] = random0N(4);
            }
        }
    }

    ++sweep;
}

void mc :: do_measurement() {
    uint N_up = 0, N_down = 0;
    for (uint s = 0; s < L; ++s) {
        ns[s] = (state[s] & 1) + ((state[s] & 2) >> 1);
        N_up += state[s] & 1;
        N_down += (state[s] & 2) >> 1;
    }
    measure.add("N_up", N_up);
    measure.add("N_down", N_down);
    measure.add("dublon_rejection_rate", dublon_rejected);
    
    // skip measurement if particle numbers are not right
    if (N_up != N_el_up || N_down != N_el_down)
        return;

    double C = (U > -abs(mu)/4) ? (U/4 + 2*abs(mu)) : (-U/4);
    double energy = -T * n + (C + epsilon)*NB + U/4*COORD*(N_up+N_down)
                    -U*NB/4 - 2*mu*(L-N_up-N_down);

    // add data to measurement
    measure.add("Energy", energy);
    measure.add("n_i", ns);
    
    // calculate and measure density-density correlation
    for (uint s = L; s > 0; --s) {
        ns[s-1] *= ns[0];
    }
    measure.add("n_1n_i", ns);
}


bool mc :: is_thermalized() {
    return (sweep>therm);
}

void mc :: init() {
    random_init();

    // initialize states randomly
    bool place_holes_up = N_el_up > L/2;
    bool place_holes_down = N_el_down > L/2;
    int initial_state = (place_holes_up ? 1 : 0)
                        | (place_holes_down ? 2 : 0);
    for (uint i = 0; i < L; i++) {
        state[i] = initial_state;
    }
    // place down spins
    uint things_to_place = place_holes_down ? (L-N_el_down) : N_el_down;
    while (things_to_place > 0) {
        int site = random0N(L);
        if ((state[site] & 2) >> 1 == place_holes_down) {
            state[site] ^= 2;
            things_to_place--;
        }
    }
    // place up spins
    things_to_place = place_holes_up ? (L-N_el_up) : N_el_up;
    while (things_to_place > 0) {
        int site = random0N(L);
        if ((state[site] & 1) == place_holes_up) {
            state[site] ^= 1;
            things_to_place--;
        }
    }

    n = 0;
    n_Hubb = 0;
    sweep=0;
    M = (uint)(a * init_n_max);
    dublon_rejected = true;
    avg_worm_len = 0;
    worm_len_sample_size = 0;
    N_loop = vtx_visited * M;
    sm.resize(M, 0);

    // add observables
    measure.add_observable("N_up");
    measure.add_observable("N_down");
    measure.add_observable("dublon_rejection_rate");
    measure.add_observable("Energy");
    measure.add_vectorobservable("n_i", L);
    measure.add_vectorobservable("n_1n_i", L);
}

void mc :: write(string dir) {
    odump d(dir + "dump");
    random_write(d);
    d.write(sweep);
    d.write(state);
    d.write(occ);
    d.write(sm);
    d.write(n);
    d.write(n_Hubb);
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
        d.read(n_Hubb);
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
    f << "SIMULATION PROPERTIES" << endl;
    f << "operator string max. length: " << M << endl;
    f << "average worm length: " << avg_worm_len << endl;
    f << "number of loops per MCS: " << N_loop << endl;
    measure.get_statistics(f);
}
