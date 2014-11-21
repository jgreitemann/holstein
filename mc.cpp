#include "mc.h"

#define N_BOND 3                // number of bond operator flavors
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
    epsilon = param.value_or_default<double>("EPSILON", -1.);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 10000);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    N_loop = param.value_or_default<double>("N_LOOP", 2.0);
    assert(N_el_up <= L && N_el_down <= L);
    assert(N_el_up % 2 == 1 && N_el_down % 2 == 1);
    assert(t > 0);

    // resize vectors
    state.resize(L);
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
    double a[][7] = {
                        {0, 0, 0, 0, 0, 0, 0}, // b_1
                        {0, 0, 0, 0, 0, 0, 0}, // b_2
                        { // a
                            0.5 * (W[2] + W[4] - W[6]),
                            0.5 * (W[0] + W[3] - W[6]),
                            0.5 * (W[1] + W[4] - W[6]),
                            0.5 * (W[1] + W[3] - W[6]),
                            0.5 * (W[4] + W[5] - W[6]),
                            0.5 * (W[3] + W[5] - W[6]),
                            t
                        },
                        { // b
                            0.5 * (W[2] - W[4] + W[6]),
                            0.5 * (W[0] - W[3] + W[6]),
                            0.5 * (W[1] - W[4] + W[6]),
                            0.5 * (W[1] - W[3] + W[6]),
                            0.5 * (W[4] - W[5] + W[6]),
                            0.5 * (W[3] - W[5] + W[6]),
                            0
                        },
                        { // c
                            0.5 * (-W[2] + W[4] + W[6]),
                            0.5 * (-W[0] + W[3] + W[6]),
                            0.5 * (-W[1] + W[4] + W[6]),
                            0.5 * (-W[1] + W[3] + W[6]),
                            0.5 * (-W[4] + W[5] + W[6]),
                            0.5 * (-W[3] + W[5] + W[6]),
                            0
                        }
                    };
    for (uint i = 0; i < 7; ++i) {
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
    for (uint i = 0; i < 7; ++i) {
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
    for (uint i = 0; i < 7; ++i) {
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
        prob[vtx] = a[j][i] / weight[vtx & 255];
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
    sm.clear();
    weight.clear();
    vtx_type.clear();
    prob.clear();
    ns.clear();
}

void mc :: do_update() {
    vector<int> current_state(state);

    // diagonal update
    for (uint i = 0; i < M; ++i) {
        if (sm[i] == 0) {   // identity operator
            int b = random0N(NB)+1;
            int vtx = current_state[LEFT_SITE(b)]
                      + (current_state[RIGHT_SITE(b)]<<2)
                      + (current_state[LEFT_SITE(b)]<<6)
                      + (current_state[RIGHT_SITE(b)]<<4);
            if (random01() < NB/T*weight[vtx]/(M-n)) {
                sm[i] = N_BOND*b;
                n++;
            }
        } else if (sm[i] % N_BOND == 0) {   // diagonal Hubbard U
            int b = sm[i] / N_BOND;
            int vtx = current_state[LEFT_SITE(b)]
                      + (current_state[RIGHT_SITE(b)]<<2)
                      + (current_state[LEFT_SITE(b)]<<6)
                      + (current_state[RIGHT_SITE(b)]<<4);
            if (random01() < (M-n+1)/(NB/T*weight[vtx])) {
                sm[i] = 0;
                n--;
            }
        } else if (sm[i] % N_BOND == 1) {   // spin up hopping
            int left_up = current_state[LEFT_SITE(sm[i]/N_BOND)] & 1;
            current_state[LEFT_SITE(sm[i]/N_BOND)] =
                (current_state[LEFT_SITE(sm[i]/N_BOND)] & 2)
                + (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 1);
            current_state[RIGHT_SITE(sm[i]/N_BOND)] =
                (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 2)
                + left_up;
        } else if (sm[i] % N_BOND == 2) {   // spin down hopping
            int left_down = (current_state[LEFT_SITE(sm[i]/N_BOND)] & 2);
            current_state[LEFT_SITE(sm[i]/N_BOND)] =
                (current_state[LEFT_SITE(sm[i]/N_BOND)] & 1)
                + (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 2);
            current_state[RIGHT_SITE(sm[i]/N_BOND)] =
                (current_state[RIGHT_SITE(sm[i]/N_BOND)] & 1)
                + left_down;
        }
    }

    // adjust M during thermalization
    if (!is_thermalized() && a*n > M) {
        M = a*n;
        sm.resize(M, 0);
    }
    assert(n <= M); // You might need to increase "a" or the
                    // thermalization time if this assertion fails.

    if (n > 0) {
        // linked list construction
        vector<int> vtx(n, -1);
        vector<int> link(4*n, 0);
        vector<int> first(L, -1);
        vector<int> last(L, -1);
        current_state = state;
        uint p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == 0)
                continue;
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
        for (uint i = 0; i < N_loop*M; ++i) {
            j0 = random0N(N_WORM*4*n);
            worm = j0 / (4*n);
            j0 %= 4*n;
            j = j0;

            // check if dublon worm can start from here
            if (worm == 2) {
                int ent_state = (vtx[j0/4] >> (2*(j0%4))) & 3;
                if ((ent_state & 1) ^ (ent_state >> 1)) {
                    continue; // not an empty or fully occupied site
                }
            }

            for (uint k = 0; ; ++k) {
                if (k == loop_term*M) {
                    do_update();
                    return;
                }
                assert(j/4 < (int)n);
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
        }
        link.clear();
        last.clear();

        // mapping back to operator sequence
        p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == 0)
                continue;
            sm[i] = N_BOND*(sm[i]/N_BOND) + vtx_type[vtx[p]] - 1;
            ++p;
        }

        // updating the state, flipping states randomly on sites not
        // affected by the bond operators
        for (uint s = 0; s < L; ++s) {
            if (first[s] == -1) {
                state[s] = random0N(4);
            } else {
                uint leg = first[s] % 4;
                state[s] = (vtx[first[s]/4] & (3<<(2*leg))) >> (2*leg);
            }
        }
        vtx.clear();
        first.clear();
    } else { // empty operator sequence
        for (uint s = 0; s < L; ++s) {
            state[s] = random0N(4);
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
    sweep=0;
    M = (uint)(a * init_n_max);
    sm.resize(M, 0);

    // add observables
    measure.add_observable("N_up");
    measure.add_observable("N_down");
    measure.add_observable("Energy");
    measure.add_vectorobservable("n_i", L);
    measure.add_vectorobservable("n_1n_i", L);
}

void mc :: write(string dir) {
    odump d(dir + "dump");
    random_write(d);
    d.write(sweep);
    d.write(state);
    d.write(sm);
    d.write(n);
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
        d.read(sm);
        M = sm.size();
        d.read(n);
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
}
