#include "mc.h"

#define N_BOND 8                // number of bond operator flavors
#define NB (L)                  // number of bonds
#define COORD 2                 // coordination number
#define LEFT_SITE(b) (b-1)      // lattice -
#define RIGHT_SITE(b) ((b)%L)   //           geometry
#define LEFT_BOND(s) ((s+L-1)%L+1)
#define RIGHT_BOND(s) (s+1)
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

struct Node {
    int i;
    int Nd;
    int m;
    double r;
};

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

    // diagonal update & subsequence construction
    vector<vector<Node> > subseq(L, vector<Node>());
    vector<int> initial_Nd(L, 0);
    vector<int> current_state(state);
    vector<int> current_occ(occ);
    for (uint i = 0; i < M; ++i) {
        // electronic diagonal update
        if (sm[i] == 0) { // identity operator
            int b = random0N(NB)+1;
            int vtx = current_state[LEFT_SITE(b)]
                      + (current_state[RIGHT_SITE(b)]<<2)
                      + (current_state[LEFT_SITE(b)]<<6)
                      + (current_state[RIGHT_SITE(b)]<<4);
            if (random01() < NB/T*weight[vtx]/(M-n)) {
                sm[i] = N_BOND*b;
                n++;
            }
        } else if (sm[i] % N_BOND == 0) { // diagonal Hubbard U
            int b = sm[i] / N_BOND;
            int vtx = current_state[LEFT_SITE(b)]
                      + (current_state[RIGHT_SITE(b)]<<2)
                      + (current_state[LEFT_SITE(b)]<<6)
                      + (current_state[RIGHT_SITE(b)]<<4);
            if (random01() < (M-n+1)/(NB/T*weight[vtx])) {
                sm[i] = 0;
                n--;
            }
        }

        // phononic diagonal update
        if (sm[i] == 0) { // identity operator
            int b = random0N(NB)+1;
            int cocc = current_occ[LEFT_SITE(b)];
            if (random01() < NB/T*omega*(Np-cocc)/(M-n)) {
                sm[i] = N_BOND*b + 7;
                n++;
            }
        } else if (sm[i] % N_BOND == 7) {
            int b = sm[i] / N_BOND;
            int cocc = current_occ[LEFT_SITE(b)];
            if (random01() < (M-n+1)/(NB/T*omega*(Np-cocc))) {
                sm[i] = 0;
                n--;
            }
        }

        if (sm[i] != 0) {
            int b = sm[i] / N_BOND;

            // append to the phonon subsequences
            if (sm[i] % N_BOND == 7) {
                if (subseq[LEFT_SITE(b)].empty()) {
                    ++initial_Nd[LEFT_SITE(b)];
                } else {
                    ++(subseq[LEFT_SITE(b)].back().Nd);
                }
            } else if (sm[i] % N_BOND != 1 && sm[i] % N_BOND != 2) {
                int site;
                if (sm[i] % N_BOND == 0) {
                    site = random0N(2) ? LEFT_SITE(b) : RIGHT_SITE(b);
                } else if (sm[i] % N_BOND == 3 || sm[i] % N_BOND == 4) {
                    site = (sm[i]%N_BOND==3) ? LEFT_SITE(b)
                                             : RIGHT_SITE(b);
                } else if (sm[i] % N_BOND == 5 || sm[i] % N_BOND == 6) {
                    site = (sm[i]%N_BOND==5) ? LEFT_SITE(b)
                                             : RIGHT_SITE(b);
                }
                int vtx = current_state[LEFT_SITE(b)]
                        + (current_state[RIGHT_SITE(b)]<<2)
                        + (current_state[LEFT_SITE(b)]<<6)
                        + (current_state[RIGHT_SITE(b)]<<4);
                int n_el = (current_state[site] & 1)
                           + (current_state[site] >> 1);
                Node newNode;
                newNode.i = i;
                newNode.Nd = 0;
                newNode.m = current_occ[site];
                newNode.r = weight[vtx]/n_el;
                subseq[site].push_back(newNode);
            }

            // propagation of state
            if (sm[i] % N_BOND == 1) {   // spin up hopping
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
        vector<Node>::iterator i1, i2;
        if (subseq[s].size() < 2)
            continue;
        for (uint i = 0; i < M/L; ++i) {
            i2 = subseq[s].begin() + random0N(subseq[s].size());
            i1 = i2++;
            if (i2 == subseq[s].end()) {
                i2 = subseq[s].begin();
            }
            if (sm[i1->i] % N_BOND == 0 && sm[i2->i] % N_BOND == 0) {
                if (random0N(2)) { // (H_1, H_1) -> (H_5, H_4)
                    double prob = g*g * i1->m / (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m+1)/(Np-i1->m), i1->Nd);
                    if (random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_5^R
                            sm[i1->i] = LEFT_BOND(s)*N_BOND + 6;
                        } else { // H_5^L
                            sm[i1->i] = RIGHT_BOND(s)*N_BOND + 5;
                        }
                        if (type % 2) { // H_4^R
                            sm[i2->i] = LEFT_BOND(s)*N_BOND + 4;
                        } else { // H_4^L
                            sm[i2->i] = RIGHT_BOND(s)*N_BOND + 3;
                        }
                        --(i2->m);
                    }
                } else { // (H_1, H_1) -> (H_4, H_5)
                    double prob = g*g * (i1->m+1) / (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m-1)/(Np-i1->m), i1->Nd);
                    if (i1->m < Np && random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_4^R
                            sm[i1->i] = LEFT_BOND(s)*N_BOND + 4;
                        } else { // H_4^L
                            sm[i1->i] = RIGHT_BOND(s)*N_BOND + 3;
                        }
                        if (type % 2) { // H_5^R
                            sm[i2->i] = LEFT_BOND(s)*N_BOND + 6;
                        } else { // H_5^L
                            sm[i2->i] = RIGHT_BOND(s)*N_BOND + 5;
                        }
                        ++(i2->m);
                    }
                }
            } else if ((sm[i1->i] % N_BOND == 3 || sm[i1->i] % N_BOND == 4)
                    && (sm[i2->i] % N_BOND == 5 || sm[i2->i] % N_BOND == 6)) {
                if (random0N(2)) { // (H_4, H_5) -> (H_5, H_4)
                    double prob = 1. * (i1->m) / (i1->m+1)
                                  * pow(1.*(Np-i1->m+1)/(Np-i1->m-1), i1->Nd);
                    if (random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_5^R
                            sm[i1->i] = LEFT_BOND(s)*N_BOND + 6;
                        } else { // H_5^L
                            sm[i1->i] = RIGHT_BOND(s)*N_BOND + 5;
                        }
                        if (type % 2) { // H_4^R
                            sm[i2->i] = LEFT_BOND(s)*N_BOND + 4;
                        } else { // H_4^L
                            sm[i2->i] = RIGHT_BOND(s)*N_BOND + 3;
                        }
                        i2->m -= 2;
                    }
                } else { // (H_4, H_5) -> (H_1, H_1)
                    double prob = 1./g/g / (i1->m+1) * (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m)/(Np-i1->m-1), i1->Nd);
                    if (random01() < prob) {
                        sm[i1->i] = (sm[i1->i]/N_BOND) * N_BOND;
                        sm[i2->i] = (sm[i2->i]/N_BOND) * N_BOND;
                        --(i2->m);
                    }
                }
            } else if ((sm[i1->i] % N_BOND == 5 || sm[i1->i] % N_BOND == 6)
                    && (sm[i2->i] % N_BOND == 3 || sm[i2->i] % N_BOND == 4)) {
                if (random0N(2)) { // (H_5, H_4) -> (H_4, H_5)
                    double prob = 1. / (i1->m) * (i1->m+1)
                                  * pow(1.*(Np-i1->m-1)/(Np-i1->m+1), i1->Nd);
                    if (i1->m < Np && random01() < prob) {
                        int type = random0N(4);
                        if (type / 2) { // H_4^R
                            sm[i1->i] = LEFT_BOND(s)*N_BOND + 4;
                        } else { // H_4^L
                            sm[i1->i] = RIGHT_BOND(s)*N_BOND + 3;
                        }
                        if (type % 2) { // H_5^R
                            sm[i2->i] = LEFT_BOND(s)*N_BOND + 6;
                        } else { // H_5^L
                            sm[i2->i] = RIGHT_BOND(s)*N_BOND + 5;
                        }
                        i2->m += 2;
                    }
                } else { // (H_5, H_4) -> (H_1, H_1)
                    double prob = 1./g/g / i1->m * (i1->r*i2->r)
                                  * pow(1.*(Np-i1->m)/(Np-i1->m+1), i1->Nd);
                    if (random01() < prob) {
                        sm[i1->i] = (sm[i1->i]/N_BOND) * N_BOND;
                        sm[i2->i] = (sm[i2->i]/N_BOND) * N_BOND;
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
        sm.resize(M, 0);
    }
    assert(n <= M); // You might need to increase "a" or the
                    // thermalization time if this assertion fails.

    // directed loops electron update
    if (n > 0) {
        // linked list construction
        vector<int> vtx(n, -1);
        vector<char> lock(n, 0);
        vector<int> link(4*n, 0);
        vector<int> first(L, -1);
        vector<int> last(L, -1);
        current_state = state;
        uint p = 0;
        for (uint i = 0; p < n; ++i) {
            assert(i < M);
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

            if (sm[i] % N_BOND == 3 || sm[i] % N_BOND == 5) {
                lock[p] = 1; // require at least one electron on left site
            } else if (sm[i] % N_BOND == 4 || sm[i] % N_BOND == 6) {
                lock[p] = 2; // require at least one electron on right site
            } else if (sm[i] % N_BOND == 7) {
                lock[p] = 3;
            }

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
            j0 = random0N(N_WORM*4*n);
            worm = j0 / (4*n);
            j0 %= 4*n;
            j = j0;

            // check if dublon worm can start from here
            int ent_state = (vtx[j0/4] >> (2*(j0%4))) & 3;
            if (worm == 2) {
                dublon_rejected = (ent_state & 1) ^ (ent_state >> 1);
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
                assert(j/4 < (int)n);
                if (lock[j/4] == 3) {
                    exit_leg = (j%4) ^ 3; // continue straight
                } else if (lock[j/4] != 0) {
                    if (lock[j/4] == (((j%4) & 1) ^ ((j%4) >> 1))+1) {
                        if ((((vtx[j/4] >> 2*(j%4)) & 3) ^ (worm+1)) == 0) {
                            exit_leg = j%4; // bounce
                        } else {
                            exit_leg = (j%4) ^ 3; // continue straight
                        }
                    } else {
                        exit_leg = (j%4) ^ 3; // continue straight
                    }
                } else {
                    ent_vtx = (worm << 12) | ((j%4) << 8) | vtx[j/4];
                    r = random01();
                    for (exit_leg = 0; exit_leg < 4; ++exit_leg)
                        if (r < prob[(exit_leg << 10) | ent_vtx])
                            break;
                    assert(exit_leg < 4); // assert that break was called
                }
                // do not count bounces into the worm length
                if (j%4 == exit_leg) {
                    --k;
                }
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
        lock.clear();

        // mapping back to operator sequence
        p = 0;
        for (uint i = 0; p < n; ++i) {
            if (sm[i] == 0)
                continue;
            if (sm[i] % N_BOND <= 2)
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
