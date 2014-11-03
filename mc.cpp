#include "mc.h"
#include <assert.h>

#define N_BOND 3

double epsilon_min(double U, double t) {
    if (t > abs(U)/4)   // region I
        return t/2 + abs(U)/8 - U/4;
    else if (U > 0)     // region II
        return 0;
    else                // region III
        return -U/2;
}

mc :: mc (string dir) {
    // initialize job parameters
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    T = param.value_or_default<double>("T", 1.);
    filling = param.value_or_default<double>("FILLING", .5);
    a = param.value_or_default<double>("A", 1.3);
    U = param.value_or_default<double>("U", 1.);
    t = param.value_or_default<double>("HOPPING", 1.);
    epsilon = param.value_or_default<double>("EPSILON", epsilon_min(U,t));
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 10000);
    loop_term = param.value_or_default<int>("LOOP_TERMINATION", 100);
    M = (uint)(a * init_n_max);
    
    // resize vectors
    state.resize(L);
    sm.resize(M, 0);
    weight.resize(256, 0);
    vtx_type.resize(256, 0);
    prob.resize(8192, 0);

    // parse vertex weights
    int vtx, j, i;
    double W[] = {epsilon, U/4+epsilon, U/2+epsilon, t};
    ifstream file1("../vertex_types.txt");
    while (file1 >> vtx >> j >> i) {
        weight[vtx] = W[i-1];
        vtx_type[vtx] = j;
    }
    file1.close();

    // calculate transition probabilities
    double b1 = (t < -U/4) ? (-U/4-t) : 0;
    double b2 = (t < +U/4) ? (+U/4-t) : 0;
    double a[] = {
                    b1,
                    b2,
                    U/8+epsilon-t/2-b1/2-b2/2,
                    -U/8+t/2-b1/2+b2/2,
                    U/8+t/2+b1/2-b2/2,
                    3*U/8+epsilon-t/2-b1/2-b2/2,
                    -U/8+t/2-b1/2+b2/2,
                    U/8+t/2+b1/2-b2/2,
                    t
                 };
    ifstream file2("../assignments.txt");
    while (file2 >> vtx >> j) {
        prob[vtx] = a[j] / weight[vtx & 255];
    }
    file2.close();

    // cumulate transition probabilities
    for (i = 0; i < 2048; ++i) {
        prob[(1<<11)+i] += prob[i];
        prob[(2<<11)+i] += prob[(1<<11)+i];
        prob[(3<<11)+i] += prob[(2<<11)+i];
    }
}

mc :: ~mc() {
    random_clear();
    state.clear();
    sm.clear();
    weight.clear();
    vtx_type.clear();
    prob.clear();
}

void mc :: do_update() {
    vector<int> current_state(state);

    // diagonal update
    for (uint i = 0; i < M; ++i) {
        if (sm[i] == 0) {   // identity operator
            int b = random0N(L-1)+1;
            int vtx = current_state[b-1] + (current_state[b]<<2)
                      + (current_state[b-1]<<6) + (current_state[b]<<4);
            if (random01() < (L-1)/T*weight[vtx]/(M-n)) {
                sm[i] = N_BOND*b;
                n++;
            }
        } else if (sm[i] % N_BOND == 0) {   // diagonal Hubbard U
            int b = sm[i] / N_BOND;
            int vtx = current_state[b-1] + (current_state[b]<<2)
                      + (current_state[b-1]<<6) + (current_state[b]<<4);
            if (random01() < (M-n+1)/((L-1)/T*weight[vtx])) {
                sm[i] = 0;
                n--;
            }
        } else if (sm[i] % N_BOND == 1) {   // spin up hopping
            int left_up = current_state[sm[i]/N_BOND-1] & 1;
            current_state[sm[i]/N_BOND-1] =
                (current_state[sm[i]/N_BOND-1] & 2)
                + (current_state[sm[i]/N_BOND] & 1);
            current_state[sm[i]/N_BOND] =
                (current_state[sm[i]/N_BOND] & 2)
                + left_up;
        } else if (sm[i] % N_BOND == 2) {   // spin down hopping
            int left_down = (current_state[sm[i]/N_BOND-1] & 2);
            current_state[sm[i]/N_BOND-1] =
                (current_state[sm[i]/N_BOND-1] & 1)
                + (current_state[sm[i]/N_BOND] & 2);
            current_state[sm[i]/N_BOND] =
                (current_state[sm[i]/N_BOND] & 1)
                + left_down;
        }
    }
    
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
        if (first[sm[i]/N_BOND-1] == -1) {
            first[sm[i]/N_BOND-1] = 4 * p;
            last[sm[i]/N_BOND-1] = 4 * p + 3;
        } else {
            link[4*p] = last[sm[i]/N_BOND-1];
            link[last[sm[i]/N_BOND-1]] = 4 * p;
            last[sm[i]/N_BOND-1] = 4 * p + 3;
        }
        if (first[sm[i]/N_BOND] == -1) {
            first[sm[i]/N_BOND] = 4 * p + 1;
            last[sm[i]/N_BOND] = 4 * p + 2;
        } else {
            link[4*p+1] = last[sm[i]/N_BOND];
            link[last[sm[i]/N_BOND]] = 4 * p + 1;
            last[sm[i]/N_BOND] = 4 * p + 2;
        }

        // determine vertex type
        vtx[p] = current_state[sm[i]/N_BOND-1]
                 + (current_state[sm[i]/N_BOND] << 2);
        if (sm[i] % N_BOND == 1) {   // spin up hopping
            int left_up = current_state[sm[i]/N_BOND-1] & 1;
            current_state[sm[i]/N_BOND-1] =
                (current_state[sm[i]/N_BOND-1] & 2)
                + (current_state[sm[i]/N_BOND] & 1);
            current_state[sm[i]/N_BOND] =
                (current_state[sm[i]/N_BOND] & 2)
                + left_up;
        } else if (sm[i] % N_BOND == 2) {   // spin down hopping
            int left_down = (current_state[sm[i]/N_BOND-1] & 2);
            current_state[sm[i]/N_BOND-1] =
                (current_state[sm[i]/N_BOND-1] & 1)
                + (current_state[sm[i]/N_BOND] & 2);
            current_state[sm[i]/N_BOND] =
                (current_state[sm[i]/N_BOND] & 1)
                + left_down;
        }
        vtx[p] += (current_state[sm[i]/N_BOND] << 4)
                  + (current_state[sm[i]/N_BOND-1] << 6);

        ++p;
    }
    for (uint s = 0; s < L; ++s) {
        link[first[s]] = last[s];
        link[last[s]] = first[s];
    }

    // directed loop construction
    int j, j0, ent_vtx, exit_leg;
    bool right_flag;
    double r;
    for (uint i = 0; i < 2*M; ++i) {
        j0 = random0N(8*n);
        right_flag = j0 & 1;
        j0 >>= 1;
        j = j0;
        for (uint k = 0; ; ++k) {
            if (k == loop_term*M) {
                do_update();
                return;
            }
            ent_vtx = ((j%4) << 9) + (right_flag << 8) + vtx[j/4];
            r = random01();
            for (exit_leg = 0; exit_leg < 4; ++exit_leg)
                if (r < prob[(exit_leg << 11) + ent_vtx])
                    break;
            assert(r < prob[(exit_leg << 11) + ent_vtx]);
            // flip the vertex:
            vtx[j/4] ^= ((2-right_flag) << 2*(j%4))
                        ^ ((2-right_flag) << 2*exit_leg);
            j += exit_leg - (j%4);   // exit leg position in linked list
            if (j == j0)    // loop closed (SS02, Fig. 4b)
                break;
            j = link[j];
            if (j == j0)    // loop closed (SS02, Fig. 4a)
                break;
        }
    }

    // mapping back to operator sequence
    p = 0;
    for (uint i = 0; p < n; ++i) {
        if (sm[i] == 0)
            continue;
        sm[i] = N_BOND*(sm[i]/N_BOND) + vtx_type[vtx[p]] - 1;
        ++p;
    }

    // updating the state, flipping states randomly on sites not affected
    // by the bond operators
    for (uint s = 0; s < L; ++s) {
        if (first[s] == -1) {
            uint leg = first[4] % 4;
            state[s] = (vtx[first[s]/4] & (3<<(2*leg))) >> (2*leg);
        } else {
            state[s] = random0N(4);
        }
    }

    ++sweep;
}

void mc :: do_measurement() {
    // add data to measurement
}


bool mc :: is_thermalized() {
    return (sweep>therm);
}

void mc :: init() {
    random_init();

    // initialize states randomly
    bool place_holes = (filling > 0.51);
    int initial_state = place_holes ? 3 : 0;
    for (uint i = 0; i < L; i++) {
        state[i] = initial_state;
    }
    // place down spins
    uint things_to_place = (uint)(filling * L);
    while (things_to_place > 0) {
        int site = random0N(L);
        if (state[site] / 2 == place_holes) {
            state[site] ^= 2;
            things_to_place--;
        }
    }
    // place up spins
    things_to_place = (uint)(filling * L);
    while (things_to_place > 0) {
        int site = random0N(L);
        if (state[site] % 2 == place_holes) {
            state[site] ^= 1;
            things_to_place--;
        }
    }

    n = 0;
    sweep=0;
    // add observables
}

void mc :: write(string dir)
{
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


bool mc :: read(string dir)
{
    idump d(dir+"dump");
    if (!d)
        return false;
    else {
        random_read(d);
        d.read(sweep);
        d.read(state);
        d.read(sm);
        d.read(n);
        d.close();
        return true;
    }
}

void mc :: write_output(string dir)
{
    // add evalables
    ofstream f;
    f.open(dir.c_str());
    f << "PARAMETERS" << endl;
    param.get_all(f);
    measure.get_statistics(f);
}

