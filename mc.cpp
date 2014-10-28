#include "mc.h"

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
    state.resize(L);
    M = (int)(a * init_n_max);
    sm.resize(M, 0);
}

mc :: ~mc() {
    random_clear();
    state.clear();
    sm.clear();
}

void mc :: do_update() {
    vector<int> current_state = state;

    // diagonal update
    for (uint p = 0; p < M; ++p) {
        if (sm[p] == 0) {   // identity operator
            int b = random0N(L-1)+1;
            int vtx = current_state[b-1] << 2 + current_state[b];
            if (random01() < (L-1)/T*weights[vtx]/(M-n)) {
                sm[p] = N_BOND*b;
                n++;
            }
        } else if (sm[p] % N_BOND == 0) {   // diagonal Hubbard U
            int b = sm[p] / N_BOND;
            int vtx = current_state[b-1] << 2 + current_state[b];
            if (random01() < (M-n+1)/((L-1)/T*weights[vtx])) {
                sm[p] = 0;
                n--;
            }
        } else if (sm[p] % N_BOND == 1) {   // spin up hopping
            int left_up = current_state[sm[p]/N_BOND-1] | 1;
            current_state[sm[p]/N_BOND-1] =
                current_state[sm[p]/N_BOND-1] | 2
                + current_state[sm[p]/N_BOND] | 1;
            current_state[sm[p]/N_BOND] =
                current_state[sm[p]/N_BOND] | 2
                + left_up;
        } else if (sm[p] % N_BOND == 2) {   // spin down hopping
            int left_down = current_state[sm[p]/N_BOND-1] | 2;
            current_state[sm[p]/N_BOND-1] =
                current_state[sm[p]/N_BOND-1] | 1
                + current_state[sm[p]/N_BOND] | 2;
            current_state[sm[p]/N_BOND] =
                current_state[sm[p]/N_BOND] | 1
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
            int left_up = current_state[sm[i]/N_BOND-1] | 1;
            current_state[sm[i]/N_BOND-1] =
                current_state[sm[i]/N_BOND-1] | 2
                + current_state[sm[i]/N_BOND] | 1;
            current_state[sm[i]/N_BOND] =
                current_state[sm[i]/N_BOND] | 2
                + left_up;
        } else if (sm[i] % N_BOND == 2) {   // spin down hopping
            int left_down = current_state[sm[i]/N_BOND-1] | 2;
            current_state[sm[i]/N_BOND-1] =
                current_state[sm[i]/N_BOND-1] | 1
                + current_state[sm[i]/N_BOND] | 2;
            current_state[sm[i]/N_BOND] =
                current_state[sm[i]/N_BOND] | 1
                + left_down;
        }
        vtx[p] += (current_state[sm[i]/N_BOND-1] << 4)
                  + (current_state[sm[i]/N_BOND] << 6);

        ++p;
    }
    for (uint s = 0; s < L; ++s) {
        link[first[s]] = last[s];
        link[last[s]] = first[s];
    }

    // directed loop construction
    

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
    uint thing_to_place = (uint)(filling * L);
    while (things_to_place > 0) {
        int site = random0N(L);
        if (state[i] / 2 == place_holes) {
            state[i] ^= 2;
            things_to_place--;
        }
    }
    // place up spins
    thing_to_place = (uint)(filling * L);
    while (things_to_place > 0) {
        int site = random0N(L);
        if (state[i] % 2 == place_holes) {
            state[i] ^= 1;
            things_to_place--;
        }
    }
    
    // calculate vertex weights
    for (int i = 0; i < 16; i++) {
        weights[i] = epsilon + U/4 * (i|1 + (i|2)>>1 + (i|4)>>2 + (i|8)>>3)
                     - U/2 * (((i|3) == 3) + ((i|12) == 12));
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
        d.read(n)
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

