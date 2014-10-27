#include "mc.h"

mc :: mc (string dir) {
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    T = param.value_or_default<double>("T", 1.);
    filling = param.value_or_default<double>("FILLING", .5);
    a = param.value_or_default<double>("A", 1.3);
    init_n_max = param.value_or_default<int>("INIT_N_MAX", 100);
    therm = param.value_or_default<int>("THERMALIZATION", 10000);    
    state.resize(L);
    sm.resize(M, 0);
}

mc :: ~mc() {
    random_clear();
    state.clear();
    sm.clear();
}

void mc :: do_update() {
    for (uint j = 0; j < M; ++j) {
        
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

