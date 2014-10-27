#include "mc.h"

mc :: mc (string dir) {
    param_init(dir);
    L = param.value_or_default<int>("L", 10);
    T = param.value_or_default<double>("T", 1.);
    filling = param.value_or_default<double>("FILLING", .5);
    therm = param.value_or_default<int>("THERMALIZATION", 10000);    
    state.resize(L);
}

mc :: ~mc() {
    random_clear();
    state.clear();    
}

void mc :: do_update() {
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

    sweep=0;
    // add observables
}

void mc :: write(string dir)
{
    odump d(dir + "dump");
    random_write(d);
    d.write(sweep);
    d.write(state);
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

