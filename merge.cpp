#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mc.h"
#include "measurements.h"
#include "dump.h"
#include "parser.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc==1) {
        cout << " usage: "
             << argv[0]
             << " jobfilename [-s number_of_bins_to_be_skipped] "
                "[-t min_task max_task] [-r min_run max_run]"
             << endl;
        exit(1);
    }
    std::string jobfile(argv[1]);
    int task_min=1;
    int task_max=-1;
    int run_min=1;
    int run_max=-1;
    int skipbins=0;
    int rarg=2;
    while (argc>rarg) {
        if (argv[rarg][0] == 's' || argv[rarg][1] == 's') {
            rarg++;
            skipbins=atoi(argv[rarg]);
            rarg++;
        } else if (argv[rarg][0] == 't' || argv[rarg][1] == 't') {
            rarg++;
            task_min=atoi(argv[rarg]);
            rarg++;
            task_max=atoi(argv[rarg]);
            rarg++;
        } else if (argv[rarg][0] == 'r' || argv[rarg][1] == 'r') {
            rarg++;
            run_min=atoi(argv[rarg]);
            rarg++;
            run_max=atoi(argv[rarg]);
            rarg++;
        } else {
            cout << "unknown option " << argv[rarg] << endl;
            exit(1);
        }
    }
    cout << "Merging " << jobfile;
    if (task_max == -1) {
        cout << " all task";
    } else {
        cout <<" task " <<task_min << " to " << task_max;
    }
    if (run_max == -1) {
        cout << " all runs";
    } else {
        cout << " run " << run_min << " to " << run_max;
    }
    if (skipbins)
        cout << " skipping the first " << skipbins << " bins";
    cout << endl;

    parser parsedfile(jobfile + ".alltasks");
    std::vector<string> taskfiles;
    taskfiles = parsedfile.return_vector<string>("@taskfiles");
    std::string masterfile = parsedfile.value_or_default<string>("masterfile",
            jobfile + ".master");
    if (task_max == -1)
        task_max = taskfiles.size();
    for (int i = task_min - 1; i < task_max; ++i) {
        std::string taskfile = taskfiles[i];
        parser cfg(taskfile);
        std::string taskdir = cfg.value_of("taskdir");
        std::stringstream rb;
        rb << taskdir << "/run" << run_min << ".";
        std::string rundir = rb.str();
        mc* sys = new mc(taskfile);
        if ((*sys).read(rundir)) {
            cout << rundir; //<< endl;
            if ((*sys).measure.merge(rundir,skipbins))
                cout << " done." << endl;
            if (1) {
                int run_counter = run_min + 1;
                bool success = true;
                while (success && ((run_max == -1) || (run_counter <= run_max))) {
                    stringstream b;
                    b << taskdir << "/run" << run_counter << ".";
                    cout << b.str(); //<<endl;
                    success = (*sys).measure.merge(b.str(), skipbins);
                    if (success)
                        cout <<" done."<<endl;
                    else
                        cout <<" not found."<<endl;
                    ++run_counter;
                }
                std::stringstream mfb;
                if (run_max == -1) {
                    mfb << taskdir << ".";
                } else {
                    if (run_min == run_max) {
                        mfb << taskdir << "." << run_min << ".";
                    } else {
                        mfb << taskdir << "." <<run_min << "." << run_max << ".";
                    }
                }
                if (skipbins)
                    mfb << "s" << skipbins << ".";
                mfb << "out";
                (*sys).write_output(mfb.str());
            }
        }
        delete sys;
    }
    return 0;
}
