#pragma once
#include <string>
#include <map>
#include <vector>
#include <GraphMol/GraphMol.h>
#include "protac.hpp"
struct Parameters {
    std::string protacs_file = "protac.sdf";
    std::string warhead_anchor_file = "e3.sdf";
    std::string warhead_flex_file = "poi.sdf";
    std::string protein_anchor_file = "e3.pdb";
    std::string protein_flex_file = "poi.pdb";
    std::string output_protac_file = "TC_protac.sdf";
    std::string output_protein_file = "";
    double output_rmsd_cutoff = 1.0;
    int n_ini = 10000;
    int n_search = 900;
    int n_keep = 900;
    int n_processes = 1;
    int verbose = 0;
    bool score_only = false;
    bool local_only = false;
    VolRegion interface;
};

void read_parameters(const std::string& filename, Parameters& params);
void run_ternify(const Parameters& params);