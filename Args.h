// This file is part of TetWild, a software for generating tetrahedral meshes.
// 
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 3/29/17.
//

#ifndef ARGS_H
#define ARGS_H

#include <string>

namespace tetwild {

struct Args{
    std::string input;
    std::string output = "";
    std::string postfix = "_";
    double i_ideal_edge_length = 20;
    double i_epsilon = 1000;
    int i_dd = -1;
    int stage = 1;
    double adaptive_scalar = 0.6;
    double filter_energy = 10;
    double delta_energy = 0.1;
    int max_pass = 80;
    int is_output_csv = true;
    std::string csv_file = "";
    std::string slz_file = "";

    int mid_result = -1;
    bool is_using_voxel = true;
    bool is_laplacian = false;

    int targeted_num_v = -1;
    std::string bg_mesh = "";

    bool is_quiet = false;
    std::string output_mesh_format = "";
};

}

#endif
