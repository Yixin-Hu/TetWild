// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 5/31/18.
//
#ifndef TETWILD_TETWILD_H
#define TETWILD_TETWILD_H

#include <array>
#include <vector>
#include <string>

#ifdef WIN32
	#if defined TETWILD_EXPORTS
		#define TETWILD_EXPORT __declspec( dllexport )
	#else
		#define TETWILD_EXPORT __declspec( dllimport )
	#endif
#else
	#define TETWILD_EXPORT
#endif

namespace tetwild {
    struct Args{
        double i_ideal_edge_length = 20;
        double i_epsilon = 1000;
        int stage = 1;
        double filter_energy = 10;
        int max_pass = 80;

        bool is_laplacian = false;
        int targeted_num_v = -1;
        std::string bg_mesh = "";
        bool is_quiet = false;
    };
    extern Args parameters;

    void TETWILD_EXPORT tetrahedralization(const std::vector<std::array<double, 3>>& V_in,
                                           const std::vector<std::array<int, 3>>& F_in,
                                           std::vector<std::array<double, 3>>& V_out,
                                           std::vector<std::array<int, 4>>& T_out);
}

#endif //TETWILD_TETWILD_H
