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

namespace tetwild {

    void tetrahedralization(const std::vector<std::array<double, 3>>& V_in,
                            const std::vector<std::array<int, 3>>& F_in,
                            std::vector<std::array<double, 3>>& V_out,
                            std::vector<std::array<int, 4>>& T_out);
}

#endif //TETWILD_TETWILD_H
