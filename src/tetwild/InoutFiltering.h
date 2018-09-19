// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by yihu on 6/13/17.
//

#ifndef NEW_GTET_INOUTFILTERING_H
#define NEW_GTET_INOUTFILTERING_H

#include <tetwild/TetmeshElements.h>
#include <Eigen/Dense>

namespace tetwild {

class InoutFiltering {
public:
    const State &state;
    std::vector<TetVertex>& tet_vertices;
    std::vector<std::array<int, 4>>& tets;
    std::vector<std::array<int, 4>>& is_surface_fs;
    std::vector<bool>& v_is_removed;
    std::vector<bool>& t_is_removed;
    std::vector<TetQuality>& tet_qualities;

    std::vector<bool> is_inside;
    InoutFiltering(std::vector<TetVertex>& t_vs, std::vector<std::array<int, 4>>& ts,
                   std::vector<std::array<int, 4>>& is_sf_fs,
                   std::vector<bool>& v_is_rm, std::vector<bool>& t_is_rm, std::vector<TetQuality>& tet_qs,
                   const State &st):
            tet_vertices(t_vs), tets(ts), is_surface_fs(is_sf_fs), v_is_removed(v_is_rm), t_is_removed(t_is_rm),
            tet_qualities(tet_qs), state(st)
    { }

    void getSurface(Eigen::MatrixXd& V_sf, Eigen::MatrixXi& F_sf);
    void filter();

    void outputWindingNumberField(const Eigen::VectorXd& W);
};

} // namespace tetwild

#endif //NEW_GTET_INOUTFILTERING_H
