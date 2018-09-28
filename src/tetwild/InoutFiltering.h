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
    const std::vector<TetVertex>& tet_vertices;
    const std::vector<std::array<int, 4>>& tets;
    const std::vector<std::array<int, 4>>& is_surface_fs;
    const std::vector<bool>& t_is_removed;

    InoutFiltering(const std::vector<TetVertex>& t_vs, const std::vector<std::array<int, 4>>& ts,
                   const std::vector<std::array<int, 4>>& is_sf_fs,
                   const std::vector<bool>& t_is_rm,
                   const State &st)
        : state(st)
        , tet_vertices(t_vs)
        , tets(ts)
        , is_surface_fs(is_sf_fs)
        , t_is_removed(t_is_rm)
    { }

    void getSurface(Eigen::MatrixXd& V_sf, Eigen::MatrixXi& F_sf);
    std::vector<bool> filter();

    void outputWindingNumberField(const Eigen::VectorXd& W);
};

} // namespace tetwild

#endif //NEW_GTET_INOUTFILTERING_H
