// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 3/31/17.
//

#ifndef NEW_GTET_SIMPLETETRAHEDRALIZATION_H
#define NEW_GTET_SIMPLETETRAHEDRALIZATION_H

#include <tetwild/MeshConformer.h>
#include <tetwild/TetmeshElements.h>

namespace tetwild {

class SimpleTetrahedralization {
public:
    const State & state;
    MeshConformer& MC;
    std::vector<Point_3> centers;

    //useful infos
//    int centroid_size;
    int m_vertices_size;
//    std::vector<bool> is_visited;

    SimpleTetrahedralization(const State &st, MeshConformer& mc) : state(st), MC(mc) { }

    void tetra(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets);
    void triangulation(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets);

    void labelSurface(const std::vector<int>& m_f_tags, const std::vector<int>& m_e_tags,
                      const std::vector<std::vector<int>>& conn_e4v,
                      std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets,
                      std::vector<std::array<int, 4>>& is_surface_fs);
    void labelBbox(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets);
    void labelBoundary(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets,
                       const std::vector<std::array<int, 4>>& is_surface_fs);

    void constructPlane(int bsp_f_id, Plane_3& pln);
};

} // namespace tetwild

#endif //NEW_GTET_SIMPLETETRAHEDRALIZATION_H
