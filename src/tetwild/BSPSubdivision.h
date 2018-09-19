// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 4/3/17.
//

#ifndef NEW_GTET_BSPSUBDIVISION_H
#define NEW_GTET_BSPSUBDIVISION_H

#include <tetwild/ForwardDecls.h>
#include <tetwild/CGALTypes.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <array>

namespace tetwild {

class BSPSubdivision {
public:
    MeshConformer& MC;
    BSPSubdivision(MeshConformer& mc): MC(mc){}
    void init();

    std::queue<int> processing_n_ids;
    void subdivideBSPNodes();

    const int V_POS=0;
    const int V_NEG=1;
    const int V_ON=2;
    void calVertexSides(const Plane_3& pln, const std::unordered_set<int>& v_ids, const std::vector<Point_3>& vs,
                        std::unordered_map<int, int>& v_sides);

    const int DIVFACE_POS=0;
    const int DIVFACE_NEG=1;
    const int DIVFACE_ON=2;
    const int DIVFACE_CROSS=3;
    int divfaceSide(const Plane_3& pln, const std::array<int, 3>& p_ids, const std::vector<Point_3>& ps);
    void getVertices(BSPFace& face);

};

} // namespace tetwild

#endif //NEW_GTET_BSPSUBDIVISION_H
