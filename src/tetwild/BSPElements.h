// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by yihu on 8/22/17.
//

#ifndef NEW_GTET_BSPELEMENTS_H
#define NEW_GTET_BSPELEMENTS_H
#include <vector>
#include <unordered_set>

namespace tetwild {

class BSPEdge{
public:
    std::vector<int> vertices;
    std::unordered_set<int> conn_faces;

    BSPEdge(){}
    BSPEdge(int v1, int v2){
        vertices={v1, v2};
    }
};

class BSPFace{
public:
    std::vector<int> vertices;
    std::vector<int> edges;
    std::unordered_set<int> conn_nodes;
    std::unordered_set<int> div_faces;

    int matched_f_id=-1;
};

class BSPtreeNode{
public:
    bool is_leaf=false;
    std::vector<int> faces;
    std::unordered_set<int> div_faces;
};

} // namespace tetwild

#endif //NEW_GTET_BSPELEMENTS_H
