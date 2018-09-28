// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Jeremie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Jeremie Dumas on 09/04/18.
//

#pragma once

#include <tetwild/ForwardDecls.h>
#include <geogram/basic/geometry.h>
#include <unordered_set>
#include <vector>

#define TIMING_BREAKDOWN true

namespace tetwild {

////////////////////////////////////////////////////////////////////////////////

void pausee();

inline bool isHaveCommonEle(const std::unordered_set<int>& v1, const std::unordered_set<int>& v2);
void setIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& s);
void setIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& s);
void sampleTriangle(const std::array<GEO::vec3, 3>& vs, std::vector<GEO::vec3>& ps, double sampling_dist);

void addRecord(const MeshRecord& record, const Args &args, const State &state);

////////////////////////////////////////////////////////////////////////////////

inline bool isHaveCommonEle(const std::unordered_set<int>& v1, const std::unordered_set<int>& v2) {
#if 0
    for (auto it = v1.begin(); it != v1.end(); it++)
        if(std::find(v2.begin(), v2.end(), *it)!=v2.end())
            return true;
#else
    if (v2.size() < v1.size()) {
        return isHaveCommonEle(v2, v1);
    }
    for (int x : v1) {
        if (v2.count(x)) {
            return true;
        }
    }
#endif
    return false;
}

} // namespace tetwild
