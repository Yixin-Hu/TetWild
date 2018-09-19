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

#ifndef NEW_GTET_TETMESHELEMENTS_H
#define NEW_GTET_TETMESHELEMENTS_H

#include <tetwild/State.h>
#include <tetwild/CGALTypes.h>
#include <unordered_set>

namespace tetwild {

//const int ON_SURFACE_FALSE = 0;//delete
//const int ON_SURFACE_TRUE_INSIDE = 1;//delete
//const int ON_SURFACE_TRUE_OUTSIDE = 2;//delete
class TetVertex {
public:
    Point_3 pos;//todo: how to remove it?

    ///for surface conforming
    int on_fixed_vertex = -1;
    std::unordered_set<int> on_edge;//fixed points can be on more than one edges
    std::unordered_set<int> on_face;
    bool is_on_surface = false;

    ///for local operations
    std::unordered_set<int> conn_tets;

    ///for hybrid rationals
    Point_3f posf;
    bool is_rounded = false;

    void round() {
        posf = Point_3f(CGAL::to_double(pos[0]), CGAL::to_double(pos[1]), CGAL::to_double(pos[2]));
    }

    ///for bbox
    bool is_on_bbox = false;

    ///for boundary
    bool is_on_boundary = false;

    //for adaptive refinement
    double adaptive_scale = 1.0;

    TetVertex() = default;

    TetVertex(const Point_3& p) {
        pos = p;
    }

    void printInfo() const;

    bool is_locked = false;
    bool is_inside = false;
};

class TetQuality {
public:
    double min_d_angle = 0;
    double max_d_angle = 0;
//    double asp_ratio_2;

    double slim_energy = 0;
    double volume = 0;

    TetQuality() = default;
    TetQuality(double d_min, double d_max, double r)
        : min_d_angle(d_min), max_d_angle(d_max)
    { }

//    bool operator < (const TetQuality& tq) {
//        if (min_d_angle < tq.min_d_angle)
//            return true;
//        if (max_d_angle > tq.max_d_angle)
//            return true;
//        if (asp_ratio_2 > tq.asp_ratio_2)
//            return true;
//        return false;
//    }

    bool isBetterThan(const TetQuality& tq, int energy_type, const State &state) {
        if (energy_type == state.ENERGY_AMIPS || energy_type == state.ENERGY_DIRICHLET) {
            return slim_energy < tq.slim_energy;
        }
        else if (energy_type == state.ENERGY_AD) {
            return min_d_angle > tq.min_d_angle && max_d_angle < tq.max_d_angle;
        }
        else
            return false;
    }

    bool isBetterOrEqualThan(const TetQuality& tq, int energy_type, const State &state) {
        if (energy_type == state.ENERGY_AMIPS || energy_type == state.ENERGY_DIRICHLET) {
            return slim_energy <= tq.slim_energy;
        }
        else if (energy_type == state.ENERGY_AD) {
            return min_d_angle >= tq.min_d_angle && max_d_angle <= tq.max_d_angle;
        }
        else
            return false;
    }
};

///for visualization
class Stage {
public:
    std::vector<TetVertex> tet_vertices;
    std::vector<std::array<int, 4>> tets;
    std::vector<std::array<int, 4>> is_surface_fs;
    std::vector<bool> t_is_removed;
    std::vector<bool> v_is_removed;
    std::vector<TetQuality> tet_qualities;

    std::vector<bool> is_shown;
    double resolution;

    Stage() = default;
    Stage(const std::vector<TetVertex>& tet_vs, const std::vector<std::array<int, 4>>& ts,
          const std::vector<std::array<int, 4>>& is_sf_fs,
          const std::vector<bool>& v_is_rd, const std::vector<bool>& t_is_rd, const std::vector<TetQuality>& tet_qs)
        : tet_vertices(tet_vs), tets(ts), is_surface_fs(is_sf_fs)
        , v_is_removed(v_is_rd), t_is_removed(t_is_rd), tet_qualities(tet_qs)
    { }

    void serialize(std::string serialize_file);
    void deserialize(std::string serialize_file);
};

} // namespace tetwild

#endif //NEW_GTET_TETMESHELEMENTS_H
