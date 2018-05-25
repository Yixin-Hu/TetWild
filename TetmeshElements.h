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

#include "heads.h"

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

    TetVertex(){}
    TetVertex(const Point_3& p) {
        pos = p;
    }

    void printInfo(){
        cout<<"is_on_surface = "<<is_on_surface<<endl;
        cout<<"is_on_bbox = "<<is_on_bbox<<endl;
        cout<<"conn_tets: ";
        for(auto it=conn_tets.begin();it!=conn_tets.end();it++){
            cout<<*it<<" ";
        }
        cout<<endl;
    }

    bool is_locked = false;
    bool is_inside = false;
};

class TetQuality {
public:
    double min_d_angle;
    double max_d_angle;
//    double asp_ratio_2;

    double slim_energy;
    double volume;

    TetQuality(){}
    TetQuality(double d_min, double d_max, double r):
            min_d_angle(d_min), max_d_angle(d_max){}

//    bool operator < (const TetQuality& tq) {
//        if (min_d_angle < tq.min_d_angle)
//            return true;
//        if (max_d_angle > tq.max_d_angle)
//            return true;
//        if (asp_ratio_2 > tq.asp_ratio_2)
//            return true;
//        return false;
//    }

    bool isBetterThan(const TetQuality& tq, int energy_type) {
        if (energy_type == ENERGY_AMIPS || energy_type == ENERGY_DIRICHLET) {
            return slim_energy < tq.slim_energy;
        }
        else if (energy_type == ENERGY_AD) {
            return min_d_angle > tq.min_d_angle && max_d_angle < tq.max_d_angle;
        }
        else
            return false;
    }
    bool isBetterOrEqualThan(const TetQuality& tq, int energy_type) {
        if (energy_type == ENERGY_AMIPS || energy_type == ENERGY_DIRICHLET) {
            return slim_energy <= tq.slim_energy;
        }
        else if (energy_type == ENERGY_AD) {
            return min_d_angle >= tq.min_d_angle && max_d_angle <= tq.max_d_angle;
        }
        else
            return false;
    }
};

///for visualization
class Stage{
public:
    std::vector<TetVertex> tet_vertices;
    std::vector<std::array<int, 4>> tets;
    std::vector<std::array<int, 4>> is_surface_fs;
    std::vector<bool> t_is_removed;
    std::vector<bool> v_is_removed;
    std::vector<TetQuality> tet_qualities;

    std::vector<bool> is_shown;
    double resolution;

    Stage(){}
    Stage(const std::vector<TetVertex>& tet_vs, const std::vector<std::array<int, 4>>& ts,
          const std::vector<std::array<int, 4>>& is_sf_fs,
          const std::vector<bool>& v_is_rd, const std::vector<bool>& t_is_rd, const std::vector<TetQuality>& tet_qs):
            tet_vertices(tet_vs), tets(ts), is_surface_fs(is_sf_fs),
            v_is_removed(v_is_rd), t_is_removed(t_is_rd), tet_qualities(tet_qs) {}

    void serialize(std::string serialize_file){
        igl::serialize(tet_vertices, "tet_vertices", serialize_file, true);
        igl::serialize(tets, "tets", serialize_file);
        igl::serialize(is_surface_fs, "tets", serialize_file);
        igl::serialize(v_is_removed, "v_is_removed", serialize_file);
        igl::serialize(t_is_removed, "t_is_removed", serialize_file);
        igl::serialize(tet_qualities, "tet_qualities", serialize_file);

        igl::serialize(is_shown, "is_shown", serialize_file);
        igl::serialize(resolution, "resolution", serialize_file);
    }

    void deserialize(std::string serialize_file){
        igl::deserialize(tet_vertices, "tet_vertices", serialize_file);
        igl::deserialize(tets, "tets", serialize_file);
        igl::deserialize(is_surface_fs, "tets", serialize_file);
        igl::deserialize(v_is_removed, "v_is_removed", serialize_file);
        igl::deserialize(t_is_removed, "t_is_removed", serialize_file);
        igl::deserialize(tet_qualities, "tet_qualities", serialize_file);

        igl::deserialize(is_shown, "is_shown", serialize_file);
        igl::deserialize(resolution, "resolution", serialize_file);
    }
};

namespace igl {
    namespace serialization {
        template<>
        inline void serialize(const TetVertex &v, std::vector<char> &buffer) {
            ::igl::serialize(v.pos, std::string("pos"), buffer);
            ::igl::serialize(v.posf, std::string("posf"), buffer);

            ::igl::serialize(v.is_rounded, std::string("is_rounded"), buffer);
            ::igl::serialize(v.is_on_surface, std::string("is_on_surface"), buffer);
            ::igl::serialize(v.is_on_bbox, std::string("is_on_bbox"), buffer);
            ::igl::serialize(v.is_on_boundary, std::string("is_on_boundary"), buffer);

            ::igl::serialize(v.adaptive_scale, std::string("adaptive_scale"), buffer);

//            ::igl::serialize(v.on_fixed_vertex, std::string("on_fixed_vertex"), buffer);
//            std::vector<int> tmp;
//            for(auto it=v.on_edge.begin();it!=v.on_edge.end();it++)
//                tmp.push_back(*it);
//            ::igl::serialize(tmp, std::string("on_edge"), buffer);
//            tmp.clear();
//            for(auto it=v.on_face.begin();it!=v.on_face.end();it++)
//                tmp.push_back(*it);
//            ::igl::serialize(tmp, std::string("on_face"), buffer);
//            tmp.clear();
//            for(auto it=v.conn_tets.begin();it!=v.conn_tets.end();it++)
//                tmp.push_back(*it);
//            ::igl::serialize(tmp, std::string("conn_tets"), buffer);
        }

        template<>
        inline void deserialize(TetVertex &v, const std::vector<char> &buffer) {
            ::igl::deserialize(v.pos, std::string("pos"), buffer);
            ::igl::deserialize(v.posf, std::string("posf"), buffer);

            ::igl::deserialize(v.is_rounded, std::string("is_rounded"), buffer);
            ::igl::deserialize(v.is_on_surface, std::string("is_on_surface"), buffer);
            ::igl::deserialize(v.is_on_bbox, std::string("is_on_bbox"), buffer);
            ::igl::deserialize(v.is_on_boundary, std::string("is_on_boundary"), buffer);

            ::igl::deserialize(v.adaptive_scale, std::string("adaptive_scale"), buffer);

//            ::igl::deserialize(v.on_fixed_vertex, std::string("on_fixed_vertex"), buffer);
//            std::vector<int> tmp;
//            ::igl::deserialize(tmp, std::string("on_edge"), buffer);
//            for(int i=0;i<tmp.size();i++)
//                v.on_edge.insert(tmp[i]);
//            tmp.clear();
//            ::igl::deserialize(tmp, std::string("on_face"), buffer);
//            for(int i=0;i<tmp.size();i++)
//                v.on_face.insert(tmp[i]);
//            ::igl::deserialize(tmp, std::string("conn_tets"), buffer);
//            for(int i=0;i<tmp.size();i++)
//                v.conn_tets.insert(tmp[i]);
        }

//        template<>
//        inline void serialize(const TetQuality &q, std::vector<char> &buffer) {
//            ::igl::serialize(q.min_d_angle, std::string("min_d_angle"), buffer);
//            ::igl::serialize(q.max_d_angle, std::string("max_d_angle"), buffer);
//            ::igl::serialize(q.asp_ratio_2, std::string("asp_ratio_2"), buffer);
//            ::igl::serialize(q.slim_energy, std::string("slim_energy"), buffer);
//            ::igl::serialize(q.volume, std::string("volume"), buffer);
//
//        }
//
//        template<>
//        inline void deserialize(TetQuality &q, const std::vector<char> &buffer) {
//            ::igl::deserialize(q.min_d_angle, std::string("min_d_angle"), buffer);
//            ::igl::deserialize(q.max_d_angle, std::string("max_d_angle"), buffer);
//            ::igl::deserialize(q.asp_ratio_2, std::string("asp_ratio_2"), buffer);
//            ::igl::deserialize(q.slim_energy, std::string("slim_energy"), buffer);
//            ::igl::deserialize(q.volume, std::string("volume"), buffer);
//
//        }
    }
}

#endif //NEW_GTET_TETMESHELEMENTS_H
