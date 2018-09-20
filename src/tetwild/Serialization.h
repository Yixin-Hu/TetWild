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

#include <tetwild/CGALTypes.h>
#include <tetwild/TetmeshElements.h>
#include <tetwild/DisableWarnings.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#include <igl/serialize.h>
#include <igl/STR.h>
#include <tetwild/EnableWarnings.h>

//for serialization
namespace igl {
    namespace serialization {
        template<>
        inline void serialize(const tetwild::Point_3 &p, std::vector<char> &buffer) {
            ::igl::serialize(STR(CGAL::exact(p[0])), std::string("x"), buffer);
            ::igl::serialize(STR(CGAL::exact(p[1])), std::string("y"), buffer);
            ::igl::serialize(STR(CGAL::exact(p[2])), std::string("z"), buffer);
        }

        template<>
        inline void deserialize(tetwild::Point_3 &p, const std::vector<char> &buffer) {
            using namespace tetwild;
            std::string s1, s2, s3;
            ::igl::deserialize(s1, std::string("x"), buffer);
            ::igl::deserialize(s2, std::string("y"), buffer);
            ::igl::deserialize(s3, std::string("z"), buffer);
            p=Point_3(CGAL_FT(s1), CGAL_FT(s2), CGAL_FT(s3));
        }

        template<>
        inline void serialize(const tetwild::Point_3f &p, std::vector<char> &buffer) {
            ::igl::serialize(p[0], std::string("x"), buffer);
            ::igl::serialize(p[1], std::string("y"), buffer);
            ::igl::serialize(p[2], std::string("z"), buffer);
        }

        template<>
        inline void deserialize(tetwild::Point_3f &p, const std::vector<char> &buffer) {
            double x, y, z;
            ::igl::deserialize(x, std::string("x"), buffer);
            ::igl::deserialize(y, std::string("y"), buffer);
            ::igl::deserialize(z, std::string("z"), buffer);
            p=tetwild::Point_3f(x, y, z);
        }

        template<>
        inline void serialize(const std::array<int, 3> &arr, std::vector<char> &buffer) {
            for(int i=0;i<3;i++)
                ::igl::serialize(arr[i], std::to_string(i), buffer);
        }

        template<>
        inline void deserialize(std::array<int, 3> &arr, const std::vector<char> &buffer) {
            for(int i=0;i<3;i++)
                ::igl::deserialize(arr[i], std::to_string(i), buffer);
        }

        template<>
        inline void serialize(const std::array<int, 4> &arr, std::vector<char> &buffer) {
            for(int i=0;i<4;i++)
                ::igl::serialize(arr[i], std::to_string(i), buffer);
        }

        template<>
        inline void deserialize(std::array<int, 4> &arr, const std::vector<char> &buffer) {
            for(int i=0;i<4;i++)
                ::igl::deserialize(arr[i], std::to_string(i), buffer);
        }
        template<>
        inline void serialize(const tetwild::TetVertex &v, std::vector<char> &buffer) {
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
        inline void deserialize(tetwild::TetVertex &v, const std::vector<char> &buffer) {
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
