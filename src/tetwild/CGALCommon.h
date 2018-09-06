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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/centroid.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/bounding_box.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#include <igl/serialize.h>
#include <igl/STR.h>

namespace tetwild {

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Triangle_2 Triangle_2;
typedef K::Intersect_2 Intersect_2;
//typedef CGAL::Polygon_2<K> Polygon_2;

typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Segment_3 Segment_3;
typedef K::Line_3 Line_3;
typedef K::Plane_3 Plane_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Intersect_3 Intersect_3;
typedef K::Tetrahedron_3 Tetrahedron_3;
typedef K::Direction_3 Direction_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kf;
typedef Kf::Point_3 Point_3f;
typedef Kf::Vector_3 Vector_3f;
typedef Kf::Plane_3 Plane_3f;
typedef Kf::Triangle_3 Triangle_3f;
typedef Kf::Segment_3 Segment_3f;
typedef Kf::Line_3 Line_3f;

typedef CGAL::Epeck::FT CGAL_FT;
//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<CGAL::Gmpq>::FT CGAL_FT;

typedef K::Iso_cuboid_3 Bbox_3;

} // namespace tetwild

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
            std::string s1, s2, s3;
            ::igl::deserialize(s1, std::string("x"), buffer);
            ::igl::deserialize(s2, std::string("y"), buffer);
            ::igl::deserialize(s3, std::string("z"), buffer);
//            p=Point_3(CGAL_FT(s1), CGAL_FT(s2), CGAL_FT(s3));
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
    }
}
