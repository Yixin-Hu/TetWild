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

#include <tetwild/DisableWarnings.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <tetwild/EnableWarnings.h>

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
