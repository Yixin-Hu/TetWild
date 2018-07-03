// This file is part of TetWild, a software for generating tetrahedral meshes.
// 
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 3/29/17.
//

#ifndef GTET_HEADS_H
#define GTET_HEADS_H

#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <vector>
#include <array>
#include <queue>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <cmath>

#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/STR.h>
#include <igl/serialize.h>
#include <igl/Timer.h>
#include <igl/unique_rows.h>

// PyMesh
#include "MshSaver.h"
#include "MshLoader.h"

using std::cerr;
using std::cout;
using std::cin;
using std::endl;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <igl/copyleft/cgal/assign_scalar.h>

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

#include <CGAL/convex_hull_2.h>
#include <CGAL/centroid.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/bounding_box.h>
typedef K::Iso_cuboid_3 Bbox_3;

#include <geogram/mesh/mesh_AABB.h>

void pausee();
bool isHaveCommonEle(const std::unordered_set<int>& v1, const std::unordered_set<int>& v2);
void setIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& s);
void setIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& s);
void sampleTriangle(const std::array<GEO::vec3, 3>& vs, std::vector<GEO::vec3>& ps);

#define TIMING_BREAKDOWN true

const int EPSILON_INFINITE=-2;
const int EPSILON_NA=-1;

const int ENERGY_NA=0;
const int ENERGY_AD=1;
const int ENERGY_AMIPS=2;
const int ENERGY_DIRICHLET=3;

const double MAX_ENERGY = 1e50;

extern int NOT_SURFACE;

//global parameters
extern std::string g_working_dir;
extern std::string g_stat_file;
extern std::string g_postfix;
extern std::string g_output_file;

struct Args{
    std::string input;
    std::string output = "";
    std::string postfix = "_";
    double i_ideal_edge_length = 20;
    double i_epsilon = 1000;
    int i_dd = -1;
    int stage = 1;
    double adaptive_scalar = 0.6;
    double filter_energy = 10;
    double delta_energy = 0.1;
    int max_pass = 80;
    int is_output_csv = true;
    std::string csv_file = "";
    std::string slz_file = "";

    int mid_result = -1;
    bool is_using_voxel = true;
    bool is_laplacian = false;

    int targeted_num_v = -1;
    std::string bg_mesh = "";

    bool is_quiet = false;
};

extern Args args;

extern double g_eps;
extern double g_eps_2;
extern double g_dd;
extern double g_ideal_l;
extern double g_diag_l;
extern bool g_is_close;

extern double g_eps_input;
extern double g_eps_delta;
extern int g_cur_stage;

//for test
extern bool is_using_energy_max;
extern bool is_using_sampling;
extern bool is_use_project;
extern int cnt_geo_aabb;

extern bool is_print_tmp;

struct MeshRecord {
    enum OpType {
        OP_INIT = 0,
        OP_PREPROCESSING,
        OP_DELAUNEY_TETRA,
        OP_DIVFACE_MATCH,
        OP_BSP,
        OP_SIMPLE_TETRA,

        OP_OPT_INIT,
        OP_SPLIT,
        OP_COLLAPSE,
        OP_SWAP,
        OP_SMOOTH,
        OP_ADAP_UPDATE,
        OP_WN,
        OP_UNROUNDED
    };
    int op;
    double timing;
    int n_v;
    int n_t;
    double min_min_d_angle = -1;
    double avg_min_d_angle = -1;
    double max_max_d_angle = -1;
    double avg_max_d_angle = -1;
    double max_energy = -1;
    double avg_energy = -1;
    MeshRecord(int op, double timing, int n_v, int n_t, double min_min_d_angle, double avg_min_d_angle,
               double max_max_d_angle, double avg_max_d_angle, double max_energy, double avg_energy) {
        this->op = op;
        this->timing = timing;
        this->n_v = n_v;
        this->n_t = n_t;
        this->min_min_d_angle = min_min_d_angle;
        this->avg_min_d_angle = avg_min_d_angle;
        this->max_max_d_angle = max_max_d_angle;
        this->avg_max_d_angle = avg_max_d_angle;
        this->max_energy = max_energy;
        this->avg_energy = avg_energy;
    }
    MeshRecord(int op, double timing, int n_v, int n_t) {
        this->op = op;
        this->timing = timing;
        this->n_v = n_v;
        this->n_t = n_t;
    }
};
//extern std::vector<MeshRecord> mesh_records;
extern bool is_app_csv;
void addRecord(const MeshRecord& record);

//for serialization
namespace igl {
    namespace serialization {
        template<>
        inline void serialize(const Point_3 &p, std::vector<char> &buffer) {
            ::igl::serialize(STR(CGAL::exact(p[0])), std::string("x"), buffer);
            ::igl::serialize(STR(CGAL::exact(p[1])), std::string("y"), buffer);
            ::igl::serialize(STR(CGAL::exact(p[2])), std::string("z"), buffer);
        }

        template<>
        inline void deserialize(Point_3 &p, const std::vector<char> &buffer) {
            std::string s1, s2, s3;
            ::igl::deserialize(s1, std::string("x"), buffer);
            ::igl::deserialize(s2, std::string("y"), buffer);
            ::igl::deserialize(s3, std::string("z"), buffer);
//            p=Point_3(CGAL_FT(s1), CGAL_FT(s2), CGAL_FT(s3));
        }

        template<>
        inline void serialize(const Point_3f &p, std::vector<char> &buffer) {
            ::igl::serialize(p[0], std::string("x"), buffer);
            ::igl::serialize(p[1], std::string("y"), buffer);
            ::igl::serialize(p[2], std::string("z"), buffer);
        }

        template<>
        inline void deserialize(Point_3f &p, const std::vector<char> &buffer) {
            double x, y, z;
            ::igl::deserialize(x, std::string("x"), buffer);
            ::igl::deserialize(y, std::string("y"), buffer);
            ::igl::deserialize(z, std::string("z"), buffer);
            p=Point_3f(x, y, z);
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
#endif //GTET_HEADS_H
