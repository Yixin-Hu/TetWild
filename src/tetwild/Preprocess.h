// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 10/12/17.
//

#ifndef NEW_GTET_PREPROCESS_H
#define NEW_GTET_PREPROCESS_H

#include <tetwild/ForwardDecls.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/geogram/mesh_AABB.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
#include <unordered_set>
#include <queue>

namespace tetwild {

class ElementInQueue_sm{
public:
    std::array<int, 2> v_ids;
    double weight;

    ElementInQueue_sm(){}
    ElementInQueue_sm(const std::array<int, 2>& ids, double w):
            v_ids(ids), weight(w){}
};

struct cmp_sm {
    bool operator()(const ElementInQueue_sm &e1, const ElementInQueue_sm &e2) {
        if (e1.weight == e2.weight)
            return e1.v_ids < e2.v_ids;
        return e1.weight > e2.weight;
    }
};

class Preprocess {
    std::priority_queue<ElementInQueue_sm, std::vector<ElementInQueue_sm>, cmp_sm> sm_queue;
    int c=0;
public:
    State &state;

    Eigen::MatrixXd V_in;
    Eigen::MatrixXi F_in;
    std::vector<bool> v_is_removed;
    std::vector<bool> f_is_removed;
    std::vector<std::unordered_set<int>> conn_fs;

    Preprocess(State &st) : state(st) { }

    bool init(const Eigen::MatrixXd& V_tmp, const Eigen::MatrixXi& F_tmp, GEO::Mesh& geo_b_mesh, GEO::Mesh& geo_sf_mesh, const Args &args);

    void getBoundaryMesh(GEO::Mesh& b_mesh);
    void process(GEO::Mesh& geo_sf_mesh, std::vector<Point_3>& m_vertices, std::vector<std::array<int, 3>>& m_faces, const Args &args);

    void simplify(const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree);
    void postProcess(const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree);
    bool removeAnEdge(int v1_id, int v2_id, const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree);

    void swap(const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree);
    double getCosAngle(int v_id, int v1_id, int v2_id);

    double getEdgeLength(const std::array<int, 2>& v_ids);
    double getEdgeLength(int v1_id, int v2_id);

    bool isEdgeValid(const std::array<int, 2>& v_ids, double old_weight);
    bool isEdgeValid(const std::array<int, 2>& v_ids);
    bool isOneRingClean(int v_id);
    bool isOutEnvelop(const std::unordered_set<int>& new_f_ids, const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree);
    bool isPointOutEnvelop(int v_id, const GEO::MeshFacetsAABBWithEps& face_aabb_tree);
    bool isEuclideanValid(int v1_id, int v2_id);
    //when call this function, the coordinate of v1 has already been changed

    int ts=0;
    std::vector<std::array<int, 2>> inf_es;
    std::vector<int> inf_e_tss;
    std::vector<int> f_tss;

    void outputSurfaceColormap(const GEO::MeshFacetsAABBWithEps& geo_face_tree, const GEO::Mesh& geo_sf_mesh);
};

} // namespace tetwild

#endif //NEW_GTET_PREPROCESS_H
