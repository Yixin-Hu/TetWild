// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 5/6/17.
//

#ifndef NEW_GTET_LOCALOPERATIONS_H
#define NEW_GTET_LOCALOPERATIONS_H

#include <tetwild/ForwardDecls.h>
#include <tetwild/TetmeshElements.h>
#include <tetwild/geogram/mesh_AABB.h>
#include <igl/grad.h>
#include <igl/Timer.h>

#ifdef TETWILD_WITH_ISPC
#include <ispc/energy.h>
#endif

namespace tetwild {

enum class EnvelopSide{
    OUTSIDE=0,
    INSIDE=1,
    UNCERTAIN=2
};

class LocalOperations {
public:
    const Args & args;
    State & state;

    std::vector<TetVertex>& tet_vertices;
    std::vector<std::array<int, 4>>& tets;
    std::vector<std::array<int, 4>>& is_surface_fs;
    std::vector<bool>& v_is_removed;
    std::vector<bool>& t_is_removed;
    std::vector<TetQuality>& tet_qualities;

    int energy_type;

    const GEO::Mesh &geo_sf_mesh;
    const GEO::MeshFacetsAABBWithEps& geo_sf_tree;
    const GEO::MeshFacetsAABBWithEps& geo_b_tree;

    int counter=0;
    int suc_counter=0;

    std::array<double, 6> cmp_d_angles = {{6/180.0*M_PI, 12/180.0*M_PI, 18/180.0*M_PI, 162/180.0*M_PI, 168/180.0*M_PI, 174/180.0*M_PI}};

    LocalOperations(std::vector<TetVertex>& t_vs, std::vector<std::array<int, 4>>& ts, std::vector<std::array<int, 4>>& is_sf_fs,
                    std::vector<bool>& v_is_rm, std::vector<bool>& t_is_rm, std::vector<TetQuality>& tet_qs,
                    int e_type, const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& geo_tree, const GEO::MeshFacetsAABBWithEps& b_t,
                    const Args &ar, State &st) :
        tet_vertices(t_vs), tets(ts), is_surface_fs(is_sf_fs), v_is_removed(v_is_rm), t_is_removed(t_is_rm),
        tet_qualities(tet_qs), energy_type(e_type),
        geo_sf_mesh(geo_mesh), geo_sf_tree(geo_tree), geo_b_tree(b_t),
        args(ar), state(st)
    { }

    void check();
    void outputInfo(int op_type, double time, bool is_log = true);

    void calTetQualities(const std::vector<std::array<int, 4>>& new_tets, std::vector<TetQuality>& tet_qs, bool all_measure = false);
    void calTetQualities(const std::vector<int>& t_ids, bool all_measure = false);

    double calEdgeLength(const std::array<int, 2>& v_ids);
    double calEdgeLength(int v1_id, int v2_id, bool is_over_refine=false);
    void calTetQuality_AD(const std::array<int, 4>& tet, TetQuality& t_quality);
    void calTetQuality_AMIPS(const std::array<int, 4>& tet, TetQuality& t_quality);

    bool isFlip(const std::vector<std::array<int, 4>>& new_tets);
    bool isTetFlip(const std::array<int, 4>& t);
    bool isTetFlip(int t_id);

//    void getWorstQuality(TetQuality& tq);
    void getAvgMaxEnergy(double& avg_tq, double& max_tq);
    double getMaxEnergy();
    double getSecondMaxEnergy(double max_energy);
    double getFilterEnergy(bool& is_clean_up);

    void getCheckQuality(const std::vector<TetQuality>& tet_qs, TetQuality& tq);
    void getCheckQuality(const std::vector<int>& t_ids, TetQuality& tq);

    bool isEdgeOnSurface(int v1_id, int v2_id);
    bool isEdgeOnBbox(int v1_id, int v2_id);
    bool isEdgeOnSurface(int v1_id, int v2_id, const std::vector<int>& t_ids);
    bool isEdgeOnBbox(int v1_id, int v2_id, const std::vector<int>& t_ids);
    bool isEdgeOnBoundary(int v1_id, int v2_id);

//    EnvelopSide getUpperLowerBounds(const Triangle_3f& tri);
    bool isFaceOutEnvelop(const Triangle_3f& tri);
    bool isPointOutEnvelop(const Point_3f& p);
    bool isFaceOutEnvelop_sampling(const Triangle_3f& tri);
    bool isPointOutBoundaryEnvelop(const Point_3f& p);
    bool isBoundarySlide(int v1_id, int v2_id, Point_3f& pf);

    bool isTetOnSurface(int t_id);
    bool isTetRounded(int t_id);
    void getFaceConnTets(int v1_id, int v2_id, int v3_id, std::vector<int>& t_ids);
    bool isIsolated(int v_id);
    bool isBoundaryPoint(int v_id);

    static double comformalAMIPSEnergy_new(const double * T);
    static void comformalAMIPSJacobian_new(const double * T, double *result_0);
    static void comformalAMIPSHessian_new(const double * T, double *result_0);

    igl::Timer igl_timer0;
    int id_sampling=0;
    int id_aabb=1;
    std::array<double, 2> breakdown_timing0;
    std::array<std::string, 2> breakdown_name0={{"Envelop_sampling", "Envelop_AABBtree"}};

    void checkUnrounded();
    int mid_id=0;
    void outputSurfaceColormap(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in, double old_eps);

    bool isLocked_ui(const std::array<int, 2>& e);
    bool isTetLocked_ui(int tid);
};

} // namespace tetwild

#endif //NEW_GTET_LOCALOPERATIONS_H
