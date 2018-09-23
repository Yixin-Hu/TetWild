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
////////////////////////////////////////////////////////////////////////////////
#include "wrapper.h"
#include "internal/utils.h"
#include <tetwild/Logger.h>
#include <tetwild/geogram/utils.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/boundary_facets.h>
#include <igl/signed_distance.h>
#include <igl/barycenter.h>
#include <igl/winding_number.h>
#include <tetwild/geogram/utils.h>
#include <geogram/mesh/mesh_io.h>
////////////////////////////////////////////////////////////////////////////////
//
// Wrapper for 3D remeshing comes from:
// https://github.com/mxncr/mmgig
//

namespace tetwild {

namespace {

bool mmgs_tri_remesh(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI,
    Eigen::MatrixXd &VO, Eigen::MatrixXi &FO, const tetwild::MmgOptions& opt)
{
    MMG5_pMesh mesh = NULL;
    MMG5_pSol met = NULL;
    bool ok = eigen_to_mmg(VI, FI, mesh, met);
    if (!ok) {
        logger().error("mmgs: failed to convert mesh to MMG5_pMesh");
        mmgs_free(mesh, met);
        return false;
    }

    // Set remeshing options
    MMGS_Set_dparameter(mesh, met, MMGS_DPARAM_angleDetection, opt.angle_value);
    if (opt.enable_anisotropy) {
        MMGS_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor);
    }
    if (opt.hsiz == 0.0) {
        MMGS_Set_dparameter(mesh, met, MMGS_DPARAM_hmin, opt.hmin);
        MMGS_Set_dparameter(mesh, met, MMGS_DPARAM_hmax, opt.hmax);
    } else {
        met->np = 0;
        MMGS_Set_dparameter(mesh, met, MMGS_DPARAM_hsiz, opt.hsiz);
    }
    MMGS_Set_dparameter(mesh, met, MMGS_DPARAM_hausd, opt.hausd);
    MMGS_Set_dparameter(mesh, met, MMGS_DPARAM_hgrad, opt.hgrad);
    MMGS_Set_iparameter(mesh, met, MMGS_IPARAM_angle, int(opt.angle_detection));
    MMGS_Set_iparameter(mesh, met, MMGS_IPARAM_noswap, int(opt.noswap));
    MMGS_Set_iparameter(mesh, met, MMGS_IPARAM_noinsert, int(opt.noinsert));
    MMGS_Set_iparameter(mesh, met, MMGS_IPARAM_nomove, int(opt.nomove));

    int ier = MMGS_mmgslib(mesh,met);
    if (ier != MMG5_SUCCESS) {
        logger().error("mmgs_remesh: failed to remesh");
        mmgs_free(mesh, met);
        return false;
    }

    ok = mmg_to_eigen(mesh, VO, FO);

    mmgs_free(mesh, met);
    return ok;
}

bool mmg3d_tet_remesh(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, const Eigen::MatrixXi &TI,
    Eigen::MatrixXd &VO, Eigen::MatrixXi &FO, Eigen::MatrixXi &TO, const tetwild::MmgOptions& opt)
{
    MMG5_pMesh mesh = NULL;
    MMG5_pSol met = NULL;
    bool ok = eigen_to_mmg(VI, FI, TI, mesh, met);
    if (!ok) {
        logger().error("mmg3d: failed to convert mesh to MMG5_pMesh");
        mmg3d_free(mesh, met);
        return false;
    }

    // Set remeshing options
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_angleDetection, opt.angle_value);
    if (opt.enable_anisotropy) {
        MMG3D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor);
    }
    if (opt.hsiz == 0.0) {
        MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hmin, opt.hmin);
        MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hmax, opt.hmax);
    } else {
        met->np = 0;
        MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hsiz, opt.hsiz);
    }
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hausd, opt.hausd);
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hgrad, opt.hgrad);
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_angle, int(opt.angle_detection));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_noswap, int(opt.noswap));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_noinsert, int(opt.noinsert));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_nomove, int(opt.nomove));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_nosurf, int(opt.nosurf));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_opnbdy, int(opt.opnbdy));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_optim, int(opt.optim));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_optimLES, int(opt.optimLES));

    int ier = MMG3D_mmg3dlib(mesh,met);
    if (ier != MMG5_SUCCESS) {
        logger().error("mmg3d_remesh: failed to remesh");
        mmg3d_free(mesh, met);
        return false;
    }

    ok = mmg_to_eigen(mesh, VO, FO, TO);

    mmg3d_free(mesh, met);
    return ok;
}

bool mmg3d_extract_iso(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, const Eigen::MatrixXi &TI, const Eigen::VectorXd &SI,
    Eigen::MatrixXd &VO, Eigen::MatrixXi &FO, Eigen::MatrixXi &TO, const tetwild::MmgOptions& opt)
{
    if (!opt.level_set) {
        logger().error("mmg3d_iso: opt.level_set is set to false, cancel");
        return false;
    }
    if (opt.angle_detection) {
        logger().warn("mmg3d_iso: angle_detection should probably be disabled because level set functions are smooth");
    }

    MMG5_pMesh mesh = NULL;
    MMG5_pSol met = NULL;
    bool ok = eigen_to_mmg(VI, FI, TI, mesh, met);
    if (!ok) {
        logger().error("mmg3d: failed to convert mesh to MMG5_pMesh");
        mmg3d_free(mesh, met);
        return false;
    }
    assert(VI.rows() == SI.size());
    for (int v = 0; v < SI.size(); ++v) {
        met->m[v+1] = SI[v];
    }

    // Flag border for future deletion
    // std::vector<bool> on_border(M.vertices.nb(), false);
    // for (index_t t = 0; t < M.cells.nb(); ++t) {
    //     for (index_t lf = 0; lf < M.cells.nb_facets(t); ++lf) {
    //         if (M.cells.adjacent(t,lf) != GEO::NO_CELL) continue;
    //         for (index_t lv = 0; lv < M.cells.facet_nb_vertices(t,lf); ++lv) {
    //             on_border[M.cells.facet_vertex(t,lf,lv)] = true;
    //         }
    //     }
    // }

    // Set remeshing options
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_iso, 1);
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_ls, opt.ls_value);
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_angleDetection, opt.angle_value);
    if (opt.hsiz == 0.0) {
        MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hmin, opt.hmin);
        MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hmax, opt.hmax);
    } else {
        logger().error("mmg3d_iso: should not use hsiz parameter for level set mode");
        mmg3d_free(mesh, met);
        return false;
    }
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hausd, opt.hausd);
    MMG3D_Set_dparameter(mesh, met, MMG3D_DPARAM_hgrad, opt.hgrad);
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_angle, int(opt.angle_detection));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_noswap, int(opt.noswap));
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_noinsert, 1);
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_nomove, 1);
    MMG3D_Set_iparameter(mesh, met, MMG3D_IPARAM_nosurf, 1);

    int ier = MMG3D_mmg3dls(mesh,met);
    if (ier != MMG5_SUCCESS) {
        logger().error("mmg3d_iso: failed to remesh isovalue");
        mmg3d_free(mesh, met);
        return false;
    }

    // Convert back
    ok = mmg_to_eigen(mesh, VO, FO, TO);
    // GEO::Attribute<double> ls_out(M_out.vertices.attributes(), opt.ls_attribute);
    // for(uint v = 0; v < M_out.vertices.nb(); ++v) {
    //     ls_out[v] = met->m[v+1];
    // }
    /* Extract only the border */
    // M_out.cells.clear(false,false);
    // M_out.vertices.remove_isolated();
    // GEO::vector<index_t> to_del(M_out.facets.nb(), 0);
    // for (index_t f = 0; f < M_out.facets.nb(); ++f) {
    //     double d = 0;
    //     bool f_on_border = true;
    //     for (index_t lv = 0; lv < M_out.facets.nb_vertices(f); ++lv) {
    //         d = geo_max(d,std::abs(ls_out[M_out.facets.vertex(f,0)] - opt.ls_value));
    //         if (M_out.facets.vertex(f,lv) < M.vertices.nb()) {
    //             if (!on_border[M_out.facets.vertex(f,lv)]) f_on_border = false;
    //         } else {
    //             f_on_border = false;
    //         }
    //     }
    //     // if (d > 1.1 * opt.hmin) {
    //     //     to_del[f] = 1;
    //     // }
    //     if (f_on_border) to_del[f] = 1;
    // }
    // M_out.facets.delete_elements(to_del, true);

    mmg3d_free(mesh, met);
    return ok;
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

void remesh_uniform_sf(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, const MmgOptions &opt)
{
    assert(V.cols() == 3);
    mmgs_tri_remesh(V, F, OV, OF, opt);
}

void remesh_uniform_3d(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT, const MmgOptions &opt)
{
    assert(V.cols() == 3);
    Eigen::MatrixXi F;
    igl::boundary_facets(T, F);
    mmg3d_tet_remesh(V, F, T, OV, OF, OT, opt);
}

void isosurface_remeshing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXd &S,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT, const MmgOptions &opt)
{
    assert(V.cols() == 3);
    MmgOptions lv_opt = opt;
    lv_opt.level_set = true;
    lv_opt.angle_detection = false;
    Eigen::MatrixXi F;
    igl::boundary_facets(T, F);
    mmg3d_extract_iso(V, F, T, S, OV, OF, OT, lv_opt);
}

void isosurface_remeshing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int num_samples,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT, const MmgOptions &opt)
{
    // Compute ambient points by sampling CVT points inside a given volume
    Eigen::MatrixXd ambient_vertices;
    Eigen::MatrixXi ambient_tets;
    sample_bbox(V, num_samples, 0.01 * igl::bounding_box_diagonal(V), ambient_vertices);
    delaunay_tetrahedralization(ambient_vertices, ambient_tets);

    // Compute signed distance function
    Eigen::VectorXd S;
    Eigen::VectorXi face_ids;
    Eigen::MatrixXd closest_pts, normals;
    igl::signed_distance(ambient_vertices, V, F,
        igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER,
        S, face_ids, closest_pts, normals);
    assert(S.size() == V.rows());
    for (int v = 0; v < S.size(); ++v) {
        if (std::isnan(S[v])) {
            S[v] = 0.0;
        }
    }

    // Extract iso-surface from level set
    MmgOptions lv_opt = opt;
    lv_opt.level_set = true;
    lv_opt.angle_detection = false;
    Eigen::MatrixXi ambient_facets;
    igl::boundary_facets(ambient_tets, ambient_facets);

    // GEO::Mesh M;
    // to_geogram_mesh(ambient_vertices, ambient_facets, ambient_tets, M);
    // GEO::Attribute<double> attr(M.vertices.attributes(), "sdf");
    // for (int v = 0; v < ambient_facets.rows(); ++v) {
    //     attr[v] = S[v];
    // }
    // mesh_save(M, "tmp.geogram");

    mmg3d_extract_iso(ambient_vertices, ambient_facets, ambient_tets, S, OV, OF, OT, lv_opt);

    // Keep tets inside using winding number
    Eigen::MatrixXd P;
    Eigen::VectorXd W;
    igl::barycenter(OV, OT, P);
    igl::winding_number(V, F, P, W);
    Eigen::MatrixXi TT(OT.rows(), OT.cols());
    int cnt = 0;
    for (int e = 0; e < OT.rows(); ++e) {
        if (W(e) > 0.5) {
            TT.row(cnt++) = OT.row(e);
        }
    }
    OT = TT.topRows(cnt);

    igl::boundary_facets(OT, OF);
}

} // namespace tetwild
