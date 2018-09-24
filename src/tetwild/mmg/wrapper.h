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

#include <Eigen/Dense>

namespace tetwild {

// See MmgTools documentation for interpreation
// https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options
//
struct MmgOptions {
    // Remeshing
    bool angle_detection = true;
    double angle_value = 45.;
    double hausd = 0.01;
    double hsiz = 0.0; // using hmin and hmax if set to 0
    double hmin = 0.01;
    double hmax = 2.0;
    double hgrad = 1.105171;
    bool enable_anisotropy = false;
    bool optim = false;
    bool optimLES = false;
    bool opnbdy = false;
    bool noinsert = false;
    bool noswap = false;
    bool nomove = false;
    bool nosurf = false;
    // Level set extraction
    bool level_set = false;
    double ls_value = 0.0;
};

///
/// Remesh a surface triangle mesh uniformly using mmgs.
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  F     { #T x 4 input mesh tetrahedra }
/// @param[out] OV    { #OV x (2|3) output mesh vertices }
/// @param[out] OF    { #OF x F output mesh triangles }
/// @param[in]  opt   { MMG parameters }
///
void remesh_uniform_sf(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
	Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, const MmgOptions &opt = MmgOptions());

///
/// Remesh a 3d tet-mesh uniformly using mmgs.
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  T     { #T x 4 input mesh tetrahedra }
/// @param[out] OV    { #OV x 3 output mesh vertices }
/// @param[out] OF    { #OF x F output mesh boundary triangles }
/// @param[out] OT    { #OT x 4 output mesh tetrahedra }
/// @param[in]  opt   { MMG parameters }
///
void remesh_uniform_3d(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T,
	Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT, const MmgOptions &opt = MmgOptions());

///
/// Mesh isosurface defined by the 0 level set of a tetrahedral mesh.
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  T     { #T x 4 input mesh tetrahedra }
/// @param[in]  S     { $V x 1 per-vertex with level set scalar values }
/// @param[out] OV    { #OV x 3 output mesh vertices }
/// @param[out] OF    { #OF x 3 output mesh boundary triangles }
/// @param[out] OT    { #OT x 4 output mesh tetrahedra }
/// @param[in]  opt   { MMG parameters }
///
void isosurface_remeshing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXd &S,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT, const MmgOptions &opt = MmgOptions());

///
/// Isosurface extraction from the 0-level-set of a signed distance field built on top of the given
/// mesh.
///
/// @param[in]  V            { #V x 3 input mesh vertices }
/// @param[in]  F            { #F x 3 input mesh triangles }
/// @param[in]  num_samples  { Number of samples for the SDF }
/// @param[out] OV           { #OV x 3 output mesh vertices }
/// @param[out] OF           { #OF x 3 output mesh boundary triangles }
/// @param[out] OT           { #OT x 4 output mesh tetrahedra }
/// @param[in]  opt          { MMG parameters }
///
void isosurface_remeshing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int num_samples,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT, const MmgOptions &opt = MmgOptions());

} // namespace tetwild
