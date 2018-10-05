// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 5/31/18.
//

#pragma once

#include <tetwild/Args.h>
#include <Eigen/Dense>

namespace tetwild {

///
/// Robust tetrahedralization of an input triangle soup, with an envelope constraint.
///
/// @param[in]  VI    { #VI x 3 input mesh vertices }
/// @param[in]  FI    { #FI x 3 input mesh triangles }
/// @param[out] VO    { #VO x 3 output mesh vertices }
/// @param[out] TO    { #TO x 4 output mesh tetrahedra }
/// @param[out] AO    { #TO x 1 array of min dihedral angle over each tet }
/// @param[in]  args  { Extra arguments controlling the behavior of TetWild }
///
void tetrahedralization(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI,
    Eigen::MatrixXd &VO, Eigen::MatrixXi &TO, Eigen::VectorXd &AO, const Args &args = Args());

///
/// Extract the boundary facets of a triangle mesh, removing unreferenced vertices
///
/// @param[in]  VI    { #VI x 3 input mesh vertices }
/// @param[in]  TI    { #TI x 4 input mesh tetrahedra }
/// @param[out] VS    { #VS x 3 output mesh vertices }
/// @param[out] FS    { #FS x 3 output mesh triangles }
///
void extractSurfaceMesh(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &TI,
    Eigen::MatrixXd &VS, Eigen::MatrixXi &FS);

} // namespace tetwild

