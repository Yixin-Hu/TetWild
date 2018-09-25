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
#pragma once

#include <Eigen/Dense>
#include <geogram/mesh/mesh.h>

namespace tetwild {

///
/// @brief      { Converts a surface mesh to a Geogram mesh }
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  F     { #F x (3|4) input mesh surface (triangles or quads) }
/// @param[out] M     { Output Geogram mesh }
///
void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, GEO::Mesh &M);

///
/// @brief      { Converts a tet-mesh to a Geogram mesh }
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  F     { #F x (3|4) input mesh surface (triangles or quads) }
/// @param[in]  T     { #F x 4 input mesh tetrahedra }
/// @param[out] M     { Output Geogram mesh }
///
void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &T, GEO::Mesh &M);

///
/// @brief      { Extract simplices from a Geogram mesh }
///
/// @param[in]  M     { Input Geogram mesh }
/// @param[out] V     { #V x 3 output mesh vertices }
/// @param[out] F     { #F x 3 output mesh faces }
/// @param[out] T     { #T x 4 output mesh tets }
///
void from_geogram_mesh(const GEO::Mesh &M, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &T);

///
/// Samples points inside the bbox of a given set of points
///
/// @param[in]  V            { #V x 3 input mesh vertices }
/// @param[in]  num_samples  { Number of desired samples }
/// @param[in]  padding      { Padding for the bbox }
/// @param[out] P            { num_samples x 3 sample points positions }
/// @param[in]  num_lloyd    { Number of Lloyd iterations }
/// @param[in]  num_newton   { Number of Newton iterations }
///
void sample_bbox(const Eigen::MatrixXd &V, int num_samples, double padding,
    Eigen::MatrixXd &P, int num_lloyd = 10, int num_newton = 10);

///
/// Sample points on surface.
///
/// @param[in]  V            { #V x 3 input mesh vertices }
/// @param[in]  F            { #F x 3 input mesh triangles }
/// @param[in]  num_samples  { Number of target samples }
/// @param      P            { num_samples x 3 output point samples }
/// @param[in]  num_lloyd    { Number of Lloyd iterations }
/// @param[in]  num_newton   { Number of Newton iterations }
///
void resample_surface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int num_samples,
    Eigen::MatrixXd &P, int num_lloyd = 10, int num_newton = 10);

///
/// Computes a Delaunay tiangulation of a point cloud in 3d
///
/// @param[in]  V     { #V x dims input point positions }
/// @param[out] T     { #T x 4 output mesh tetrahedra }
///
void delaunay_tetrahedralization(const Eigen::MatrixXd &V, Eigen::MatrixXi &T);

} // namespace tetwild
