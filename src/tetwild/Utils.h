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

#include <tetwild/ForwardDecls.h>
#include <tetwild/TetmeshElements.h>
#include <Eigen/Dense>

namespace tetwild {

///
/// @brief      { Extract rounded vertices positions in matrix form }
///
/// @param[in]  verts  { Input mesh vertices }
/// @param[out] V      { #V x 3 matrix of output vertices positions }
///
void extractVertices(const std::vector<TetVertex> &verts, Eigen::MatrixXd &V);

///
/// @brief      { Extract boundary triangles in matrix form. }
///
/// @param[in]  verts           { Input mesh vertices }
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
/// @param[in]  is_surface_fs   { How many (oriented) boundary facets are covered (duplicated) on this facet. }
/// @param[out] F               { #F x 3 matrix of output facets }
/// @param[in]  state           { Global state of the program }
///
void extractTriangles(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    const std::vector<std::array<int, 4>> &is_surface_fs,
    Eigen::MatrixXi &F,
    const State &state);

///
/// @brief      { Extract tetrahedra in matrix form, skipping over removed tets }
///
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
/// @param[out] T               { #T x 4 matrix of output tets }
///
void extractTetrahedra(const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    Eigen::MatrixXi &T);

///
/// @brief      { Extract boundary triangle mesh in matrix form, remapping vertex indices }
///
/// @param[in]  verts           { Input mesh vertices }
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
/// @param[in]  is_surface_fs   { How many (oriented) boundary facets are covered (duplicated) on this facet. }
/// @param[out] V               { #V x 3 matrix of output vertices positions }
/// @param[out] F               { #F x 3 matrix of output facets }
/// @param[in]  state           { Global state of the program }
///
void extractSurfaceMesh(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    const std::vector<std::array<int, 4>> &is_surface_fs,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    const State &state);

///
/// @brief      { Extract ambient tet mesh in matrix form, skipping over removed tets }
///
/// @param[in]  verts           { Input mesh vertices }
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
/// @param[out] V               { #V x 3 matrix of output vertices positions }
/// @param[out] T               { #T x 4 matrix of output tets }
///
void extractVolumeMesh(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &T);

///
/// @brief      { Filter tet mesh corresponding to a given region }
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  T     { #T x 4 input mesh tetrahedra }
/// @param[in]  R     { #T x 1 input region tags }
/// @param[in]  id    { Region tag to keep }
/// @param[out] OV    { #OV x 3 output mesh vertices }
/// @param[out] OT    { #OT x 4 output mesh tetrahedra }
///
void filterRegion(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXi &R, int id,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OT);

} // namespace tetwild
