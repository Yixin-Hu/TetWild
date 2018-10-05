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
/// @brief      Determines if the mesh quality is ok before passing on to mmg.
///
/// @param[in]  verts           { Input mesh vertices }
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
///
/// @return     true if mesh quality ok, false otherwise.
///
bool isMeshQualityOk(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed);

// Same as above
bool isMeshQualityOk(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);

///
/// @brief      Check that tetrahedra have non-zero volumes.
///
/// @param[in]  verts           { Input mesh vertices }
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
///
/// @return     true if all tets have non-zero unoriented volumes, false otherwise.
///
bool checkVolume(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed);

// Same as above
bool checkVolume(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);

///
/// @brief      Determines if the tet mesh has no slivers.
///
/// @param[in]  verts           { Input mesh vertices }
/// @param[in]  tets            { List of input tets }
/// @param[in]  tet_is_removed  { Whether a tet has been removed }
/// @param[in]  angle_thres     { Angle threshold in degree }
///
/// @return     true if there are no slivers below the given threshold, false otherwise.
///
bool hasNoSlivers(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    double angle_thres);

} // namespace tetwild
