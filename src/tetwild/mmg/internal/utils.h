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
#include <mmg/libmmg.h>

namespace tetwild {

bool mmg_to_eigen(const MMG5_pMesh mmg, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
bool mmg_to_eigen(const MMG5_pMesh mmg, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &T, Eigen::VectorXi *R = nullptr);

bool eigen_to_mmg(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, MMG5_pMesh& mmg, MMG5_pSol& sol);
bool eigen_to_mmg(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &T, MMG5_pMesh& mmg, MMG5_pSol& sol);

void mmgs_free(MMG5_pMesh mmg, MMG5_pSol sol);
void mmg3d_free(MMG5_pMesh mmg, MMG5_pSol sol);

} // namespace tetwild
