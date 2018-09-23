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
#include "utils.h"
#include <tetwild/Logger.h>
////////////////////////////////////////////////////////////////////////////////

namespace tetwild {

bool mmg_to_eigen(const MMG5_pMesh mmg, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &T) {
    logger().debug("converting MMG5_pMesh to Eigen matrices ...");

    // Note: indexing seems to start at 1 in MMG
    assert(mmg->dim == 3);
    V.resize(mmg->np, 3);
    F.resize(mmg->nt, 3);
    T.resize(mmg->ne, 4);

    for (int v = 0; v < V.rows(); ++v) {
        for (int d = 0; d < mmg->dim; ++d) {
            V(v, d) = mmg->point[v+1].c[d];
        }
    }
    for (int t = 0; t < F.rows(); ++t) {
        F(t,0) = mmg->tria[t+1].v[0] - 1;
        F(t,1) = mmg->tria[t+1].v[1] - 1;
        F(t,2) = mmg->tria[t+1].v[2] - 1;
    }
    for (int c = 0; c < T.rows(); ++c) {
        T(c,0) = mmg->tetra[c+1].v[0] - 1;
        T(c,1) = mmg->tetra[c+1].v[1] - 1;
        T(c,2) = mmg->tetra[c+1].v[2] - 1;
        T(c,3) = mmg->tetra[c+1].v[3] - 1;
    }

    return true;
}

bool eigen_to_mmg(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &T, MMG5_pMesh& mmg, MMG5_pSol& sol) {
    logger().debug("converting Eigen matrices to MMG5_pMesh ...");
    assert(V.cols() == 3);
    bool volume_mesh = (T.rows() > 0);

    if (volume_mesh) {
        MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh,&mmg,MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
    } else {
        MMGS_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh,&mmg,MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
    }

    if (volume_mesh && MMG3D_Set_meshSize(
                mmg,
                V.rows(),
                T.rows(),
                0, /* nb prisms */
                F.rows(),
                0, /* nb quad */
                0 /* nb edges */
                ) != 1 ) {
        logger().error("failed to MMG3D_Set_meshSize");
        return false;
    } else if (!volume_mesh && MMGS_Set_meshSize(
                mmg,
                V.rows(),
                F.rows(),
                0  /* nb edges */
                ) != 1 ) {
        logger().error("failed to MMGS_Set_meshSize");
        return false;
    }

    for (int v = 0; v < mmg->np; ++v) {
        for (int d = 0; d < V.cols(); ++d) {
            mmg->point[v+1].c[d] = V(v,d);
        }
    }
    for (int t = 0; t < mmg->nt; ++t) {
        mmg->tria[t+1].v[0] = F(t,0) + 1;
        mmg->tria[t+1].v[1] = F(t,1) + 1;
        mmg->tria[t+1].v[2] = F(t,2) + 1;
    }
    if (volume_mesh) {
        for (int c = 0; c < mmg->ne; ++c) {
            MMG3D_Set_tetrahedron(mmg, T(c,0)+1, T(c,1)+1, T(c,2)+1, T(c,3)+1, 0, c+1);
            // mmg->tetra[c+1].v[0] = T(c,0) + 1;
            // mmg->tetra[c+1].v[1] = T(c,1) + 1;
            // mmg->tetra[c+1].v[2] = T(c,2) + 1;
            // mmg->tetra[c+1].v[3] = T(c,3) + 1;
        }
    }

    if (volume_mesh && MMG3D_Set_solSize(mmg,sol,MMG5_Vertex,V.rows(),MMG5_Scalar) != 1 ) {
        logger().error("failed to MMG3D_Set_solSize");
        return false;
    } else if (!volume_mesh && MMGS_Set_solSize(mmg,sol,MMG5_Vertex,V.rows(),MMG5_Scalar) != 1 ) {
        logger().error("failed to MMGS_Set_solSize");
        return false;
    }
    for (int v = 0; v < V.rows(); ++v) {
        sol->m[v+1] = 1.0;
    }
    if (volume_mesh && MMG3D_Chk_meshData(mmg,sol) != 1) {
        logger().error("error in mmg: inconsistent mesh and sol");
        return false;
    } else if (!volume_mesh && MMGS_Chk_meshData(mmg,sol) != 1) {
        logger().error("error in mmg: inconsistent mesh and sol");
        return false;
    }

    if (volume_mesh) {
        // Orient each tetrahedra to have positive volume.
        // TODO: Use the API function MMG3D_Set_tetrahedron instead
        // MMG3D_Set_handGivenMesh(mmg);
    }

    return true;
}

bool mmg_to_eigen(const MMG5_pMesh mmg, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    Eigen::MatrixXi T;
    return mmg_to_eigen(mmg, V, F, T);
}

bool eigen_to_mmg(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, MMG5_pMesh& mmg, MMG5_pSol& sol) {
    Eigen::MatrixXi T;
    return eigen_to_mmg(V, F, T, mmg, sol);
}

void mmg3d_free(MMG5_pMesh mmg, MMG5_pSol sol) {
    MMG3D_Free_all(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmg,MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
}

void mmgs_free(MMG5_pMesh mmg, MMG5_pSol sol) {
    MMGS_Free_all(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmg,MMG5_ARG_ppMet,&sol, MMG5_ARG_end);
}

// bool mmg_wrapper_test_geo2mmg2geo(const Mesh& M_in, Mesh& M_out) {
//     MMG5_pMesh mmg = NULL;
//     MMG5_pSol sol = NULL;
//     bool ok = geo_to_mmg(M_in, mmg, sol);
//     if (!ok) return false;
//     ok = mmg_to_geo(mmg, M_out);
//     mmg3d_free(mmg, sol);
//     return ok;
// }

} // namespace tetwild
