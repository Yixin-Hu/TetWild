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

#include <tetwild/Utils.h>
#include <igl/remove_unreferenced.h>

namespace tetwild {

void extractVertices(const std::vector<TetVertex> &verts, Eigen::MatrixXd &V) {
    V.resize(verts.size(), 3);
    for (int v = 0; v < V.rows(); ++v) {
        V.row(v) << verts[v].posf[0], verts[v].posf[1], verts[v].posf[2];
    }
}

void extractTriangles(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    const std::vector<std::array<int, 4>> &is_surface_fs,
    Eigen::MatrixXi &F,
    const State &state)
{
    std::vector<std::array<int, 3>> fs;
    for (int i = 0; i < tets.size(); ++i) {
        if (tet_is_removed[i]) {
            continue;
        }
        for (int j = 0; j < 4; j++) {
            if (is_surface_fs[i][j] != state.NOT_SURFACE && is_surface_fs[i][j] > 0) {//outside
                std::array<int, 3> v_ids = {{tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4]}};
                if (CGAL::orientation(verts[v_ids[0]].pos, verts[v_ids[1]].pos,
                                      verts[v_ids[2]].pos, verts[tets[i][j]].pos) != CGAL::POSITIVE) {
                    std::swap(v_ids[0], v_ids[2]);
                }
                // push back duplicated faces as many times as needed
                for (int k = 0; k < is_surface_fs[i][j]; k++) {
                    fs.push_back(v_ids);
                }
            }
        }
    }
    F.resize(fs.size(), 3);
    for (int f = 0; f < F.rows(); ++f) {
        F.row(f) << fs[f][0], fs[f][1], fs[f][2];
    }
}

void extractTetrahedra(const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    Eigen::MatrixXi &T)
{
    int num_tets = std::count(tet_is_removed.begin(), tet_is_removed.end(), false);
    T.resize(num_tets, 4);
    for (int t = 0, cnt = 0; t < tets.size(); ++t) {
        if (!tet_is_removed[t]) {
            T.row(cnt++) << tets[t][0], tets[t][1], tets[t][2], tets[t][3];
        }
    }
}

void extractSurfaceMesh(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    const std::vector<std::array<int, 4>> &is_surface_fs,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    const State &state)
{
    Eigen::MatrixXd V0;
    Eigen::MatrixXi F0;
    extractVertices(verts, V0);
    extractTriangles(verts, tets, tet_is_removed, is_surface_fs, F0, state);
    Eigen::VectorXi I;
    igl::remove_unreferenced(V0, F0, V, F, I);
}

void extractVolumeMesh(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &T)
{
    Eigen::MatrixXd V0;
    Eigen::MatrixXi T0;
    extractVertices(verts, V0);
    extractTetrahedra(tets, tet_is_removed, T0);
    Eigen::VectorXi I;
    igl::remove_unreferenced(V0, T0, V, T, I);
}

void filterRegion(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXi &R, int id,
    Eigen::MatrixXd &OV, Eigen::MatrixXi &OT)
{
    int num_in_region = (R.array() == id).count();
    OT.resize(num_in_region, T.cols());
    for (int e = 0, cnt = 0; e < T.rows(); ++e) {
        if (R(e) == id) {
            OT.row(cnt++) = T.row(e);
        }
    }
    Eigen::VectorXi I;
    igl::remove_unreferenced(V, OT, OV, OT, I);
}

} // namespace tetwild
