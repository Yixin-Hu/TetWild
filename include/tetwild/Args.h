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

#include <string>

namespace tetwild {

// Global arguments controlling the behavior of TetWild
struct Args {
    // Initial target edge-length at every vertex (in % of the bbox diagonal)
    double initial_edge_len_rel = 1/20.0;

    // Initial absolute target edge-length at every vertex. Only used if -a is specified.
    double initial_edge_len_abs = 0.0;

    // convenience function to get the correct absolute edge length depending
    // on what was set by the CLI reader
    double getAbsoluteEdgeLength(const double bbox_diag) const
    {
        return initial_edge_len_abs != 0.0 ? initial_edge_len_abs
          : initial_edge_len_rel*bbox_diag;
    }

    // convenience function to get the correct relative edge length depending
    // on what was set by the CLI reader
    double getRelativeEdgeLength(const double bbox_diag) const
    {
        return initial_edge_len_abs != 0.0 ? initial_edge_len_abs/bbox_diag
          : initial_edge_len_rel;
    }

    // Target epsilon (in % of the bbox diagonal)
    double eps_rel = 1e-3;

    //////////////////////
    // Advanced options //
    //////////////////////

    // Explicitly specify a sampling distance for triangles (in % of the bbox diagonal)
    int sampling_dist_rel = -1;

    // Run the algorithm in stage (as explain in p.8 of the paper)
    // If the first stage didn't succeed, call again with `stage = 2`,  etc.
    int stage = 1;

    // Multiplier for resizing the target-edge length around bad-quality vertices
    // See MeshRefinement::updateScalarField() for more details
    double adaptive_scalar = 0.6;

    // Energy threshold
    // If the max tet energy is below this threshold, the mesh optimization process is stopped.
    // Also used to determine where to resize the scalar field (if a tet incident to a vertex has larger energy than this threshold, then resize around this vertex).
    double filter_energy_thres = 10;

    // Threshold on the energy delta (avg and max) below which to rescale the target edge length scalar field
    double delta_energy_thres = 0.1;

    // Maximum number of mesh optimization iterations
    int max_num_passes = 80;

    // Sample points at voxel centers for initial Delaunay triangulation
    bool not_use_voxel_stuffing = false;

    // Use Laplacian smoothing on the faces/vertices covering an open boundary after the mesh optimization step (post-processing)
    bool smooth_open_boundary = false;

    // Target number of vertices (minimum), within 5% of tolerance
    int target_num_vertices = -1;

    // Background mesh for the edge length sizing field
    std::string background_mesh = "";

    // [debug] logging
    bool write_csv_file = true;
    std::string working_dir = "";
    std::string postfix = "_";
    std::string csv_file = "";
    int save_mid_result = -1; // save intermediate result

    bool is_quiet = false;
};

} // namespace tetwild
