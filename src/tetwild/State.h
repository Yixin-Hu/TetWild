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

#include <string>

namespace tetwild {

struct State {
    const int EPSILON_INFINITE=-2;
    const int EPSILON_NA=-1;
    const int ENERGY_NA=0;
    const int ENERGY_AD=1;
    const int ENERGY_AMIPS=2;
    const int ENERGY_DIRICHLET=3;
    const double MAX_ENERGY = 1e50;
    int NOT_SURFACE = 0;

    // paths used for i/o
    std::string g_working_dir;
    std::string g_stat_file;
    std::string g_postfix;
    std::string g_output_file;

    double g_eps = 0; // effective epsilon at the current stage (see \hat{\epsilon} in the paper)
    double g_eps_2 = 0;
    double g_dd = 0; // sampling distance for triangles (see d_k p.8 of the paper)
    double g_ideal_l = 0; // target edge-length
    double g_diag_l = 0; // bbox diagonal
    bool g_is_close = 0; // open mesh or closed mesh?

    double g_eps_input = 0; // target epsilon desired by the user
    double g_eps_delta = 0; // ??
    int g_cur_stage = 1; // k

    ///////////////
    // [testing] //
    ///////////////

    // Whether to use the max or the total energy when checking improvements in local operations
    bool is_using_energy_max = true;

    // Use sampling to determine whether a face lies outside the envelope during mesh optimization
    // (if false, then only its vertices are tested)
    bool is_using_sampling = true;

    // Project vertices to the plane of their one-ring instead of the original surface during vertex smoothing
    bool is_use_project = false;

    // [debug]
    bool is_print_tmp = false;

    static State & state() {
        static State st;
        return st;
    }

private:
    State() = default;
};


struct GArgs {
    // [I/O] Filename
    std::string input;
    std::string output = "";
    std::string postfix = "_";

    // [input] User-defined arguments
    double i_ideal_edge_length = 20; // target edge-length = bbox diagonal / i_ideal_edge_length
    double i_epsilon = 1000; // target epsilon = bbox_diagonal / i_epsilon

    // Advanced options
    int i_dd = -1; // explicitly specify a sampling distance for triangles
    int stage = 1; // max_num_stage [p.8]

    // Multiplier for resizing the target-edge length around bad-quality vertices
    // See MeshRefinement::updateScalarField() for more details
    double adaptive_scalar = 0.6;

    // Energy threshold
    // If the max tet energy is below this threshold, the mesh optimization process is stopped.
    // Also used to determine where to resize the scalar field (if a tet incident to a vertex has larger energy than this threshold, then resize around this vertex).
    double filter_energy = 10;

    // Threshold on the energy delta (avg and max) below which to rescale the target edge length scalar field
    double delta_energy = 0.1;

    // Maximum number of mesh optimization iterations
    int max_pass = 80;

    // [debug] log files
    int is_output_csv = true;
    std::string csv_file = "";
    std::string slz_file = "";
    int mid_result = -1; // save intermediate result

    // Sample points at voxel centers for initial Delaunay triangulation
    bool is_using_voxel = true;

    // Use Laplacian smoothing on the faces/vertices covering an open boundary after the mesh optimization step (post-processing)
    bool is_laplacian = false;

    // Target number of vertices (minimum), within 5% of tolerance
    int targeted_num_v = -1;

    // Background mesh for the edge length sizing field
    std::string bg_mesh = "";

    bool is_quiet = false;

    static GArgs & args() {
        static GArgs ag;
        return ag;
    }

private:
    GArgs() = default;
};


struct MeshRecord {
    enum OpType {
        OP_INIT = 0,
        OP_PREPROCESSING,
        OP_DELAUNEY_TETRA,
        OP_DIVFACE_MATCH,
        OP_BSP,
        OP_SIMPLE_TETRA,

        OP_OPT_INIT,
        OP_SPLIT,
        OP_COLLAPSE,
        OP_SWAP,
        OP_SMOOTH,
        OP_ADAP_UPDATE,
        OP_WN,
        OP_UNROUNDED
    };

    int op;
    double timing;
    int n_v;
    int n_t;
    double min_min_d_angle = -1;
    double avg_min_d_angle = -1;
    double max_max_d_angle = -1;
    double avg_max_d_angle = -1;
    double max_energy = -1;
    double avg_energy = -1;

    MeshRecord(int op_, double timing_, int n_v_, int n_t_, double min_min_d_angle_, double avg_min_d_angle_,
               double max_max_d_angle_, double avg_max_d_angle_, double max_energy_, double avg_energy_) {
        this->op = op_;
        this->timing = timing_;
        this->n_v = n_v_;
        this->n_t = n_t_;
        this->min_min_d_angle = min_min_d_angle_;
        this->avg_min_d_angle = avg_min_d_angle_;
        this->max_max_d_angle = max_max_d_angle_;
        this->avg_max_d_angle = avg_max_d_angle_;
        this->max_energy = max_energy_;
        this->avg_energy = avg_energy_;
    }

    MeshRecord(int op_, double timing_, int n_v_, int n_t_) {
        this->op = op_;
        this->timing = timing_;
        this->n_v = n_v_;
        this->n_t = n_t_;
    }
};

} // namespace tetwild
