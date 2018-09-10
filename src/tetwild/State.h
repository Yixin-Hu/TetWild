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

    //global parameters
    std::string g_working_dir;
    std::string g_stat_file;
    std::string g_postfix;
    std::string g_output_file;

    double g_eps = 0;
    double g_eps_2 = 0;
    double g_dd = 0;
    double g_ideal_l = 0;
    double g_diag_l = 0;
    bool g_is_close = 0;

    double g_eps_input = 0;
    double g_eps_delta = 0;
    int g_cur_stage = 1;

    //for test
    bool is_using_energy_max = true;
    bool is_using_sampling = true;
    bool is_use_project = false;
    bool is_print_tmp = false;

    static State & state() {
        static State st;
        return st;
    }

private:
    State() = default;
};


struct GArgs {
    std::string input;
    std::string output = "";
    std::string postfix = "_";
    double i_ideal_edge_length = 20;
    double i_epsilon = 1000;
    int i_dd = -1;
    int stage = 1;
    double adaptive_scalar = 0.6;
    double filter_energy = 10;
    double delta_energy = 0.1;
    int max_pass = 80;
    int is_output_csv = true;
    std::string csv_file = "";
    std::string slz_file = "";

    int mid_result = -1;
    bool is_using_voxel = true;
    bool is_laplacian = false;

    int targeted_num_v = -1;
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
