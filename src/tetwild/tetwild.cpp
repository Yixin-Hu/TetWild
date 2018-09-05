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

#include <tetwild/tetwild.h>
#include <tetwild/Preprocess.h>
#include <tetwild/DelaunayTetrahedralization.h>
#include <tetwild/BSPSubdivision.h>
#include <tetwild/SimpleTetrahedralization.h>
#include <tetwild/MeshRefinement.h>
#include <tetwild/InoutFiltering.h>
#include <pymesh/MshSaver.h>

namespace tetwild {

//    MeshRefinement MR;

void outputFinalQuality(double time, const std::vector<TetVertex>& tet_vertices, const std::vector<std::array<int, 4>>& tets,
                        const std::vector<bool> &t_is_removed, const std::vector<TetQuality>& tet_qualities,
                        const std::vector<int>& v_ids) {
    logger().debug("final quality:");
    double min = 10, max = 0;
    double min_avg = 0, max_avg = 0;
//    double max_asp_ratio = 0, avg_asp_ratio = 0;
    double max_slim_energy = 0, avg_slim_energy = 0;
    std::array<double, 6> cmp_cnt = {0, 0, 0, 0, 0, 0};
    std::array<double, 6> cmp_d_angles = {6 / 180.0 * M_PI, 12 / 180.0 * M_PI, 18 / 180.0 * M_PI,
                                          162 / 180.0 * M_PI, 168 / 180.0 * M_PI, 174 / 180.0 * M_PI};
    int cnt = 0;
    for (int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        cnt++;
        if (tet_qualities[i].min_d_angle < min)
            min = tet_qualities[i].min_d_angle;
        if (tet_qualities[i].max_d_angle > max)
            max = tet_qualities[i].max_d_angle;
//        if (tet_qualities[i].asp_ratio_2 > max_asp_ratio)
//            max_asp_ratio = tet_qualities[i].asp_ratio_2;
        if (tet_qualities[i].slim_energy > max_slim_energy)
            max_slim_energy = tet_qualities[i].slim_energy;
        min_avg += tet_qualities[i].min_d_angle;
        max_avg += tet_qualities[i].max_d_angle;
//        avg_asp_ratio += tet_qualities[i].asp_ratio_2;
        avg_slim_energy += tet_qualities[i].slim_energy;

        for (int j = 0; j < 3; j++) {
            if (tet_qualities[i].min_d_angle < cmp_d_angles[j])
                cmp_cnt[j]++;
        }
        for (int j = 0; j < 3; j++) {
            if (tet_qualities[i].max_d_angle > cmp_d_angles[j + 3])
                cmp_cnt[j + 3]++;
        }
    }
    logger().debug("min_d_angle = {}, max_d_angle = {}, max_slim_energy = {}", min, max, max_slim_energy);
    logger().debug("avg_min_d_angle = {}, avg_max_d_angle = {}, avg_slim_energy = {}", min_avg / cnt, max_avg / cnt, avg_slim_energy / cnt);
    logger().debug("min_d_angle: <6 {};   <12 {};  <18 {}", cmp_cnt[0] / cnt, cmp_cnt[1] / cnt, cmp_cnt[2] / cnt);
    logger().debug("max_d_angle: >174 {}; >168 {}; >162 {}", cmp_cnt[5] / cnt, cmp_cnt[4] / cnt, cmp_cnt[3] / cnt);

    addRecord(MeshRecord(MeshRecord::OpType::OP_WN, time, v_ids.size(), cnt,
                         min, min_avg / cnt, max, max_avg / cnt, max_slim_energy, avg_slim_energy / cnt));

    ////output unrounded vertices:
    cnt = 0;
    for (int v_id: v_ids) {
        if (!tet_vertices[v_id].is_rounded)
            cnt++;
    }
    logger().debug("{}/{} vertices are unrounded!!!", cnt, v_ids.size());
    addRecord(MeshRecord(MeshRecord::OpType::OP_UNROUNDED, -1, cnt, -1));
}

void outputFinalTetmesh(MeshRefinement& MR,
                        std::vector<std::array<double, 3>>& out_vertices,
                        std::vector<std::array<int, 4>>& out_tets) {
    std::vector<TetVertex> &tet_vertices = MR.tet_vertices;
    std::vector<std::array<int, 4>> &tets = MR.tets;
    std::vector<bool> &v_is_removed = MR.v_is_removed;
    std::vector<bool> &t_is_removed = MR.t_is_removed;
    std::vector<TetQuality> &tet_qualities = MR.tet_qualities;
    int t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
    double tmp_time = 0;
    if (!GArgs::args().is_laplacian) {
        InoutFiltering IOF(tet_vertices, tets, MR.is_surface_fs, v_is_removed, t_is_removed, tet_qualities);
#ifndef MUTE_COUT
        igl::Timer igl_timer;
        igl_timer.start();
#endif
        IOF.filter();
#ifndef MUTE_COUT
        tmp_time = igl_timer.getElapsedTime();
        logger().debug("time = {}s", tmp_time);
        t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
        logger().debug("{} tets inside!", t_cnt);
#endif
    }

    //output result
    std::vector<int> v_ids;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            v_ids.push_back(tets[i][j]);
    }
    std::sort(v_ids.begin(), v_ids.end());
    v_ids.erase(std::unique(v_ids.begin(), v_ids.end()), v_ids.end());
    std::unordered_map<int, int> map_ids;
    for (int i = 0; i < v_ids.size(); i++)
        map_ids[v_ids[i]] = i;

    out_vertices.reserve(v_ids.size());
    for (int i = 0; i < v_ids.size(); i++) {
        out_vertices.push_back(std::array<double, 3>({tet_vertices[v_ids[i]].posf[0],
                                                      tet_vertices[v_ids[i]].posf[1],
                                                      tet_vertices[v_ids[i]].posf[2]}));
    }
    out_tets.reserve(std::count(t_is_removed.begin(), t_is_removed.end(), false));
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        out_tets.push_back(std::array<int, 4>({map_ids[tets[i][0]], map_ids[tets[i][1]], map_ids[tets[i][2]],
                                               map_ids[tets[i][3]]}));
    }

    if(GArgs::args().is_quiet)
        return;

#ifndef MUTE_COUT
    outputFinalQuality(tmp_time, tet_vertices, tets, t_is_removed, tet_qualities, v_ids);
#endif
}

void gtet_new(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in,
              std::vector<std::array<double, 3>>& out_vertices,
              std::vector<std::array<int, 4>>& out_tets) {
    State::state().is_using_energy_max = true;
    State::state().is_use_project = false;
    State::state().is_using_sampling = true;

    int energy_type = State::state().ENERGY_AMIPS;
    bool is_sm_single = true;
    bool is_preprocess = true;
    bool is_check_correctness = false;
    bool is_ec_check_quality = true;

#ifndef MUTE_COUT
    igl::Timer igl_timer;
    double tmp_time = 0;
    double sum_time = 0;
#endif

    ////pipeline
    MeshRefinement MR;
    {/// STAGE 1
        //preprocess
#ifndef MUTE_COUT
        igl_timer.start();
        logger().debug("Preprocessing...");
#endif
        Preprocess pp;
        if (!pp.init(V_in, F_in, MR.geo_b_mesh, MR.geo_sf_mesh)) {
            logger().debug("Empty!");
            //todo: output a empty tetmesh
            PyMesh::MshSaver mSaver(State::state().g_working_dir + State::state().g_postfix + ".msh", true);
            Eigen::VectorXd oV;
            Eigen::VectorXi oT;
            oV.resize(0);
            oT.resize(0);
            mSaver.save_mesh(oV, oT, 3, mSaver.TET);
            exit(250);
        }
#ifndef MUTE_COUT
        addRecord(MeshRecord(MeshRecord::OpType::OP_INIT, 0, MR.geo_sf_mesh.vertices.nb(), MR.geo_sf_mesh.facets.nb()));
#endif

        std::vector<Point_3> m_vertices;
        std::vector<std::array<int, 3>> m_faces;
        pp.process(MR.geo_sf_mesh, m_vertices, m_faces);
#ifndef MUTE_COUT
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_PREPROCESSING, tmp_time, m_vertices.size(), m_faces.size()));
        sum_time += tmp_time;
        logger().debug("time = {}s", tmp_time);
#endif

        //delaunay tetrahedralization
#ifndef MUTE_COUT
        igl_timer.start();
        logger().debug("Delaunay tetrahedralizing...");
#endif
        DelaunayTetrahedralization DT;
        std::vector<int> raw_e_tags;
        std::vector<std::vector<int>> raw_conn_e4v;
        std::vector<int> m_f_tags;//need to use it in ST
        DT.init(m_vertices, m_faces, m_f_tags, raw_e_tags, raw_conn_e4v);
        std::vector<Point_3> bsp_vertices;
        std::vector<BSPEdge> bsp_edges;
        std::vector<BSPFace> bsp_faces;
        std::vector<BSPtreeNode> bsp_nodes;
        DT.tetra(m_vertices, MR.geo_sf_mesh, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes);
#ifndef MUTE_COUT
        logger().debug("# bsp_vertices = {}", bsp_vertices.size());
        logger().debug("# bsp_edges = {}", bsp_edges.size());
        logger().debug("# bsp_faces = {}", bsp_faces.size());
        logger().debug("# bsp_nodes = {}", bsp_nodes.size());
        logger().debug("Delaunay tetrahedralization done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_DELAUNEY_TETRA, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        sum_time += tmp_time;
        logger().debug("time = {}s", tmp_time);
#endif

        //mesh conforming
#ifndef MUTE_COUT
        igl_timer.start();
        logger().debug("Divfaces matching...");
#endif
        MeshConformer MC(m_vertices, m_faces, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes);
        MC.match();
#ifndef MUTE_COUT
        logger().debug("Divfaces matching done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_DIVFACE_MATCH, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        logger().debug("time = {}s", tmp_time);
#endif

        //bsp subdivision
#ifndef MUTE_COUT
        igl_timer.start();
        logger().debug("BSP subdivision ...");
#endif
        BSPSubdivision BS(MC);
        BS.init();
        BS.subdivideBSPNodes();
#ifndef MUTE_COUT
        logger().debug("Output: ");
        logger().debug("# node = {}", MC.bsp_nodes.size());
        logger().debug("# face = {}", MC.bsp_faces.size());
        logger().debug("# edge = {}", MC.bsp_edges.size());
        logger().debug("# vertex = {}", MC.bsp_vertices.size());
        logger().debug("BSP subdivision done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_BSP, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        sum_time += tmp_time;
        logger().debug("time = {}s", tmp_time);
#endif

        //simple tetrahedralization
#ifndef MUTE_COUT
        igl_timer.start();
        logger().debug("Tetrehedralizing ...");
#endif
        SimpleTetrahedralization ST(MC);
        ST.tetra(MR.tet_vertices, MR.tets);
        ST.labelSurface(m_f_tags, raw_e_tags, raw_conn_e4v, MR.tet_vertices, MR.tets, MR.is_surface_fs);
        ST.labelBbox(MR.tet_vertices, MR.tets);
        if (!State::state().g_is_close)//if input is an open mesh
            ST.labelBoundary(MR.tet_vertices, MR.tets, MR.is_surface_fs);
#ifndef MUTE_COUT
        logger().debug("# tet_vertices = {}", MR.tet_vertices.size());
        logger().debug("# tets = {}", MR.tets.size());
        logger().debug("Tetrahedralization done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_SIMPLE_TETRA, tmp_time, MR.tet_vertices.size(), MR.tets.size()));
        sum_time += tmp_time;
        logger().debug("time = {}s", tmp_time);

        logger().debug("Total time for the first stage = {}", sum_time);
#endif
    }

    /// STAGE 2
    //init
#ifndef MUTE_COUT
    logger().debug("Refinement initializing...");
#endif
    MR.prepareData();
#ifndef MUTE_COUT
    logger().debug("Refinement intialization done!");
#endif

    //improvement
    MR.refine(energy_type);

    outputFinalTetmesh(MR, out_vertices, out_tets); //do winding number and output the tetmesh
}

void tetrahedralization(const std::vector<std::array<double, 3>>& in_vertices,
                        const std::vector<std::array<int, 3>>& in_faces,
                        std::vector<std::array<double, 3>>& out_vertices,
                        std::vector<std::array<int, 4>>& out_tets,
                        Args parameters)
{
    out_vertices.clear();
    out_tets.clear();

    GArgs::args().i_epsilon = parameters.i_epsilon;
    GArgs::args().bg_mesh = parameters.bg_mesh;
    GArgs::args().filter_energy = parameters.filter_energy;
    GArgs::args().i_ideal_edge_length = parameters.i_ideal_edge_length;
    GArgs::args().is_laplacian = parameters.is_laplacian;
    GArgs::args().is_quiet = parameters.is_quiet;
    GArgs::args().max_pass = parameters.max_pass;
    GArgs::args().stage = parameters.stage;
    GArgs::args().targeted_num_v = parameters.targeted_num_v;

    //initalization
    GEO::initialize();
    State::state().g_postfix = GArgs::args().postfix;
    State::state().g_working_dir = GArgs::args().input.substr(0, GArgs::args().input.size() - 4);

    if (GArgs::args().csv_file == "")
        State::state().g_stat_file = State::state().g_working_dir + State::state().g_postfix + ".csv";
    else
        State::state().g_stat_file = GArgs::args().csv_file;

    if (GArgs::args().output == "")
        State::state().g_output_file = State::state().g_working_dir + State::state().g_postfix + ".msh";
    else
        State::state().g_output_file = GArgs::args().output;

    if (GArgs::args().is_quiet) {
        GArgs::args().is_output_csv = false;
        std::cout.setstate(std::ios_base::failbit);//use std::cout.clear() to get it back
    }

    //do tetrahedralization
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    V.resize(in_vertices.size(), 3);
    F.resize(in_faces.size(), 3);
    for (int i = 0; i < in_vertices.size(); i++)
        for (int j = 0; j < 3; j++)
            V(i, j) = in_vertices[i][j];
    for (int i = 0; i < in_faces.size(); i++)
        for (int j = 0; j < 3; j++)
            F(i, j) = in_faces[i][j];
    gtet_new(V, F, out_vertices, out_tets);

    if (GArgs::args().is_quiet) {
        std::cout.clear();
    }
}

} // namespace tetwild

