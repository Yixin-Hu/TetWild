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
#include "tetwild.h"
#include "Preprocess.h"
#include "DelaunayTetrahedralization.h"
#include "BSPSubdivision.h"
#include "SimpleTetrahedralization.h"
#include "MeshRefinement.h"
#include "InoutFiltering.h"
#include "CLI11.hpp"

namespace tetwild {
//    MeshRefinement MR;
    Args parameters;

    void outputFinalQuality(double time, const std::vector<TetVertex>& tet_vertices, const std::vector<std::array<int, 4>>& tets,
                            const std::vector<bool> &t_is_removed, const std::vector<TetQuality>& tet_qualities,
                            const std::vector<int>& v_ids) {
        cout << "final quality:" << endl;
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
        cout << "min_d_angle = " << min
             << ", max_d_angle = " << max
             //         << ", max_aspect_ratio = " << max_asp_ratio
             << ", max_slim_energy = " << max_slim_energy
             << endl;
        cout << "avg_min_d_angle = " << min_avg / cnt
             << ", avg_max_d_angle = " << max_avg / cnt
             //         << ", avg_aspect_ratio = " << avg_asp_ratio / cnt
             << ", avg_slim_energy = " << avg_slim_energy / cnt
             << endl;
        cout << "min_d_angle: <6 " << cmp_cnt[0] / cnt << ";   <12 " << cmp_cnt[1] / cnt << ";  <18 " << cmp_cnt[2] / cnt
             << endl;
        cout << "max_d_angle: >174 " << cmp_cnt[5] / cnt << "; >168 " << cmp_cnt[4] / cnt << "; >162 " << cmp_cnt[3] / cnt
             << endl;

        addRecord(MeshRecord(MeshRecord::OpType::OP_WN, time, v_ids.size(), cnt,
                             min, min_avg / cnt, max, max_avg / cnt, max_slim_energy, avg_slim_energy / cnt));

        ////output unrounded vertices:
        cnt = 0;
        for (int v_id: v_ids) {
            if (!tet_vertices[v_id].is_rounded)
                cnt++;
        }
        cout << cnt << "/" << v_ids.size() << " vertices are unrounded!!!" << endl;
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
        if (!args.is_laplacian) {
            InoutFiltering IOF(tet_vertices, tets, MR.is_surface_fs, v_is_removed, t_is_removed, tet_qualities);
            igl::Timer igl_timer;
            igl_timer.start();
            IOF.filter();
            tmp_time = igl_timer.getElapsedTime();
            cout << "time = " << tmp_time << "s" << endl;
            t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
            cout << t_cnt << " tets inside!" << endl;
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

        if(args.is_quiet)
            return;

        outputFinalQuality(tmp_time, tet_vertices, tets, t_is_removed, tet_qualities, v_ids);
    }

    void gtet_new(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in,
                  std::vector<std::array<double, 3>>& out_vertices,
                  std::vector<std::array<int, 4>>& out_tets) {
        is_using_energy_max = true;
        is_use_project = false;
        is_using_sampling = true;

        int energy_type = ENERGY_AMIPS;
        bool is_sm_single = true;
        bool is_preprocess = true;
        bool is_check_correctness = false;
        bool is_ec_check_quality = true;

        igl::Timer igl_timer;
        double tmp_time = 0;
        double sum_time = 0;

        ////pipeline
        MeshRefinement MR;
        {/// STAGE 1
            //preprocess
            igl_timer.start();
            cout << "Preprocessing..." << endl;
            Preprocess pp;
            if (!pp.init(V_in, F_in, MR.geo_b_mesh, MR.geo_sf_mesh)) {
                cout << "Empty!" << endl;
                //todo: output a empty tetmesh
                PyMesh::MshSaver mSaver(g_working_dir + g_postfix + ".msh", true);
                Eigen::VectorXd oV;
                Eigen::VectorXi oT;
                oV.resize(0);
                oT.resize(0);
                mSaver.save_mesh(oV, oT, 3, mSaver.TET);
                exit(250);
            }
            addRecord(MeshRecord(MeshRecord::OpType::OP_INIT, 0, MR.geo_sf_mesh.vertices.nb(), MR.geo_sf_mesh.facets.nb()));

            std::vector<Point_3> m_vertices;
            std::vector<std::array<int, 3>> m_faces;
            pp.process(MR.geo_sf_mesh, m_vertices, m_faces);
            tmp_time = igl_timer.getElapsedTime();
            addRecord(MeshRecord(MeshRecord::OpType::OP_PREPROCESSING, tmp_time, m_vertices.size(), m_faces.size()));
            sum_time += tmp_time;
            cout << "time = " << tmp_time << "s" << endl;
            cout << endl;

            //delaunay tetrahedralization
            igl_timer.start();
            cout << "Delaunay tetrahedralizing..." << endl;
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
            cout << "# bsp_vertices = " << bsp_vertices.size() << endl;
            cout << "# bsp_edges = " << bsp_edges.size() << endl;
            cout << "# bsp_faces = " << bsp_faces.size() << endl;
            cout << "# bsp_nodes = " << bsp_nodes.size() << endl;
            cout << "Delaunay tetrahedralization done!" << endl;
            tmp_time = igl_timer.getElapsedTime();
            addRecord(MeshRecord(MeshRecord::OpType::OP_DELAUNEY_TETRA, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
            sum_time += tmp_time;
            cout << "time = " << tmp_time << "s" << endl;
            cout << endl;

            //mesh conforming
            igl_timer.start();
            cout << "Divfaces matching..." << endl;
            MeshConformer MC(m_vertices, m_faces, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes);
            MC.match();
            cout << "Divfaces matching done!" << endl;
            tmp_time = igl_timer.getElapsedTime();
            addRecord(MeshRecord(MeshRecord::OpType::OP_DIVFACE_MATCH, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
            cout << "time = " << tmp_time << "s" << endl;
            cout << endl;

            //bsp subdivision
            igl_timer.start();
            BSPSubdivision BS(MC);
            BS.init();
            BS.subdivideBSPNodes();
            cout << "Output: " << endl;
            cout << "# node = " << MC.bsp_nodes.size() << endl;
            cout << "# face = " << MC.bsp_faces.size() << endl;
            cout << "# edge = " << MC.bsp_edges.size() << endl;
            cout << "# vertex = " << MC.bsp_vertices.size() << endl;
            cout << "BSP subdivision done!" << endl;
            tmp_time = igl_timer.getElapsedTime();
            addRecord(MeshRecord(MeshRecord::OpType::OP_BSP, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
            sum_time += tmp_time;
            cout << "time = " << tmp_time << "s" << endl;
            cout << endl;

            //simple tetrahedralization
            igl_timer.start();
            cout << "Tetrehedralizing ..." << endl;
            SimpleTetrahedralization ST(MC);
            ST.tetra(MR.tet_vertices, MR.tets);
            ST.labelSurface(m_f_tags, raw_e_tags, raw_conn_e4v, MR.tet_vertices, MR.tets, MR.is_surface_fs);
            ST.labelBbox(MR.tet_vertices, MR.tets);
            if (!g_is_close)//if input is an open mesh
                ST.labelBoundary(MR.tet_vertices, MR.tets, MR.is_surface_fs);
            cout << "# tet_vertices = " << MR.tet_vertices.size() << endl;
            cout << "# tets = " << MR.tets.size() << endl;
            cout << "Tetrahedralization done!" << endl;
            tmp_time = igl_timer.getElapsedTime();
            addRecord(MeshRecord(MeshRecord::OpType::OP_SIMPLE_TETRA, tmp_time, MR.tet_vertices.size(), MR.tets.size()));
            sum_time += tmp_time;
            cout << "time = " << tmp_time << "s" << endl;
            cout << endl;

            cout << "Total time for the first stage = " << sum_time << endl;
        }

        /// STAGE 2
        //init
        cout << "Refinement initializing..." << endl;
        MR.prepareData();
        cout << "Refinement intialization done!" << endl;
        cout << endl;

        //improvement
        MR.refine(energy_type);

        outputFinalTetmesh(MR, out_vertices, out_tets); //do winding number and output the tetmesh
    }

    void tetrahedralization(const std::vector<std::array<double, 3>>& in_vertices,
                            const std::vector<std::array<int, 3>>& in_faces,
                            std::vector<std::array<double, 3>>& out_vertices,
                            std::vector<std::array<int, 4>>& out_tets) {
        out_vertices.clear();
        out_tets.clear();

        args.i_epsilon = parameters.i_epsilon;
        args.bg_mesh = parameters.bg_mesh;
        args.filter_energy = parameters.filter_energy;
        args.i_ideal_edge_length = parameters.i_ideal_edge_length;
        args.is_laplacian = parameters.is_laplacian;
        args.is_quiet = parameters.is_quiet;
        args.max_pass = parameters.max_pass;
        args.stage = parameters.stage;
        args.targeted_num_v = parameters.targeted_num_v;

        //initalization
        GEO::initialize();
        g_postfix = args.postfix;
        g_working_dir = args.input.substr(0, args.input.size() - 4);

        if (args.csv_file == "")
            g_stat_file = g_working_dir + g_postfix + ".csv";
        else
            g_stat_file = args.csv_file;

        if (args.output == "")
            g_output_file = g_working_dir + g_postfix + ".msh";
        else
            g_output_file = args.output;

        if (args.is_quiet) {
            args.is_output_csv = false;
            cout.setstate(std::ios_base::failbit);//use std::cout.clear() to get it back
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

        if (args.is_quiet) {
            std::cout.clear();
        }
    }
}
