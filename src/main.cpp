// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "heads.h"
#include "Preprocess.h"
#include "DelaunayTetrahedralization.h"
#include "BSPSubdivision.h"
#include "SimpleTetrahedralization.h"
#include "MeshRefinement.h"
#include "InoutFiltering.h"
#include "CLI11.hpp"

MeshRefinement MR;

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

void outputFinalSurface(MeshRefinement& MR){
    std::vector<TetVertex> &tet_vertices = MR.tet_vertices;
    std::vector<std::array<int, 4>> &tets = MR.tets;
    std::vector<bool> &v_is_removed = MR.v_is_removed;
    std::vector<bool> &t_is_removed = MR.t_is_removed;
    int t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);

    std::vector<std::array<int, 4>> tet_faces;
    tet_faces.reserve(t_cnt * 4);
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {tets[i][j], tets[i][(j + 1) % 4], tets[i][(j + 2) % 4]};
            std::sort(f.begin(), f.end());
            tet_faces.push_back(std::array<int, 4>({f[0], f[1], f[2], i}));
        }
    }
    std::sort(tet_faces.begin(), tet_faces.end());//we can sort all 4 digits!!!

    std::vector<int> sf_tets;
    std::vector<int> sf_vs;
    for (int i = 1; i < tet_faces.size(); i++) {
        if (tet_faces[i][0] == tet_faces[i - 1][0] && tet_faces[i][1] == tet_faces[i - 1][1]
            && tet_faces[i][2] == tet_faces[i - 1][2]) {
            i++;
            continue;
        }
        sf_tets.push_back(i - 1);
        for(int j=0;j<3;j++)
            sf_vs.push_back(tet_faces[i - 1][j]);
    }
    std::sort(sf_vs.begin(), sf_vs.end());
    sf_vs.erase(std::unique(sf_vs.begin(), sf_vs.end()), sf_vs.end());
    std::unordered_map<int, int> map_vs;
    Eigen::MatrixXd V_sf(sf_vs.size(), 3);
    for(int i=0;i<sf_vs.size();i++) {
        map_vs[sf_vs[i]] = i;
        V_sf.row(i) = Eigen::Vector3d(tet_vertices[sf_vs[i]].posf[0], tet_vertices[sf_vs[i]].posf[1],
                                      tet_vertices[sf_vs[i]].posf[2]);
    }

    Eigen::MatrixXi F_sf(sf_tets.size(), 3);
    int I=0;
    for (int i:sf_tets) {
        int t_id = tet_faces[i][3];
        for (int j = 0; j < 4; j++)
            if (tets[t_id][j] != tet_faces[i][0] && tets[t_id][j] != tet_faces[i][1] &&
                tets[t_id][j] != tet_faces[i][2]) {
                if (CGAL::orientation(tet_vertices[tets[t_id][j]].pos, tet_vertices[tet_faces[i][0]].pos,
                                      tet_vertices[tet_faces[i][1]].pos, tet_vertices[tet_faces[i][2]].pos)
                    != CGAL::POSITIVE)
                    std::swap(tet_faces[i][0], tet_faces[i][2]);
                F_sf.row(I) = Eigen::Vector3i(map_vs[tet_faces[i][0]], map_vs[tet_faces[i][1]], map_vs[tet_faces[i][2]]);
                I++;
                break;
            }
    }
    igl::writeOBJ(g_working_dir+g_postfix+"_sf.obj", V_sf, F_sf);
}

void outputFinalTetmesh(MeshRefinement& MR) {
    std::vector<TetVertex> &tet_vertices = MR.tet_vertices;
    std::vector<std::array<int, 4>> &tets = MR.tets;
    std::vector<bool> &v_is_removed = MR.v_is_removed;
    std::vector<bool> &t_is_removed = MR.t_is_removed;
    std::vector<TetQuality> &tet_qualities = MR.tet_qualities;
    int t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
    double tmp_time = 0;
    if (!args.is_laplacian) {
        InoutFiltering IOF(tet_vertices, tets, MR.is_surface_fs, v_is_removed, t_is_removed, tet_qualities);
#ifndef MUTE_COUT
        igl::Timer igl_timer;
        igl_timer.start();
#endif
        IOF.filter();
        t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
#ifndef MUTE_COUT
        tmp_time = igl_timer.getElapsedTime();
        cout << "time = " << tmp_time << "s" << endl;
        cout << t_cnt << " tets inside!" << endl;
#endif
    }

    //output result
    cout<<"Writing mesh to "<<g_output_file<<"..."<<endl;
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

    Eigen::VectorXd oV(v_ids.size() * 3);
    Eigen::VectorXi oT(t_cnt * 4);
    for (int i = 0; i < v_ids.size(); i++) {
        for (int j = 0; j < 3; j++)
            oV(i * 3 + j) = tet_vertices[v_ids[i]].posf[j];
    }
    int cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            oT(cnt * 4 + j) = map_ids[tets[i][j]];
        cnt++;
    }
    cout << "#v = " << oV.rows() / 3 << endl;
    cout << "#t = " << oT.rows() / 4 << endl;

    std::string output_format = g_output_file.substr(g_output_file.size() - 4, 4);
    if (output_format == "mesh") {
        std::fstream f(g_output_file, std::ios::out);
        f.precision(std::numeric_limits<double>::digits10 + 1);
        f << "MeshVersionFormatted 1" << std::endl;
        f << "Dimension 3" << std::endl;

        f << "Vertices" << std::endl << oV.rows() / 3 << std::endl;
        for (int i = 0; i < oV.rows() / 3; i++)
            f << oV(i * 3) << " " << oV(i * 3 + 1) << " " << oV(i * 3 + 2) << " " << 0 << std::endl;
        f << "Triangles" << endl << 0 <<endl;
        f << "Tetrahedra" << endl;
        f << oT.rows() / 4 << std::endl;
        for (int i = 0; i < oT.rows() / 4; i++) {
            for (int j = 0; j < 4; j++)
                f << oT(i * 4 + j) + 1 << " ";
            f << 0 << std::endl;
        }

        f << "End";
        f.close();
    } else {
        PyMesh::MshSaver mSaver(g_output_file, true);
        mSaver.save_mesh(oV, oT, 3, mSaver.TET);
#ifndef MUTE_COUT
        Eigen::VectorXd angle(t_cnt);
        cnt = 0;
        for (int i = 0; i < tet_qualities.size(); i++) {
            if (t_is_removed[i])
                continue;
            angle(cnt) = tet_qualities[i].min_d_angle;
            cnt++;
        }
        mSaver.save_elem_scalar_field("min_dihedral_angle", angle);
#endif
    }

#ifndef MUTE_COUT
    if (args.is_quiet)
        return;
    outputFinalQuality(tmp_time, tet_vertices, tets, t_is_removed, tet_qualities, v_ids);
    outputFinalSurface(MR);
#endif
}

void gtet_new() {
    is_using_energy_max = true;
    is_use_project = false;
    is_using_sampling = true;

    int energy_type = ENERGY_AMIPS;
//    int energy_type = ENERGY_DIRICHLET;
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
//    MeshRefinement MR;
    {/// STAGE 1
        //preprocess
#ifndef MUTE_COUT
        igl_timer.start();
        cout << "Preprocessing..." << endl;
#endif
        Preprocess pp;
        if (!pp.init(MR.geo_b_mesh, MR.geo_sf_mesh)) {
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
        cout << "time = " << tmp_time << "s" << endl;
        cout << endl;
#endif

        //delaunay tetrahedralization
#ifndef MUTE_COUT
        igl_timer.start();
        cout << "Delaunay tetrahedralizing..." << endl;
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
#endif

        //mesh conforming
#ifndef MUTE_COUT
        igl_timer.start();
        cout << "Divfaces matching..." << endl;
#endif
        MeshConformer MC(m_vertices, m_faces, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes);
        MC.match();
#ifndef MUTE_COUT
        cout << "Divfaces matching done!" << endl;
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_DIVFACE_MATCH, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        cout << "time = " << tmp_time << "s" << endl;
        cout << endl;
#endif

        //bsp subdivision
#ifndef MUTE_COUT
        igl_timer.start();
        cout << "BSP subdivision ..." << endl;
#endif
        BSPSubdivision BS(MC);
        BS.init();
        BS.subdivideBSPNodes();
#ifndef MUTE_COUT
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
#endif

        //simple tetrahedralization
#ifndef MUTE_COUT
        igl_timer.start();
        cout << "Tetrehedralizing ..." << endl;
#endif
        SimpleTetrahedralization ST(MC);
        ST.tetra(MR.tet_vertices, MR.tets);
        ST.labelSurface(m_f_tags, raw_e_tags, raw_conn_e4v, MR.tet_vertices, MR.tets, MR.is_surface_fs);
        ST.labelBbox(MR.tet_vertices, MR.tets);
        if (!g_is_close)//if input is an open mesh
            ST.labelBoundary(MR.tet_vertices, MR.tets, MR.is_surface_fs);
#ifndef MUTE_COUT
        cout << "# tet_vertices = " << MR.tet_vertices.size() << endl;
        cout << "# tets = " << MR.tets.size() << endl;
        cout << "Tetrahedralization done!" << endl;
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_SIMPLE_TETRA, tmp_time, MR.tet_vertices.size(), MR.tets.size()));
        sum_time += tmp_time;
        cout << "time = " << tmp_time << "s" << endl;
        cout << endl;
        cout << "Total time for the first stage = " << sum_time << endl;
#endif
    }

    /// STAGE 2
    //init
#ifndef MUTE_COUT
    cout << "Refinement initializing..." << endl;
#endif
    MR.prepareData();
#ifndef MUTE_COUT
    cout << "Refinement intialization done!" << endl;
    cout << endl;
#endif

    //improvement
    MR.refine(energy_type);

    outputFinalTetmesh(MR); //do winding number and output the tetmesh
}

void gtet_new_slz(const std::string& sf_file, const std::string& slz_file, int max_pass,
                  const std::array<bool, 4>& ops){

    MeshRefinement MR;
    MR.deserialization(sf_file, slz_file);

//    MR.is_dealing_unrounded = true;
    MR.refine(ENERGY_AMIPS, ops, false, true);

    outputFinalTetmesh(MR);
}

int main(int argc, char *argv[]) {
#ifdef MUTE_COUT
    cout<<"Unnecessary checks are muted."<<endl;
#endif
    CLI::App app{"RobustTetMeshing"};
    app.add_option("--input", args.input, "--input INPUT. Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")->required();
    app.add_option("--postfix", args.postfix, "--postfix P. Postfix P for output files. (string, optinal, default: '_')");
    app.add_option("--output", args.output, "--output OUTPUT. Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')");
    app.add_option("--ideal-edge-length", args.i_ideal_edge_length, "--ideal-edge-length L. ideal_edge_length = diag_of_bbox / L. (double, optional, default: 20)");
    app.add_option("--epsilon", args.i_epsilon, "--epsilon EPS. epsilon = diag_of_bbox / EPS. (double, optional, default: 1000)");
    app.add_option("--stage", args.stage, "--stage STAGE. Run pipeline in stage STAGE. (integer, optional, default: 1)");
    app.add_option("--filter-energy", args.filter_energy, "--filter-energy ENERGY. Stop mesh improvement when the maximum energy is smaller than ENERGY. (double, optional, default: 10)");
    app.add_option("--max-pass", args.max_pass, "--max-pass PASS. Do PASS mesh improvement passes in maximum. (integer, optional, default: 80)");

    app.add_option("--is-laplacian", args.is_laplacian, "--is-laplacian ISLAP. Do Laplacian smoothing for the surface of output on the holes of input, if ISLAP = 1. Otherwise, ISLAP = 0. (integer, optinal, default: 0)");
    app.add_option("--targeted-num-v", args.targeted_num_v, "--targeted-num-v TV. Output tetmesh that contains TV vertices. (integer, optinal, tolerance: 5%)");
    app.add_option("--bg-mesh", args.bg_mesh, "--bg-mesh BGMESH. Background tetmesh BGMESH in .msh format for applying sizing field. (string, optional)");
    app.add_option("--is-quiet", args.is_quiet, "--is-quiet Q. Mute log info and only output tetmesh if Q = 1. (integer, optional, default: 0)");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    //initalization
    GEO::initialize();
    g_postfix = args.postfix;
    if(args.slz_file != "")
        g_working_dir = args.input.substr(0, args.slz_file.size() - 4);
    else
        g_working_dir = args.input.substr(0, args.input.size() - 4);

    if(args.csv_file == "")
        g_stat_file = g_working_dir + g_postfix + ".csv";
    else
        g_stat_file = args.csv_file;

    if(args.output == "")
        g_output_file = g_working_dir + g_postfix + ".msh";
    else
        g_output_file = args.output;

    if(args.is_quiet) {
        args.is_output_csv = false;
        std::streambuf* orig_buf = cout.rdbuf();
        cout.rdbuf(NULL);
//        cout.setstate(std::ios_base::failbit);//use std::cout.clear() to get it back
    }

    //do tetrahedralization
    if(args.slz_file != "")
        gtet_new_slz(args.input, args.slz_file, args.max_pass, std::array<bool, 4>({true, false, true, true}));
    else
        gtet_new();

    return 0;
}
