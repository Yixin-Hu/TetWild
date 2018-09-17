// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <tetwild/Common.h>
#include <tetwild/Preprocess.h>
#include <tetwild/DelaunayTetrahedralization.h>
#include <tetwild/BSPSubdivision.h>
#include <tetwild/SimpleTetrahedralization.h>
#include <tetwild/MeshRefinement.h>
#include <tetwild/InoutFiltering.h>
#include <igl/writeOBJ.h>
#include <pymesh/MshSaver.h>
#include <tetwild/DisableWarnings.h>
#include <CLI/CLI.hpp>
#include <tetwild/EnableWarnings.h>

using namespace tetwild;

MeshRefinement MR;

void outputFinalQuality(double time, const std::vector<TetVertex>& tet_vertices, const std::vector<std::array<int, 4>>& tets,
                        const std::vector<bool> &t_is_removed, const std::vector<TetQuality>& tet_qualities,
                        const std::vector<int>& v_ids) {
    logger().debug("final quality:");
    double min = 10, max = 0;
    double min_avg = 0, max_avg = 0;
//    double max_asp_ratio = 0, avg_asp_ratio = 0;
    double max_slim_energy = 0, avg_slim_energy = 0;
    std::array<double, 6> cmp_cnt = {{0, 0, 0, 0, 0, 0}};
    std::array<double, 6> cmp_d_angles = {{6 / 180.0 * M_PI, 12 / 180.0 * M_PI, 18 / 180.0 * M_PI,
                                           162 / 180.0 * M_PI, 168 / 180.0 * M_PI, 174 / 180.0 * M_PI}};
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
            std::array<int, 3> f = {{tets[i][j], tets[i][(j + 1) % 4], tets[i][(j + 2) % 4]}};
            std::sort(f.begin(), f.end());
            tet_faces.push_back(std::array<int, 4>({{f[0], f[1], f[2], i}}));
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
    igl::writeOBJ(State::state().g_working_dir+State::state().g_postfix+"_sf.obj", V_sf, F_sf);
}

void outputFinalTetmesh(MeshRefinement& MR) {
    std::vector<TetVertex> &tet_vertices = MR.tet_vertices;
    std::vector<std::array<int, 4>> &tets = MR.tets;
    std::vector<bool> &v_is_removed = MR.v_is_removed;
    std::vector<bool> &t_is_removed = MR.t_is_removed;
    std::vector<TetQuality> &tet_qualities = MR.tet_qualities;
    int t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
    double tmp_time = 0;
    if (!GArgs::args().is_laplacian) {
        InoutFiltering IOF(tet_vertices, tets, MR.is_surface_fs, v_is_removed, t_is_removed, tet_qualities);
        igl::Timer igl_timer;
        igl_timer.start();
        IOF.filter();
        t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
        tmp_time = igl_timer.getElapsedTime();
        logger().info("time = {}s", tmp_time);
        logger().debug("{} tets inside!", t_cnt);
    }

    //output result
    logger().debug("Writing mesh to {}...", State::state().g_output_file);
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
    logger().debug("#v = {}", oV.rows() / 3);
    logger().debug("#t = {}", oT.rows() / 4);

    std::string output_format = State::state().g_output_file.substr(State::state().g_output_file.size() - 4, 4);
    if (output_format == "mesh") {
        std::fstream f(State::state().g_output_file, std::ios::out);
        f.precision(std::numeric_limits<double>::digits10 + 1);
        f << "MeshVersionFormatted 1" << std::endl;
        f << "Dimension 3" << std::endl;

        f << "Vertices" << std::endl << oV.rows() / 3 << std::endl;
        for (int i = 0; i < oV.rows() / 3; i++)
            f << oV(i * 3) << " " << oV(i * 3 + 1) << " " << oV(i * 3 + 2) << " " << 0 << std::endl;
        f << "Triangles" << std::endl << 0 <<std::endl;
        f << "Tetrahedra" << std::endl;
        f << oT.rows() / 4 << std::endl;
        for (int i = 0; i < oT.rows() / 4; i++) {
            for (int j = 0; j < 4; j++)
                f << oT(i * 4 + j) + 1 << " ";
            f << 0 << std::endl;
        }

        f << "End";
        f.close();
    } else {
        PyMesh::MshSaver mSaver(State::state().g_output_file, true);
        mSaver.save_mesh(oV, oT, 3, mSaver.TET);
        Eigen::VectorXd angle(t_cnt);
        cnt = 0;
        for (int i = 0; i < tet_qualities.size(); i++) {
            if (t_is_removed[i])
                continue;
            angle(cnt) = tet_qualities[i].min_d_angle;
            cnt++;
        }
        mSaver.save_elem_scalar_field("min_dihedral_angle", angle);
    }

    if (GArgs::args().is_quiet)
        return;
    outputFinalQuality(tmp_time, tet_vertices, tets, t_is_removed, tet_qualities, v_ids);
    outputFinalSurface(MR);
}

void gtet_new() {
    State::state().is_using_energy_max = true;
    State::state().is_use_project = false;
    State::state().is_using_sampling = true;

    int energy_type = State::state().ENERGY_AMIPS;
//    int energy_type = State::state().ENERGY_DIRICHLET;
    bool is_sm_single = true;
    bool is_preprocess = true;
    bool is_check_correctness = false;
    bool is_ec_check_quality = true;

    igl::Timer igl_timer;
    double tmp_time = 0;
    double sum_time = 0;

    ////pipeline
//    MeshRefinement MR;
    {/// STAGE 1
        //preprocess
        igl_timer.start();
        logger().info("Preprocessing...");
        Preprocess pp;
        if (!pp.init(MR.geo_b_mesh, MR.geo_sf_mesh)) {
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
        addRecord(MeshRecord(MeshRecord::OpType::OP_INIT, 0, MR.geo_sf_mesh.vertices.nb(), MR.geo_sf_mesh.facets.nb()));
        std::vector<Point_3> m_vertices;
        std::vector<std::array<int, 3>> m_faces;
        pp.process(MR.geo_sf_mesh, m_vertices, m_faces);
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_PREPROCESSING, tmp_time, m_vertices.size(), m_faces.size()));
        sum_time += tmp_time;
        logger().info("time = {}s", tmp_time);

        //delaunay tetrahedralization
        igl_timer.start();
        logger().info("Delaunay tetrahedralizing...");
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
        logger().debug("# bsp_vertices = {}", bsp_vertices.size());
        logger().debug("# bsp_edges = {}", bsp_edges.size());
        logger().debug("# bsp_faces = {}", bsp_faces.size());
        logger().debug("# bsp_nodes = {}", bsp_nodes.size());
        logger().info("Delaunay tetrahedralization done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_DELAUNEY_TETRA, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        sum_time += tmp_time;
        logger().info("time = {}s", tmp_time);

        //mesh conforming
        igl_timer.start();
        logger().info("Divfaces matching...");
        MeshConformer MC(m_vertices, m_faces, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes);
        MC.match();
        logger().info("Divfaces matching done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_DIVFACE_MATCH, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        logger().info("time = {}s", tmp_time);

        //bsp subdivision
        igl_timer.start();
        logger().info("BSP subdivision ...");
        BSPSubdivision BS(MC);
        BS.init();
        BS.subdivideBSPNodes();
        logger().debug("Output: ");
        logger().debug("# node = {}", MC.bsp_nodes.size());
        logger().debug("# face = {}", MC.bsp_faces.size());
        logger().debug("# edge = {}", MC.bsp_edges.size());
        logger().debug("# vertex = {}", MC.bsp_vertices.size());
        logger().info("BSP subdivision done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_BSP, tmp_time, bsp_vertices.size(), bsp_nodes.size()));
        sum_time += tmp_time;
        logger().info("time = {}s", tmp_time);

        //simple tetrahedralization
        igl_timer.start();
        logger().info("Tetrehedralizing ...");
        SimpleTetrahedralization ST(MC);
        ST.tetra(MR.tet_vertices, MR.tets);
        ST.labelSurface(m_f_tags, raw_e_tags, raw_conn_e4v, MR.tet_vertices, MR.tets, MR.is_surface_fs);
        ST.labelBbox(MR.tet_vertices, MR.tets);
        if (!State::state().g_is_close)//if input is an open mesh
            ST.labelBoundary(MR.tet_vertices, MR.tets, MR.is_surface_fs);
        logger().debug("# tet_vertices = {}", MR.tet_vertices.size());
        logger().debug("# tets = {}", MR.tets.size());
        logger().info("Tetrahedralization done!");
        tmp_time = igl_timer.getElapsedTime();
        addRecord(MeshRecord(MeshRecord::OpType::OP_SIMPLE_TETRA, tmp_time, MR.tet_vertices.size(), MR.tets.size()));
        sum_time += tmp_time;
        logger().info("time = {}s", tmp_time);
        logger().info("Total time for the first stage = {}s", sum_time);
    }

    /// STAGE 2
    //init
    logger().info("Refinement initializing...");
    MR.prepareData();
    logger().info("Refinement intialization done!");

    //improvement
    MR.refine(energy_type);

    outputFinalTetmesh(MR); //do winding number and output the tetmesh
}

void gtet_new_slz(const std::string& sf_file, const std::string& slz_file, int max_pass,
                  const std::array<bool, 4>& ops){

    MeshRefinement MR;
    MR.deserialization(sf_file, slz_file);

//    MR.is_dealing_unrounded = true;
    MR.refine(State::state().ENERGY_AMIPS, ops, false, true);

    outputFinalTetmesh(MR);
}

int main(int argc, char *argv[]) {
#ifdef MUTE_COUT
    logger().debug("Unnecessary checks are muted.");
#endif
    int log_level = 1; // debug
    std::string log_filename = "";

    CLI::App app{"RobustTetMeshing"};
    app.add_option("input,--input", GArgs::args().input, "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")->required();
    app.add_option("output,--output", GArgs::args().output, "Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')");
    app.add_option("--postfix", GArgs::args().postfix, "Postfix P for output files. (string, optional, default: '_')");
    app.add_option("-l,--ideal-edge-length", GArgs::args().i_ideal_edge_length, "ideal_edge_length = diag_of_bbox / L. (double, optional, default: 20)");
    app.add_option("-e,--epsilon", GArgs::args().i_epsilon, "epsilon = diag_of_bbox / EPS. (double, optional, default: 1000)");
    app.add_option("--stage", GArgs::args().stage, "Run pipeline in stage STAGE. (integer, optional, default: 1)");
    app.add_option("--filter-energy", GArgs::args().filter_energy, "Stop mesh improvement when the maximum energy is smaller than ENERGY. (double, optional, default: 10)");
    app.add_option("--max-pass", GArgs::args().max_pass, "Do PASS mesh improvement passes in maximum. (integer, optional, default: 80)");

    app.add_flag("--is-laplacian", GArgs::args().is_laplacian, "Do Laplacian smoothing for the surface of output on the holes of input (optional)");
    app.add_option("--targeted-num-v", GArgs::args().targeted_num_v, "Output tetmesh that contains TV vertices. (integer, optional, tolerance: 5%)");
    app.add_option("--bg-mesh", GArgs::args().bg_mesh, "Background tetmesh BGMESH in .msh format for applying sizing field. (string, optional)");
    app.add_flag("-q,--is-quiet", GArgs::args().is_quiet, "Mute console output. (optional)");
    app.add_option("--log", log_filename, "Log info to given file.");
    app.add_option("--level", log_level, "Log level (0 = most verbose, 6 = off).");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    Logger::init(!GArgs::args().is_quiet, log_filename);
    log_level = std::max(0, std::min(6, log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    // logger().info("this is a test");
    // logger().debug("debug stuff");

    //initalization
    GEO::initialize();
    State::state().g_postfix = GArgs::args().postfix;
    if(GArgs::args().slz_file != "")
        State::state().g_working_dir = GArgs::args().input.substr(0, GArgs::args().slz_file.size() - 4);
    else
        State::state().g_working_dir = GArgs::args().input.substr(0, GArgs::args().input.size() - 4);

    if(GArgs::args().csv_file == "")
        State::state().g_stat_file = State::state().g_working_dir + State::state().g_postfix + ".csv";
    else
        State::state().g_stat_file = GArgs::args().csv_file;

    if(GArgs::args().output == "")
        State::state().g_output_file = State::state().g_working_dir + State::state().g_postfix + ".msh";
    else
        State::state().g_output_file = GArgs::args().output;

    if(GArgs::args().is_quiet) {
        GArgs::args().is_output_csv = false;
        std::streambuf* orig_buf = std::cout.rdbuf();
        std::cout.rdbuf(NULL);
//        std::cout.setstate(std::ios_base::failbit);//use std::std::cout.clear() to get it back
    }

    //do tetrahedralization
    if(GArgs::args().slz_file != "")
        gtet_new_slz(GArgs::args().input, GArgs::args().slz_file, GArgs::args().max_pass, std::array<bool, 4>({{true, false, true, true}}));
    else
        gtet_new();

    spdlog::shutdown();

    return 0;
}
