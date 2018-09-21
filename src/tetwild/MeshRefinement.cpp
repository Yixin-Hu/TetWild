// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 4/11/17.
//

#include <tetwild/MeshRefinement.h>
#include <tetwild/Common.h>
#include <tetwild/Args.h>
#include <tetwild/Logger.h>
#include <tetwild/Serialization.h>
#include <tetwild/EdgeCollapser.h>
#include <tetwild/EdgeSplitter.h>
#include <tetwild/EdgeRemover.h>
#include <tetwild/VertexSmoother.h>
#include <tetwild/DisableWarnings.h>
#include <tetwild/geogram/mesh_AABB.h>
#include <CGAL/centroid.h>
#include <tetwild/EnableWarnings.h>
#include <pymesh/MshLoader.h>
#include <pymesh/MshSaver.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/points/kd_tree.h>
#include <igl/winding_number.h>

namespace tetwild {

void MeshRefinement::prepareData(bool is_init) {
    igl_timer.start();
    if (is_init) {
        t_is_removed = std::vector<bool>(tets.size(), false);//have to
        v_is_removed = std::vector<bool>(tet_vertices.size(), false);
        for (int i = 0; i < tet_vertices.size(); i++) {
            if (tet_vertices[i].is_rounded)
                continue;
            tet_vertices[i].round();
        }
        round();
    }

    GEO::Mesh simple_mesh;
    getSimpleMesh(simple_mesh);
    GEO::MeshFacetsAABBWithEps simple_tree(simple_mesh);
    LocalOperations localOperation(tet_vertices, tets, is_surface_fs, v_is_removed, t_is_removed, tet_qualities,
                                   state.ENERGY_AMIPS, simple_mesh, simple_tree, simple_tree, args, state);
    localOperation.calTetQualities(tets, tet_qualities, true);//cal all measure
    double tmp_time = igl_timer.getElapsedTime();
    logger().debug("{}s", tmp_time);
    localOperation.outputInfo(MeshRecord::OpType::OP_OPT_INIT, tmp_time);
}

void MeshRefinement::round() {
    int cnt = 0;
    int sub_cnt = 0;
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        if (tet_vertices[i].is_rounded)
            continue;
        tet_vertices[i].is_rounded = true;
        Point_3 old_p = tet_vertices[i].pos;
        tet_vertices[i].pos = Point_3(tet_vertices[i].posf[0], tet_vertices[i].posf[1], tet_vertices[i].posf[2]);

        for (auto it = tet_vertices[i].conn_tets.begin(); it != tet_vertices[i].conn_tets.end(); it++) {
            CGAL::Orientation ori;
//            bool is_rounded = true;
//            for (int j = 0; j < 4; j++)
//                if (!tet_vertices[tets[*it][j]].is_rounded) {
//                    is_rounded = false;
//                    break;
//                }
//            if (is_rounded)
//                ori = CGAL::orientation(tet_vertices[tets[*it][0]].posf, tet_vertices[tets[*it][1]].posf,
//                                        tet_vertices[tets[*it][2]].posf, tet_vertices[tets[*it][3]].posf);
//            else
                ori = CGAL::orientation(tet_vertices[tets[*it][0]].pos, tet_vertices[tets[*it][1]].pos,
                                        tet_vertices[tets[*it][2]].pos, tet_vertices[tets[*it][3]].pos);

            if (ori != CGAL::POSITIVE) {
                tet_vertices[i].is_rounded = false;
                break;
            }
        }
        if (!tet_vertices[i].is_rounded)
            tet_vertices[i].pos = old_p;
        else {
            cnt++;
            sub_cnt++;
        }
    }
    logger().debug("round: {}({})", cnt, tet_vertices.size());

    //for check
//    for (int i = 0; i < tets.size(); i++) {
//        if (t_is_removed[i])
//            continue;
//        CGAL::Orientation ori = CGAL::orientation(tet_vertices[tets[i][0]].pos, tet_vertices[tets[i][1]].pos,
//                                                  tet_vertices[tets[i][2]].pos, tet_vertices[tets[i][3]].pos);
//        if (ori != CGAL::POSITIVE) {
//            logger().debug("round hehe");
//        }
//    }
}

void MeshRefinement::clear() {
    tet_vertices.clear();
    tets.clear();

    t_is_removed.clear();
    v_is_removed.clear();
    is_surface_fs.clear();
    tet_qualities.clear();
}

int MeshRefinement::doOperations(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                 VertexSmoother& smoother, const std::array<bool, 4>& ops){
    int cnt0=0;
    for(int i=0;i<tet_vertices.size();i++){
        if(v_is_removed[i] || tet_vertices[i].is_locked || tet_vertices[i].is_rounded)
            continue;
        cnt0++;
    }
    bool is_log = true;
    double tmp_time;

    if (ops[0]) {
        igl_timer.start();
        logger().info("edge splitting...");
        splitter.init();
        splitter.split();
        tmp_time = igl_timer.getElapsedTime();
        splitter.outputInfo(MeshRecord::OpType::OP_SPLIT, tmp_time, is_log);
        logger().info("edge splitting done!");
        logger().info("time = {}s", tmp_time);
    }

    if (ops[1]) {
        igl_timer.start();
        logger().info("edge collapsing...");
        collapser.init();
        collapser.collapse();
        tmp_time = igl_timer.getElapsedTime();
        collapser.outputInfo(MeshRecord::OpType::OP_COLLAPSE, tmp_time, is_log);
        logger().info("edge collasing done!");
        logger().info("time = {}s", tmp_time);
    }

    if (ops[2]) {
        igl_timer.start();
        logger().info("edge removing...");
        edge_remover.init();
        edge_remover.swap();
        tmp_time = igl_timer.getElapsedTime();
        edge_remover.outputInfo(MeshRecord::OpType::OP_SWAP, tmp_time, is_log);
        logger().info("edge removal done!");
        logger().info("time = {}s", tmp_time);
    }

    if (ops[3]) {
        igl_timer.start();
        logger().info("vertex smoothing...");
        smoother.smooth();
        tmp_time = igl_timer.getElapsedTime();
        smoother.outputInfo(MeshRecord::OpType::OP_SMOOTH, tmp_time, is_log);
        logger().info("vertex smooth done!");
        logger().info("time = {}s", tmp_time);
    }

    round();

    int cnt1=0;
    for(int i=0;i<tet_vertices.size();i++){
        if(v_is_removed[i] || tet_vertices[i].is_locked || tet_vertices[i].is_rounded)
            continue;
        cnt1++;
    }
    return cnt0-cnt1;
}

int MeshRefinement::doOperationLoops(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
    VertexSmoother& smoother, int max_pass, const std::array<bool, 4>& ops)
{
    double avg_energy, max_energy;
    splitter.getAvgMaxEnergy(avg_energy, max_energy);

    int loop_cnt = 0;
    for (int i = 0; i < max_pass; i++) {
        doOperations(splitter, collapser, edge_remover, smoother, ops);
        loop_cnt++;

        double tmp_avg_energy, tmp_max_energy;
        splitter.getAvgMaxEnergy(tmp_avg_energy, tmp_max_energy);
        if (std::abs(tmp_avg_energy - avg_energy) < args.delta_energy_thres
            && std::abs(tmp_max_energy - max_energy) < args.delta_energy_thres)
            break;
        avg_energy = tmp_avg_energy;
        max_energy = tmp_max_energy;
    }

    return loop_cnt;
}

void MeshRefinement::refine(int energy_type, const std::array<bool, 4>& ops, bool is_pre, bool is_post, int scalar_update) {
    GEO::MeshFacetsAABBWithEps geo_sf_tree(geo_sf_mesh);
    if (geo_b_mesh.vertices.nb() == 0) {
        getSimpleMesh(geo_b_mesh);//for constructing aabb tree, the mesh cannot be empty
    }
    GEO::MeshFacetsAABBWithEps geo_b_tree(geo_b_mesh);

    if (is_dealing_unrounded)
        min_adaptive_scale = state.eps / state.initial_edge_len * 0.5; //min to eps/2
    else
//        min_adaptive_scale = state.eps_input / state.initial_edge_len; // state.eps_input / state.initial_edge_len * 0.5 is too small
        min_adaptive_scale = (state.bbox_diag / 1000) / state.initial_edge_len; // set min_edge_length to diag / 1000 would be better

    LocalOperations localOperation(tet_vertices, tets, is_surface_fs, v_is_removed, t_is_removed, tet_qualities,
                                   energy_type, geo_sf_mesh, geo_sf_tree, geo_b_tree, args, state);
    EdgeSplitter splitter(localOperation, state.initial_edge_len * (4.0 / 3.0) * state.initial_edge_len * (4.0 / 3.0));
    EdgeCollapser collapser(localOperation, state.initial_edge_len * (4.0 / 5.0) * state.initial_edge_len * (4.0 / 5.0));
    EdgeRemover edge_remover(localOperation, state.initial_edge_len * (4.0 / 3.0) * state.initial_edge_len * (4.0 / 3.0));
    VertexSmoother smoother(localOperation);

    collapser.is_check_quality = true;

    if (args.save_mid_result == 1)
        outputMidResult(false, 1);

//    double old_state.eps = state.eps;
//    state.eps = 0.5 * old_state.eps;
//    state.eps_2 = state.eps * state.eps;

    if (is_pre)
        refine_pre(splitter, collapser, edge_remover, smoother);

    /// apply the local operations
    if (is_dealing_unrounded) {
        for (int i = 0; i < tet_vertices.size(); i++) {
            if (v_is_removed[i] || tet_vertices[i].is_rounded)
                continue;
            smoother.outputOneRing(i, "");
        }
    }

    double avg_energy0, max_energy0;
    localOperation.getAvgMaxEnergy(avg_energy0, max_energy0);
    double target_energy0 = 1e6;
    int update_buget = 2;
    int update_cnt = 0;
    int is_output = true;
//    const double eps_s = 0.8;
//    state.eps *= eps_s;
//    state.eps_2 *= eps_s*eps_s;
    bool is_split = true;
    for (int pass = old_pass; pass < old_pass + args.max_num_passes; pass++) {
        if (is_dealing_unrounded && pass == old_pass) {
            updateScalarField(false, false, args.filter_energy_thres);
        }

        logger().info("//////////////// Pass {} ////////////////", pass);
        if (is_dealing_unrounded)
            collapser.is_limit_length = false;
        doOperations(splitter, collapser, edge_remover, smoother,
                     std::array<bool, 4>({{is_split, ops[1], ops[2], ops[3]}}));
        update_cnt++;

        if (is_dealing_unrounded) {
            bool is_finished = true;
            for (int i = 0; i < tet_vertices.size(); i++) {
                if (v_is_removed[i])
                    continue;
                if (!tet_vertices[i].is_rounded)
                    is_finished = false;
            }
            if (is_finished) {
                logger().debug("all vertices rounded!!");
//                break;
            }
        }

        if (localOperation.getMaxEnergy() < args.filter_energy_thres)
            break;

        //check and mark is_bad_element
        double avg_energy, max_energy;
        localOperation.getAvgMaxEnergy(avg_energy, max_energy);
        if (pass > 0 && pass < old_pass + args.max_num_passes - 1
            && avg_energy0 - avg_energy < args.delta_energy_thres && max_energy0 - max_energy < args.delta_energy_thres) {

//            if (args.target_num_vertices > 0 && getInsideVertexSize() > 1.05 * args.target_num_vertices && isRegionFullyRounded()) {
//                if (state.sub_stage < args.stage) {
//                    state.eps += state.eps_delta;
//                    state.eps_2 = state.eps * state.eps;
//                    state.sub_stage++;
//                    avg_energy0 = avg_energy;
//                    max_energy0 = max_energy;
//                    continue;
//                } else {
//                    is_split = false;
//                    continue;
////                    break;
//                }
//            }
//            is_split = true;

            if (update_cnt == 1) {
                if (is_hit_min) {
                    update_buget--;
                } else
                    continue;
            }
            if (update_buget == 0) {
                if (state.sub_stage > 1 && state.sub_stage < args.stage) {
                    state.eps += state.eps_delta;
                    state.eps_2 = state.eps * state.eps;
                    state.sub_stage++;
                    update_buget = 2;
//                    logger().debug("[[[[[[[[[[[[[[UPDATE EPSILON {}]]]]]]]]]]]]]]]]", state.eps);
                } else
                    break;
            }
            update_cnt = 0;

//            if(is_start_adaptive && !isRegionFullyRounded()) {
//            if(is_hit_min && !isRegionFullyRounded()) {
//                refine_unrounded(splitter, collapser, edge_remover, smoother);
//            }

            //get target energy
            double target_energy = localOperation.getMaxEnergy() / 100;
            target_energy = std::min(target_energy, target_energy0 / 10);
            target_energy = std::max(target_energy, args.filter_energy_thres * 0.8);
            target_energy0 = target_energy;
            updateScalarField(false, false, target_energy);

            if (state.sub_stage == 1 && state.sub_stage < args.stage
                && target_energy < args.filter_energy_thres) {
                state.eps += state.eps_delta;
                state.eps_2 = state.eps * state.eps;
                state.sub_stage++;
//                logger().debug("[[[[[[[[[[[[[[UPDATE EPSILON {}]]]]]]]]]]]]]]]]", state.eps);
            }

            if (is_output && args.save_mid_result == 1) {
                outputMidResult(false, 1.5);
                is_output = false;
            }

//            collapser.is_soft = true;

//            if(is_hit_min) {
//                Eigen::MatrixXd V_tmp;
//                Eigen::MatrixXi F_tmp;
//                getSurface(V_tmp, F_tmp);
//                localOperation.outputSurfaceColormap(V_tmp, F_tmp, old_state.eps);
//            }
        }
        avg_energy0 = avg_energy;
        max_energy0 = max_energy;
    }

    old_pass = old_pass + args.max_num_passes;

//    if (!isRegionFullyRounded()) {
//        refine_unrounded(splitter, collapser, edge_remover, smoother);
//    }
//    if (max_energy0 > 1e3) {
//        refine_local(splitter, collapser, edge_remover, smoother, args.filter_energy_thres);
//    }

//    if (!isRegionFullyRounded() || max_energy0 > 1e3)
//        serialization(state.working_dir + state.postfix_str + ".slz");

    if (!args.is_quiet) {
        double max_e = localOperation.getMaxEnergy();
        if (max_e > 100) {
            bool is_print = false;
            std::ofstream f;
            f.open(state.working_dir + args.postfix + ".tmp");
            for (int i = 0; i < tet_qualities.size(); i++) {
                if (t_is_removed[i])
                    continue;
                if (tet_qualities[i].slim_energy > max_e * 0.9) {
                    is_print = true;
                    f << "tet " << i << ": energy = " << tet_qualities[i].slim_energy << "; ";
                    std::array<double, 6> l;
                    for (int j = 0; j < 3; j++) {
                        l[j * 2] = CGAL::squared_distance(tet_vertices[tets[i][0]].posf,
                                                          tet_vertices[tets[i][j + 1]].posf);
                        l[j * 2 + 1] = CGAL::squared_distance(tet_vertices[tets[i][j + 1]].posf,
                                                              tet_vertices[tets[i][(j + 1) % 3 + 1]].posf);
                    }
                    auto it = std::min_element(l.begin(), l.end());
                    f << "min_el = " << std::sqrt(*it) << "; ";
                    int n = it - l.begin();
                    int v1_id, v2_id;
                    if (n % 2 == 0) {
                        v1_id = 0;
                        v2_id = n / 2 + 1;
                    } else {
                        v1_id = (n - 1) / 2 + 1;
                        v2_id = ((n - 1) / 2 + 1) % 3 + 1;
                    }
                    f << "v1 " << tets[i][v1_id] << " " << tet_vertices[tets[i][v1_id]].is_on_surface << " "
                      << tet_vertices[tets[i][v1_id]].is_on_boundary << " "
                      << localOperation.isPointOutEnvelop(tet_vertices[tets[i][v1_id]].posf) << " "
                      << localOperation.isPointOutBoundaryEnvelop(tet_vertices[tets[i][v1_id]].posf) << "; "

                      << "v2 " << tets[i][v2_id] << " " << tet_vertices[tets[i][v2_id]].is_on_surface << " "
                      << tet_vertices[tets[i][v2_id]].is_on_boundary << " "
                      << localOperation.isPointOutEnvelop(tet_vertices[tets[i][v2_id]].posf) << " "
                      << localOperation.isPointOutBoundaryEnvelop(tet_vertices[tets[i][v2_id]].posf) << std::endl;
                }
            }
            if (is_print)
                f << state.eps << std::endl;
            f.close();
        }
    }

    if (is_post) {
        if (args.target_num_vertices > 0) {
            double n = getInsideVertexSize();
            if (n > args.target_num_vertices) {
                collapser.is_limit_length = false;
                collapser.is_soft = true;
                collapser.soft_energy = localOperation.getMaxEnergy();
                collapser.budget =
                        (n - args.target_num_vertices) * std::count(v_is_removed.begin(), v_is_removed.end(), false) / n *
                        1.5;
            }
        }
        refine_post(splitter, collapser, edge_remover, smoother);
    }


    if (args.target_num_vertices > 0)
        applyTargetedVertexNum(splitter, collapser, edge_remover, smoother);

    if (args.background_mesh != "") {
        applySizingField(splitter, collapser, edge_remover, smoother);
    }

    if (args.save_mid_result == 2)
        outputMidResult(true, 2);//mark in/out


//    if (!args.is_quiet) {
////        Eigen::MatrixXd V_tmp;
////        Eigen::MatrixXi F_tmp;
////        getTrackedSurface(V_tmp, F_tmp);
////        igl::writeOBJ(state.g_working_dir + state.postfix + "_tracked_sf1.obj", V_tmp, F_tmp);
////        getSurface(V_tmp, F_tmp);
////        igl::writeOBJ(state.g_working_dir + state.postfix + "_tracked_sf2.obj", V_tmp, F_tmp);
//////        localOperation.outputSurfaceColormap(V_tmp, F_tmp, state.g_eps_input);//compared with user input eps
////        localOperation.checkUnrounded();
//    }

    if (args.smooth_open_boundary)
        postProcess(smoother);
}

void MeshRefinement::refine_pre(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                VertexSmoother& smoother){
    logger().info("////////////////// Pre-processing //////////////////");
    collapser.is_limit_length = false;
    doOperations(splitter, collapser, edge_remover, smoother, std::array<bool, 4>{{false, true, false, false}});
    collapser.is_limit_length = true;
}

void MeshRefinement::refine_post(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                 VertexSmoother& smoother){
    logger().info("////////////////// Post-processing //////////////////");
    collapser.is_limit_length = true;
    for (int i = 0; i < tet_vertices.size(); i++) {
        tet_vertices[i].adaptive_scale = 1;
    }

    doOperations(splitter, collapser, edge_remover, smoother, std::array<bool, 4>{{false, true, false, false}});
}

void MeshRefinement::refine_local(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                  VertexSmoother& smoother, double target_energy) {
    EdgeSplitter &localOperation = splitter;
    double old_min_adaptive_scale = min_adaptive_scale;
    min_adaptive_scale = state.eps / state.initial_edge_len * 0.5;

    double avg_energy0, max_energy0;
    localOperation.getAvgMaxEnergy(avg_energy0, max_energy0);
    if(target_energy<0) {
        target_energy = max_energy0 / 100;
        target_energy = std::max(target_energy, args.filter_energy_thres);
    }
    updateScalarField(false, true, target_energy * 0.8, true);
    for (int pass = 0; pass < 20; pass++) {
        logger().info("////////////////// Local Pass {} //////////////////", pass);
        doOperations(splitter, collapser, edge_remover, smoother);

        double avg_energy, max_energy;
        localOperation.getAvgMaxEnergy(avg_energy, max_energy);
        if (max_energy < target_energy)
            break;
        avg_energy0 = avg_energy;
        max_energy0 = max_energy;

        if (pass > 0 && pass < args.max_num_passes - 1
            && avg_energy0 - avg_energy < args.delta_energy_thres && max_energy - max_energy0 < args.delta_energy_thres) {
            updateScalarField(false, true, target_energy);
        }
    }
    min_adaptive_scale = old_min_adaptive_scale;
    refine_revert(splitter, collapser, edge_remover, smoother);

    for(int i=0;i<tet_vertices.size();i++)
        tet_vertices[i].is_locked = false;
}

bool MeshRefinement::refine_unrounded(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                      VertexSmoother& smoother) {
    EdgeSplitter &localOperation = splitter;
    int scalar_update = 3;
    double old_min_adaptive_scale = min_adaptive_scale;
    min_adaptive_scale = state.eps / state.initial_edge_len * 0.5;

    collapser.is_limit_length = false;
    updateScalarField(true, false, -1, true);
    for (int pass = 0; pass < 5 * scalar_update; pass++) {
        logger().info("////////////////// Local Pass {} //////////////////", pass);
        doOperations(splitter, collapser, edge_remover, smoother);

        if (isRegionFullyRounded())
            break;

        if (scalar_update > 0 && pass % scalar_update == scalar_update - 1 && pass < args.max_num_passes * scalar_update - 1) {
            updateScalarField(true, false, -1);
        }
    }
    collapser.is_limit_length = true;
    min_adaptive_scale = old_min_adaptive_scale;
    refine_revert(splitter, collapser, edge_remover, smoother);

    for(int i=0;i<tet_vertices.size();i++)
        tet_vertices[i].is_locked = false;

    return false;
}

void MeshRefinement::refine_revert(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                   VertexSmoother& smoother){
    EdgeSplitter &localOperation = splitter;
    collapser.is_limit_length = false;
    collapser.is_soft = true;

    for(int i=0;i<tet_vertices.size();i++) {
        if(!v_is_removed[i] && !tet_vertices[i].is_locked)
            tet_vertices[i].adaptive_scale = 1;
    }

    int n_v0 = std::count(v_is_removed.begin(), v_is_removed.end(), false);
    for (int pass = 0; pass < 10; pass++) {
        logger().info("////////////////// Local (revert) Pass {} //////////////////", pass);
        doOperations(splitter, collapser, edge_remover, smoother, std::array<bool, 4>({{false, true, true, true}}));
//        doOperations(splitter, collapser, edge_remover, smoother);

        int n_v = std::count(v_is_removed.begin(), v_is_removed.end(), false);
        if (n_v0 - n_v < 1) //when number of vertices becomes stable
            break;
        n_v0 = n_v;
    }

    collapser.is_limit_length = true;
    collapser.is_soft = false;
}

int MeshRefinement::getInsideVertexSize(){
    std::vector<bool> tmp_t_is_removed;
    markInOut(tmp_t_is_removed);
    std::unordered_set<int> inside_vs;
    for (int i = 0; i < tets.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            inside_vs.insert(tets[i][j]);
    }
    return inside_vs.size();
}

void MeshRefinement::markInOut(std::vector<bool>& tmp_t_is_removed){
    tmp_t_is_removed = t_is_removed;
    Eigen::MatrixXd C(std::count(tmp_t_is_removed.begin(), tmp_t_is_removed.end(), false), 3);
    int cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        std::vector<Point_3f> vs;
        vs.reserve(4);
        for (int j = 0; j < 4; j++)
            vs.push_back(tet_vertices[tets[i][j]].posf);
        Point_3f p = CGAL::centroid(vs.begin(), vs.end(), CGAL::Dimension_tag<0>());
        for (int j = 0; j < 3; j++)
            C(cnt, j) = p[j];
        cnt++;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    getSurface(V, F);
    Eigen::VectorXd W;
    logger().debug("winding number...");
    igl::winding_number(V, F, C, W);
    logger().debug("winding number done");

    cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        tmp_t_is_removed[i] = !(W(cnt) > 0.5);
        cnt++;
    }
}

void MeshRefinement::applySizingField(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                      VertexSmoother& smoother) {
    PyMesh::MshLoader mshLoader(args.background_mesh);
    Eigen::VectorXd V_in = mshLoader.get_nodes();
    Eigen::VectorXi T_in = mshLoader.get_elements();
    Eigen::VectorXd values = mshLoader.get_node_field("values");
    if (V_in.rows() == 0 || T_in.rows() == 0 || values.rows() == 0)
        return;

    logger().debug("Applying sizing field...");

    GEO::Mesh bg_mesh;
    bg_mesh.vertices.clear();
    bg_mesh.vertices.create_vertices((int) V_in.rows() / 3);
    for (int i = 0; i < V_in.rows() / 3; i++) {
        GEO::vec3 &p = bg_mesh.vertices.point(i);
		for (int j = 0; j < 3; j++)
			p[j] = V_in(i * 3 + j);
    }
    bg_mesh.cells.clear();
	bg_mesh.cells.create_tets((int) T_in.rows() / 4);
    for (int i = 0; i < T_in.rows() / 4; i++) {
        for (int j = 0; j < 4; j++)
            bg_mesh.cells.set_vertex(i, j, T_in(i * 4 + j));
    }

    // background_mesh.cells.compute_borders();
    // background_mesh.cells.connect();

    GEO::MeshCellsAABB bg_aabb(bg_mesh, false);
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        GEO::vec3 p(tet_vertices[i].posf[0], tet_vertices[i].posf[1], tet_vertices[i].posf[2]);
        int bg_t_id = bg_aabb.containing_tet(p);
        if (bg_t_id == GEO::MeshCellsAABB::NO_TET)
            continue;

        //compute barycenter
        double value = 0;
        std::array<Point_3f, 4> vs;
        for (int j = 0; j < 4; j++)
            vs[j] = Point_3f(V_in(T_in(bg_t_id * 4 + j) * 3), V_in(T_in(bg_t_id * 4 + j) * 3 + 1),
                             V_in(T_in(bg_t_id * 4 + j) * 3 + 2));

        std::array<double, 4> weights;
        for (int j = 0; j < 4; j++) {
            Plane_3f pln(vs[j], vs[(j + 1) % 4], vs[(j + 2) % 4]);
            double weight = std::sqrt(
                    CGAL::squared_distance(tet_vertices[i].posf, pln) / CGAL::squared_distance(vs[(j + 3) % 4], pln));
            weights[j] = weight;
            value += weight * values(T_in(bg_t_id * 4 + (j + 3) % 4));
        }

        tet_vertices[i].adaptive_scale = value / state.initial_edge_len; //we allow .adaptive_scale > 1
    }

//     for debugging
    outputMidResult(true, -1);

    //do more refinement
    collapser.is_limit_length = true;
    collapser.is_soft = true;
    collapser.soft_energy = splitter.getMaxEnergy();
//    state.is_print_tmp = true;//debugging splitting
    doOperationLoops(splitter, collapser, edge_remover, smoother, 20);
}

void MeshRefinement::applyTargetedVertexNum(EdgeSplitter& splitter, EdgeCollapser& collapser, EdgeRemover& edge_remover,
                                            VertexSmoother& smoother) {
    if (args.target_num_vertices < 0)
        return;
    if (args.target_num_vertices == 0)
        for (int i = 0; i < t_is_removed.size(); i++)
            t_is_removed[i] = true;

    double N = args.target_num_vertices; //targeted #v

    //marking in/out
    std::vector<bool> tmp_t_is_removed;
    markInOut(tmp_t_is_removed);

    for (int i = 0; i < tet_vertices.size(); i++)
        tet_vertices[i].is_locked = true;

    for (int i = 0; i < tets.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            tet_vertices[tets[i][j]].is_locked = false;
    }

    int cnt = 0;
    for (int i = 0; i < tet_vertices.size(); i++)
        if (!v_is_removed[i] && !tet_vertices[i].is_locked)
            cnt++;

    const double size_threshold = 0.05;
    if (std::abs(cnt - N) / N < size_threshold)
        return;

    logger().debug("{} -> target {}", cnt, N);

    if (cnt > N) {//reduce vertices
        double max_energy = splitter.getMaxEnergy();
        for (int i = 0; i < tet_vertices.size(); i++)
            tet_vertices[i].adaptive_scale = 10;

        collapser.is_soft = true;
        collapser.soft_energy = max_energy;
//        state.eps *= 1.5;
//        state.eps_2 *= 1.5 * 1.5;
//        for (int i = 0; i < tet_vertices.size(); i++)
//            tet_vertices[i].is_locked = false;

        collapser.budget = cnt - N;
        for (int pass = 0; pass < 10; pass++) {
            doOperations(splitter, collapser, edge_remover, smoother, std::array<bool, 4>({{false, true, false, false}}));
            doOperationLoops(splitter, collapser, edge_remover, smoother, 5, std::array<bool, 4>({{false, false, true, true}}));
            if (collapser.budget / N < size_threshold)
                break;
//            collapser.soft_energy *= 1.5;
        }
    } else {//increase vertices
        for (int i = 0; i < tet_vertices.size(); i++)
            tet_vertices[i].adaptive_scale = 0;

        splitter.budget = N - cnt;
        while(splitter.budget / N >= size_threshold) {
            doOperations(splitter, collapser, edge_remover, smoother, std::array<bool, 4>({{true, false, false, false}}));
            doOperationLoops(splitter, collapser, edge_remover, smoother, 5, std::array<bool, 4>({{false, false, true, true}}));
            splitter.budget = N - getInsideVertexSize();
        }
    }
}

bool MeshRefinement::isRegionFullyRounded(){
    for(int i=0;i<tet_vertices.size();i++){
        if(v_is_removed[i] || tet_vertices[i].is_locked)
            continue;
        if(!tet_vertices[i].is_rounded)
            return false;
    }
    return true;
}

void MeshRefinement::updateScalarField(bool is_clean_up_unrounded, bool is_clean_up_local, double filter_energy, bool is_lock)
{
    // Whenever the mesh energy cannot be optimized too much (delta of avg and
    // max energy is < `delta_energy_thres`), we update the scalar field of the
    // target edge length. The update is performed as follows:
    // - Every vertex whose incident tets have a total energy below a given
    //   threshold is selected.
    // - For each selected vertex, place a ball around it (see code below for the
    //   radius).
    // - The scalar field is * `adaptive_scalar` (0.6 by default) at the center
    //   of the ball, left untouched at its boundary, and linearly interpolated
    //   in-between.

    igl_timer.start();
    logger().debug("marking adaptive scales...");
    double tmp_time = 0;

    double radius0 = state.initial_edge_len * 1.8;//increasing the radius would increase the #v in output
    if(is_hit_min)
        radius0 *= 2;
    if (is_clean_up_local)
        radius0 = state.initial_edge_len;
    if (is_clean_up_unrounded)
        radius0 *= 2;

    logger().debug("filter_energy_thres = {}", filter_energy);
    std::vector<double> adap_tmp(tet_vertices.size(), 1.5);
    double dynamic_adaptive_scale = args.adaptive_scalar;

    const int N = -int(std::log2(min_adaptive_scale) - 1);
    std::vector<std::vector<int>> v_ids(N, std::vector<int>());
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i] || tet_vertices[i].is_locked)
            continue;

        if (is_clean_up_unrounded) {
            if (tet_vertices[i].is_rounded)
                continue;
        } else {
            bool is_refine = false;
            for (int t_id: tet_vertices[i].conn_tets) {
                if (tet_qualities[t_id].slim_energy > filter_energy)
                    is_refine = true;
            }
            if (!is_refine)
                continue;
        }

        int n = -int(std::log2(tet_vertices[i].adaptive_scale) - 0.5);
        if (n >= N)
            n = N - 1;
        v_ids[n].push_back(i);
    }

    for (int n = 0; n < N; n++) {
        if(v_ids[n].size() == 0)
            continue;

        double radius = radius0 / std::pow(2, n);
//        double radius = radius0 / 1.5;

        std::unordered_set<int> is_visited;
        std::queue<int> v_queue;

        std::vector<double> pts;
        pts.reserve(v_ids[n].size() * 3);
        for (int i = 0; i < v_ids[n].size(); i++) {
            for (int j = 0; j < 3; j++)
                pts.push_back(tet_vertices[v_ids[n][i]].posf[j]);

            v_queue.push(v_ids[n][i]);
            is_visited.insert(v_ids[n][i]);
            adap_tmp[v_ids[n][i]] = dynamic_adaptive_scale;
        }
        // construct the kdtree
        GEO::NearestNeighborSearch_var nnsearch = GEO::NearestNeighborSearch::create(3, "BNN");
        nnsearch->set_points(int(v_ids[n].size()), pts.data());

        while (!v_queue.empty()) {
            int v_id = v_queue.front();
            v_queue.pop();

            for (int t_id:tet_vertices[v_id].conn_tets) {
                for (int k = 0; k < 4; k++) {
                    if (is_visited.find(tets[t_id][k]) != is_visited.end())
                        continue;
                    GEO::index_t _;
                    double sq_dist;
                    const double p[3] = {tet_vertices[tets[t_id][k]].posf[0], tet_vertices[tets[t_id][k]].posf[1],
                                         tet_vertices[tets[t_id][k]].posf[2]};
                    nnsearch->get_nearest_neighbors(1, p, &_, &sq_dist);
                    double dis = sqrt(sq_dist);

                    if (dis < radius && !tet_vertices[tets[t_id][k]].is_locked) {
                        v_queue.push(tets[t_id][k]);
                        double new_ss =
                                (dis / radius) * (1 - dynamic_adaptive_scale) + dynamic_adaptive_scale;
                        if (new_ss < adap_tmp[tets[t_id][k]])
                            adap_tmp[tets[t_id][k]] = new_ss;
                    }
                    is_visited.insert(tets[t_id][k]);
                }
            }
        }
    }

    // update scalars
    int cnt = 0;
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        if (is_clean_up_unrounded && is_lock && adap_tmp[i] > 1) {
            tet_vertices[i].is_locked = true;
            cnt++;
        }
        double new_scale = tet_vertices[i].adaptive_scale * adap_tmp[i];
        if (new_scale > 1)
            tet_vertices[i].adaptive_scale = 1;
        else if (new_scale < min_adaptive_scale) {
            if (!is_clean_up_unrounded)
                is_hit_min = true;
            tet_vertices[i].adaptive_scale = min_adaptive_scale;
        } else
            tet_vertices[i].adaptive_scale = new_scale;
    }
    if (is_clean_up_unrounded && is_lock)
        logger().debug("{} vertices locked", cnt);

    logger().debug("marked!");
    tmp_time = igl_timer.getElapsedTime();
    logger().debug("time = {}s", tmp_time);
    addRecord(MeshRecord(MeshRecord::OpType::OP_ADAP_UPDATE, tmp_time, -1, -1), args, state);
//    outputMidResult(true);
}

void MeshRefinement::getSimpleMesh(GEO::Mesh& mesh){
    mesh.vertices.clear();
    mesh.vertices.create_vertices(1);
    mesh.vertices.point(0)=GEO::vec3(0,0,0);

    mesh.facets.clear();
    mesh.facets.create_triangles(1);
    mesh.facets.set_vertex(0, 0, 0);
    mesh.facets.set_vertex(0, 1, 0);
    mesh.facets.set_vertex(0, 2, 0);

    mesh.facets.compute_borders();//for what??
}

void MeshRefinement::postProcess(VertexSmoother& smoother) {
    igl_timer.start();

    std::vector<bool> tmp_t_is_removed;
    markInOut(tmp_t_is_removed);

    //get final surface and do smoothing
    std::vector<int> b_v_ids;
    std::vector<bool> tmp_is_on_surface(tet_vertices.size(), false);
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        if (tet_vertices[i].is_on_bbox)
            continue;
        bool has_removed = false;
        bool has_unremoved = false;
        for (int t_id:tet_vertices[i].conn_tets) {
            if (tmp_t_is_removed[t_id]) {
                has_removed = true;
                if (has_unremoved)
                    break;
            }
            if (!tmp_t_is_removed[t_id]) {
                has_unremoved = true;
                if (has_removed)
                    break;
            }
        }
        if (!has_removed || !has_unremoved)
            continue;
        tmp_is_on_surface[i] = true;
//        if(tet_vertices[i].is_on_boundary && !tet_vertices[i].is_on_surface)
//            pausee();
        if (!tet_vertices[i].is_on_surface)
            b_v_ids.push_back(i);
    }
    logger().debug("tmp_is_on_surface.size = {}", tmp_is_on_surface.size());
    logger().debug("b_v_ids.size = {}", b_v_ids.size());
    for (int i = 0; i < 20; i++) {
        if (smoother.laplacianBoundary(b_v_ids, tmp_is_on_surface, tmp_t_is_removed) == 0) {
            break;
        }
    }
    smoother.outputInfo(MeshRecord::OpType::OP_SMOOTH, igl_timer.getElapsedTime());

    t_is_removed = tmp_t_is_removed;
}

/////just for check

void MeshRefinement::check() {
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        for (auto it = tet_vertices[i].conn_tets.begin(); it != tet_vertices[i].conn_tets.end(); it++)
            if (t_is_removed[*it])
                logger().debug("v {} should have conn_tet t{}", i, *it);
        if (tet_vertices[i].conn_tets.size() == 0) {
            logger().debug("empty conn_tets: v {}", i);
        }
    }

    std::vector<std::array<int, 3>> tet_faces;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;

        CGAL::Orientation ori = CGAL::orientation(tet_vertices[tets[i][0]].pos,
                                                  tet_vertices[tets[i][1]].pos,
                                                  tet_vertices[tets[i][2]].pos,
                                                  tet_vertices[tets[i][3]].pos);

        if (ori == CGAL::COPLANAR) {
            logger().debug("tet {} is degenerate!", i);
        } else if (ori == CGAL::NEGATIVE) {
            logger().debug("tet {} is flipped!", i);
        }

        for (int j = 0; j < 4; j++) {
            bool is_found = false;
            for (auto it = tet_vertices[tets[i][j]].conn_tets.begin();
                 it != tet_vertices[tets[i][j]].conn_tets.end(); it++) {
                if (*it == i) {
                    is_found = true;
                }
                if (t_is_removed[*it])
                    logger().debug("tet {} is removed!", *it);
            }
            if (!is_found) {
                logger().debug("{} {} {} {}", tets[i][0], tets[i][1], tets[i][2], tets[i][3]);
                logger().debug("tet {} should be conn to v {}", i, tets[i][j]);
            }

            std::array<int, 3> f = {{tets[i][j], tets[i][(j + 1) % 4], tets[i][(j + 2) % 4]}};
            std::sort(f.begin(), f.end());
            tet_faces.push_back(f);
        }
    }
    std::sort(tet_faces.begin(), tet_faces.end());
    tet_faces.erase(std::unique(tet_faces.begin(), tet_faces.end()), tet_faces.end());

    for (int i = 0; i < tet_faces.size(); i++) {
        std::unordered_set<int> tmp;
        setIntersection(tet_vertices[tet_faces[i][0]].conn_tets, tet_vertices[tet_faces[i][1]].conn_tets, tmp);
        setIntersection(tet_vertices[tet_faces[i][2]].conn_tets, tmp, tmp);

        if (tmp.size() != 1 && tmp.size() != 2)
            logger().debug("{}", tmp.size());
    }
}

void MeshRefinement::outputMidResult(bool is_with_bbox, double id) {
    std::vector<bool> tmp_t_is_removed = t_is_removed;
    Eigen::VectorXd in_out(std::count(t_is_removed.begin(), t_is_removed.end(), false));
//    if (!is_with_bbox) {
        Eigen::MatrixXd C(std::count(tmp_t_is_removed.begin(), tmp_t_is_removed.end(), false), 3);
        int cnt = 0;
        for (int i = 0; i < tets.size(); i++) {
            if (tmp_t_is_removed[i])
                continue;
            std::vector<Point_3f> vs;
            vs.reserve(4);
            for (int j = 0; j < 4; j++)
                vs.push_back(tet_vertices[tets[i][j]].posf);
            Point_3f p = CGAL::centroid(vs.begin(), vs.end(), CGAL::Dimension_tag<0>());
            for (int j = 0; j < 3; j++)
                C(cnt, j) = p[j];
            cnt++;
        }

        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        getSurface(V, F);
        Eigen::VectorXd W;
        logger().debug("winding number...");
        igl::winding_number(V, F, C, W);
        logger().debug("winding number done");

        cnt = 0;
        for (int i = 0; i < tets.size(); i++) {
            if (tmp_t_is_removed[i])
                continue;
            if(!is_with_bbox)
                tmp_t_is_removed[i] = !(W(cnt) > 0.5);
            else
                in_out(cnt) = !(W(cnt) > 0.5);
            cnt++;
        }
//    }

    int t_cnt = 0;
    std::vector<int> v_ids;
    for (int i = 0; i < tets.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            v_ids.push_back(tets[i][j]);
        t_cnt++;
    }
    std::sort(v_ids.begin(), v_ids.end());
    v_ids.erase(std::unique(v_ids.begin(), v_ids.end()), v_ids.end());
    std::unordered_map<int, int> map_ids;
    for (int i = 0; i < v_ids.size(); i++)
        map_ids[v_ids[i]] = i;

    PyMesh::MshSaver mSaver(state.working_dir + state.postfix + "_mid" + std::to_string(id) + ".msh", true);
    Eigen::VectorXd oV(v_ids.size() * 3);
    Eigen::VectorXi oT(t_cnt * 4);
    for (int i = 0; i < v_ids.size(); i++) {
        for (int j = 0; j < 3; j++)
            oV(i * 3 + j) = tet_vertices[v_ids[i]].posf[j];
    }
//    int cnt = 0;
    cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            oT(cnt * 4 + j) = map_ids[tets[i][j]];
        cnt++;
    }
    mSaver.save_mesh(oV, oT, 3, mSaver.TET);
    logger().debug("#v = {}", oV.rows() / 3);
    logger().debug("#t = {}", oT.rows() / 4);

    Eigen::VectorXd angle(t_cnt);
    Eigen::VectorXd energy(t_cnt);
    cnt = 0;
    for (int i = 0; i < tet_qualities.size(); i++) {
        if (tmp_t_is_removed[i])
            continue;
        angle(cnt) = tet_qualities[i].min_d_angle;
        energy(cnt) = tet_qualities[i].slim_energy;
        cnt++;
    }
    mSaver.save_elem_scalar_field("min_dihedral_angle", angle);
    mSaver.save_elem_scalar_field("energy", energy);
    if(is_with_bbox)
        mSaver.save_elem_scalar_field("in/out", in_out);

    // for debugging
    Eigen::VectorXd scalar(v_ids.size());
    cnt = 0;
    for (int i = 0; i < v_ids.size(); i++) {
        scalar(cnt) = tet_vertices[v_ids[i]].adaptive_scale;
        cnt++;
    }
    mSaver.save_scalar_field("scalar field", scalar);
}

void MeshRefinement::getSurface(Eigen::MatrixXd& V, Eigen::MatrixXi& F){
    std::vector<std::array<int, 3>> fs;
    std::vector<int> vs;
    for(int i=0;i<tets.size();i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++) {
            if (is_surface_fs[i][j] != state.NOT_SURFACE && is_surface_fs[i][j] > 0) {//outside
                std::array<int, 3> v_ids = {{tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4]}};
                if (CGAL::orientation(tet_vertices[v_ids[0]].pos, tet_vertices[v_ids[1]].pos,
                                      tet_vertices[v_ids[2]].pos, tet_vertices[tets[i][j]].pos) != CGAL::POSITIVE) {
                    int tmp = v_ids[0];
                    v_ids[0] = v_ids[2];
                    v_ids[2] = tmp;
                }
                for (int k = 0; k < is_surface_fs[i][j]; k++)
                    fs.push_back(v_ids);
                for (int k = 0; k < 3; k++)
                    vs.push_back(v_ids[k]);
            }
        }
    }
    std::sort(vs.begin(), vs.end());
    vs.erase(std::unique(vs.begin(), vs.end()), vs.end());

    V.resize(vs.size(), 3);
    std::map<int, int> map_ids;
    for(int i=0;i<vs.size();i++){
        map_ids[vs[i]]=i;
        for(int j=0;j<3;j++)
            V(i, j)=tet_vertices[vs[i]].posf[j];
    }

    F.resize(fs.size(), 3);
    for(int i=0;i<fs.size();i++){
        for(int j=0;j<3;j++)
            F(i, j)=map_ids[fs[i][j]];
    }
}

void MeshRefinement::getTrackedSurface(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    std::vector<std::array<int, 6>> fs;
    std::vector<int> vs;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++) {
            if (is_surface_fs[i][j] != state.NOT_SURFACE && is_surface_fs[i][j] >= 0) {//outside
                std::array<int, 3> v_ids = {{tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4]}};
                if (CGAL::orientation(tet_vertices[v_ids[0]].pos, tet_vertices[v_ids[1]].pos,
                                      tet_vertices[v_ids[2]].pos, tet_vertices[tets[i][j]].pos) != CGAL::POSITIVE) {
                    int tmp = v_ids[0];
                    v_ids[0] = v_ids[2];
                    v_ids[2] = tmp;
                }
                std::array<int, 3> v_ids1 = v_ids;
                std::sort(v_ids1.begin(), v_ids1.end());
                fs.push_back(std::array<int, 6>({{v_ids1[0], v_ids1[1], v_ids1[2], v_ids[0], v_ids[1], v_ids[2]}}));
                for (int k = 0; k < 3; k++)
                    vs.push_back(v_ids[k]);
            }
        }
    }
    std::sort(vs.begin(), vs.end());
    vs.erase(std::unique(vs.begin(), vs.end()), vs.end());

    V.resize(vs.size(), 3);
    std::unordered_map<int, int> map_ids;
    for (int i = 0; i < vs.size(); i++) {
        map_ids[vs[i]] = i;
        for (int j = 0; j < 3; j++)
            V(i, j) = tet_vertices[vs[i]].posf[j];
    }

    F.resize(fs.size(), 3);
    for (int i = 0; i < fs.size(); i++)
        for (int j = 0; j < 3; j++) {
            F(i, j) = map_ids[fs[i][j + 3]];
    }

    return;

    std::sort(fs.begin(), fs.end());
    int nf = 0;
    for (int i = 0; i < fs.size(); i++) {
        if (i > 0 && fs[i][0] == fs[i - 1][0] && fs[i][1] == fs[i - 1][1] && fs[i][2] == fs[i - 1][2])
            continue;
        nf++;
    }

    F.resize(nf, 3);
    int cnt = 0;
    for (int i = 0; i < fs.size(); i++) {
        if (i > 0 && fs[i][0] == fs[i - 1][0] && fs[i][1] == fs[i - 1][1] && fs[i][2] == fs[i - 1][2])
            continue;
        for (int j = 0; j < 3; j++) {
            F(cnt, j) = map_ids[fs[i][j + 3]];
            cnt++;
        }
    }
}

namespace {

bool getSurfaceMesh(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in, GEO::Mesh& geo_sf_mesh){
    geo_sf_mesh.vertices.clear();
    geo_sf_mesh.vertices.create_vertices((int) V_in.rows());
    for (int i = 0; i < V_in.rows(); i++) {
        GEO::vec3 &p = geo_sf_mesh.vertices.point(i);
        for (int j = 0; j < 3; j++) {
            p[j] = V_in(i, j);
        }
    }
    geo_sf_mesh.facets.clear();
    geo_sf_mesh.facets.create_triangles((int) F_in.rows());
    for (int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            geo_sf_mesh.facets.set_vertex(i, j, F_in(i, j));
        }
    }
    geo_sf_mesh.facets.compute_borders();
    return true;
}

void getBoundaryMesh(const Eigen::MatrixXd& V_sf, const Eigen::MatrixXi& F_sf, GEO::Mesh& b_mesh){
    std::vector<std::vector<int>> conn_f4v(V_sf.rows(), std::vector<int>());
    for (int i = 0; i < F_sf.rows(); i++) {
        for (int j = 0; j < 3; j++)
            conn_f4v[F_sf(i, j)].push_back(i);
    }

    std::vector<std::array<int, 2>> b_edges;
    for (int i = 0; i < F_sf.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            std::vector<int> tmp;
            std::set_intersection(conn_f4v[F_sf(i, j)].begin(), conn_f4v[F_sf(i, j)].end(),
                                  conn_f4v[F_sf(i, (j + 1) % 3)].begin(), conn_f4v[F_sf(i, (j + 1) % 3)].end(),
                                  std::back_inserter(tmp));
            if (tmp.size() == 1)
                b_edges.push_back(std::array<int, 2>({{F_sf(i, j), F_sf(i, (j + 1) % 3)}}));
        }
    }

    if(b_edges.size()==0){
        b_mesh.vertices.clear();
        return;
    }

    std::unordered_set<int> v_ids;
    for (int i = 0; i < b_edges.size(); i++) {
        v_ids.insert(b_edges[i][0]);
        v_ids.insert(b_edges[i][1]);
    }
    std::unordered_map<int, int> v_ids_map;
    int cnt = 0;
    for (int v_id : v_ids) {
        v_ids_map[v_id] = cnt;
        cnt++;
    }

    b_mesh.vertices.clear();
    b_mesh.vertices.create_vertices((int) v_ids.size());
    for (int v_id : v_ids) {
        GEO::vec3 &p = b_mesh.vertices.point(v_ids_map[v_id]);
        for (int j = 0; j < 3; j++)
            p[j] = V_sf(v_id, j);
    }
    b_mesh.facets.clear();
    b_mesh.facets.create_triangles((int) b_edges.size());
    for (int i = 0; i < b_edges.size(); i++) {
        b_mesh.facets.set_vertex(i, 0, v_ids_map[b_edges[i][0]]);
        b_mesh.facets.set_vertex(i, 1, v_ids_map[b_edges[i][1]]);
        b_mesh.facets.set_vertex(i, 2, v_ids_map[b_edges[i][1]]);
    }
}

} // anonymous namespace

bool MeshRefinement::deserialization(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in,
    const std::string& slz_file)
{
    logger().debug("deserializing ...");

    //process sf_file
    if(!getSurfaceMesh(V_in, F_in, geo_sf_mesh))
        return false;
    getBoundaryMesh(V_in, F_in, geo_b_mesh);
    state.is_mesh_closed = (geo_b_mesh.vertices.nb() == 0);

    //deserialization
    igl::deserialize(state.bbox_diag, "bbox_diag", slz_file);
    igl::deserialize(state.eps, "eps", slz_file);
    igl::deserialize(state.eps_2, "eps_2", slz_file);
    igl::deserialize(state.sampling_dist, "sampling_dist", slz_file);
    igl::deserialize(state.initial_edge_len, "initial_edge_len", slz_file);
    // igl::deserialize(state.NOT_SURFACE, "NOT_SURFACE", slz_file);
    igl::deserialize(old_pass, "old_pass", slz_file);

    igl::deserialize(tet_vertices, "tet_vertices", slz_file);
    igl::deserialize(tets, "tets", slz_file);
    igl::deserialize(is_surface_fs, "is_surface_fs", slz_file);

    t_is_removed = std::vector<bool>(tets.size(), false);
    v_is_removed = std::vector<bool>(tet_vertices.size(), false);
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            tet_vertices[tets[i][j]].conn_tets.insert(i);
        }
    }

    prepareData(false);

//    for (int i = 0; i < tets.size(); i++) {
//        for (int j = 0; j < 4; j++)
//            tet_vertices[tets[i][j]].conn_tets.insert(i);
//    }
    logger().debug("deserialization done");

    return true;
}

void MeshRefinement::serialization(const std::string& slz_file) {
    logger().debug("serializing ...");
    igl::serialize(state.bbox_diag, "bbox_diag", slz_file, true);
    igl::serialize(state.eps, "eps", slz_file);
    igl::serialize(state.eps_2, "eps_2", slz_file);
    igl::serialize(state.sampling_dist, "sampling_dist", slz_file);
    igl::serialize(state.initial_edge_len, "initial_edge_len", slz_file);
    // igl::serialize(state.NOT_SURFACE, "NOT_SURFACE", slz_file);
    igl::serialize(old_pass, "old_pass", slz_file);

    //relabel
    std::vector<int> new_v_ids(tet_vertices.size(), -1);
    int cnt = 0;
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        new_v_ids[i] = cnt++;
    }

    // serialize
    std::vector <TetVertex> slz_tet_vertices;
    slz_tet_vertices.reserve(std::count(v_is_removed.begin(), v_is_removed.end(), false));
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        slz_tet_vertices.push_back(tet_vertices[i]);
        slz_tet_vertices.back().conn_tets.clear();
        slz_tet_vertices.back().on_face.clear();
        slz_tet_vertices.back().on_edge.clear();
    }
    igl::serialize(slz_tet_vertices, "tet_vertices", slz_file);
    slz_tet_vertices.clear();

    std::vector <std::array<int, 4>> slz_tets;
    slz_tets.reserve(std::count(t_is_removed.begin(), t_is_removed.end(), false));
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        slz_tets.push_back(tets[i]);
        for (int j = 0; j < 4; j++)
            slz_tets.back()[j] = new_v_ids[tets[i][j]];
    }
    igl::serialize(slz_tets, "tets", slz_file);
    slz_tets.clear();

    std::vector<std::array<int, 4>> slz_is_surface_fs;
    slz_is_surface_fs.reserve(std::count(t_is_removed.begin(), t_is_removed.end(), false));
    for (int i = 0; i < is_surface_fs.size(); i++) {
        if (t_is_removed[i])
            continue;
        slz_is_surface_fs.push_back(is_surface_fs[i]);
    }
    igl::serialize(slz_is_surface_fs, "is_surface_fs", slz_file);
    slz_is_surface_fs.clear();

    logger().debug("serialization done");
}

} // namespace tetwild
