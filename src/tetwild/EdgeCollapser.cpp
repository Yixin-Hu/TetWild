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

#include <tetwild/EdgeCollapser.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>
#include <igl/Timer.h>

namespace tetwild {

void EdgeCollapser::init() {
    energy_time = 0;

    ////cal dir_edge
    //find all edges
    //check if collapsable 1
    //if yes, insert it into queue
    const unsigned int tets_size = tets.size();
    std::vector<std::array<int, 2>> edges;
    edges.reserve(tets_size*6);
//    for (unsigned int i = 0; i < tets_size; i++) {
//        if (t_is_removed[i])
//            continue;
//        for (int j = 0; j < 4; j++) {
//            std::array<int, 2> e = {{tets[i][j], tets[i][(j + 1) % 4]}};
//            if (e[0] > e[1]) e = {{e[1], e[0]}};
//            if(!isLocked_ui(e))
//                edges.push_back(e);
//        }
//    }
    for (unsigned int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{tets[i][0], tets[i][j + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                edges.push_back(e);
            e = {{tets[i][j + 1], tets[i][(j + 1) % 3 + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    const unsigned int edges_size = edges.size();
    for (unsigned int i = 0; i < edges_size; i++) {
        double weight = -1;
//        if (isCollapsable_cd1(edges[i][0], edges[i][1]) && isCollapsable_cd2(edges[i][0], edges[i][1])) {
        if (isCollapsable_cd1(edges[i][0], edges[i][1])) {
            weight = calEdgeLength(edges[i][0], edges[i][1]);
            if (isCollapsable_cd3(edges[i][0], edges[i][1], weight)) {
                ElementInQueue_ec ele(edges[i], weight);
                ec_queue.push(ele);
            }
        }
//        if (isCollapsable_cd1(edges[i][1], edges[i][0]) && isCollapsable_cd2(edges[i][0], edges[i][1])) {
        if (isCollapsable_cd1(edges[i][1], edges[i][0])) {
            weight = weight == -1 ? calEdgeLength(edges[i][0], edges[i][1]) : weight;
            if (isCollapsable_cd3(edges[i][0], edges[i][1], weight)) {
                ElementInQueue_ec ele({{edges[i][1], edges[i][0]}}, weight);
                ec_queue.push(ele);
            }
        }
    }

    counter = 0;
    suc_counter = 0;
    breakdown_timing = {{0, 0, 0, 0, 0}};
    breakdown_timing0 = {{0, 0}};
}

void EdgeCollapser::collapse() {
    tet_tss.assign(tets.size(), 0);
    int cnt = 0;
    logger().debug("edge queue size = {}", ec_queue.size());
    while (!ec_queue.empty()) {
        std::array<int, 2> v_ids = ec_queue.top().v_ids;
        double old_weight = ec_queue.top().weight;
        ec_queue.pop();

        if (!isEdgeValid(v_ids)) {
            continue;
        }

        //during operations, the length of edges in the queue may be changed
        //also, we need to eliminate the old edges, that is, the edges have an wrong/old weight in the queue
        double weight = calEdgeLength(v_ids);
        if (weight != old_weight || !isCollapsable_cd3(v_ids[0], v_ids[1], weight)) {
            continue;
        }
//        if(!isCollapsable_cd2(v_ids[0], v_ids[1])){
//            continue;
//        }

        while (!ec_queue.empty()) {
            std::array<int, 2> tmp_v_ids = ec_queue.top().v_ids;
            if (tmp_v_ids == v_ids)
                ec_queue.pop();
            else
                break;
        }

#if TIMING_BREAKDOWN
        igl_timer.start();
#endif
        int return_code = collapseAnEdge(v_ids[0], v_ids[1]);
        if (return_code == SUCCESS) {
#if TIMING_BREAKDOWN
            breakdown_timing[id_success] += igl_timer.getElapsedTime();
#endif
            suc_counter++;
            if (budget > 0) {
                budget--;
                if(budget == 0)
                    return;
            }
        } else if (return_code == ENVELOP_SUC) {
#if TIMING_BREAKDOWN
            breakdown_timing[id_env_success] += igl_timer.getElapsedTime();
#endif
            suc_counter++;
            cnt++;
            if (budget > 0) {
                budget--;
                if(budget == 0)
                    return;
            }
        } else {
            if (return_code == ENVELOP) {
#if TIMING_BREAKDOWN
                breakdown_timing[id_env_fail] += igl_timer.getElapsedTime();
#endif
            } else if (return_code == FLIP) {
#if TIMING_BREAKDOWN
                breakdown_timing[id_flip_fail] += igl_timer.getElapsedTime();
#endif
            } else {
#if TIMING_BREAKDOWN
                breakdown_timing[id_energy_fail] += igl_timer.getElapsedTime();
#endif
            }

            inf_es.push_back(v_ids);
            inf_e_tss.push_back(ts);
        }

        counter++;
    }
    logger().debug("{} {} {}", suc_counter, counter, inf_es.size());
    logger().debug("envelop accept = {}", envelop_accept_cnt);

    if (suc_counter == 0 || inf_es.size() == 0) {
//        logger().debug("checking.......");
//        init();
//        logger().debug("{}", ec_queue.size());
//        int cnt_flip=0, cnt_quality=0, cnt_envelop=0, cnt_suc=0;
//        while (!ec_queue.empty()) {
//            std::array<int, 2> v_ids = ec_queue.top().v_ids;
//            double old_weight=ec_queue.top().weight;
//            ec_queue.pop();
//
//            if (!isEdgeValid(v_ids)) {
//                continue;
//            }
//
//            double weight = calEdgeLength(v_ids);
//            if (weight != old_weight || !isCollapsable_cd3(v_ids[0], v_ids[1], weight)) {
//                continue;
//            }
//            if(!isCollapsable_cd2(v_ids[0], v_ids[1])){
//                continue;
//            }
//
//            while (!ec_queue.empty()) {
//                std::array<int, 2> tmp_v_ids = ec_queue.top().v_ids;
//                if (tmp_v_ids == v_ids)
//                    ec_queue.pop();
//                else
//                    break;
//            }
//
//            int return_code=collapseAnEdge(v_ids[0], v_ids[1]);
//            if (return_code == FLIP)
//                cnt_flip++;
//            else if(return_code==QUALITY)
//                cnt_quality++;
//            else if(return_code==ENVELOP)
//                cnt_envelop++;
//            else
//                cnt_suc++;
//        }
//        logger().debug("{} {} {} {}", cnt_flip, cnt_quality, cnt_envelop, cnt_suc);

        logger().debug("{}: {}s", breakdown_name0[id_sampling], breakdown_timing0[id_sampling]);
        logger().debug("{}: {}s", breakdown_name0[id_aabb], breakdown_timing0[id_aabb]);
        logger().debug("----");
        for (int i = 0; i < breakdown_timing.size(); i++)
            logger().debug("{}: {}s", breakdown_name[i], breakdown_timing[i]);

//        std::ofstream of(timing_log_file_name, std::fstream::app);
//        if (of.is_open()) {
//            of<<breakdown_name0[id_sampling]<<": "<<breakdown_timing0[id_sampling]<<"s"<<std::endl;
//            of<<breakdown_name0[id_aabb]<<": "<<breakdown_timing0[id_aabb]<<"s"<<std::endl;
//            of<<"----"<<std::endl;
//            for(int i=0;i<breakdown_timing.size();i++)
//                of<<breakdown_name[i]<<": "<<breakdown_timing[i]<<"s"<<std::endl;
//            of.close();
//        }

        logger().debug("energy_time = {}", energy_time);

        return;
    }

    postProcess();
}

void EdgeCollapser::postProcess() {
    logger().debug("postProcess!");
    counter = 0;
    suc_counter = 0;
    envelop_accept_cnt = 0;

    //you CANNOT sort it here!! every inf_es has a time stamp!!!
//    std::sort(inf_es.begin(), inf_es.end());
//    inf_es.erase(std::unique(inf_es.begin(), inf_es.end()), inf_es.end());

#if TIMING_BREAKDOWN
    igl_timer.start();
#endif
    std::vector<std::array<int, 2>> tmp_inf_es;
    const unsigned int inf_es_size = inf_es.size();
    tmp_inf_es.reserve(inf_es_size/4.0+1);
    for (unsigned int i = 0; i < inf_es_size; i++) {
        if (!isEdgeValid(inf_es[i]))
            continue;
        double weight = calEdgeLength(inf_es[i][0], inf_es[i][1]);
        if (!isCollapsable_cd3(inf_es[i][0], inf_es[i][1], weight))
            continue;

        bool is_recal = false;
        for (auto it = tet_vertices[inf_es[i][0]].conn_tets.begin(); it != tet_vertices[inf_es[i][0]].conn_tets.end();
             it++) {
            if (tet_tss[*it] > inf_e_tss[i]) {
                is_recal = true;
                break;
            }
        }

//        if (is_recal && isCollapsable_cd1(inf_es[i][0], inf_es[i][1]) && isCollapsable_cd2(inf_es[i][0], inf_es[i][1])) {
        if (is_recal && isCollapsable_cd1(inf_es[i][0], inf_es[i][1])) {
            if(!isLocked_ui(inf_es[i])) {
                ElementInQueue_ec ele(inf_es[i], weight);
                ec_queue.push(ele);
            }
        } else
            tmp_inf_es.push_back(inf_es[i]);
    }

    std::sort(tmp_inf_es.begin(), tmp_inf_es.end());
    tmp_inf_es.erase(std::unique(tmp_inf_es.begin(), tmp_inf_es.end()), tmp_inf_es.end());//it's better
    inf_es = tmp_inf_es;
    ts++;
    inf_e_tss = std::vector<int>(inf_es.size(), ts);

#if TIMING_BREAKDOWN
    breakdown_timing[id_postprocessing]+=igl_timer.getElapsedTime();
#endif

    collapse();
}

int EdgeCollapser::collapseAnEdge(int v1_id, int v2_id) {
    bool is_edge_too_short = false;
    bool is_edge_degenerate = false;
    double length = sqrt(CGAL::squared_distance(tet_vertices[v1_id].posf, tet_vertices[v2_id].posf));
    if(length == 0) {
        is_edge_degenerate = true;
    }
//    else if(length < 1e-30) {
//        logger().debug("{} {} {}{}{} {} {}{}{}", v1_id, tet_vertices[v1_id].is_on_surface, tet_vertices[v1_id].is_on_boundary, ", "
//, v2_id, tet_vertices[v2_id].is_on_surface, tet_vertices[v2_id].is_on_boundary, ": "
//, length);
//        logger().debug("{} {}", tet_vertices[v1_id].is_rounded, tet_vertices[v2_id].is_rounded);
//        is_edge_too_short = true;
//    }

    //check isolated
    if(tet_vertices[v1_id].is_on_surface && isIsolated(v1_id)) {
        tet_vertices[v1_id].is_on_surface = false;
        tet_vertices[v1_id].is_on_boundary = false;
        tet_vertices[v1_id].on_fixed_vertex = -1;
        tet_vertices[v1_id].on_face.clear();
        tet_vertices[v1_id].on_edge.clear();
    }
    if(!isBoundaryPoint(v1_id))
        tet_vertices[v1_id].is_on_boundary = false;

    //check boundary
    if(tet_vertices[v1_id].is_on_boundary && !tet_vertices[v2_id].is_on_boundary)
        if(!is_edge_degenerate && isPointOutBoundaryEnvelop(tet_vertices[v2_id].posf)) {
//            if(is_edge_too_short) {
//                logger().debug("v2 bonndary");
//                logger().debug("v1 boundary = {}", isPointOutBoundaryEnvelop(tet_vertices[v1_id].posf));
//            }
            return ENVELOP;
        }

    //check envelop
    if(tet_vertices[v1_id].is_on_surface && !tet_vertices[v2_id].is_on_surface){
        if(!is_edge_degenerate && isPointOutEnvelop(tet_vertices[v2_id].posf)) {
//            if(is_edge_too_short) {
//                logger().debug("v2 envelop");
//                logger().debug("v1 envelop = {}", isPointOutEnvelop(tet_vertices[v1_id].posf));
//            }
            return ENVELOP;
        }
    }

    //old_t_ids
    std::vector<int> old_t_ids;
    old_t_ids.reserve(tet_vertices[v1_id].conn_tets.size());
    for (auto it = tet_vertices[v1_id].conn_tets.begin(); it != tet_vertices[v1_id].conn_tets.end(); it++)
        old_t_ids.push_back(*it);
    std::vector<bool> is_removed(old_t_ids.size(), false);

    //new_tets
    std::vector<std::array<int, 4>> new_tets;
    new_tets.reserve(old_t_ids.size());
    std::unordered_set<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    for (int i = 0; i < old_t_ids.size(); i++) {
        auto it = std::find(tets[old_t_ids[i]].begin(), tets[old_t_ids[i]].end(), v2_id);
        if (it == tets[old_t_ids[i]].end()) {
            std::array<int, 4> t = tets[old_t_ids[i]];
            auto jt = std::find(t.begin(), t.end(), v1_id);
            *jt = v2_id;
            new_tets.push_back(t);
        } else {
            is_removed[i] = true;
            for (int j = 0; j < 4; j++)
                if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id)
                    n12_v_ids.insert(tets[old_t_ids[i]][j]);
            n12_t_ids.push_back(old_t_ids[i]);
        }
    }

    //check is_valid
    //check 1 //todo: look in details later
//    for (auto it = n12_v_ids.begin(); it != n12_v_ids.end(); it++) {
//        bool is_degenerate = true;
//        for (auto jt = tet_vertices[*it].conn_tets.begin(); jt != tet_vertices[*it].conn_tets.end(); jt++) {
//            auto kt = std::find(n12_t_ids.begin(), n12_t_ids.end(), *jt);
//            if (kt == n12_t_ids.end()) {
//                is_degenerate = false;
//                break;
//            }
//        }
//        if (is_degenerate)
//            return FLIP;
//
//        is_degenerate = true;
//        for (auto jt = tet_vertices[v2_id].conn_tets.begin(); jt != tet_vertices[v2_id].conn_tets.end(); jt++) {
//            auto kt = std::find(n12_t_ids.begin(), n12_t_ids.end(), *jt);
//            if (kt == n12_t_ids.end()) {
//                is_degenerate = false;
//                break;
//            }
//        }
//        if (is_degenerate)
//            return FLIP;
//    }

    //check 2
    if (isFlip(new_tets)) {
//        if(is_edge_too_short)
//            logger().debug("flip");
        return FLIP;
    }
    std::vector<TetQuality> tet_qs;
    igl::Timer tmp_timer;
    tmp_timer.start();
    calTetQualities(new_tets, tet_qs);
    energy_time+=tmp_timer.getElapsedTime();

    if (energy_type != state.ENERGY_NA && is_check_quality) {
        TetQuality old_tq, new_tq;
        getCheckQuality(old_t_ids, old_tq);
        getCheckQuality(tet_qs, new_tq);
//        if (is_edge_too_short)
//            logger().debug("old {} new {}", old_tq.slim_energy, new_tq.slim_energy);
//        if (is_soft && old_tq.slim_energy < soft_energy) {
//            old_tq.slim_energy = Args::args().filter_energy_thres;
//        }
        if(is_soft)
            old_tq.slim_energy = soft_energy;
        if (!tet_vertices[v1_id].is_rounded) //remove an unroundable vertex anyway
            new_tq.slim_energy = 0;
        if (!is_edge_degenerate && !new_tq.isBetterOrEqualThan(old_tq, energy_type, state)) {
//            if (is_edge_too_short)
//                logger().debug("quality");
            return QUALITY;
        }
    }

    //check 2.5
    if (tet_vertices[v1_id].is_on_boundary) {
        Point_3 old_p = tet_vertices[v1_id].pos;
        Point_3f old_pf = tet_vertices[v1_id].posf;
        tet_vertices[v1_id].posf = tet_vertices[v2_id].posf;
        tet_vertices[v1_id].pos = tet_vertices[v2_id].pos;
        if (!is_edge_degenerate && isBoundarySlide(v1_id, v2_id, old_pf)) {
            tet_vertices[v1_id].posf = old_pf;
            tet_vertices[v1_id].pos = old_p;
//            if (is_edge_too_short)
//                logger().debug("boundary");
            return ENVELOP;
        }
        tet_vertices[v1_id].posf = old_pf;
        tet_vertices[v1_id].pos = old_p;
    }

    //check 3
    bool is_envelop_suc = false;
    if (state.eps != state.EPSILON_NA && state.eps != state.EPSILON_INFINITE && tet_vertices[v1_id].is_on_surface) {
        if (!is_edge_degenerate && !isCollapsable_epsilon(v1_id, v2_id)) {
//            if (is_edge_too_short)
//                logger().debug("envelop");
            return ENVELOP;
        }
        is_envelop_suc = true;
        envelop_accept_cnt++;
        if (envelop_accept_cnt % 1000 == 0)
            logger().debug("1000 accepted!");
    }


    //real update
//    if(is_edge_too_short)
//        logger().debug("success");
    if(tet_vertices[v1_id].is_on_boundary)
        tet_vertices[v2_id].is_on_boundary=true;

    std::vector<std::array<int, 2>> update_sf_t_ids(n12_t_ids.size(), std::array<int, 2>());
    if (tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface) {
        for (int i = 0; i < n12_t_ids.size(); i++) {
            for (int j = 0; j < 4; j++) {
                if (tets[n12_t_ids[i]][j] == v1_id || tets[n12_t_ids[i]][j] == v2_id) {
                    std::vector<int> ts;
                    getFaceConnTets(tets[n12_t_ids[i]][(j + 1) % 4], tets[n12_t_ids[i]][(j + 2) % 4],
                                    tets[n12_t_ids[i]][(j + 3) % 4], ts);

//                    if(ts.size() != 2) {
//                        logger().debug("ts.size() != 2 but = {}", ts.size());
//                        logger().debug("v1 info:");
//                        tet_vertices[v1_id].printInfo();
//                        logger().debug("v2 info:");
//                        tet_vertices[v2_id].printInfo();
//
//                        tet_vertices[tets[n12_t_ids[i]][(j + 1) % 4]].printInfo();
//                        tet_vertices[tets[n12_t_ids[i]][(j + 2) % 4]].printInfo();
//                        tet_vertices[tets[n12_t_ids[i]][(j + 3) % 4]].printInfo();
//                        pausee();
//                    }

                    if(tets[n12_t_ids[i]][j] == v1_id)
                        update_sf_t_ids[i][1] = ts[0] != n12_t_ids[i] ? ts[0] : ts[1];
                    else
                        update_sf_t_ids[i][0] = ts[0] != n12_t_ids[i] ? ts[0] : ts[1];
                }
            }
        }
    }

    std::unordered_set<int> n1_v_ids;
    int cnt = 0;
    for (int i = 0; i < old_t_ids.size(); i++) {
        if (is_removed[i]) {
            t_is_removed[old_t_ids[i]] = true;
            for (int j = 0; j < 4; j++)
                if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id) {
                    tet_vertices[tets[old_t_ids[i]][j]].conn_tets.erase(
                            std::find(tet_vertices[tets[old_t_ids[i]][j]].conn_tets.begin(),
                                      tet_vertices[tets[old_t_ids[i]][j]].conn_tets.end(), old_t_ids[i]));
                }
            tet_vertices[v2_id].conn_tets.erase(std::find(tet_vertices[v2_id].conn_tets.begin(),
                                                          tet_vertices[v2_id].conn_tets.end(), old_t_ids[i]));
        } else {
            tet_vertices[v2_id].conn_tets.insert(old_t_ids[i]);
            tet_qualities[old_t_ids[i]] = tet_qs[cnt];
            for (int j = 0; j < 4; j++) {
                if (tets[old_t_ids[i]][j] != v1_id)
                    n1_v_ids.insert(tets[old_t_ids[i]][j]);//n12_v_ids would still be inserted
            }
            tets[old_t_ids[i]] = new_tets[cnt];
            cnt++;
        }
    }


    if (tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface) {
        tet_vertices[v2_id].is_on_surface = true;

        bool is_check_isolated = false;
        for (int i = 0; i < n12_t_ids.size(); i++) {
            std::array<int, 2> is_sf_fs;
            std::vector<int> es;
            for (int j = 0; j < 4; j++) {
                if (tets[n12_t_ids[i]][j] != v1_id && tets[n12_t_ids[i]][j] != v2_id)
                    es.push_back(tets[n12_t_ids[i]][j]);
                else if (tets[n12_t_ids[i]][j] == v1_id)
                    is_sf_fs[0] = is_surface_fs[n12_t_ids[i]][j];
                else
                    is_sf_fs[1] = is_surface_fs[n12_t_ids[i]][j];
            }
            //be careful about the order!!

//            if (is_sf_fs[0] == is_sf_fs[1]) {
//                if (is_sf_fs[0] != ON_SURFACE_FALSE) {
//                    is_sf_fs = {ON_SURFACE_FALSE, ON_SURFACE_FALSE};
//                    is_check_isolated = true;
//                } else
//                    continue;
//            } else if (is_sf_fs[0] != ON_SURFACE_FALSE && is_sf_fs[1] != ON_SURFACE_FALSE)
//                continue;
//            else if (is_sf_fs[0] != ON_SURFACE_FALSE)
//                is_sf_fs[1] = is_sf_fs[0] == ON_SURFACE_TRUE_INSIDE ? ON_SURFACE_TRUE_OUTSIDE : ON_SURFACE_TRUE_INSIDE;
//            else
//                is_sf_fs[0] = is_sf_fs[1] == ON_SURFACE_TRUE_INSIDE ? ON_SURFACE_TRUE_OUTSIDE : ON_SURFACE_TRUE_INSIDE;

            if (is_sf_fs[0] == is_sf_fs[1] && is_sf_fs[0] == state.NOT_SURFACE)
                continue;
            if(is_sf_fs[0] == state.NOT_SURFACE)
                is_sf_fs[0] = 0;
            if(is_sf_fs[1] == state.NOT_SURFACE)
                is_sf_fs[1] = 0;

            int tmp0 = is_sf_fs[0];
            int tmp1 = is_sf_fs[1];
            is_sf_fs[0] += -tmp1;
            is_sf_fs[1] += -tmp0;

            for (int j = 0; j < 4; j++) {
                int v_id0 = tets[update_sf_t_ids[i][0]][j];
                if (v_id0 != v2_id && v_id0 != es[0] && v_id0 != es[1])
                    is_surface_fs[update_sf_t_ids[i][0]][j] = is_sf_fs[0];

                int v_id1 = tets[update_sf_t_ids[i][1]][j];
                if (v_id1 != v2_id && v_id1 != es[0] && v_id1 != es[1])
                    is_surface_fs[update_sf_t_ids[i][1]][j] = is_sf_fs[1];
            }
        }
    }

    //update boundary points //todo: Pls figure out a more efficient way
//    if(tet_vertices[v2_id].is_on_boundary && !isBoundaryPoint(v2_id)) {
//        tet_vertices[v2_id].is_on_boundary = false;
////        logger().debug("a boundary vertex is removed");
//    }
//    for(int v_id:n12_v_ids) {
//        if (tet_vertices[v_id].is_on_boundary && !isBoundaryPoint(v_id)) {
//            tet_vertices[v_id].is_on_boundary = false;
////            logger().debug("a boundary vertex is removed");
//        }
//    }

    v_is_removed[v1_id] = true;

    //update time stamps
    ts++;
    for (int i = 0; i < old_t_ids.size(); i++) {
        tet_tss[old_t_ids[i]] = ts;
    }

    //add new elements
//    std::vector<std::array<int, 2>> es;
//    for (int i = 0; i < new_tets.size(); i++) {
//        for (int j = 0; j < 3; j++) {
//            std::array<int, 2> e = {new_tets[i][0], new_tets[i][j + 1]};
//            std::sort(e.begin(), e.end());
//            es.push_back(e);
//            e = {new_tets[i][j + 1], new_tets[i][(j + 1) % 3 + 1]};
//            es.push_back(e);
//        }
//    }
//    std::sort(es.begin(), es.end());
//    es.erase(std::unique(es.begin(), es.end()), es.end());
//    for (int i = 0; i < es.size(); i++) {
//        addNewEdge(es[i]);
//    }

//    logger().debug("{}{}jt==tri.end()", n1_v_ids.size(), "->";
    std::vector<int> n1_v_ids_vec, n12_v_ids_vec;
    n1_v_ids_vec.reserve(n1_v_ids.size());
    n12_v_ids_vec.reserve(n12_v_ids.size());
    for(auto it = n1_v_ids.begin(); it != n1_v_ids.end(); it++)
        n1_v_ids_vec.push_back(*it);
    for(auto it = n12_v_ids.begin(); it != n12_v_ids.end(); it++)
        n12_v_ids_vec.push_back(*it);
    std::sort(n1_v_ids_vec.begin(), n1_v_ids_vec.end());
    std::sort(n12_v_ids_vec.begin(), n12_v_ids_vec.end());
    n1_v_ids.clear();
    std::set_difference(n1_v_ids_vec.begin(), n1_v_ids_vec.end(),n12_v_ids_vec.begin(), n12_v_ids_vec.end(),
                        std::inserter(n1_v_ids, n1_v_ids.end()));

    for (auto it = n1_v_ids.begin(); it != n1_v_ids.end(); it++) {
        double weight = -1;
//        if (isCollapsable_cd1(v2_id, *it) && isCollapsable_cd2(v2_id, *it)) {
        if (isCollapsable_cd1(v2_id, *it)) {
            weight = calEdgeLength(v2_id, *it);
            if (isCollapsable_cd3(v2_id, *it, weight)) {
                std::array<int, 2> e={{v2_id, *it}};
                if(!isLocked_ui(e)) {
                    ElementInQueue_ec ele(e, weight);
                    ec_queue.push(ele);
                }
            }
        }
//        if (isCollapsable_cd1(*it, v2_id) && isCollapsable_cd2(v2_id, *it)) {
        if (isCollapsable_cd1(*it, v2_id)) {
            weight = weight == -1 ? calEdgeLength(*it, v2_id) : weight;
            if (isCollapsable_cd3(*it, v2_id, weight)) {
                std::array<int, 2> e={{*it, v2_id}};
                if(!isLocked_ui(e)) {
                    ElementInQueue_ec ele(e, weight);
                    ec_queue.push(ele);
                }
            }
        }
    }

    if(is_envelop_suc)
        return ENVELOP_SUC;
    return SUCCESS;
}

//bool EdgeCollapser::isCollapsable_cd2(int v1_id, int v2_id) {
//    return true;
//
//    //check envelop
//    if (tet_vertices[v1_id].is_on_surface && !tet_vertices[v2_id].is_on_surface) {
//        if (isPointOutEnvelop(tet_vertices[v2_id].posf))
//            return false;
//    }
//    return true;
//}

bool EdgeCollapser::isCollapsable_cd1(int v1_id, int v2_id) {
    //check the bbox tags //if the moved vertex is on the bbox
    bool is_movable = false;
    if (tet_vertices[v1_id].on_fixed_vertex < -1)
        return false;
    if (tet_vertices[v1_id].is_on_bbox && !tet_vertices[v2_id].is_on_bbox)
        return false;
    else if (tet_vertices[v1_id].is_on_bbox && tet_vertices[v2_id].is_on_bbox) {
        if (tet_vertices[v1_id].on_edge.size() == 0) {//inside the face
            is_movable = isHaveCommonEle(tet_vertices[v1_id].on_face, tet_vertices[v2_id].on_face);
        } else {//on the edge
            is_movable = isHaveCommonEle(tet_vertices[v1_id].on_edge, tet_vertices[v2_id].on_edge);
        }
        return is_movable;
    }

    //check the surface tags //if the vertex is on the surface
//    if (state.eps != state.EPSILON_NA) {
//        return true;
//    }
    return true;

    ////////////////////////////////////
//    //without envelop
//    if (tet_vertices[v1_id].on_fixed_vertex >= 0)
//        return false;
//
//    is_movable = false;
//    if (tet_vertices[v1_id].on_face.size() == 0) //inside the volumn
//        is_movable = true;
//    else {
//        if (tet_vertices[v1_id].on_edge.size() == 0) {//inside the face
//            if (isHaveCommonEle(tet_vertices[v1_id].on_face, tet_vertices[v2_id].on_face))
//                is_movable = true;
//        } else {//on the edge
//            if (isHaveCommonEle(tet_vertices[v1_id].on_edge, tet_vertices[v2_id].on_edge))
//                is_movable = true;
//        }
//    }
//    return is_movable;
}

//bool EdgeCollapser::isCollapsable_cd3(double weight) {
//    if (!is_limit_length)
//        return true;
//
//    if (weight < ideal_weight)
//        return true;
//    return false;
//}

bool EdgeCollapser::isCollapsable_cd3(int v1_id, int v2_id, double weight) {
    if (!is_limit_length)
        return true;

    double adaptive_scale = (tet_vertices[v1_id].adaptive_scale + tet_vertices[v2_id].adaptive_scale) / 2;
    if (weight < ideal_weight * adaptive_scale * adaptive_scale)
        return true;
//    if (tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface) {
//        if (weight < ideal_weight * adaptive_scale * adaptive_scale)
//            return true;
//    } else {
//        if (weight < ideal_weight)
//            return true;
//    }
    return false;
}

bool EdgeCollapser::isCollapsable_epsilon(int v1_id, int v2_id) {
//    std::vector<Triangle_3f> tris;
//    for (auto it = tet_vertices[v1_id].conn_tets.begin(); it != tet_vertices[v1_id].conn_tets.end(); it++) {
//        for (int j = 0; j < 4; j++) {
//            if (tets[*it][j] == v2_id && is_surface_fs[*it][j] != state.NOT_SURFACE) {
//                std::array<int, 3> tri = {tets[*it][(j + 1) % 4], tets[*it][(j + 2) % 4], tets[*it][(j + 3) % 4]};
//                auto jt = std::find(tri.begin(), tri.end(), v1_id);
//                if(jt==tri.end()){
//                    std::cout);
//                    throw TetWildError("");
//                }
//                *jt = v2_id;
//                Triangle_3f tr(tet_vertices[tri[0]].posf, tet_vertices[tri[1]].posf, tet_vertices[tri[2]].posf);
//                tris.push_back(tr);
//            }
//        }
//    }

    std::vector<std::array<int, 3>> tri_ids;
    for (auto it = tet_vertices[v1_id].conn_tets.begin(); it != tet_vertices[v1_id].conn_tets.end(); it++) {
        for (int j = 0; j < 4; j++) {
            if (tets[*it][j] != v1_id && is_surface_fs[*it][j] != state.NOT_SURFACE) {
                std::array<int, 3> tri = {{tets[*it][(j + 1) % 4], tets[*it][(j + 2) % 4], tets[*it][(j + 3) % 4]}};
                std::sort(tri.begin(), tri.end());
                tri_ids.push_back(tri);
            }
        }
    }
    std::sort(tri_ids.begin(), tri_ids.end());
    tri_ids.erase(std::unique(tri_ids.begin(), tri_ids.end()), tri_ids.end());

    std::vector<Triangle_3f> tris;
    for (int i = 0; i < tri_ids.size(); i++) {
        if (std::find(tri_ids[i].begin(), tri_ids[i].end(), v2_id) != tri_ids[i].end())
            continue;
        auto jt = std::find(tri_ids[i].begin(), tri_ids[i].end(), v1_id);
        *jt = v2_id;
        Triangle_3f tri(tet_vertices[tri_ids[i][0]].posf, tet_vertices[tri_ids[i][1]].posf, tet_vertices[tri_ids[i][2]].posf);
        tris.push_back(tri);
    }

    ///note that tris.size() can be 0 when v1 is on the boundary of the surface!!!
    for (int i = 0; i < tris.size(); i++) {
        if (isFaceOutEnvelop(tris[i]))
            return false;
    }

    return true;
}

bool EdgeCollapser::isEdgeValid(const std::array<int, 2>& e){
    if(v_is_removed[e[0]] || v_is_removed[e[1]])
        return false;
    return isHaveCommonEle(tet_vertices[e[0]].conn_tets, tet_vertices[e[1]].conn_tets);

//    if(!isHaveCommonEle(tet_vertices[e[0]].conn_tets, tet_vertices[e[1]].conn_tets))
//        return false;
//    return true;
}

//void EdgeCollapser::addNewEdge(const std::array<int, 2>& e){
//    double weight = -1;
//    if (isCollapsable_cd1(e[0], e[1])) {
//        weight = calEdgeLength(e[0], e[1]);
//        ElementInQueue_ec ele(e, weight);
//        if (isCollapsable_cd3(weight))
//            ec_queue.push(ele);
//    }
//    if (isCollapsable_cd1(e[1], e[0])) {
//        ElementInQueue_ec ele(e, weight == -1 ? calEdgeLength(e[0], e[1]) : weight);
//        if (isCollapsable_cd3(weight))
//            ec_queue.push(ele);
//    }
//}

} // namespace tetwild
