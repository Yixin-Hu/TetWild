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

#include <tetwild/VertexSmoother.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>
#include <pymesh/MshSaver.h>

namespace tetwild {

void VertexSmoother::smooth() {
    tets_tss = std::vector<int>(tets.size(), 1);
    tet_vertices_tss = std::vector<int>(tet_vertices.size(), 0);
    ts = 1;

    igl::Timer tmp_timer0;
    int max_pass = 1;
    double v_cnt = std::count(v_is_removed.begin(), v_is_removed.end(), false);
    for (int i = 0; i < max_pass; i++) {
        double suc_in = 0;
        double suc_surface = 0;
        smoothSingle();
        suc_in = suc_counter;
        if (state.eps >= 0) {
            smoothSurface();
            suc_surface = suc_counter;
        }
        logger().debug("{}", (suc_in + suc_surface) / v_cnt);
        if (suc_in + suc_surface < v_cnt * 0.1) {
            logger().debug("{}", i);
            break;
        }
    }
    for (int i = 0; i < breakdown_timing.size(); i++) {
        logger().debug("{}: {}s", breakdown_name[i], breakdown_timing[i]);
        breakdown_timing[i] = 0;//reset
    }
}

bool VertexSmoother::smoothSingleVertex(int v_id, bool is_cal_energy){
    std::vector<std::array<int, 4>> new_tets;
    std::vector<int> t_ids;
    for (int t_id:tet_vertices[v_id].conn_tets) {
        new_tets.push_back(tets[t_id]);
        t_ids.push_back(t_id);
    }

    ///try to round the vertex
    if(!tet_vertices[v_id].is_rounded) {
        Point_3 old_p = tet_vertices[v_id].pos;
        tet_vertices[v_id].pos = Point_3(tet_vertices[v_id].posf[0], tet_vertices[v_id].posf[1],
                                         tet_vertices[v_id].posf[2]);
        if(isFlip(new_tets))
            tet_vertices[v_id].pos = old_p;
        else
            tet_vertices[v_id].is_rounded = true;
    }

    ///check if should use exact smoothing
    bool is_valid = true;
    for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
        CGAL::Orientation ori = CGAL::orientation(tet_vertices[tets[*it][0]].posf, tet_vertices[tets[*it][1]].posf,
                                                  tet_vertices[tets[*it][2]].posf, tet_vertices[tets[*it][3]].posf);
        if (ori != CGAL::POSITIVE) {
            is_valid = false;
            break;
        }
    }
    if (!is_valid) {
        return false;
    } else {
        Point_3f pf;
        if(energy_type == state.ENERGY_AMIPS) {
            if (!NewtonsMethod(t_ids, new_tets, v_id, pf))
                return false;
        }

        //assign new coordinate and try to round it
        Point_3 old_p = tet_vertices[v_id].pos;
        Point_3f old_pf = tet_vertices[v_id].posf;
        bool old_is_rounded = tet_vertices[v_id].is_rounded;
        Point_3 p = Point_3(pf[0], pf[1], pf[2]);
        tet_vertices[v_id].pos = p;
        tet_vertices[v_id].posf = pf;
        tet_vertices[v_id].is_rounded = true;
        if (isFlip(new_tets)) {//TODO: why it happens?
            logger().debug("flip in the end");
            tet_vertices[v_id].pos = old_p;
            tet_vertices[v_id].posf = old_pf;
            tet_vertices[v_id].is_rounded = old_is_rounded;
        }
    }

    if(is_cal_energy){
        std::vector<TetQuality> tet_qs;
        calTetQualities(new_tets, tet_qs);
        int cnt = 0;
        for (int t_id:tet_vertices[v_id].conn_tets) {
            tet_qualities[t_id] = tet_qs[cnt++];
        }
    }

    return true;
}

void VertexSmoother::smoothSingle() {
    double old_ts = ts;
    counter = 0;
    suc_counter = 0;
    for (int v_id = 0; v_id < tet_vertices.size(); v_id++) {
        if (v_is_removed[v_id])
            continue;
        if (tet_vertices[v_id].is_on_bbox)
            continue;
        if (state.eps != state.EPSILON_INFINITE && tet_vertices[v_id].is_on_surface)
            continue;

        if (tet_vertices[v_id].is_locked)
            continue;

        ///check if its one-ring is changed
//        bool is_changed=false;
//        for(auto it=tet_vertices[v_id].conn_tets.begin();it!=tet_vertices[v_id].conn_tets.end();it++){
//            if(tets_tss[*it]>tet_vertices_tss[v_id]){
//                is_changed=true;
//                break;
//            }
//        }
//        if(!is_changed)
//            continue;

        counter++;

#if TIMING_BREAKDOWN
        igl_timer.start();
#endif
        std::vector<std::array<int, 4>> new_tets;
        std::vector<int> t_ids;
        for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
            new_tets.push_back(tets[*it]);
            t_ids.push_back(*it);
        }

        ///try to round the vertex
        if (!tet_vertices[v_id].is_rounded) {
            Point_3 old_p = tet_vertices[v_id].pos;
            tet_vertices[v_id].pos = Point_3(tet_vertices[v_id].posf[0], tet_vertices[v_id].posf[1],
                                             tet_vertices[v_id].posf[2]);
            if (isFlip(new_tets))
                tet_vertices[v_id].pos = old_p;
            else
                tet_vertices[v_id].is_rounded = true;
        }

        ///check if should use exact smoothing
        bool is_valid = true;
        for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
            CGAL::Orientation ori = CGAL::orientation(tet_vertices[tets[*it][0]].posf, tet_vertices[tets[*it][1]].posf,
                                                      tet_vertices[tets[*it][2]].posf, tet_vertices[tets[*it][3]].posf);
            if (ori != CGAL::POSITIVE) {
                is_valid = false;
                break;
            }
        }
#if TIMING_BREAKDOWN
        breakdown_timing[id_round] += igl_timer.getElapsedTime();
#endif

        if (!is_valid) {
            continue;
        } else {
            Point_3f pf;
            if (energy_type == state.ENERGY_AMIPS) {
                if (!NewtonsMethod(t_ids, new_tets, v_id, pf))
                    continue;
            }
#if TIMING_BREAKDOWN
            igl_timer.start();
#endif
            //assign new coordinate and try to round it
            Point_3 old_p = tet_vertices[v_id].pos;
            Point_3f old_pf = tet_vertices[v_id].posf;
            bool old_is_rounded = tet_vertices[v_id].is_rounded;
            Point_3 p = Point_3(pf[0], pf[1], pf[2]);
            tet_vertices[v_id].pos = p;
            tet_vertices[v_id].posf = pf;
            tet_vertices[v_id].is_rounded = true;
            if (isFlip(new_tets)) {//TODO: why it happens?
                logger().debug("flip in the end");
                tet_vertices[v_id].pos = old_p;
                tet_vertices[v_id].posf = old_pf;
                tet_vertices[v_id].is_rounded = old_is_rounded;
            }
#if TIMING_BREAKDOWN
            breakdown_timing[id_round] += igl_timer.getElapsedTime();
#endif
        }

        ///update timestamps
        ts++;
        for(auto it=tet_vertices[v_id].conn_tets.begin();it!=tet_vertices[v_id].conn_tets.end();it++)
            tets_tss[*it]=ts;
        tet_vertices_tss[v_id]=ts;

        suc_counter++;
    }

    //calculate the quality for all tets
    std::vector<std::array<int, 4>> new_tets;//todo: can be improve
    new_tets.reserve(std::count(t_is_removed.begin(), t_is_removed.end(), false));
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
//        if(tets_tss[i]<=old_ts)
//            continue;
        new_tets.push_back(tets[i]);
    }
    std::vector<TetQuality> tet_qs;
#if TIMING_BREAKDOWN
    igl_timer.start();
#endif
    calTetQualities(new_tets, tet_qs);
#if TIMING_BREAKDOWN
    breakdown_timing[id_value_e] += igl_timer.getElapsedTime();
#endif
    int cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
//        if(tets_tss[i]<=old_ts)
//            continue;
        tet_qualities[i] = tet_qs[cnt++];
    }
}

void VertexSmoother::smoothSurface() {//smoothing surface using two methods
//    suc_counter = 0;
//    counter = 0;
    int sf_suc_counter = 0;
    int sf_counter = 0;

    for (int v_id = 0; v_id < tet_vertices.size(); v_id++) {
        if (v_is_removed[v_id])
            continue;
        if (!tet_vertices[v_id].is_on_surface)
            continue;

        if (tet_vertices[v_id].is_locked)
            continue;

        if (isIsolated(v_id)) {
            tet_vertices[v_id].is_on_surface = false;
            tet_vertices[v_id].is_on_boundary = false;
            tet_vertices[v_id].on_fixed_vertex = -1;
            tet_vertices[v_id].on_face.clear();
            tet_vertices[v_id].on_edge.clear();
            continue;
        }
        if (!isBoundaryPoint(v_id))
            tet_vertices[v_id].is_on_boundary = false;

        counter++;
        sf_counter++;

        std::vector<std::array<int, 4>> new_tets;
        std::vector<int> old_t_ids;
        for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
            new_tets.push_back(tets[*it]);
            old_t_ids.push_back(*it);
        }

        if (!tet_vertices[v_id].is_rounded) {
            Point_3 old_p = tet_vertices[v_id].pos;
            tet_vertices[v_id].pos = Point_3(tet_vertices[v_id].posf[0], tet_vertices[v_id].posf[1],
                                             tet_vertices[v_id].posf[2]);
            if (isFlip(new_tets))
                tet_vertices[v_id].pos = old_p;
            else
                tet_vertices[v_id].is_rounded = true;
        }

        bool is_valid = true;
        for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
            CGAL::Orientation ori = CGAL::orientation(tet_vertices[tets[*it][0]].posf, tet_vertices[tets[*it][1]].posf,
                                                      tet_vertices[tets[*it][2]].posf, tet_vertices[tets[*it][3]].posf);
            if (ori != CGAL::POSITIVE) {
                is_valid = false;
                break;
            }
        }

        Point_3 p_out;
        Point_3f pf_out;
        if (!is_valid) {
            continue;
        } else {
            if (energy_type == state.ENERGY_AMIPS) {
                if (!NewtonsMethod(old_t_ids, new_tets, v_id, pf_out))
                    continue;
            }
            p_out = Point_3(pf_out[0], pf_out[1], pf_out[2]);
        }

        ///find one-ring surface faces
#if TIMING_BREAKDOWN
        igl_timer.start();
#endif
        std::vector<std::array<int, 3>> tri_ids;
        for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
            for (int j = 0; j < 4; j++) {
                if (tets[*it][j] != v_id && is_surface_fs[*it][j] != state.NOT_SURFACE) {
                    std::array<int, 3> tri = {{tets[*it][(j + 1) % 4], tets[*it][(j + 2) % 4], tets[*it][(j + 3) % 4]}};
                    std::sort(tri.begin(), tri.end());
                    tri_ids.push_back(tri);
                }
            }
        }
        std::sort(tri_ids.begin(), tri_ids.end());
        tri_ids.erase(std::unique(tri_ids.begin(), tri_ids.end()), tri_ids.end());

        Point_3f pf;
        Point_3 p;
        if (state.use_onering_projection) {//we have to use exact construction here. Or the projecting points may be not exactly on the plane.
            std::vector<Triangle_3> tris;
            for (int i = 0; i < tri_ids.size(); i++) {
                tris.push_back(Triangle_3(tet_vertices[tri_ids[i][0]].pos, tet_vertices[tri_ids[i][1]].pos,
                                          tet_vertices[tri_ids[i][2]].pos));
            }

            is_valid = false;
            for (int i = 0; i < tris.size(); i++) {
                if (tris[i].is_degenerate())
                    continue;
                Plane_3 pln = tris[i].supporting_plane();
                p = pln.projection(p_out);
                if (tris[i].has_on(p)) {
                    is_valid = true;
                    break;
                }
            }
            if (!is_valid)
                continue;
            pf = Point_3f(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
            p = Point_3(pf[0], pf[1], pf[2]);
        } else {
            GEO::vec3 geo_pf(pf_out[0], pf_out[1], pf_out[2]);
            GEO::vec3 nearest_pf;
            double _;
            if (tet_vertices[v_id].is_on_boundary)
                geo_b_tree.nearest_facet(geo_pf, nearest_pf, _);
            else
                geo_sf_tree.nearest_facet(geo_pf, nearest_pf, _);
            pf = Point_3f(nearest_pf[0], nearest_pf[1], nearest_pf[2]);
            p = Point_3(nearest_pf[0], nearest_pf[1], nearest_pf[2]);
        }
#if TIMING_BREAKDOWN
        breakdown_timing[id_project] += igl_timer.getElapsedTime();
#endif

        Point_3 old_p = tet_vertices[v_id].pos;
        Point_3f old_pf = tet_vertices[v_id].posf;
        std::vector<TetQuality> tet_qs;
        bool is_found = false;

        tet_vertices[v_id].posf = pf;
        tet_vertices[v_id].pos = p;
        if (isFlip(new_tets)) {
            tet_vertices[v_id].pos = old_p;
            tet_vertices[v_id].posf = old_pf;
            continue;
        }
        TetQuality old_tq, new_tq;
        getCheckQuality(old_t_ids, old_tq);
        calTetQualities(new_tets, tet_qs);
        getCheckQuality(tet_qs, new_tq);
        if (!new_tq.isBetterThan(old_tq, energy_type, state)) {
            tet_vertices[v_id].pos = old_p;
            tet_vertices[v_id].posf = old_pf;
            continue;
        }
        is_found = true;

        if (!is_found) {
            tet_vertices[v_id].pos = old_p;
            tet_vertices[v_id].posf = old_pf;
            continue;
        }

#if TIMING_BREAKDOWN
        igl_timer.start();
#endif
        ///check if the boundary is sliding
        if (tet_vertices[v_id].is_on_boundary) {
            if (isBoundarySlide(v_id, -1, old_pf)) {
                tet_vertices[v_id].pos = old_p;
                tet_vertices[v_id].posf = old_pf;
#if TIMING_BREAKDOWN
                breakdown_timing[id_aabb] += igl_timer.getElapsedTime();
#endif
                continue;
            }
        }

        ///check if tris outside the envelop
        std::vector<Triangle_3f> trisf;
        for (int i = 0; i < tri_ids.size(); i++) {
            auto jt = std::find(tri_ids[i].begin(), tri_ids[i].end(), v_id);
            int k = jt - tri_ids[i].begin();
            Triangle_3f tri(Point_3f(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                            tet_vertices[tri_ids[i][(k + 1) % 3]].posf, tet_vertices[tri_ids[i][(k + 2) % 3]].posf);
            if (!tri.is_degenerate())
                trisf.push_back(tri);
        }

        is_valid = true;
        for (int i = 0; i < trisf.size(); i++) {
            if (isFaceOutEnvelop(trisf[i])) {
                is_valid = false;
                break;
            }
        }
#if TIMING_BREAKDOWN
        breakdown_timing[id_aabb] += igl_timer.getElapsedTime();
#endif
        if (!is_valid) {
            tet_vertices[v_id].pos = old_p;
            tet_vertices[v_id].posf = old_pf;
            continue;
        }

        ///real update
        ///update timestamps
        ts++;
        for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++)
            tets_tss[*it] = ts;
        tet_vertices_tss[v_id] = ts;

        if (!tet_vertices[v_id].is_rounded) {
            tet_vertices[v_id].pos = Point_3(pf[0], pf[1], pf[2]);
            if (isFlip(new_tets)) {
                tet_vertices[v_id].pos = old_p;
                tet_vertices[v_id].is_rounded = false;
            } else
                tet_vertices[v_id].is_rounded = true;
        }
        for (int i = 0; i < old_t_ids.size(); i++)
            tet_qualities[old_t_ids[i]] = tet_qs[i];

        suc_counter++;
        sf_suc_counter++;
        if (sf_suc_counter % 1000 == 0)
            logger().debug("1000 accepted!");
    }
    logger().debug("Totally {}({}) vertices on surface are smoothed.", sf_suc_counter, sf_counter);
}

bool VertexSmoother::NewtonsMethod(const std::vector<int>& t_ids, const std::vector<std::array<int, 4>>& new_tets,
                                   int v_id, Point_3f& p) {
//    bool is_moved=true;
    bool is_moved = false;
    const int MAX_STEP = 15;
    const int MAX_IT = 20;
    Point_3f pf0 = tet_vertices[v_id].posf;
    Point_3 p0 = tet_vertices[v_id].pos;

    double old_energy = 0;
    Eigen::Vector3d J;
    Eigen::Matrix3d H;
    Eigen::Vector3d X0;
    for (int step = 0; step < MAX_STEP; step++) {
        if (NewtonsUpdate(t_ids, v_id, old_energy, J, H, X0) == false)
            break;
        Point_3f old_pf = tet_vertices[v_id].posf;
        Point_3 old_p = tet_vertices[v_id].pos;
        double a = 1;
        bool step_taken = false;
        double new_energy;

        for (int it = 0; it < MAX_IT; it++) {
            //solve linear system
            //check flip
            //check energy
            igl_timer.start();
            Eigen::Vector3d X = H.colPivHouseholderQr().solve(H * X0 - a * J);
            breakdown_timing[id_solve] += igl_timer.getElapsedTime();
            if (!X.allFinite()) {
                a /= 2.0;
                continue;
            }

            tet_vertices[v_id].posf = Point_3f(X(0), X(1), X(2));
            tet_vertices[v_id].pos = Point_3(X(0), X(1), X(2));
//            tet_vertices[v_id].is_rounded=true;//need to remember old value?

            //check flipping
            if (isFlip(new_tets)) {
                tet_vertices[v_id].posf = old_pf;
                tet_vertices[v_id].pos = old_p;
                a /= 2.0;
                continue;
            }

            //check quality
            igl_timer.start();
            new_energy = getNewEnergy(t_ids);
            breakdown_timing[id_value_e] += igl_timer.getElapsedTime();
            if (new_energy >= old_energy || std::isinf(new_energy) || std::isnan(new_energy)) {
                tet_vertices[v_id].posf = old_pf;
                tet_vertices[v_id].pos = old_p;
                a /= 2.0;
                continue;
            }

            step_taken = true;
            break;
        }
        if (std::abs(new_energy - old_energy) < 1e-5)
            step_taken = false;

        if (!step_taken) {
            if (step == 0)
                is_moved = false;
            else
                is_moved = true;
            break;
        } else
            is_moved = true;
    }
    p = tet_vertices[v_id].posf;
    tet_vertices[v_id].posf = pf0;
    tet_vertices[v_id].pos = p0;

    return is_moved;
}

double VertexSmoother::getNewEnergy(const std::vector<int>& t_ids) {
    double s_energy = 0;

#ifdef TETWILD_WITH_ISPC
    int n = t_ids.size();

    static thread_local std::vector<double> T0;
    static thread_local std::vector<double> T1;
    static thread_local std::vector<double> T2;
    static thread_local std::vector<double> T3;
    static thread_local std::vector<double> T4;
    static thread_local std::vector<double> T5;
    static thread_local std::vector<double> T6;
    static thread_local std::vector<double> T7;
    static thread_local std::vector<double> T8;
    static thread_local std::vector<double> T9;
    static thread_local std::vector<double> T10;
    static thread_local std::vector<double> T11;
    static thread_local std::vector<double> energy;

    if (T0.empty()) {
        // logger().trace("Initial ISPC allocation: n = {}", n);
    } else if (T0.size() != n) {
        // logger().trace("ISPC reallocation: n = {}", n);
    }

    T0.resize(n);
    T1.resize(n);
    T2.resize(n);
    T3.resize(n);
    T4.resize(n);
    T5.resize(n);
    T6.resize(n);
    T7.resize(n);
    T8.resize(n);
    T9.resize(n);
    T10.resize(n);
    T11.resize(n);
    energy.resize(n);

    for (int i = 0; i < n; i++) {
        T0[i] = tet_vertices[tets[t_ids[i]][0]].posf[0];
        T1[i] = tet_vertices[tets[t_ids[i]][0]].posf[1];
        T2[i] = tet_vertices[tets[t_ids[i]][0]].posf[2];
        T3[i] = tet_vertices[tets[t_ids[i]][1]].posf[0];
        T4[i] = tet_vertices[tets[t_ids[i]][1]].posf[1];
        T5[i] = tet_vertices[tets[t_ids[i]][1]].posf[2];
        T6[i] = tet_vertices[tets[t_ids[i]][2]].posf[0];
        T7[i] = tet_vertices[tets[t_ids[i]][2]].posf[1];
        T8[i] = tet_vertices[tets[t_ids[i]][2]].posf[2];
        T9[i] = tet_vertices[tets[t_ids[i]][3]].posf[0];
        T10[i] = tet_vertices[tets[t_ids[i]][3]].posf[1];
        T11[i] = tet_vertices[tets[t_ids[i]][3]].posf[2];
    }

    ispc::energy_ispc(T0.data(), T1.data(), T2.data(), T3.data(), T4.data(),
        T5.data(), T6.data(), T7.data(), T8.data(),
        T9.data(), T10.data(), T11.data(), energy.data(), n);

    for (int i = 0; i < n; i++) {
        s_energy += energy[i]; //s_energy intialized in the beginning
    }
#else
    for (int i = 0; i < t_ids.size(); i++) {
        std::array<double, 12> t;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 3; k++) {
                t[j*3 + k] = tet_vertices[tets[t_ids[i]][j]].posf[k];
            }
        }
        if (energy_type == state.ENERGY_AMIPS) {
            s_energy += comformalAMIPSEnergy_new(t.data());
        }
    }
#endif
    if (std::isinf(s_energy) || std::isnan(s_energy) || s_energy <= 0 || s_energy > state.MAX_ENERGY) {
        logger().debug("new E inf");
        s_energy = state.MAX_ENERGY;
    }

    return s_energy;
}

bool VertexSmoother::NewtonsUpdate(const std::vector<int>& t_ids, int v_id,
                                   double& energy, Eigen::Vector3d& J, Eigen::Matrix3d& H, Eigen::Vector3d& X0) {
    energy = 0;
    for (int i = 0; i < 3; i++) {
        J(i) = 0;
        for (int j = 0; j < 3; j++) {
            H(i, j) = 0;
        }
        X0(i) = tet_vertices[v_id].posf[i];
    }

    for (int i = 0; i < t_ids.size(); i++) {
        std::array<double, 12> t;
        int start = 0;
        for (int j = 0; j < 4; j++) {
            if (tets[t_ids[i]][j] == v_id) {
                start = j;
                break;
            }
        }
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 3; k++) {
                t[j*3+k] = tet_vertices[tets[t_ids[i]][(start + j) % 4]].posf[k];
            }
        }
#ifndef TETWILD_WITH_ISPC
        igl_timer.start();
        energy += comformalAMIPSEnergy_new(t.data());
        breakdown_timing[id_value_e] += igl_timer.getElapsedTime();
#endif

        double J_1[3];
        double H_1[9];
        igl_timer.start();
        comformalAMIPSJacobian_new(t.data(), J_1);
        breakdown_timing[id_value_j] += igl_timer.getElapsedTime();
        igl_timer.start();
        comformalAMIPSHessian_new(t.data(), H_1);
        breakdown_timing[id_value_h] += igl_timer.getElapsedTime();

        for (int j = 0; j < 3; j++) {
            J(j) += J_1[j];
            H(j, 0) += H_1[j * 3 + 0];
            H(j, 1) += H_1[j * 3 + 1];
            H(j, 2) += H_1[j * 3 + 2];
        }
    }
#ifdef TETWILD_WITH_ISPC
    igl_timer.start();
    energy = getNewEnergy(t_ids);
    breakdown_timing[id_value_e] += igl_timer.getElapsedTime();
#endif

    if (std::isinf(energy)) {
        logger().debug("{} E inf", v_id);
        energy = state.MAX_ENERGY;
    }
    if (std::isnan(energy)) {
        logger().debug("{} E nan", v_id);
        return false;
    }
    if (energy <= 0) {
        logger().debug("{} E < 0", v_id);
        return false;
    }
    if (!J.allFinite()) {
        logger().debug("{} J inf/nan", v_id);
        return false;
    }
    if (!H.allFinite()) {
        logger().debug("{} H inf/nan", v_id);
        return false;
    }

    return true;
}

int VertexSmoother::laplacianBoundary(const std::vector<int>& b_v_ids, const std::vector<bool>& tmp_is_on_surface,
                                      const std::vector<bool>& tmp_t_is_removed){
    int cnt_suc = 0;
    double max_slim_evergy = 0;
    for(unsigned int i=0;i<tet_qualities.size();i++) {
        if (tmp_t_is_removed[i])
            continue;
        if (tet_qualities[i].slim_energy > max_slim_evergy)
            max_slim_evergy = tet_qualities[i].slim_energy;
    }

    for(int v_id:b_v_ids){
        // do laplacian on v_id
        std::vector<std::array<int, 4>> new_tets;
        std::unordered_set<int> n_v_ids;
//        std::unordered_set<int> n_v_ids2;
        std::unordered_set<int> tmp_n_sf_v_ids;
        std::unordered_set<int> n_sf_v_ids;
        for(int t_id:tet_vertices[v_id].conn_tets){
            for(int j=0;j<4;j++){
                if(tmp_is_on_surface[tets[t_id][j]]) {
                    tmp_n_sf_v_ids.insert(tets[t_id][j]);
                    continue;
                }
                if(!tet_vertices[tets[t_id][j]].is_on_surface)
                    n_v_ids.insert(tets[t_id][j]);
            }
            new_tets.push_back(tets[t_id]);
        }
        for(int n_sf_v_id:tmp_n_sf_v_ids){
            std::vector<int> t_ids;
            setIntersection(tet_vertices[v_id].conn_tets, tet_vertices[n_sf_v_id].conn_tets, t_ids);
            bool has_removed = false;
            bool has_unremoved = false;
            for(int t_id:t_ids){
                if(tmp_t_is_removed[t_id])
                    has_removed=true;
                if(!tmp_t_is_removed[t_id])
                    has_unremoved=true;
            }
            if(has_removed && has_unremoved)
                n_sf_v_ids.insert(n_sf_v_id);
        }
//        for(int n_v_id:n_v_ids){
//            for(int t_id:tet_vertices[n_v_id].conn_tets){
//                for(int j=0;j<4;j++)
//                    if(!tmp_is_on_surface[tets[t_id][j]] && !tet_vertices[tets[t_id][j]].is_on_surface)
//                        n_v_ids2.insert(tets[t_id][j]);
//            }
//        }
        std::array<double, 3> vec ={{0, 0, 0}};
        for(int n_sf_v_id:n_sf_v_ids) {
            for (int j = 0; j < 3; j++)
                vec[j] += tet_vertices[n_sf_v_id].posf[j];
        }
        for(int j=0;j<3;j++) {
            vec[j] = (vec[j] / n_sf_v_ids.size()) - tet_vertices[v_id].posf[j];
        }

        // do bisection and check flipping
        Point_3 old_p = tet_vertices[v_id].pos;
        Point_3f old_pf = tet_vertices[v_id].posf;
        double a = 1;
        bool is_suc = false;
        while(true) {
            //give stop condition
            bool is_stop = true;
            for (int j = 0; j < 3; j++)
                if (vec[j] * a > state.eps)
                    is_stop = false;
            if (is_stop)
                break;
            tet_vertices[v_id].pos = Point_3(old_pf[0] + vec[0] * a, old_pf[1] + vec[1] * a, old_pf[2] + vec[2] * a);
            tet_vertices[v_id].posf = Point_3f(old_pf[0] + vec[0] * a, old_pf[1] + vec[1] * a, old_pf[2] + vec[2] * a);
            if (isFlip(new_tets)) {
                a /= 2;
                continue;
            }
            //check quality
            std::vector<TetQuality> tet_qs;
            calTetQualities(new_tets, tet_qs);
            bool is_valid=true;
            for (int i = 0; i < tet_qs.size(); i++) {
                if (tet_qs[i].slim_energy > max_slim_evergy)
                    is_valid=false;
            }
            if(!is_valid) {
                a /= 2;
                continue;
            }

            int cnt = 0;
            for (int t_id:tet_vertices[v_id].conn_tets) {
                tet_qualities[t_id] = tet_qs[cnt++];
            }

            is_suc = true;
            cnt_suc++;
            break;
        }
        if(!is_suc) {
            tet_vertices[v_id].pos = old_p;
            tet_vertices[v_id].posf = old_pf;
            continue;
        }

        std::vector<TetQuality> tet_qs;
        calTetQualities(new_tets, tet_qs);
        int cnt = 0;
        for (int t_id:tet_vertices[v_id].conn_tets) {
            tet_qualities[t_id] = tet_qs[cnt++];
        }

        // do normal smoothing on neighbor vertices
//        logger().debug("n_v_ids.size = {}", n_v_ids.size());
//        logger().debug("n_v_ids2.size = {}", n_v_ids2.size());
        for(int n_v_id:n_v_ids){
            smoothSingleVertex(n_v_id, true);
        }
//        for(int n_v_id:n_v_ids2){
//            smoothSingleVertex(n_v_id, true);
//        }
//        for(int n_v_id:n_v_ids){
//            smoothSingleVertex(n_v_id, true);
//        }
    }

    logger().debug("suc.size = {}", cnt_suc);
    return cnt_suc;
}

void VertexSmoother::outputOneRing(int v_id, std::string s){
    PyMesh::MshSaver mSaver(state.working_dir+state.postfix+"_smooth_"+std::to_string(v_id)+s+".msh", true);
    std::vector<int> v_ids;
    std::vector<int> new_ids(tet_vertices.size(), -1);
    for(int t_id: tet_vertices[v_id].conn_tets){
        for(int j=0;j<4;j++)
            v_ids.push_back(tets[t_id][j]);
    }
    std::sort(v_ids.begin(), v_ids.end());
    v_ids.erase(std::unique(v_ids.begin(), v_ids.end()), v_ids.end());
    int cnt=0;
    for(int id:v_ids){
        new_ids[id] = cnt;
        cnt++;
    }

    Eigen::VectorXd oV(v_ids.size() * 3);
    Eigen::VectorXi oT(tet_vertices[v_id].conn_tets.size() * 4);
    for (int i = 0; i < v_ids.size(); i++) {
        for (int j = 0; j < 3; j++)
            oV(i * 3 + j) = tet_vertices[v_ids[i]].posf[j];
    }
    cnt = 0;
    for (int t_id: tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            oT(cnt * 4 + j) = new_ids[tets[t_id][j]];
        }
        cnt++;
    }

    mSaver.save_mesh(oV, oT, 3, mSaver.TET);

    Eigen::VectorXd cv(v_ids.size());
    for(int i=0;i<v_ids.size();i++){
        cv[i] = (v_ids[i] == v_id);
    }
    mSaver.save_scalar_field("center vertex", cv);

    Eigen::VectorXd q(tet_vertices[v_id].conn_tets.size());
    cnt = 0;
    for (int t_id: tet_vertices[v_id].conn_tets) {
        q[cnt] = tet_qualities[t_id].slim_energy;
        cnt++;
    }
    mSaver.save_elem_scalar_field("quality", q);
}

} // namespace tetwild
