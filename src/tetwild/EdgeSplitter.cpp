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

#include <tetwild/EdgeSplitter.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>

namespace tetwild {

void EdgeSplitter::getMesh_ui(const std::vector<std::array<int, 4>>& tets, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    ///get V, F, C
    V.resize(tets.size() * 4, 3);
    F.resize(tets.size() * 4, 3);
    Eigen::VectorXd Z(F.rows());
    int i = 0;
    for (unsigned j = 0; j < tets.size(); ++j) {
        for (int k = 0; k < 4; k++) {
            for (int r = 0; r < 3; r++)
                V(i * 4 + k, r) = tet_vertices[tets[j][k]].posf[r];
        }
        F.row(i * 4 + 0) << (i * 4) + 0, (i * 4) + 1, (i * 4) + 3;
        F.row(i * 4 + 1) << (i * 4) + 0, (i * 4) + 2, (i * 4) + 1;
        F.row(i * 4 + 2) << (i * 4) + 3, (i * 4) + 2, (i * 4) + 0;
        F.row(i * 4 + 3) << (i * 4) + 1, (i * 4) + 2, (i * 4) + 3;
        i++;
    }
}

void EdgeSplitter::init() {
    std::vector<std::array<int, 2>> edges;
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

    for (unsigned int i = 0; i < edges.size(); i++) {
        double weight = calEdgeLength(edges[i][0], edges[i][1]);
        if (isSplittable_cd1(edges[i][0], edges[i][1], weight)) {
            ElementInQueue_es ele(edges[i], weight);
            es_queue.push(ele);
        }
    }

    t_empty_start = 0;
    v_empty_start = 0;

//    if(budget > 0)
//        is_cal_quality_end = true;

    counter=0;
    suc_counter=0;

}

void EdgeSplitter::split() {

    if(budget >0) {
        int v_slots = std::count(v_is_removed.begin(), v_is_removed.end(), true);
        v_slots = budget - v_slots;
        if (v_slots > 0) {
            tet_vertices.reserve(tet_vertices.size() + v_slots);
            v_is_removed.reserve(tet_vertices.size() + v_slots);
        }
        int t_slots = std::count(t_is_removed.begin(), t_is_removed.end(), true);
        t_slots = budget * 6 - t_slots;
        if (t_slots > 0) {
            tet_vertices.reserve(tet_vertices.size() + t_slots);
            v_is_removed.reserve(tet_vertices.size() + t_slots);
        }
    } else {
        // reserve space
        int v_slot_size = std::count(v_is_removed.begin(), v_is_removed.end(), true);
        int t_slot_size = std::count(t_is_removed.begin(), t_is_removed.end(), true);
        if (v_slot_size < es_queue.size() * 2)
            tet_vertices.reserve(es_queue.size() * 2 - v_slot_size);
        if (t_slot_size < es_queue.size() * 6 * 2)
            tets.reserve(es_queue.size() * 6 * 2 - t_slot_size + 1);
    }
    logger().debug("{}", es_queue.size());
    logger().debug("ideal_weight = {}", ideal_weight);

    while (!es_queue.empty()) {
        const ElementInQueue_es &ele = es_queue.top();

        std::array<int, 2> v_ids = ele.v_ids;
//        if (state.is_print_tmp)
//            logger().debug("{}{}{} {} {} {} {}", v_ids[0], ' ', v_ids[1]
//, std::sqrt(calEdgeLength(v_ids))
//, std::sqrt(ideal_weight) *
//                           (tet_vertices[v_ids[0]].adaptive_scale + tet_vertices[v_ids[1]].adaptive_scale) / 2.0
//, tet_vertices[v_ids[0]].adaptive_scale
//, tet_vertices[v_ids[1]].adaptive_scale
//);
        es_queue.pop();
        if (splitAnEdge(v_ids))
            suc_counter++;
        counter++;

        if (budget > 0) {
            budget--;
            if(budget == 0)
                break;
        }
    }

    //cal the qualities in the very end
    //not really, in the last few pass, only some of the edges are splitted
    //maybe, we can marked all the splitted tets in this pass and the update their quality in the end? -- Try it.
    //todo ...
    if (is_cal_quality_end) {
        const int tets_size = tets.size();
        std::vector<std::array<int, 4>> tmp_tets;
        for (int i = 0; i < tets_size; i++) {
            if (t_is_removed[i])
                continue;
            tmp_tets.push_back(tets[i]);
        }

        std::vector<TetQuality> tet_qs;
        calTetQualities(tmp_tets, tet_qs);
        int cnt = 0;
        for (int i = 0; i < tets_size; i++) {
            if (t_is_removed[i])
                continue;
            tet_qualities[i] = tet_qs[cnt++];
        }
    }

}

bool EdgeSplitter::splitAnEdge(const std::array<int, 2>& edge) {
    int v1_id = edge[0];
    int v2_id = edge[1];

    //add new vertex
    TetVertex v;//tet_vertices[v_id] is actually be reset
    bool is_found = false;
    for(int i=v_empty_start;i<v_is_removed.size();i++){
        v_empty_start = i;
        if(v_is_removed[i]) {
            is_found = true;
            break;
        }
    }
    if(!is_found)
        v_empty_start = v_is_removed.size();

    int v_id = v_empty_start;
    if (v_empty_start < v_is_removed.size()) {
        tet_vertices[v_id] = v;
        v_is_removed[v_id] = false;
    } else {
        tet_vertices.push_back(v);
        v_is_removed.push_back(false);
    }

    //    int v_id = -1;
//    auto empty_slot = std::find(v_is_removed.begin(), v_is_removed.end(), true);//can be improved
//    if (empty_slot != v_is_removed.end()) {
//        v_id = empty_slot - v_is_removed.begin();
//        tet_vertices[v_id] = v;
//        v_is_removed[v_id] = false;
//    } else {
//        tet_vertices.push_back(v);
//        v_is_removed.push_back(false);
//        v_id = v_is_removed.size() - 1;
//    }

    //old_t_ids
    std::vector<int> old_t_ids;
    setIntersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, old_t_ids);

    //new_tets
    std::vector<int> new_t_ids;
    std::vector<int> n12_v_ids;
    std::vector<std::array<int, 4>> new_tets;
    new_tets.reserve(old_t_ids.size() * 2);
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id)
                n12_v_ids.push_back(tets[old_t_ids[i]][j]);
        }

        std::array<int, 4> tet1 = tets[old_t_ids[i]], tet2 = tets[old_t_ids[i]];
        auto it = std::find(tet1.begin(), tet1.end(), v2_id);
        *it = v_id;
        it = std::find(tet2.begin(), tet2.end(), v1_id);
        *it = v_id;
        new_tets.push_back(tet1);
        new_tets.push_back(tet2);
    }

    //check is_valid
    tet_vertices[v_id].adaptive_scale = (tet_vertices[v1_id].adaptive_scale + tet_vertices[v2_id].adaptive_scale) / 2;
    if(tet_vertices[v1_id].is_locked && tet_vertices[v2_id].is_locked)
        tet_vertices[v_id].is_locked = true;

    tet_vertices[v_id].posf = CGAL::midpoint(tet_vertices[v1_id].posf, tet_vertices[v2_id].posf);
    tet_vertices[v_id].pos = Point_3(tet_vertices[v_id].posf[0], tet_vertices[v_id].posf[1], tet_vertices[v_id].posf[2]);
    std::vector<TetQuality> tet_qs;
    if(!is_cal_quality_end) {
        calTetQualities(new_tets, tet_qs);
    }

    if (isFlip(new_tets)) {
        tet_vertices[v_id].pos = CGAL::midpoint(tet_vertices[v1_id].pos, tet_vertices[v2_id].pos);
        tet_vertices[v_id].posf = Point_3f(CGAL::to_double(tet_vertices[v_id].pos[0]), CGAL::to_double(tet_vertices[v_id].pos[1]),
                                           CGAL::to_double(tet_vertices[v_id].pos[2]));
        tet_vertices[v_id].is_rounded = false;
    } else {
        tet_vertices[v_id].is_rounded = true;
    }

//    if(!is_cal_quality_end)
//          calTetQualities(new_tets, tet_qs);

    ////real update//
    //update boundary tags
    if(isEdgeOnBoundary(v1_id, v2_id)) {
        tet_vertices[v_id].is_on_boundary = true;
    }

    //update surface tags
    if (state.eps != state.EPSILON_INFINITE) {
        if (isEdgeOnSurface(v1_id, v2_id)) {
            tet_vertices[v_id].is_on_surface = true;
            if (state.eps == state.EPSILON_NA) {
                setIntersection(tet_vertices[v1_id].on_edge, tet_vertices[v2_id].on_edge, tet_vertices[v_id].on_edge);
                setIntersection(tet_vertices[v1_id].on_face, tet_vertices[v2_id].on_face, tet_vertices[v_id].on_face);
            }
        } else
            tet_vertices[v_id].is_on_surface = false;
    }

    //get new tet ids
    getNewTetSlots(old_t_ids.size(), new_t_ids);
    for (int i = 0; i < old_t_ids.size(); i++) {
        tets[old_t_ids[i]] = new_tets[i * 2];
        tets[new_t_ids[i]] = new_tets[i * 2 + 1];
        if(!is_cal_quality_end) {
            tet_qualities[old_t_ids[i]] = tet_qs[i * 2];
            tet_qualities[new_t_ids[i]] = tet_qs[i * 2 + 1];
        }
        t_is_removed[new_t_ids[i]] = false;
        is_surface_fs[new_t_ids[i]] = is_surface_fs[old_t_ids[i]];
    }

    //track surface
    for (int i = 0; i < new_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {//v1->v
            if (tets[new_t_ids[i]][j] == v2_id)
                is_surface_fs[new_t_ids[i]][j] = state.NOT_SURFACE;
//                else if(tets[new_t_ids[i]][j]==v_id)//no need to change
//                    is_surface_fs[new_t_ids[i]][j]=is_surface_fs[old_t_ids[i]][j];
        }
    }
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {//v2->v
            if (tets[old_t_ids[i]][j] == v1_id)
                is_surface_fs[old_t_ids[i]][j] = state.NOT_SURFACE;
        }
    }

    //update bbox tags //Note that no matter what the epsilon is, the bbox has to be preserved anyway
    if (tet_vertices[v1_id].is_on_bbox && tet_vertices[v2_id].is_on_bbox) {
        setIntersection(tet_vertices[v1_id].on_face, tet_vertices[v2_id].on_face, tet_vertices[v_id].on_face);
        if (tet_vertices[v_id].on_face.size() == 0)
            tet_vertices[v_id].is_on_bbox = false;
        else {
            tet_vertices[v_id].is_on_bbox = true;
            setIntersection(tet_vertices[v1_id].on_edge, tet_vertices[v2_id].on_edge, tet_vertices[v_id].on_edge);
        }
    }

    //update the connection
    for (int i = 0; i < old_t_ids.size(); i++) {
        tet_vertices[v2_id].conn_tets.erase(old_t_ids[i]);
        tet_vertices[v2_id].conn_tets.insert(new_t_ids[i]);
    }
    for (int i = 0; i < n12_v_ids.size(); i++) {
        tet_vertices[n12_v_ids[i]].conn_tets.insert(new_t_ids[i / 2]);
    }
    tet_vertices[v_id].conn_tets.insert(old_t_ids.begin(), old_t_ids.end());
    tet_vertices[v_id].conn_tets.insert(new_t_ids.begin(), new_t_ids.end());

    //push new ele into queue
    double weight = calEdgeLength(v1_id, v_id);
    if (isSplittable_cd1(v1_id, v_id, weight)) {
        std::array<int, 2> e={{v1_id, v_id}};
        if(!isLocked_ui(e)) {
            ElementInQueue_es ele(e, weight);
            es_queue.push(ele);
        }
    }

    weight = calEdgeLength(v2_id, v_id);
    if (isSplittable_cd1(v2_id, v_id, weight)) {
        std::array<int, 2> e={{v2_id, v_id}};
        if(!isLocked_ui(e)) {
            ElementInQueue_es ele(e, weight);
            es_queue.push(ele);
        }
    }

    std::sort(n12_v_ids.begin(), n12_v_ids.end());
    n12_v_ids.erase(std::unique(n12_v_ids.begin(), n12_v_ids.end()), n12_v_ids.end());
    for (auto it = n12_v_ids.begin(); it != n12_v_ids.end(); it++) {
        weight = calEdgeLength(*it, v_id);
        if (isSplittable_cd1(*it, v_id, weight)) {
            std::array<int, 2> e = {{*it, v_id}};
            if(!isLocked_ui(e)) {
                ElementInQueue_es ele(e, weight);
                es_queue.push(ele);
            }
        }
    }

    return true;
}

int EdgeSplitter::getOverRefineScale(int v1_id, int v2_id){
    return 1;

    if(is_over_refine) {
        std::vector<int> n12_t_ids;
        setIntersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, n12_t_ids);
        for(int i=0;i<n12_t_ids.size();i++) {
            if (energy_type == state.ENERGY_AMIPS &&
                tet_qualities[n12_t_ids[i]].slim_energy > 500) {//todo: add || for other types of energy
                int scale = 1;
                scale = (tet_qualities[n12_t_ids[i]].slim_energy - 500) / 500.0;
                if (scale < 1)
                    scale = 1;
                else if (scale > 5)
                    scale = 5;
                return 1 + scale;
            }
        }
    }
    return 1;
}

bool EdgeSplitter::isSplittable_cd1(double weight) {
    if(is_check_quality)
        return true;

    if (weight > ideal_weight)
        return true;
    return false;
}

bool EdgeSplitter::isSplittable_cd1(int v1_id, int v2_id, double weight) {
    double adaptive_scale = (tet_vertices[v1_id].adaptive_scale + tet_vertices[v2_id].adaptive_scale) / 2.0;
//    if(adaptive_scale==0){
//        logger().debug("adaptive_scale==0!!!");
//    }
    if (weight > ideal_weight * adaptive_scale * adaptive_scale)
        return true;
//    if (tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface) {
//        if (weight > ideal_weight * adaptive_scale * adaptive_scale)
//            return true;
//    } else {
//        if (weight > ideal_weight)
//            return true;
//    }
    return false;
}

void EdgeSplitter::getNewTetSlots(int n, std::vector<int>& new_conn_tets) {
    int cnt = 0;
    for (int i = t_empty_start; i < t_is_removed.size(); i++) {
        if (t_is_removed[i]) {
            new_conn_tets.push_back(i);
            cnt++;
            if (cnt == n) {
                t_empty_start = i + 1;
                break;
            }
        }
    }
    if (cnt < n) {
        for (int i = 0; i < n - cnt; i++)
            new_conn_tets.push_back(tets.size() + i);

        tets.resize(tets.size() + n - cnt);
        t_is_removed.resize(t_is_removed.size() + n - cnt);
        tet_qualities.resize(tet_qualities.size() + n - cnt);
        is_surface_fs.resize(is_surface_fs.size() + n - cnt);
        t_empty_start = tets.size();
    }
}

} // namespace tetwild
