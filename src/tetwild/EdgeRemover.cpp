// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 4/17/17.
//

#include <tetwild/EdgeRemover.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>
#include <unordered_map>

namespace tetwild {

void EdgeRemover::init() {
    energy_time = 0;

    const unsigned int tets_size = tets.size();
    std::vector<std::array<int, 2>> edges;
    edges.reserve(tets_size * 6);
    for (unsigned int i = 0; i < tets_size; i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{tets[i][0], tets[i][j + 1]}};
            if (e[0] > e[1]) e = {{e[1], e[0]}};
            if (!isLocked_ui(e))
                edges.push_back(e);
            e = {{tets[i][j + 1], tets[i][(j + 1) % 3 + 1]}};
            if (e[0] > e[1]) e = {{e[1], e[0]}};
            if (!isLocked_ui(e))
                edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    for (unsigned int i = 0; i < edges.size(); i++) {
        std::vector<int> t_ids;
        setIntersection(tet_vertices[edges[i][0]].conn_tets, tet_vertices[edges[i][1]].conn_tets, t_ids);
        addNewEdge(edges[i]);
//        if (isSwappable_cd1(edges[i])) {
//            double weight = calEdgeLength(edges[i]);
//            if (isSwappable_cd2(weight)) {
//                ElementInQueue_er ele(edges[i], weight);
//                er_queue.push(ele);
//            }
//        }
    }

    counter = 0;
    suc_counter = 0;
    t_empty_start = 0;
    v_empty_start = 0;

    equal_buget = 100;
}

void EdgeRemover::swap(){
    tmp_cnt3=0;
    tmp_cnt4=0;
    tmp_cnt5=0;
    tmp_cnt6=0;
    int cnt5=0;

    while(!er_queue.empty()){
        const ElementInQueue_er& ele=er_queue.top();

        if(!isEdgeValid(ele.v_ids)){
            er_queue.pop();
            continue;
        }

        std::vector<int> t_ids;
        if(!isSwappable_cd1(ele.v_ids, t_ids, true)){
            er_queue.pop();
            continue;
        }

        std::array<int, 2> v_ids=ele.v_ids;
        er_queue.pop();

//        logger().debug("{} {} {} ", v_ids[0], v_ids[1], t_ids.size());

        while(!er_queue.empty()){
            std::array<int, 2> tmp_v_ids = er_queue.top().v_ids;
            if(tmp_v_ids==v_ids)
                er_queue.pop();
            else
                break;
        }

        bool is_fail=false;
        if(removeAnEdge_32(v_ids[0], v_ids[1], t_ids))
            suc_counter++;
        else if(removeAnEdge_44(v_ids[0], v_ids[1], t_ids))
            suc_counter++;
        else if(removeAnEdge_56(v_ids[0], v_ids[1], t_ids)) {
            suc_counter++;
            cnt5++;
        } else{
            is_fail=true;
        }

//        if(is_fail){
//            logger().debug("f");
//        } else
//            logger().debug("s");

        counter++;
    }
    logger().debug("tmp_cnt3 = {}", tmp_cnt3);
    logger().debug("tmp_cnt4 = {}", tmp_cnt4);
    logger().debug("tmp_cnt5 = {}", tmp_cnt5);
    logger().debug("tmp_cnt6 = {}", tmp_cnt6);
    logger().debug("{}", cnt5);

    logger().debug("energy_time = {}", energy_time);
}

bool EdgeRemover::removeAnEdge_32(int v1_id, int v2_id, const std::vector<int>& old_t_ids) {
    if(old_t_ids.size() >= 6) tmp_cnt6++;
    if(old_t_ids.size() == 5) tmp_cnt5++;
    if(old_t_ids.size() == 4) tmp_cnt4++;
    if(old_t_ids.size() == 3) tmp_cnt3++;

    if (old_t_ids.size() != 3)
        return false;

    //new_tets
    std::array<int, 2> v_ids;
    std::vector<std::array<int, 4>> new_tets;
    std::array<int, 2> t_ids;
    int cnt = 0;
    for (int i = 0; i < 4; i++) {
        if (tets[old_t_ids[0]][i] != v1_id && tets[old_t_ids[0]][i] != v2_id) {
            v_ids[cnt++] = tets[old_t_ids[0]][i];
        }
    }
    auto it = std::find(tets[old_t_ids[1]].begin(), tets[old_t_ids[1]].end(), v_ids[0]);
    if (it != tets[old_t_ids[1]].end()) {
        new_tets.push_back(tets[old_t_ids[1]]);
        new_tets.push_back(tets[old_t_ids[2]]);
        t_ids = {{old_t_ids[1], old_t_ids[2]}};
    } else {
        new_tets.push_back(tets[old_t_ids[2]]);
        new_tets.push_back(tets[old_t_ids[1]]);
        t_ids = {{old_t_ids[2], old_t_ids[1]}};
    }
    it = std::find(new_tets[0].begin(), new_tets[0].end(), v1_id);
    *it = v_ids[1];
    it = std::find(new_tets[1].begin(), new_tets[1].end(), v2_id);
    *it = v_ids[0];

    //check is_valid
    std::vector<TetQuality> tet_qs;
    if(isFlip(new_tets))
        return false;
    TetQuality old_tq, new_tq;
    getCheckQuality(old_t_ids, old_tq);
    tmp_timer.start();
    calTetQualities(new_tets, tet_qs);
    energy_time+=tmp_timer.getElapsedTime();
    getCheckQuality(tet_qs, new_tq);
    if(equal_buget>0) {
        equal_buget--;
        if (!new_tq.isBetterOrEqualThan(old_tq, energy_type, state))
            return false;
    } else {
        if (!new_tq.isBetterThan(old_tq, energy_type, state))
            return false;
    }

    //real update
    std::vector<std::array<int, 3>> fs;
    std::vector<int> is_sf_fs;
    for(int i=0;i<old_t_ids.size();i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] == v1_id || tets[old_t_ids[i]][j] == v2_id) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][(j + 1) % 4], tets[old_t_ids[i]][(j + 2) % 4],
                                           tets[old_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                fs.push_back(tmp);
                is_sf_fs.push_back(is_surface_fs[old_t_ids[i]][j]);
            }
        }
    }

    t_is_removed[old_t_ids[0]] = true;
    tets[t_ids[0]] = new_tets[0];//v2
    tets[t_ids[1]] = new_tets[1];//v1

    for(int i=0;i<4;i++) {
        if (tets[t_ids[0]][i] != v2_id) {
            std::array<int, 3> tmp = {{tets[t_ids[0]][(i + 1) % 4], tets[t_ids[0]][(i + 2) % 4],
                                       tets[t_ids[0]][(i + 3) % 4]}};
            std::sort(tmp.begin(), tmp.end());
            auto it = std::find(fs.begin(), fs.end(), tmp);
            is_surface_fs[t_ids[0]][i] = is_sf_fs[it - fs.begin()];
        } else
            is_surface_fs[t_ids[0]][i] = state.NOT_SURFACE;

        if (tets[t_ids[1]][i] != v1_id) {
            std::array<int, 3> tmp = {{tets[t_ids[1]][(i + 1) % 4], tets[t_ids[1]][(i + 2) % 4],
                                       tets[t_ids[1]][(i + 3) % 4]}};
            std::sort(tmp.begin(), tmp.end());
            auto it = std::find(fs.begin(), fs.end(), tmp);
            is_surface_fs[t_ids[1]][i] = is_sf_fs[it - fs.begin()];
        } else
            is_surface_fs[t_ids[1]][i] = state.NOT_SURFACE;
    }

    tet_vertices[v_ids[0]].conn_tets.erase(std::find(tet_vertices[v_ids[0]].conn_tets.begin(),
                                                     tet_vertices[v_ids[0]].conn_tets.end(), old_t_ids[0]));
    tet_vertices[v_ids[1]].conn_tets.erase(std::find(tet_vertices[v_ids[1]].conn_tets.begin(),
                                                     tet_vertices[v_ids[1]].conn_tets.end(), old_t_ids[0]));

    tet_vertices[v_ids[0]].conn_tets.insert(t_ids[1]);
    tet_vertices[v_ids[1]].conn_tets.insert(t_ids[0]);

    tet_vertices[v1_id].conn_tets.erase(std::find(tet_vertices[v1_id].conn_tets.begin(),
                                                  tet_vertices[v1_id].conn_tets.end(), old_t_ids[0]));
    tet_vertices[v2_id].conn_tets.erase(std::find(tet_vertices[v2_id].conn_tets.begin(),
                                                  tet_vertices[v2_id].conn_tets.end(), old_t_ids[0]));

    tet_vertices[v1_id].conn_tets.erase(std::find(tet_vertices[v1_id].conn_tets.begin(),
                                                  tet_vertices[v1_id].conn_tets.end(), t_ids[0]));
    tet_vertices[v2_id].conn_tets.erase(std::find(tet_vertices[v2_id].conn_tets.begin(),
                                                  tet_vertices[v2_id].conn_tets.end(), t_ids[1]));

    for (int i = 0; i < 2; i++) {
        tet_qualities[t_ids[i]] = tet_qs[i];
    }

    //repush new edges
    //Note that you need to pop out the current element first!!
    std::unordered_set<int> n12_v_ids;
    for(int i=0;i<new_tets.size();i++){
        for(int j=0;j<4;j++){
            if(new_tets[i][j]!=v1_id && new_tets[i][j]!=v2_id)
                n12_v_ids.insert(new_tets[i][j]);
        }
    }

//    for(auto it=n12_v_ids.begin();it!=n12_v_ids.end();it++) {
//        addNewEdge(std::array<int, 2>({{*it, v1_id}}));
//        addNewEdge(std::array<int, 2>({{*it, v2_id}}));
//    }

    std::vector<std::array<int, 2>> es;
    es.reserve(new_tets.size()*6);
    for(int i=0;i<new_tets.size();i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{new_tets[i][0], new_tets[i][j + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                es.push_back(e);
            e = {{new_tets[i][j + 1], new_tets[i][(j + 1) % 3 + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                es.push_back(e);
        }
    }
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    for(int i=0;i<es.size();i++){
        addNewEdge(es[i]);
    }

    return true;
}

bool EdgeRemover::removeAnEdge_44(int v1_id, int v2_id, const std::vector<int>& old_t_ids) {
    const int N = 4;
    if (old_t_ids.size() != N)
        return false;

    std::vector<std::array<int, 3>> n12_es;
    n12_es.reserve(old_t_ids.size());
    for (int i = 0; i < old_t_ids.size(); i++) {
        std::array<int, 3> e;
        int cnt = 0;
        for (int j = 0; j < 4; j++)
            if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id) {
                e[cnt++] = tets[old_t_ids[i]][j];
            }
        e[cnt] = old_t_ids[i];
        n12_es.push_back(e);
    }

    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    n12_v_ids.push_back(n12_es[0][0]);
    n12_v_ids.push_back(n12_es[0][1]);
    n12_t_ids.push_back(n12_es[0][2]);
    std::vector<bool> is_visited(N, false);
    is_visited[0] = true;
    for (int i = 0; i < N - 2; i++) {
        for (int j = 0; j < N; j++) {
            if (!is_visited[j]) {
                if (n12_es[j][0] == n12_v_ids.back()) {
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][1]);
                } else if (n12_es[j][1] == n12_v_ids.back()) {//else if!!!!!!!!!!
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][0]);
                }
                if (is_visited[j]) {
                    n12_t_ids.push_back(n12_es[j][2]);
                    break;
                }
            }
        }
    }
    n12_t_ids.push_back(n12_es[std::find(is_visited.begin(), is_visited.end(), false) - is_visited.begin()][2]);

    bool is_valid = false;
    std::vector<std::array<int, 4>> new_tets;
    new_tets.reserve(4);
    std::vector<int> tags;
    std::vector<TetQuality> tet_qs;
    std::array<int, 2> v_ids;
    TetQuality old_tq, new_tq;
    getCheckQuality(old_t_ids, old_tq);
    for (int i = 0; i < 2; i++) {
        std::vector<std::array<int, 4>> tmp_new_tets;
        std::vector<int> tmp_tags;
        std::vector<TetQuality> tmp_tet_qs;
        std::array<int, 2> tmp_v_ids;
        tmp_v_ids = {{n12_v_ids[0 + i], n12_v_ids[2 + i]}};
        for (int j = 0; j < old_t_ids.size(); j++) {
            std::array<int, 4> t = tets[old_t_ids[j]];
            auto it = std::find(t.begin(), t.end(), tmp_v_ids[0]);
            if (it != t.end()) {
                auto jt = std::find(t.begin(), t.end(), v2_id);
                *jt = tmp_v_ids[1];
                tmp_tags.push_back(1);
            } else {
                auto jt = std::find(t.begin(), t.end(), v1_id);
                *jt = tmp_v_ids[0];
                tmp_tags.push_back(0);
            }
            tmp_new_tets.push_back(t);
        }

        if (isFlip(tmp_new_tets))
            continue;
        tmp_timer.start();
        calTetQualities(tmp_new_tets, tmp_tet_qs);
        energy_time+=tmp_timer.getElapsedTime();
        getCheckQuality(tmp_tet_qs, new_tq);
        if(equal_buget>0) {
            equal_buget--;
            if (!new_tq.isBetterOrEqualThan(old_tq, energy_type, state))
                return false;
        } else {
            if (!new_tq.isBetterThan(old_tq, energy_type, state))
                return false;
        }

        is_valid = true;
        old_tq = new_tq;
        new_tets = tmp_new_tets;
        tags = tmp_tags;
        tet_qs = tmp_tet_qs;
        v_ids = tmp_v_ids;
    }
    if (!is_valid)
        return false;

    //real update
    std::vector<std::array<int, 3>> fs;
    std::vector<int> is_sf_fs;
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] == v1_id || tets[old_t_ids[i]][j] == v2_id) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][(j + 1) % 4], tets[old_t_ids[i]][(j + 2) % 4],
                                           tets[old_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                fs.push_back(tmp);
                is_sf_fs.push_back(is_surface_fs[old_t_ids[i]][j]);
            }
        }
    }

    for (int j = 0; j < new_tets.size(); j++) {
        if (tags[j] == 0) {
            tet_vertices[v1_id].conn_tets.erase(
                    std::find(tet_vertices[v1_id].conn_tets.begin(),
                              tet_vertices[v1_id].conn_tets.end(), old_t_ids[j]));
            tet_vertices[v_ids[0]].conn_tets.insert(old_t_ids[j]);
        } else {
            tet_vertices[v2_id].conn_tets.erase(
                    std::find(tet_vertices[v2_id].conn_tets.begin(),
                              tet_vertices[v2_id].conn_tets.end(), old_t_ids[j]));
            tet_vertices[v_ids[1]].conn_tets.insert(old_t_ids[j]);
        }
        tets[old_t_ids[j]] = new_tets[j];
        tet_qualities[old_t_ids[j]] = tet_qs[j];
    }

    for (int i = 0; i < old_t_ids.size(); i++) {//old_t_ids contains new tets
        for (int j = 0; j < 4; j++) {
            is_surface_fs[old_t_ids[i]][j] = state.NOT_SURFACE;
            if (tets[old_t_ids[i]][j] == v_ids[0] || tets[old_t_ids[i]][j] == v_ids[1]) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][(j + 1) % 4], tets[old_t_ids[i]][(j + 2) % 4],
                                           tets[old_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                auto it = std::find(fs.begin(), fs.end(), tmp);
                is_surface_fs[old_t_ids[i]][j] = is_sf_fs[it - fs.begin()];
            }
        }
    }

    //repush
    std::vector<std::array<int, 2>> es;
    es.reserve(new_tets.size()*6);
    for (int i = 0; i < new_tets.size(); i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{new_tets[i][0], new_tets[i][j + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                es.push_back(e);
            e = {{new_tets[i][j + 1], new_tets[i][(j + 1) % 3 + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                es.push_back(e);
        }
    }
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    for (int i = 0; i < es.size(); i++) {
//        if (es[i][0] != v_ids[0] && es[i][1] != v_ids[0] && es[i][0] != v_ids[1] && es[i][1] != v_ids[1])
            addNewEdge(es[i]);
    }

    return true;
}

//void getMinQuality(const std::vector<TetQuality>& qs, double& d_min, double& d_max, double& r_max) {
//    d_min = 10;
//    d_max = 0;
//    r_max = 0;
//    for (int i = 0; i < qs.size(); i++) {
//        if (qs[i].min_d_angle < d_min)
//            d_min = qs[i].min_d_angle;
//        if (qs[i].max_d_angle > d_max)
//            d_max = qs[i].max_d_angle;
//        if (qs[i].asp_ratio_2 > r_max)
//            r_max = qs[i].asp_ratio_2;
//    }
//}

bool EdgeRemover::removeAnEdge_56(int v1_id, int v2_id, const std::vector<int>& old_t_ids) {
    if (old_t_ids.size() != 5)
        return false;

    //oriented the n12_v_ids
    std::vector<std::array<int, 3>> n12_es;
    n12_es.reserve(old_t_ids.size());
    for (int i = 0; i < old_t_ids.size(); i++) {
        std::array<int, 3> e;
        int cnt = 0;
        for (int j = 0; j < 4; j++)
            if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id) {
                e[cnt++] = tets[old_t_ids[i]][j];
            }
        e[cnt] = old_t_ids[i];
        n12_es.push_back(e);
    }

    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    n12_v_ids.push_back(n12_es[0][0]);
    n12_v_ids.push_back(n12_es[0][1]);
    n12_t_ids.push_back(n12_es[0][2]);
    std::vector<bool> is_visited(5, false);
    is_visited[0] = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 5; j++) {
            if (!is_visited[j]) {
                if (n12_es[j][0] == n12_v_ids.back()) {
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][1]);
                } else if (n12_es[j][1] == n12_v_ids.back()) {//else if!!!!!!!!!!
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][0]);
                }
                if (is_visited[j]) {
                    n12_t_ids.push_back(n12_es[j][2]);
                    break;
                }
            }
        }
    }
    n12_t_ids.push_back(n12_es[std::find(is_visited.begin(), is_visited.end(), false) - is_visited.begin()][2]);


    //check valid
    TetQuality old_tq, new_tq;
    getCheckQuality(old_t_ids, old_tq);
    std::unordered_map<int, std::array<TetQuality, 2>> tet_qs;
    std::unordered_map<int, std::array<std::array<int, 4>, 2>> new_tets;
    std::vector<bool> is_v_valid(5, true);
    for (int i = 0; i < n12_v_ids.size(); i++) {
        if (!is_v_valid[(i + 1) % 5] && !is_v_valid[(i - 1 + 5) % 5])
            continue;

        std::vector<std::array<int, 4>> new_ts;
        new_ts.reserve(6);
        std::array<int, 4> t = tets[n12_t_ids[i]];
        auto it = std::find(t.begin(), t.end(), v1_id);
        *it = n12_v_ids[(i - 1 + 5) % 5];
        new_ts.push_back(t);

        t = tets[n12_t_ids[i]];
        it = std::find(t.begin(), t.end(), v2_id);
        *it = n12_v_ids[(i - 1 + 5) % 5];
        new_ts.push_back(t);
        if (isFlip(new_ts)) {
            is_v_valid[(i + 1) % 5] = false;
            is_v_valid[(i - 1 + 5) % 5] = false;
            continue;
        }

        std::vector<TetQuality> qs;
        tmp_timer.start();
        calTetQualities(new_ts, qs);
        energy_time+=tmp_timer.getElapsedTime();
        tet_qs[i] = std::array<TetQuality, 2>({{qs[0], qs[1]}});
        new_tets[i] = std::array<std::array<int, 4>, 2>({{new_ts[0], new_ts[1]}});

//        std::vector<TetQuality> qs;
//        calTetQualities(new_ts, qs);
//        getCheckQuality(qs, new_tq);
//        if(new_tq.isBetterThan(old_tq, energy_type)){
//            tet_qs[i] = std::array<TetQuality, 2>({{qs[0], qs[1]}});
//            new_tets[i] = std::array<std::array<int, 4>, 2>({{new_ts[0], new_ts[1]}});
//        } else {
//            is_v_valid[(i + 1) % 5] = false;
//            is_v_valid[(i - 1 + 5) % 5] = false;
//        }
    }
    if (std::count(is_v_valid.begin(), is_v_valid.end(), true) == 0)
        return false;


    int selected_id = -1;
    for (int i = 0; i < is_v_valid.size(); i++) {
        if (!is_v_valid[i])
            continue;

        std::vector<std::array<int, 4>> new_ts;
        new_ts.reserve(6);
        std::array<int, 4> t = tets[n12_t_ids[(i + 2) % 5]];
        auto it = std::find(t.begin(), t.end(), v1_id);
        *it = n12_v_ids[i];
        new_ts.push_back(t);
        t = tets[n12_t_ids[(i + 2) % 5]];
        it = std::find(t.begin(), t.end(), v2_id);
        *it = n12_v_ids[i];
        new_ts.push_back(t);
        if (isFlip(new_ts))
            continue;

        std::vector<TetQuality> qs;
        tmp_timer.start();
        calTetQualities(new_ts, qs);
        energy_time+=tmp_timer.getElapsedTime();
        for (int j = 0; j < 2; j++) {
            qs.push_back(tet_qs[(i + 1) % 5][j]);
            qs.push_back(tet_qs[(i - 1 + 5) % 5][j]);
        }
        if(qs.size() != 6){
            log_and_throw("qs.size() != 6");
        }
        getCheckQuality(qs, new_tq);
        if(equal_buget>0) {
            equal_buget--;
            if (!new_tq.isBetterOrEqualThan(old_tq, energy_type, state))
                continue;
        } else {
            if (!new_tq.isBetterThan(old_tq, energy_type, state))
                continue;
        }

        old_tq = new_tq;
        selected_id = i;
        tet_qs[i + 5] = std::array<TetQuality, 2>({{qs[0], qs[1]}});
        new_tets[i + 5] = std::array<std::array<int, 4>, 2>({{new_ts[0], new_ts[1]}});
    }
    if (selected_id < 0)
        return false;

    //real update
    //update on surface -- 1
    std::vector<std::array<int, 3>> fs;
    std::vector<int> is_sf_fs;
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] == v1_id || tets[old_t_ids[i]][j] == v2_id) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][(j + 1) % 4], tets[old_t_ids[i]][(j + 2) % 4],
                                           tets[old_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                fs.push_back(tmp);
                is_sf_fs.push_back(is_surface_fs[old_t_ids[i]][j]);
            }
        }
    }

    std::vector<int> new_t_ids = old_t_ids;
    getNewTetSlots(1, new_t_ids);
    t_is_removed[new_t_ids.back()] = false;
    for (int i = 0; i < 2; i++) {
        tets[new_t_ids[i]] = new_tets[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]] = new_tets[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]] = new_tets[selected_id + 5][i];

        tet_qualities[new_t_ids[i]] = tet_qs[(selected_id + 1) % 5][i];
        tet_qualities[new_t_ids[i + 2]] = tet_qs[(selected_id - 1 + 5) % 5][i];
        tet_qualities[new_t_ids[i + 4]] = tet_qs[selected_id + 5][i];
    }

    //update on_surface -- 2
    for (int i = 0; i < new_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            is_surface_fs[new_t_ids[i]][j] = state.NOT_SURFACE;
            if (tets[new_t_ids[i]][j] != v1_id && tets[new_t_ids[i]][j] != v2_id
                && tets[new_t_ids[i]][j] != n12_v_ids[(selected_id + 1) % 5]
                && tets[new_t_ids[i]][j] != n12_v_ids[(selected_id - 1 + 5) % 5]) {
                std::array<int, 3> tmp = {{tets[new_t_ids[i]][(j + 1) % 4], tets[new_t_ids[i]][(j + 2) % 4],
                                           tets[new_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                auto it = std::find(fs.begin(), fs.end(), tmp);
                if (it != fs.end())
                    is_surface_fs[new_t_ids[i]][j] = is_sf_fs[it - fs.begin()];
            }
        }
    }

    //update conn_tets
    for (int i = 0; i < n12_v_ids.size(); i++) {
        tet_vertices[n12_v_ids[i]].conn_tets.erase(n12_t_ids[i]);
        tet_vertices[n12_v_ids[i]].conn_tets.erase(n12_t_ids[(i - 1 + 5) % 5]);
    }
    for (int i = 0; i < n12_t_ids.size(); i++) {
        tet_vertices[v1_id].conn_tets.erase(n12_t_ids[i]);
        tet_vertices[v2_id].conn_tets.erase(n12_t_ids[i]);
    }

    //add
    for (int i = 0; i < new_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++)
            tet_vertices[tets[new_t_ids[i]][j]].conn_tets.insert(new_t_ids[i]);
    }

    //repush
//    addNewEdge(std::array<int, 2>({{n12_v_ids[selected_id], n12_v_ids[(selected_id + 2) % 5]}}));
//    addNewEdge(std::array<int, 2>({{n12_v_ids[selected_id], n12_v_ids[(selected_id - 2 + 5) % 5]}}));
//
//    addNewEdge(std::array<int, 2>({{v1_id, n12_v_ids[(selected_id + 1) % 5]}}));
//    addNewEdge(std::array<int, 2>({{v1_id, n12_v_ids[(selected_id - 1 + 5) % 5]}}));
//    addNewEdge(std::array<int, 2>({{v2_id, n12_v_ids[(selected_id + 1) % 5]}}));
//    addNewEdge(std::array<int, 2>({{v2_id, n12_v_ids[(selected_id - 1 + 5) % 5]}}));

    std::vector<std::array<int, 2>> es;
    es.reserve(new_t_ids.size()*6);
    for(int i=0;i<new_t_ids.size();i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{tets[new_t_ids[i]][0], tets[new_t_ids[i]][j + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                es.push_back(e);
            e = {{tets[new_t_ids[i]][j + 1], tets[new_t_ids[i]][(j + 1) % 3 + 1]}};
            if(e[0]>e[1]) e={{e[1], e[0]}};
            if(!isLocked_ui(e))
                es.push_back(e);
        }
    }
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    for(int i=0;i<es.size();i++){
        addNewEdge(es[i]);
    }

    return true;
}

bool EdgeRemover::isSwappable_cd1(const std::array<int, 2>& v_ids){
    std::vector<int> t_ids;
    setIntersection(tet_vertices[v_ids[0]].conn_tets, tet_vertices[v_ids[1]].conn_tets, t_ids);

    if(isEdgeOnSurface(v_ids[0], v_ids[1], t_ids))
        return false;
    if(isEdgeOnBbox(v_ids[0], v_ids[1], t_ids))
        return false;

    return true;
}

bool EdgeRemover::isSwappable_cd1(const std::array<int, 2>& v_ids, std::vector<int>& t_ids, bool is_check_conn_tet_num){
//    std::vector<int> t_ids;
    setIntersection(tet_vertices[v_ids[0]].conn_tets, tet_vertices[v_ids[1]].conn_tets, t_ids);

    if(is_check_conn_tet_num)
        if(t_ids.size()<3 || t_ids.size()>5)
            return false;
    if(isEdgeOnSurface(v_ids[0], v_ids[1], t_ids))
        return false;
    if(isEdgeOnBbox(v_ids[0], v_ids[1], t_ids))
        return false;

    return true;
}

bool EdgeRemover::isSwappable_cd2(double weight){
    return true;

    if(weight>ideal_weight)
        return true;
    return false;
}

bool EdgeRemover::isEdgeValid(const std::array<int, 2>& v_ids){
    if(v_is_removed[v_ids[0]] || v_is_removed[v_ids[1]])
        return false;
    if(!isHaveCommonEle(tet_vertices[v_ids[0]].conn_tets, tet_vertices[v_ids[1]].conn_tets))
        return false;
    return true;
}

void EdgeRemover::getNewTetSlots(int n, std::vector<int>& new_conn_tets) {
    unsigned int cnt = 0;
    for (unsigned int i = t_empty_start; i < t_is_removed.size(); i++) {
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
        for (unsigned int i = 0; i < n - cnt; i++)
            new_conn_tets.push_back(tets.size() + i);

        tets.resize(tets.size() + n - cnt);
        t_is_removed.resize(t_is_removed.size() + n - cnt);
        tet_qualities.resize(tet_qualities.size() + n - cnt);
        is_surface_fs.resize(is_surface_fs.size() + n - cnt);
        t_empty_start = tets.size();
    }
}

void EdgeRemover::addNewEdge(const std::array<int, 2>& e){
    if (isSwappable_cd1(e)) {
        double weight = calEdgeLength(e);
        if (isSwappable_cd2(weight)) {
            if (e[0] > e[1]) {
                ElementInQueue_er ele(std::array<int, 2>({{e[1], e[0]}}), weight);
                er_queue.push(ele);
            } else {
                ElementInQueue_er ele(e, weight);
                er_queue.push(ele);
            }
        }
    }
}

} // namespace tetwild
