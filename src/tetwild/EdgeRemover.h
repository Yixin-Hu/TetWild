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

#ifndef NEW_GTET_EDGEREMOVER_H
#define NEW_GTET_EDGEREMOVER_H

#include <tetwild/LocalOperations.h>
#include <queue>

namespace tetwild {

class ElementInQueue_er{
public:
    std::array<int, 2> v_ids;
    double weight;

    ElementInQueue_er(){}
    ElementInQueue_er(const std::array<int, 2>& ids, double w): v_ids(ids), weight(w){}
};

struct cmp_er {
    bool operator()(const ElementInQueue_er &e1, const ElementInQueue_er &e2) {
        if (e1.weight == e2.weight)
            return e1.v_ids < e2.v_ids;
        return e1.weight < e2.weight;///choose larger edge for removal
    }
};

class EdgeRemover:public LocalOperations {
public:
    std::priority_queue<ElementInQueue_er, std::vector<ElementInQueue_er>, cmp_er> er_queue;

    double ideal_weight;

    int v_empty_start=0;
    int t_empty_start=0;

    int flag_cnt=0;

    int tmp_cnt3=0;
    int tmp_cnt4=0;
    int tmp_cnt5=0;
    int tmp_cnt6=0;

    int equal_buget = 100;

    EdgeRemover(LocalOperations lo, double i_weight): LocalOperations(lo), ideal_weight(i_weight){}

    void init();
    void swap();
    bool removeAnEdge_32(int v1_id, int v2_id, const std::vector<int>& old_t_ids);
    bool removeAnEdge_44(int v1_id, int v2_id, const std::vector<int>& old_t_ids);
    bool removeAnEdge_56(int v1_id, int v2_id, const std::vector<int>& old_t_ids);

    bool isSwappable_cd1(const std::array<int, 2>& v_ids, std::vector<int>& t_ids, bool is_check_conn_tet_num=false);
    bool isSwappable_cd1(const std::array<int, 2>& v_ids);

    bool isSwappable_cd2(double weight);
    bool isEdgeValid(const std::array<int, 2>& v_ids);
    void getNewTetSlots(int n, std::vector<int>& new_conn_tets);

    void addNewEdge(const std::array<int, 2>& e);

    igl::Timer tmp_timer;
    double energy_time = 0;
};

} // namespace tetwild

#endif //NEW_GTET_EDGEREMOVER_H
