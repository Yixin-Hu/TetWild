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

#ifndef NEW_GTET_EDGECOLLAPSER_H
#define NEW_GTET_EDGECOLLAPSER_H

#include <tetwild/LocalOperations.h>
#include <queue>

namespace tetwild {

class ElementInQueue_ec{
public:
    std::array<int, 2> v_ids;
    double weight;

    ElementInQueue_ec(){}
    ElementInQueue_ec(const std::array<int, 2>& ids, double w):
            v_ids(ids), weight(w){}
};

struct cmp_ec {
    bool operator()(const ElementInQueue_ec &e1, const ElementInQueue_ec &e2) {
        if (e1.weight == e2.weight)
            return e1.v_ids < e2.v_ids;
        return e1.weight > e2.weight;
    }
};

class EdgeCollapser: public LocalOperations {
public:
    std::priority_queue<ElementInQueue_ec, std::vector<ElementInQueue_ec>, cmp_ec> ec_queue;
    double ideal_weight=0;

    bool is_limit_length=true;
    bool is_check_quality=true;

    int envelop_accept_cnt=0;
    EdgeCollapser(LocalOperations& lo, double ideal_w): LocalOperations(lo), ideal_weight(ideal_w){}

    void init();
    void collapse();

    const int SUCCESS=0;
    const int FLIP=1;
    const int QUALITY=2;
    const int ENVELOP=3;
    const int ENVELOP_SUC=4;
    int collapseAnEdge(int v1_id, int v2_id);

    bool is_soft = false;
    double soft_energy = 6;
    int budget = 0;

    int ts = 0;
    std::vector<std::array<int, 2>> inf_es;
    std::vector<int> inf_e_tss;
    std::vector<int> tet_tss;
    void postProcess();

    bool isCollapsable_cd1(int v1_id, int v2_id);
//    bool isCollapsable_cd2(int v1_id, int v2_id);//check if a vertex is outside the envelop
//    bool isCollapsable_cd3(double weight);
    bool isCollapsable_cd3(int v1_id, int v2_id, double weight);
    bool isCollapsable_epsilon(int v1_id, int v2_id);

    bool isEdgeValid(const std::array<int, 2>& e);
//    void addNewEdge(const std::array<int, 2>& e);

    int tmp=0;
    int tmp0=0;

    double energy_time = 0;

    //for timing
    int id_postprocessing=0;
    int id_flip_fail=1;
    int id_env_fail=2;
    int id_success=3;
    int id_env_success=4;
    int id_energy_fail = 5;
    std::array<double, 6> breakdown_timing;
    std::array<std::string, 6> breakdown_name={{"Postprocessing",
                                               "Failed (flip)",
                                               "Failed (envelop)",
                                               "Successful (non-surface)",
                                               "Successful (surface)",
                                               "Failed (energy)"}};
    igl::Timer igl_timer;
};

} // namespace tetwild

#endif //NEW_GTET_EDGECOLLAPSER_H
