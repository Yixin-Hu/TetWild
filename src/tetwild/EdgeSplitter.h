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

#ifndef NEW_GTET_EDGESPLITTER_H
#define NEW_GTET_EDGESPLITTER_H

#include <tetwild/LocalOperations.h>
#include <queue>

namespace tetwild {

class ElementInQueue_es{
public:
    std::array<int, 2> v_ids;
    double weight;

    ElementInQueue_es(){}
    ElementInQueue_es(const std::array<int, 2>& ids, double w):
            v_ids(ids), weight(w){}
};

struct cmp_es {
    bool operator()(const ElementInQueue_es &e1, const ElementInQueue_es &e2) {
        return e1.weight < e2.weight;
    }
};

class EdgeSplitter:public LocalOperations {
public:
    bool is_check_quality = false;
    bool is_cal_quality_end = false;

    std::priority_queue<ElementInQueue_es, std::vector<ElementInQueue_es>, cmp_es> es_queue;

    int t_empty_start=0;
    int v_empty_start=0;

    double max_weight=0;
    double ideal_weight=0;

    EdgeSplitter(LocalOperations lo, double ideal_w): LocalOperations(lo), ideal_weight(ideal_w){}

    void init();
    void split();

    bool is_over_refine=false;
    int getOverRefineScale(int v1_id, int v2_id);
    bool splitAnEdge(const std::array<int, 2>& edge);

    bool isSplittable_cd1(double weight);
    bool isSplittable_cd1(int v1_id, int v2_id, double weight);
    void getNewTetSlots(int n, std::vector<int>& new_conn_tets);
//    igl::viewer::Viewer viewer;
    void getMesh_ui(const std::vector<std::array<int, 4>>& tets, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

    unsigned int budget = 0;
};

} // namespace tetwild

#endif //NEW_GTET_EDGESPLITTER_H
