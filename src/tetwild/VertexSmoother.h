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

#ifndef NEW_GTET_VERTEXSMOOTHER_H
#define NEW_GTET_VERTEXSMOOTHER_H

#include <tetwild/LocalOperations.h>

namespace tetwild {

class VertexSmoother:public LocalOperations {
public:
    VertexSmoother(LocalOperations lo): LocalOperations(lo){}

    void smooth();
    void smoothSingle();
    bool smoothSingleVertex(int v_id, bool is_cal_energy);
    void smoothSurface();

    bool NewtonsMethod(const std::vector<int>& t_ids, const std::vector<std::array<int, 4>>& new_tets, int v_id, Point_3f& p);
    bool NewtonsUpdate(const std::vector<int>& t_ids, int v_id, double& energy, Eigen::Vector3d& J, Eigen::Matrix3d& H, Eigen::Vector3d& X0);
    double getNewEnergy(const std::vector<int>& t_ids);

    int ts;
    std::vector<int> tets_tss;
    std::vector<int> tet_vertices_tss;

    void outputOneRing(int v_i, std::string s);
    //for postprocessing
    int laplacianBoundary(const std::vector<int>& b_v_ids, const std::vector<bool>& tmp_is_on_surface,
                          const std::vector<bool>& tmp_t_is_removed);

    int id_value_e=0;
    int id_value_j=1;
    int id_value_h=2;
    int id_solve=3;
    int id_aabb=4;
    int id_project = 5;
    int id_round = 6;
    std::array<double, 7> breakdown_timing={{0,0,0,0,0,0,0}};
    std::array<std::string, 7> breakdown_name={{"Computing E", "Computing J", "Computing H", "Solving linear system", "AABB", "Project", "Rounding"}};
    igl::Timer igl_timer;
};

} // namespace tetwild

#endif //NEW_GTET_VERTEXSMOOTHER_H
