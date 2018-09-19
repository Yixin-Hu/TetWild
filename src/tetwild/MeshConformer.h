// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 3/29/17.
//

#ifndef GTET_MESHCONFORMER_H
#define GTET_MESHCONFORMER_H

#include <tetwild/CGALTypes.h>
#include <tetwild/BSPElements.h>

namespace tetwild {

class MeshConformer {
public:
    const std::vector<Point_3> &m_vertices;
    const std::vector<std::array<int, 3>> &m_faces;
    std::vector<bool> is_matched;

    std::vector<Point_3>& bsp_vertices;
    std::vector<BSPEdge>& bsp_edges;
    std::vector<BSPFace>& bsp_faces;
    std::vector<BSPtreeNode>& bsp_nodes;

    MeshConformer(const std::vector<Point_3> &m_vs, const std::vector<std::array<int, 3>> &m_fs, std::vector<Point_3>& bsp_vs,
                  std::vector<BSPEdge>& bsp_es, std::vector<BSPFace>& bsp_fs, std::vector<BSPtreeNode>& bsp_ns) :
            m_vertices(m_vs), m_faces(m_fs), bsp_vertices(bsp_vs), bsp_edges(bsp_es), bsp_faces(bsp_fs), bsp_nodes(bsp_ns){}

    void match();
    void matchVertexIndices(int x, const std::vector<std::array<int, 2>>& seed_v_list, std::vector<int>& f_list);
    void matchDivFaces();
    void getOrientedVertices(int bsp_f_id);

    const int COPLANAR_INT=0;
    const int CROSS_INT=1;
    const int POINT_INT=2;
    const int NONE_INT=3;
    int triangleIntersection3d(const std::array<Point_3, 3>& tri1, const std::array<Point_3, 3>& tri2,
                               bool intersect_known=true);

    int t = 0;
    void initT(const Vector_3 &nv);
    Point_2 to2d(const Point_3 &p);
    Point_3 to3d(const Point_2 &p, const Plane_3 &pln);
};

} // namespace tetwild

#endif //GTET_MESHCONFORMER_H
