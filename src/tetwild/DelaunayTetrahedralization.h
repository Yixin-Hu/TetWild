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

#ifndef GTET_DELAUNAYTETRAHEDRALIZATION_H
#define GTET_DELAUNAYTETRAHEDRALIZATION_H

#include <tetwild/ForwardDecls.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/BSPElements.h>
#include <geogram/mesh/mesh.h>
#include <tetwild/DisableWarnings.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <tetwild/EnableWarnings.h>
#include <Eigen/Dense>

namespace tetwild {

// More CGAL types
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point_d;

class DelaunayTetrahedralization {
public:
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;

    void init(const std::vector<Point_3>& m_vertices, const std::vector<std::array<int, 3>>& m_faces,
                  std::vector<int>& m_f_tags, std::vector<int>& raw_e_tags, std::vector<std::vector<int>>& raw_conn_e4v);

    void getVoxelPoints(const Point_3& p_min, const Point_3& p_max, GEO::Mesh& geo_surface_mesh,
                        std::vector<Point_d>& voxel_points, const Args &args, const State &state);
    void tetra(const std::vector<Point_3>& m_vertices, GEO::Mesh& geo_surface_mesh,
               std::vector<Point_3>& bsp_vertices, std::vector<BSPEdge>& bsp_edges,
               std::vector<BSPFace>& bsp_faces, std::vector<BSPtreeNode>& bsp_nodes,
               const Args &args, const State &state);
    void outputTetmesh(const std::vector<Point_3>& m_vertices, std::vector<std::array<int, 4>>& cells,
                       const std::string& output_file);
};

} // namespace tetwild

#endif //GTET_DELAUNAYTETRAHEDRALIZATION_H
