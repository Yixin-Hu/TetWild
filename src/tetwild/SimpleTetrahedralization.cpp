// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 3/31/17.
//

#include <tetwild/SimpleTetrahedralization.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>
#include <igl/winding_number.h>
#include <igl/Timer.h>
#include <bitset>

//centroid
#include <tetwild/DisableWarnings.h>
#include <CGAL/centroid.h>
#include <tetwild/EnableWarnings.h>

//triangulation
#include <tetwild/DisableWarnings.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <tetwild/EnableWarnings.h>
typedef CGAL::Constrained_Delaunay_triangulation_2<tetwild::K, CGAL::Default, CGAL::Exact_predicates_tag> CDT;
typedef CDT::Point Point_cdt_2;
typedef CGAL::Polygon_2<tetwild::K> Polygon_2;

//arrangement
#include <tetwild/DisableWarnings.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <tetwild/EnableWarnings.h>
typedef CGAL::Arr_segment_traits_2<tetwild::K> Traits_2;
typedef Traits_2::Point_2 Point_arr_2;
typedef Traits_2::X_monotone_curve_2 Segment_arr_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

namespace tetwild {

void SimpleTetrahedralization::tetra(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets) {
    std::vector<BSPFace> &faces = MC.bsp_faces;
    std::vector<Point_3> &vertices = MC.bsp_vertices;

    ///cal arrangement & tetrahedralization
    triangulation(tet_vertices, tets);
    logger().debug("#v = {} #t = {}", tet_vertices.size(), tets.size());

    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            tet_vertices[tets[i][j]].conn_tets.insert(i);
        }
    }
}

void SimpleTetrahedralization::triangulation(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets) {
    std::vector<BSPtreeNode> &bsp_nodes = MC.bsp_nodes;
    std::vector<BSPFace> &bsp_faces = MC.bsp_faces;
    std::vector<BSPEdge> &bsp_edges = MC.bsp_edges;
    std::vector<Point_3> &bsp_vertices = MC.bsp_vertices;
    const std::vector<Point_3>& m_vertices=MC.m_vertices;
    const std::vector<std::array<int, 3>>& m_faces=MC.m_faces;

    igl::Timer tmp_timer;

    tmp_timer.start();
    ///init
    m_vertices_size = MC.m_vertices.size();

    int original_v_size = bsp_vertices.size();

    ///cal arrangement
    for (int i = 0; i < bsp_faces.size(); i++) {
        if (bsp_faces[i].div_faces.size() == 0)
            continue;

        ///cal arrangement
        Plane_3 pln;
        constructPlane(i, pln);
        MC.initT(pln.orthogonal_vector());
        MC.getOrientedVertices(i);
        Polygon_2 poly;
        for (int j = 0; j < bsp_faces[i].vertices.size(); j++) {
            poly.push_back(MC.to2d(bsp_vertices[bsp_faces[i].vertices[j]]));
        }
        assert(poly.is_simple());

        Arrangement_2 arr;
        std::vector<Segment_arr_2> arr_segs;
        for (int j = 0; j < poly.size(); j++)
            arr_segs.push_back(Segment_arr_2(poly[j], poly[(j + 1) % poly.size()]));
        for (auto it = bsp_faces[i].div_faces.begin(); it != bsp_faces[i].div_faces.end(); it++) {
            std::array<Point_2, 3> df_vs;
            for (int j = 0; j < 3; j++)
                df_vs[j] = MC.to2d(MC.m_vertices[MC.m_faces[*it][j]]);
            for (int j = 0; j < 3; j++) {
                Line_2 l(df_vs[j], df_vs[(j + 1) % 3]);
                int cnt_pos = 0, cnt_neg = 0;
                for (int k = 0; k < poly.size(); k++) {
                    CGAL::Oriented_side side=l.oriented_side(poly[k]);
                    if(side==CGAL::ON_POSITIVE_SIDE)
                        cnt_pos++;
                    if(side==CGAL::ON_NEGATIVE_SIDE)
                        cnt_neg++;
                }
                if(cnt_pos>0 && cnt_neg>0)
                    arr_segs.push_back(Segment_arr_2(df_vs[j], df_vs[(j + 1) % 3]));
            }
        }
        CGAL::insert(arr, arr_segs.begin(), arr_segs.end());

        std::map<Point_2, int> vs_arr2bsp;
        std::vector<std::vector<int>> new_v_ids(bsp_faces[i].edges.size(), std::vector<int>());
        std::vector<Segment_3> segs(bsp_faces[i].edges.size(), Segment_3());
        for (int j = 0; j < bsp_faces[i].edges.size(); j++) {
            segs[j] = Segment_3(bsp_vertices[bsp_edges[bsp_faces[i].edges[j]].vertices[0]],
                                bsp_vertices[bsp_edges[bsp_faces[i].edges[j]].vertices[1]]);
        }
        for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
            if (poly.has_on_unbounded_side(it->point())) {
                vs_arr2bsp[it->point()] = -1;
                continue;
            }
            auto vit = std::find(poly.vertices_begin(), poly.vertices_end(), it->point());//todo
            if (vit != poly.vertices_end()) {
                int n = vit - poly.vertices_begin();
                vs_arr2bsp[it->point()] = bsp_faces[i].vertices[n];
                continue;
            }
            Point_3 p = MC.to3d(it->point(), pln);
            int on_e_local_id = -1;
            for (int j = 0; j < bsp_faces[i].edges.size(); j++) {
                if (segs[j].has_on(p)) {
                    on_e_local_id = j;
                    break;
                }
            }
            if (on_e_local_id >= 0) {//if the vertex is on the edges of this bsp_face
                bsp_vertices.push_back(p);
                vs_arr2bsp[it->point()] = bsp_vertices.size() - 1;
                new_v_ids[on_e_local_id].push_back(bsp_vertices.size() - 1);
                bsp_faces[i].vertices.push_back(bsp_vertices.size() - 1);
                continue;
            }
            //else
            bsp_vertices.push_back(p);
            vs_arr2bsp[it->point()] = bsp_vertices.size() - 1;
            bsp_faces[i].vertices.push_back(bsp_vertices.size() - 1);
        }

        std::vector<std::array<int, 2>> es;
        for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
            Point_2 &p1 = it->source()->point();
            Point_2 &p2 = it->target()->point();
            if (vs_arr2bsp[p1] < 0 || vs_arr2bsp[p2] < 0)
                continue;
            std::array<int, 2> e = {{vs_arr2bsp[p1], vs_arr2bsp[p2]}};
            std::sort(e.begin(), e.end());
            es.push_back(e);
        }
        std::sort(es.begin(), es.end());

        std::vector<std::array<int, 2>> tmp_es;
        int old_e_size = bsp_faces[i].edges.size();
        for (int j = 0; j < old_e_size; j++) {
            if (new_v_ids[j].size() == 0)
                continue;
            std::vector<int> new_es = bsp_edges[bsp_faces[i].edges[j]].vertices;
            new_es.insert(new_es.end(), new_v_ids[j].begin(), new_v_ids[j].end());
            std::sort(new_es.begin(), new_es.end(), [&](int a, int b) {
                return bsp_vertices[a] < bsp_vertices[b];
            });

            std::vector<int> new_e_ids;
            for (int k = 0; k < new_es.size() - 1; k++) {
                BSPEdge new_bsp_e(new_es[k], new_es[k + 1]);
                new_bsp_e.conn_faces = bsp_edges[bsp_faces[i].edges[j]].conn_faces;
                if (k == 0)
                    bsp_edges[bsp_faces[i].edges[j]] = new_bsp_e;
                else {
                    bsp_edges.push_back(new_bsp_e);
                    new_e_ids.push_back(bsp_edges.size() - 1);
                }

                std::array<int, 2> e = {{new_es[k], new_es[k + 1]}};
                std::sort(e.begin(), e.end());
                tmp_es.push_back(e);
            }

            for (auto it = bsp_edges[bsp_faces[i].edges[j]].conn_faces.begin();
                 it != bsp_edges[bsp_faces[i].edges[j]].conn_faces.end(); it++) {
                bsp_faces[*it].edges.insert(bsp_faces[*it].edges.end(), new_e_ids.begin(), new_e_ids.end());
                if (*it == i)
                    continue;
                bsp_faces[*it].vertices.insert(bsp_faces[*it].vertices.end(), new_v_ids[j].begin(), new_v_ids[j].end());
            }
        }
        std::sort(tmp_es.begin(), tmp_es.end());

        std::vector<std::array<int, 2>> diff_es;
        std::set_difference(es.begin(), es.end(), tmp_es.begin(), tmp_es.end(),
                            std::back_inserter(diff_es));
        for (int j = 0; j < diff_es.size(); j++) {
            BSPEdge e(diff_es[j][0], diff_es[j][1]);
            bsp_edges.push_back(e);
            bsp_faces[i].edges.push_back(bsp_edges.size() - 1);
        }
    }
    logger().debug("2D arr {}", tmp_timer.getElapsedTime());
    tmp_timer.start();

    tet_vertices.reserve(bsp_vertices.size() + bsp_nodes.size());
    for (unsigned int i = 0; i < bsp_vertices.size(); i++) {
        TetVertex v(bsp_vertices[i]);
        if (i < m_vertices_size)
            v.is_on_surface=true;
        tet_vertices.push_back(v);
    }

    ///improvement
    std::vector<bool> is_tets(bsp_nodes.size(), false);
    std::unordered_map<int, int> centroids_for_nodes;
    tets.reserve(bsp_nodes.size()*6);//approx
    for(unsigned int i = 0; i < bsp_nodes.size(); i++) {
        BSPtreeNode &node = bsp_nodes[i];
        std::vector<int> v_ids;
        for (int j = 0; j < node.faces.size(); j++) {
            for (int k = 0; k < bsp_faces[node.faces[j]].vertices.size(); k++)
                v_ids.push_back(bsp_faces[node.faces[j]].vertices[k]);
        }
        std::sort(v_ids.begin(), v_ids.end());
        v_ids.erase(std::unique(v_ids.begin(), v_ids.end()), v_ids.end());

        bool is_tet = false;
        if (bsp_nodes[i].faces.size() == 4) {
            is_tet = true;
            for (int j = 0; j < bsp_nodes[i].faces.size(); j++) {
                if (bsp_faces[bsp_nodes[i].faces[j]].vertices.size() != 3) {
                    is_tet = false;
                    break;
                }
            }
        }

        if (is_tet) {
            is_tets[i] = true;
            std::array<int, 4> t = {{v_ids[0], v_ids[1], v_ids[2], v_ids[3]}};
            if (CGAL::orientation(tet_vertices[t[0]].pos, tet_vertices[t[1]].pos, tet_vertices[t[2]].pos,
                                  tet_vertices[t[3]].pos) != CGAL::POSITIVE) {
                int tmp = t[1];
                t[1] = t[3];
                t[3] = tmp;
            }
            tets.push_back(t);
        } else {
            TetVertex v;
            tet_vertices.push_back(v);
            centroids_for_nodes[i] = tet_vertices.size() - 1;
        }
    }

    logger().debug("improvement {}", tmp_timer.getElapsedTime());
    tmp_timer.start();
    ///cal CDT & insert tets
    std::vector<std::vector<std::array<int, 3>>> cdt_faces(bsp_faces.size(), std::vector<std::array<int, 3>>());
    CDT cdt;
    for (unsigned int i = 0; i < bsp_faces.size(); i++) {
        if (bsp_faces[i].vertices.size() == 3) {
            cdt_faces[i].push_back(std::array<int, 3>({{bsp_faces[i].vertices[0], bsp_faces[i].vertices[1],
                                                        bsp_faces[i].vertices[2]}}));
            if (bsp_faces[i].conn_nodes.size() == 1) {
                for (int j = 0; j < bsp_faces[i].vertices.size(); j++)
                    tet_vertices[bsp_faces[i].vertices[j]].is_on_bbox = true;
            }
            continue;
        }

        cdt.clear();
        Plane_3 pln;
        constructPlane(i, pln);
        MC.initT(pln.orthogonal_vector());
        for (int j = 0; j < bsp_faces[i].edges.size(); j++) {
            cdt.insert_constraint(MC.to2d(bsp_vertices[bsp_edges[bsp_faces[i].edges[j]].vertices[0]]),
                                  MC.to2d(bsp_vertices[bsp_edges[bsp_faces[i].edges[j]].vertices[1]]));
        }
        if(cdt.number_of_vertices() != bsp_faces[i].vertices.size()){
            logger().debug("error: cdt.number_of_vertices() != bsp_faces[i].vertices.size()");
        }
        std::map<Point_2, int> vs_cdt2bsp;
        for (int j = 0; j < bsp_faces[i].vertices.size(); j++) {
            vs_cdt2bsp[MC.to2d(bsp_vertices[bsp_faces[i].vertices[j]])] = bsp_faces[i].vertices[j];//todo: improve to2d
        }
        for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            cdt_faces[i].push_back(std::array<int, 3>({{vs_cdt2bsp[fit->vertex(0)->point()],
                                                        vs_cdt2bsp[fit->vertex(1)->point()],
                                                        vs_cdt2bsp[fit->vertex(2)->point()]}}));
            if (bsp_faces[i].conn_nodes.size() == 1) {
                for (int j = 0; j < bsp_faces[i].vertices.size(); j++)
                    tet_vertices[bsp_faces[i].vertices[j]].is_on_bbox = true;
            }
        }
    }

    int rounded_cnt = 0;
    int all_cnt = 0;
    for(unsigned int i=0;i<bsp_nodes.size();i++) {
        if (is_tets[i])
            continue;

        all_cnt++;

        //calculate exact centroid
        int c_id = centroids_for_nodes[i];
        std::unordered_set<int> v_ids;
        for (int j = 0; j < bsp_nodes[i].faces.size(); j++) {
            for (int v_id:bsp_faces[bsp_nodes[i].faces[j]].vertices) {
                if (v_id < original_v_size)
                    v_ids.insert(v_id);
            }
        }

        std::vector<Point_3> vs;
        vs.reserve(v_ids.size());
        for (int v_id:v_ids)
            vs.push_back(bsp_vertices[v_id]);
        tet_vertices[c_id].pos = CGAL::centroid(vs.begin(), vs.end(), CGAL::Dimension_tag<0>());

        //insert new tets
        int t_cnt = 0;
        for (int j = 0; j < bsp_nodes[i].faces.size(); j++) {
            for (const std::array<int, 3> &f_ids:cdt_faces[bsp_nodes[i].faces[j]]) {
                std::array<int, 4> t = {{c_id, f_ids[0], f_ids[1], f_ids[2]}};
                if (CGAL::orientation(tet_vertices[t[0]].pos, tet_vertices[t[1]].pos, tet_vertices[t[2]].pos,
                                      tet_vertices[t[3]].pos) != CGAL::POSITIVE) {
                    int tmp = t[1];
                    t[1] = t[3];
                    t[3] = tmp;
                }
                tets.push_back(t);
                t_cnt++;
            }
        }

        //round into float
        Point_3 old_p = tet_vertices[c_id].pos;
        tet_vertices[c_id].posf = Point_3f(CGAL::to_double(old_p[0]), CGAL::to_double(old_p[1]),
                                           CGAL::to_double(old_p[2]));
        tet_vertices[c_id].pos = Point_3(tet_vertices[c_id].posf[0], tet_vertices[c_id].posf[1],
                                         tet_vertices[c_id].posf[2]);
        int tets_size = tets.size();
        bool is_rounded = true;
        for (int j = 0; j < t_cnt; j++) {
            if (CGAL::orientation(tet_vertices[tets[tets_size - 1 - j][0]].pos,
                                  tet_vertices[tets[tets_size - 1 - j][1]].pos,
                                  tet_vertices[tets[tets_size - 1 - j][2]].pos,
                                  tet_vertices[tets[tets_size - 1 - j][3]].pos) != CGAL::POSITIVE) {
                is_rounded = false;
                break;
            }
        }

        if (is_rounded) {
            tet_vertices[c_id].is_rounded = true;
            rounded_cnt++;
        } else {
            tet_vertices[c_id].pos = old_p;
            //todo: calculate a new position
        }
    }
    logger().debug("all_cnt = {}", all_cnt);
    logger().debug("rounded_cnt = {}", rounded_cnt);

    logger().debug("CDT {}", tmp_timer.getElapsedTime());
}

void SimpleTetrahedralization::labelSurface(const std::vector<int>& m_f_tags, const std::vector<int>& m_e_tags,
                                            const std::vector<std::vector<int>>& conn_e4v,
                                            std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets,
                                            std::vector<std::array<int, 4>>& is_surface_fs) {
    std::vector<BSPFace> &bsp_faces = MC.bsp_faces;
    std::vector<Point_3> &bsp_vertices = MC.bsp_vertices;
    const std::vector<Point_3> &m_vertices = MC.m_vertices;
    const std::vector<std::array<int, 3>> &m_faces = MC.m_faces;

    std::vector<Triangle_3> m_triangles;
    m_triangles.reserve(m_faces.size());
    for (unsigned int i = 0; i < m_faces.size(); i++)
        m_triangles.push_back(Triangle_3(m_vertices[m_faces[i][0]], m_vertices[m_faces[i][1]],
                                         m_vertices[m_faces[i][2]]));

    std::vector<std::vector<int>> track_on_faces;
    track_on_faces.resize(bsp_vertices.size());
    for (unsigned int i = 0; i < m_faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
//            tet_vertices[centroid_size+m_faces[i][j]].on_face.insert(i);
            tet_vertices[m_faces[i][j]].on_face.insert(i);
        }
    }

    std::vector<std::vector<int>> track_on_edges = conn_e4v;
    track_on_edges.resize(bsp_vertices.size());
    for (unsigned int i = 0; i < m_vertices.size(); i++) {
        for (int j = 0; j < track_on_edges[i].size(); j++) {
            if (m_e_tags[track_on_edges[i][j]] >= 0) {
//                tet_vertices[centroid_size+i].on_edge.insert(m_e_tags[track_on_edges[i][j]]);
                tet_vertices[i].on_edge.insert(m_e_tags[track_on_edges[i][j]]);
            }
        }
    }

    for (unsigned int i = 0; i < bsp_faces.size(); i++) {
        for (int j = 0; j < bsp_faces[i].vertices.size(); j++) {
            int v_id = bsp_faces[i].vertices[j];
            for (auto it = bsp_faces[i].div_faces.begin(); it != bsp_faces[i].div_faces.end(); it++)
                track_on_faces[v_id].push_back(*it);
            if (bsp_faces[i].matched_f_id >= 0)
                track_on_faces[v_id].push_back(bsp_faces[i].matched_f_id);
        }
    }
    for (unsigned int i = 0; i < track_on_faces.size(); i++) {
        std::sort(track_on_faces[i].begin(), track_on_faces[i].end());
        track_on_faces[i].erase(std::unique(track_on_faces[i].begin(), track_on_faces[i].end()),
                                track_on_faces[i].end());
    }

    for (unsigned int i = 0; i < bsp_vertices.size(); i++) {
        for (int j = 0; j < track_on_faces[i].size(); j++) {
            if (i < m_vertices.size() &&
//                std::find(tet_vertices[centroid_size+i].on_face.begin(), tet_vertices[centroid_size+i].on_face.end(),
//                          track_on_faces[i][j]) != tet_vertices[centroid_size+i].on_face.end()) {
                std::find(tet_vertices[i].on_face.begin(), tet_vertices[i].on_face.end(), track_on_faces[i][j]) != tet_vertices[i].on_face.end()) {
                    continue;
            }
            if (m_triangles[track_on_faces[i][j]].has_on(bsp_vertices[i])) {
                //if bsp_vertices[i] is not a vertex of the triangle
//                tet_vertices[centroid_size+i].on_face.insert(track_on_faces[i][j]);
                tet_vertices[i].on_face.insert(track_on_faces[i][j]);

                ///check on_edge
                std::array<int, 3> v_ids = {{m_faces[track_on_faces[i][j]][0], m_faces[track_on_faces[i][j]][1],
                                             m_faces[track_on_faces[i][j]][2]}};
                std::vector<int> e_ids;
                bool is_already_on_edge = false;
                for (int k = 0; k < 3; k++) {
                    std::set_intersection(conn_e4v[v_ids[k]].begin(), conn_e4v[v_ids[k]].end(),
                                          conn_e4v[v_ids[(k + 1) % 3]].begin(), conn_e4v[v_ids[(k + 1) % 3]].end(),
                                          std::back_inserter(e_ids));
//                    if (std::find(tet_vertices[centroid_size+i].on_edge.begin(), tet_vertices[centroid_size+i].on_edge.end(),
//                                  m_e_tags[e_ids.back()]) != tet_vertices[centroid_size+i].on_edge.end()) {
                    if (std::find(tet_vertices[i].on_edge.begin(), tet_vertices[i].on_edge.end(), m_e_tags[e_ids.back()]) != tet_vertices[i].on_edge.end()) {
                        is_already_on_edge = true;
                        break;
                    }
                }
                if (is_already_on_edge)
                    continue;

                assert(e_ids.size() == 3);
                for (int k = 0; k < 3; k++) {
                    if (m_e_tags[e_ids[k]] < 0)
                        continue;
                    Segment_3 seg(m_vertices[v_ids[k]], m_vertices[v_ids[(k + 1) % 3]]);
                    if (seg.has_on(bsp_vertices[i])) {
//                        tet_vertices[centroid_size+i].on_edge.insert(m_e_tags[e_ids[k]]);
                        tet_vertices[i].on_edge.insert(m_e_tags[e_ids[k]]);
                        break;
                    }
                }
            }
        }
    }

    for(unsigned int i=0;i<tet_vertices.size();i++){
        if(tet_vertices[i].on_face.size()>0)
            tet_vertices[i].is_on_surface=true;
    }

    ////is face on surface////
    // state.NOT_SURFACE = m_faces.size()+1;
    is_surface_fs=std::vector<std::array<int, 4>>(tets.size(),
                                                  std::array<int, 4>({{state.NOT_SURFACE, state.NOT_SURFACE, state.NOT_SURFACE, state.NOT_SURFACE}}));
//    std::vector<std::array<bool, 4>> is_visited(tets.size(), std::array<bool, 4>({{false, false, false, false}}));

    for(unsigned int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
//            if (is_visited[i][j])
//                continue;

            ///mark visited
//            int opp_i = getFaceOppoTets(tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4],
//                                        i, tet_vertices);
//            int opp_i = -1;
//
//            int opp_j = 0;
//            if (opp_i >= 0) {
//                for (int k = 0; k < 4; k++) {
//                    if (tets[opp_i][k] != tets[i][(j + 1) % 4] && tets[opp_i][k] != tets[i][(j + 2) % 4]
//                        && tets[opp_i][k] != tets[i][(j + 3) % 4]) {
//                        opp_j = k;
//                        break;
//                    }
//                }
//                is_visited[opp_i][opp_j] = true;
//            }
//            is_visited[i][j] = true;

            if (!tet_vertices[tets[i][(j + 1) % 4]].is_on_surface || !tet_vertices[tets[i][(j + 2) % 4]].is_on_surface
                || !tet_vertices[tets[i][(j + 3) % 4]].is_on_surface) {
                is_surface_fs[i][j] = state.NOT_SURFACE;
//                if (opp_i >= 0)
//                    is_visited[opp_i][opp_j] = state.NOT_SURFACE;
                continue;
            }
            std::unordered_set<int> sf_faces_tmp;
            setIntersection(tet_vertices[tets[i][(j + 1) % 4]].on_face, tet_vertices[tets[i][(j + 2) % 4]].on_face,
                            sf_faces_tmp);
            if (sf_faces_tmp.size() == 0) {
                is_surface_fs[i][j] = state.NOT_SURFACE;
//                if (opp_i >= 0)
//                    is_visited[opp_i][opp_j] = state.NOT_SURFACE;
                continue;
            }
            std::vector<int> sf_faces;
            setIntersection(sf_faces_tmp, tet_vertices[tets[i][(j + 3) % 4]].on_face, sf_faces);
            if (sf_faces.size() == 0) {
                is_surface_fs[i][j] = state.NOT_SURFACE;
//                if (opp_i >= 0)
//                    is_visited[opp_i][opp_j] = state.NOT_SURFACE;
                continue;
            }

//            if (tmp.size() > 1) {
//                std::array<int, 3> f = {{tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4]}};
//                std::sort(f.begin(), f.end());
//                folding_fs.push_back(std::array<int, 4>({{f[0], f[1], f[2], i}}));
//                continue;
//            }

            ////get the first ori
            is_surface_fs[i][j] = 0;
            Plane_3 pln(m_vertices[m_faces[sf_faces[0]][0]], m_vertices[m_faces[sf_faces[0]][1]],
                        m_vertices[m_faces[sf_faces[0]][2]]);
            CGAL::Oriented_side side = pln.oriented_side(tet_vertices[tets[i][j]].pos);

            if (side == CGAL::ON_ORIENTED_BOUNDARY) {
                log_and_throw("ERROR: side == CGAL::ON_ORIENTED_BOUNDARY!!");
            }
            if (side == CGAL::ON_POSITIVE_SIDE)//outside
                is_surface_fs[i][j]++;
            else//inside
                is_surface_fs[i][j]--;

            ////if there are more than one sf_faces
            int delta = is_surface_fs[i][j];
            if (sf_faces.size() > 1) {
                //cal normal vec for [0]
//                Vector_3 nv = CGAL::cross_product(
//                        m_vertices[m_faces[sf_faces[0]][0]] - m_vertices[m_faces[sf_faces[0]][1]],
//                        m_vertices[m_faces[sf_faces[0]][1]] - m_vertices[m_faces[sf_faces[0]][2]]);
//                Direction_3 dir = nv.direction();
                Direction_3 dir = pln.orthogonal_direction();

                for (int f_id = 1; f_id < sf_faces.size(); f_id++) {
                    Vector_3 nv1 = CGAL::cross_product(
                            m_vertices[m_faces[sf_faces[f_id]][0]] - m_vertices[m_faces[sf_faces[f_id]][1]],
                            m_vertices[m_faces[sf_faces[f_id]][1]] - m_vertices[m_faces[sf_faces[f_id]][2]]);
                    Direction_3 dir1 = nv1.direction();
                    if (dir1 == dir)
                        is_surface_fs[i][j] += delta;
                    else
                        is_surface_fs[i][j] -= delta;
//                    else {
//                        logger().debug("wrong direction!!");
//                        pausee();
//                    }
                }
            }

//            for (int tri_id : sf_faces) {
//                Plane_3 pln(m_vertices[m_faces[tri_id][0]], m_vertices[m_faces[tri_id][1]],
//                            m_vertices[m_faces[tri_id][2]]);
//                CGAL::Oriented_side side = pln.oriented_side(tet_vertices[tets[i][j]].pos);
//
//                if (side == CGAL::ON_ORIENTED_BOUNDARY) {
//                    log_and_throw("ERROR: side == CGAL::ON_ORIENTED_BOUNDARY!!");
//                }
//                if (side == CGAL::ON_POSITIVE_SIDE)//outside
//                    is_surface_fs[i][j]++;
//                else//inside
//                    is_surface_fs[i][j]--;
//            }

//            if (opp_i >= 0)
//                is_surface_fs[opp_i][opp_j] = -is_surface_fs[i][j];
        }
    }

    //tag the surface
    for(unsigned int i=0;i<tet_vertices.size();i++){
        std::unordered_set<int> tmp;
        for(auto it=tet_vertices[i].on_face.begin();it!=tet_vertices[i].on_face.end();it++)
            tmp.insert(m_f_tags[*it]);
        tet_vertices[i].on_face=tmp;
    }
}

void SimpleTetrahedralization::labelBbox(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets){
    std::vector<Point_3> &bsp_vertices = MC.bsp_vertices;

    //label bbox
    std::vector<int> bbx_f_tags;
    std::vector<int> bbx_e_tags;
    for (int i = 0; i < 6; i++)
        bbx_f_tags.push_back(-(i + 1));
    for (int i = 0; i < 12; i++)
        bbx_e_tags.push_back(-(i + 1));
    std::vector<int> bbx_faces={{0, 1,
                                 2, 3,
                                 4, 5}};
    std::vector<int> bbx_edges={{0, 1, 2, 3,
                                 4, 5, 6, 7,
                                 8, 9, 10, 11}};

    int i=0;
    int i0=0, i7=0;
    for (int I = 0; I < tet_vertices.size(); I++) {
        if (!tet_vertices[I].is_on_bbox)
            continue;
        if (i < 8) {
            tet_vertices[I].on_fixed_vertex = -2 - i;

            std::bitset<3> a(i);
            for (int j = 0; j < 3; j++) {
                tet_vertices[I].on_face.insert(bbx_f_tags[bbx_faces[j * 2 + a[j]]]);
                tet_vertices[I].on_edge.insert(bbx_e_tags[bbx_edges[j * 4 + a[j] * 2 + a[(j + 1) % 3]]]);
            }
        } else {
            //和bbx_v_ids[0], bbx_v_ids[7]比较
            std::array<int, 3> a = {{-1, -1, -1}};
            for (int j = 0; j < 3; j++) {
//                if (bsp_vertices[I - centroid_size][j] == bsp_vertices[i0 - centroid_size][j])
                if (bsp_vertices[I][j] == bsp_vertices[i0][j])
                    a[j] = 0;
//                else if (bsp_vertices[I - centroid_size][j] == bsp_vertices[i7 - centroid_size][j])
                else if (bsp_vertices[I][j] == bsp_vertices[i7][j])
                    a[j] = 1;
            }
            for (int j = 0; j < 3; j++) {
                if (a[j] < 0)
                    continue;
                tet_vertices[I].on_face.insert(bbx_f_tags[bbx_faces[j * 2 + a[j]]]);
                if (a[(j + 1) % 3] >= 0) {
                    tet_vertices[I].on_edge.insert(bbx_e_tags[bbx_edges[j * 4 + a[j] * 2 + a[(j + 1) % 3]]]);
                }
            }
        }
        if (i == 0)
            i0 = I;
        else if (i == 7)
            i7 = I;
        i++;
    }
    logger().debug("#v on bbox = {}", i);
}

void SimpleTetrahedralization::labelBoundary(std::vector<TetVertex>& tet_vertices, std::vector<std::array<int, 4>>& tets,
                                             const std::vector<std::array<int, 4>>& is_surface_fs) {

    std::vector<std::array<int, 2>> edges_tmp;
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 1; j < 4; j++) {
            if (tet_vertices[tets[i][0]].is_on_surface && tet_vertices[tets[i][j]].is_on_surface) {
                std::array<int, 2> e = {{tets[i][0], tets[i][j]}};
                if (e[1] < e[0])
                    e = {{e[1], e[0]}};
                edges_tmp.push_back(e);
            }
            if (tet_vertices[tets[i][j]].is_on_surface && tet_vertices[tets[i][j % 3 + 1]].is_on_surface) {
                std::array<int, 2> e = {{tets[i][j], tets[i][j % 3 + 1]}};
                if (e[1] < e[0])
                    e = {{e[1], e[0]}};
                edges_tmp.push_back(e);
            }
        }
    }
    std::sort(edges_tmp.begin(), edges_tmp.end());
    edges_tmp.erase(std::unique(edges_tmp.begin(), edges_tmp.end()), edges_tmp.end());

    for (int i = 0; i < edges_tmp.size(); i++) {
        int cnt = 0;
        for (int t_id: tet_vertices[edges_tmp[i][0]].conn_tets) {
            std::vector<int> opp_js;
            for (int j = 0; j < 4; j++) {
                if (tets[t_id][j] == edges_tmp[i][0] || tets[t_id][j] == edges_tmp[i][1])
                    continue;
                opp_js.push_back(j);
            }
            if (opp_js.size() == 2) {
                if (is_surface_fs[t_id][opp_js[0]] != state.NOT_SURFACE)
                    cnt++;
                if (is_surface_fs[t_id][opp_js[1]] != state.NOT_SURFACE)
                    cnt++;
                if (cnt > 2)
                    break;
            }
        }
        if (cnt == 2) {//is boundary edge
            tet_vertices[edges_tmp[i][0]].is_on_boundary = true;
            tet_vertices[edges_tmp[i][1]].is_on_boundary = true;
        }
    }

    int cnt_boundary = 0, cnt_surface = 0;
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (tet_vertices[i].is_on_boundary)
            cnt_boundary++;
        if (tet_vertices[i].is_on_surface)
            cnt_surface++;
    }
    logger().debug("{} vertices on boundary", cnt_boundary);
    logger().debug("{} vertices on surface", cnt_surface);
}

void SimpleTetrahedralization::constructPlane(int bsp_f_id, Plane_3& pln) {
    pln = Plane_3(MC.bsp_vertices[MC.bsp_faces[bsp_f_id].vertices[0]],
                  MC.bsp_vertices[MC.bsp_faces[bsp_f_id].vertices[1]],
                  MC.bsp_vertices[MC.bsp_faces[bsp_f_id].vertices[2]]);
    int i = 3;
    while (pln.is_degenerate() && i < MC.bsp_faces[bsp_f_id].vertices.size()) {
        pln = Plane_3(MC.bsp_vertices[MC.bsp_faces[bsp_f_id].vertices[0]],
                      MC.bsp_vertices[MC.bsp_faces[bsp_f_id].vertices[1]],
                      MC.bsp_vertices[MC.bsp_faces[bsp_f_id].vertices[i++]]);
    }
    assert(!(pln.is_degenerate()));
}

} // namespace tetwild
