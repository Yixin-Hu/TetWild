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

#include <tetwild/MeshConformer.h>
#include <tetwild/Logger.h>

namespace tetwild {

void MeshConformer::match() {
    is_matched.assign(m_faces.size(), false);
    matchDivFaces();
    return;
}

void MeshConformer::matchVertexIndices(int x, const std::vector<std::array<int, 2>>& seed_v_list,
                                       std::vector<int>& f_list) {
    int seed_id;
    int l = 0, r = seed_v_list.size() - 1;
    while (l < r) {
        int mid = (l + r) / 2;
        if (seed_v_list[mid][0] < x)
            l = mid + 1;
        else
            r = mid;
    }
    if (r >= l && seed_v_list[l][0]==x) {
        f_list.push_back(seed_v_list[l][1]);
        int s = l - 1;
        while (s >= 0) {
            if (seed_v_list[s][0] == x)
                f_list.push_back(seed_v_list[s][1]);
            else
                break;
            s--;
        }
        s = l + 1;
        while (s < seed_v_list.size()) {
            if (seed_v_list[s][0] == x)
                f_list.push_back(seed_v_list[s][1]);
            else
                break;
            s++;
        }
    }
}


void MeshConformer::matchDivFaces() {
    std::vector<std::array<int, 2>> seed_v_list;
    seed_v_list.reserve(bsp_faces.size() * 3);
    for (int i = 0; i < bsp_faces.size(); i++) {
        for (int j = 0; j < 3; j++)
            seed_v_list.push_back(std::array<int, 2>({{bsp_faces[i].vertices[j], i}}));
    }
    std::sort(seed_v_list.begin(), seed_v_list.end(), [](const std::array<int, 2> &a, const std::array<int, 2> &b) {
        return a[0] < b[0];//todo: there should be a O(n) algorithm
    });

    const int m_faces_size = m_faces.size();
    for (int i = 0; i < m_faces_size; i++) {
        int m_f_id = i;
        std::array<Point_3, 3> tri1 = {{m_vertices[m_faces[i][0]], m_vertices[m_faces[i][1]], m_vertices[m_faces[i][2]]}};

        ///find seed info
        std::unordered_set<int> seed_fids, seed_nids;
        std::vector<std::vector<int>> f_lists;
        is_matched[i] = true;
        for (int j = 0; j < 3; j++) {
            std::vector<int> f_list;
            matchVertexIndices(m_faces[i][j], seed_v_list, f_list);
            if (f_list.size() == 0)
                is_matched[i] = false;
            for (int k = 0; k < f_list.size(); k++)
                seed_fids.insert(f_list[k]);
            if (is_matched[i])//possibly matched
                f_lists.push_back(f_list);
        }
        if (is_matched[i]) {
            is_matched[i] = false;
            for (int j = 0; j < f_lists.size(); j++)
                std::sort(f_lists[j].begin(), f_lists[j].end());
            std::vector<int> tmp, tmp1;
            std::set_intersection(f_lists[0].begin(), f_lists[0].end(), f_lists[1].begin(), f_lists[1].end(),
                                  std::back_inserter(tmp));
            std::set_intersection(tmp.begin(), tmp.end(), f_lists[2].begin(), f_lists[2].end(),
                                  std::back_inserter(tmp1));//todo: find a better way to cal intersection
            //todo: you can use unordered_set and set a *good* hash function to find identical triangles on bsp tree
            if (tmp1.size() == 1) {
                std::array<int, 3> tmp_bsp = {{bsp_faces[tmp1[0]].vertices[0], bsp_faces[tmp1[0]].vertices[1],
                                               bsp_faces[tmp1[0]].vertices[2]}};
                std::array<int, 3> tmp_m = m_faces[i];
                std::sort(tmp_bsp.begin(), tmp_bsp.end());
                std::sort(tmp_m.begin(), tmp_m.end());
                if (tmp_bsp == tmp_m) {
                    is_matched[i] = true;
                    bsp_faces[tmp1[0]].matched_f_id = i;
                    continue;
                }
            }
        }
        for (auto it = seed_fids.begin(); it != seed_fids.end(); it++)
            for (auto jt = bsp_faces[*it].conn_nodes.begin(); jt != bsp_faces[*it].conn_nodes.end(); jt++)
                seed_nids.insert(*jt);

        for (auto it = seed_fids.begin(); it != seed_fids.end(); it++) {
            ////cal intersection type
            std::array<Point_3, 3> tri2 = {{bsp_vertices[bsp_faces[*it].vertices[0]],
                                            bsp_vertices[bsp_faces[*it].vertices[1]],
                                            bsp_vertices[bsp_faces[*it].vertices[2]]}};
            int int_type = triangleIntersection3d(tri1, tri2, true);
            if (int_type == COPLANAR_INT) {
                bsp_faces[*it].div_faces.insert(i);
            } else if (int_type == CROSS_INT) {
                for (auto nit = bsp_faces[*it].conn_nodes.begin(); nit != bsp_faces[*it].conn_nodes.end(); nit++)
                    bsp_nodes[*nit].div_faces.insert(i);
            }
        }

        ///dfs all the info
        std::unordered_set<int> new_fids;
        std::unordered_set<int> new_nids = seed_nids;
        while (true) {
            new_fids.clear();
            for (auto it = new_nids.begin(); it != new_nids.end(); it++) {
                for (int j = 0; j < bsp_nodes[*it].faces.size(); j++) {
                    int bsp_f_id = bsp_nodes[*it].faces[j];
                    auto fit = std::find(seed_fids.begin(), seed_fids.end(), bsp_f_id);
                    if (fit == seed_fids.end()) {
                        ///check if the plane-coplanar or plane-crossing
                        ///check if intersecting (coplanar -> do_intersection / crossing -> sort 4 interseting points)
                        ///if intersected -> insert into new_fids
                        /////if coplanar-intersecting -> divface for face
                        /////else if crossing-intersecting -> divface for node
                        std::array<Point_3, 3> tri2 = {{bsp_vertices[bsp_faces[bsp_f_id].vertices[0]],
                                                        bsp_vertices[bsp_faces[bsp_f_id].vertices[1]],
                                                        bsp_vertices[bsp_faces[bsp_f_id].vertices[2]]}};
                        int int_type = triangleIntersection3d(tri1, tri2, false);
                        if (int_type != NONE_INT)
                            new_fids.insert(bsp_f_id);

                        if (int_type == COPLANAR_INT) {
                            bsp_faces[bsp_f_id].div_faces.insert(i);
                        } else if (int_type == CROSS_INT) {
                            for (auto nit = bsp_faces[bsp_f_id].conn_nodes.begin();
                                 nit != bsp_faces[bsp_f_id].conn_nodes.end(); nit++)
                                bsp_nodes[*nit].div_faces.insert(i);
                        }
                    }
                }
            }

            if (new_fids.size() == 0)
                break;

            new_nids.clear();
            for (auto it = new_fids.begin(); it != new_fids.end(); it++) {
                for (auto jt = bsp_faces[*it].conn_nodes.begin(); jt != bsp_faces[*it].conn_nodes.end(); jt++) {
                    auto nid = std::find(seed_nids.begin(), seed_nids.end(), *jt);
                    if (nid == seed_nids.end()) {
                        new_nids.insert(*jt);
                    }
                }
            }
            if (new_nids.size() == 0)
                break;

            seed_fids.insert(new_fids.begin(), new_fids.end());
            seed_nids.insert(new_nids.begin(), new_nids.end());//c++11
        }
    }
    logger().debug("{} faces matched!", std::count(is_matched.begin(), is_matched.end(), true));
}

void MeshConformer::getOrientedVertices(int bsp_f_id){
    if(bsp_faces[bsp_f_id].vertices.size()==3)
        return;

    std::vector<int> vertices;
    int begin=bsp_edges[bsp_faces[bsp_f_id].edges[0]].vertices[0];
    int end=bsp_edges[bsp_faces[bsp_f_id].edges[0]].vertices[1];
    vertices.push_back(begin);

    std::vector<bool> is_visited(bsp_faces[bsp_f_id].edges.size(), false);
    is_visited[0]=true;
    for(int i=0;i<bsp_faces[bsp_f_id].vertices.size()-1;i++) {
        for (int j = 0; j < bsp_faces[bsp_f_id].edges.size(); j++) {
            if (is_visited[j])
                continue;
            if (bsp_edges[bsp_faces[bsp_f_id].edges[j]].vertices[0] == end) {
                vertices.push_back(end);
                end = bsp_edges[bsp_faces[bsp_f_id].edges[j]].vertices[1];
                is_visited[j]=true;
            } else if (bsp_edges[bsp_faces[bsp_f_id].edges[j]].vertices[1] == end) {
                vertices.push_back(end);
                end = bsp_edges[bsp_faces[bsp_f_id].edges[j]].vertices[0];
                is_visited[j]=true;
            }
        }
    }

    bsp_faces[bsp_f_id].vertices=vertices;

    if(vertices.size()!=bsp_faces[bsp_f_id].vertices.size()){
        logger().error("{}, {}", bsp_faces[bsp_f_id].vertices.size(), bsp_faces[bsp_f_id].edges.size());
        throw TetWildError("MeshConformer::getOrientedVertices");
    }
}

void triangleSideofPlane(const std::array<Point_3, 3>& tri, const Plane_3& pln,
                         std::vector<Point_3>& pos_vs, std::vector<Point_3>& neg_vs, std::vector<Point_3>& on_vs) {
    for (int i = 0; i < 3; i++) {
        CGAL::Oriented_side side = pln.oriented_side(tri[i]);
        if (side == CGAL::ON_POSITIVE_SIDE) {
            pos_vs.push_back(tri[i]);
        } else if (side == CGAL::ON_NEGATIVE_SIDE) {
            neg_vs.push_back(tri[i]);
        } else {
            on_vs.push_back(tri[i]);
        }
    }
}

int MeshConformer::triangleIntersection3d(const std::array<Point_3, 3>& tri1, const std::array<Point_3, 3>& tri2,
                                          bool intersect_known) {
    if (!intersect_known) {
        Triangle_3 t1(tri1[0], tri1[1], tri1[2]);
        Triangle_3 t2(tri2[0], tri2[1], tri2[2]);
        if (!do_intersect(t1, t2))
            return NONE_INT;
    }

    Plane_3 pln1(tri1[0], tri1[1], tri1[2]);
    Plane_3 pln2(tri2[0], tri2[1], tri2[2]);

    std::vector<Point_3> pos_vs1, neg_vs1, on_vs1;
    triangleSideofPlane(tri1, pln2, pos_vs1, neg_vs1, on_vs1);

    ///coplanar
    if (on_vs1.size() == 3) {
//        if (!intersect_known) {
//            Triangle_3 t1(tri1[0], tri1[1], tri1[2]);
//            Triangle_3 t2(tri2[0], tri2[1], tri2[2]);
//            if (do_intersect(t1, t2))
//                return COPLANAR_INT;
//            else
//                return NONE_INT;
//        } else
            return COPLANAR_INT;
    }

    ///on one side
    if (pos_vs1.size() == 0 || neg_vs1.size() == 0) {
//        if(intersect_known)
            return POINT_INT;
//
//        Triangle_3 t1(tri1[0], tri1[1], tri1[2]);
//        Triangle_3 t2(tri2[0], tri2[1], tri2[2]);
//        if (do_intersect(t1, t2))
//            return POINT_INT;
//        else
//            return NONE_INT;
    }

    ///cross
    std::vector<Point_3> pos_vs2, neg_vs2, on_vs2;
    triangleSideofPlane(tri2, pln1, pos_vs2, neg_vs2, on_vs2);
    if (pos_vs2.size() == 0 || neg_vs2.size() == 0) {
//        if(intersect_known)
            return POINT_INT;
//
//        Triangle_3 t1(tri1[0], tri1[1], tri1[2]);
//        Triangle_3 t2(tri2[0], tri2[1], tri2[2]);
//        if (do_intersect(t1, t2))
//            return POINT_INT;
//        else
//            return NONE_INT;
    }

    std::vector<std::pair<Point_3, int>> sorted_vs;
    ///cal intersecting points for tri1
    for (int i = 0; i < pos_vs1.size(); i++) {
        for (int j = 0; j < neg_vs1.size(); j++) {
            Segment_3 seg1(pos_vs1[i], neg_vs1[j]);
            auto result = intersection(seg1, pln2);
            assert(!(!result));
            if (result) {
                if (const Point_3 *p = boost::get<Point_3>(&*result))
                    on_vs1.push_back(*p);
                else
                    throw TetWildError("MeshConformer::triangleIntersection3d");
            }
        }
    }
    assert(!(on_vs1.size() != 2));
    for (int i = 0; i < 2; i++)
        sorted_vs.push_back(std::make_pair(on_vs1[i], 1));

    ///cal intersecting points for tri2
    for (int i = 0; i < pos_vs2.size(); i++) {
        for (int j = 0; j < neg_vs2.size(); j++) {
            Segment_3 seg2(pos_vs2[i], neg_vs2[j]);
            auto result = intersection(seg2, pln1);
            assert(!(!result));
            if (result) {
                if (const Point_3 *p = boost::get<Point_3>(&*result))
                    on_vs2.push_back(*p);
                else
                    throw TetWildError("MeshConformer::triangleIntersection3d");
            }
        }
    }
    assert(!(on_vs2.size() != 2));
    for (int i = 0; i < 2; i++)
        sorted_vs.push_back(std::make_pair(on_vs2[i], 2));
    std::sort(sorted_vs.begin(), sorted_vs.end(), [](const std::pair<Point_3, int> &p1,
                                                     const std::pair<Point_3, int> &p2) {
        return p1.first < p2.first;
    });

//    if (!intersect_known) {
//        if (sorted_vs[1].first != sorted_vs[2].first){
//            if(sorted_vs[0].second != sorted_vs[1].second)
//                return CROSS_INT;
//            else
//                return NONE_INT;
//        }
//        else
//            return POINT_INT;
//    } else {
        if (sorted_vs[1].first != sorted_vs[2].first)
            return CROSS_INT;
        else
            return POINT_INT;
//    }
}

void MeshConformer::initT(const Vector_3& nv) {
    std::vector<Vector_3> vs {Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1)};
    std::vector<bool> is_ppd(3, false);

    for (int i = 0; i < 3; i++) {
        if (nv * vs[i] == 0)
            is_ppd[i] = true;
    }
    int i = std::find(is_ppd.begin(), is_ppd.end(), false) - is_ppd.begin();
    t = i;
}

Point_2 MeshConformer::to2d(const Point_3& p){
    int x=(t+1)%3;
    int y=(t+2)%3;
    return Point_2(p[x], p[y]);
}

Point_3 MeshConformer::to3d(const Point_2& p, const Plane_3& pln) {
    Line_3 l;
    switch (t) {
        case 0:
            l = Line_3(Point_3(0, p[0], p[1]), Direction_3(1, 0, 0));
            break;
        case 1:
            l = Line_3(Point_3(p[1], 0, p[0]), Direction_3(0, 1, 0));
            break;
        case 2:
            l = Line_3(Point_3(p[0], p[1], 0), Direction_3(0, 0, 1));
    }

    auto result = intersection(l, pln);
    if (result) {
        const Point_3 *p = boost::get<Point_3>(&*result);
        return *p;
    } else {
        log_and_throw("error to3d!");
    }
}

} // namespace tetwild
