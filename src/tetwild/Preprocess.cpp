// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 10/12/17.
//

#include <tetwild/Preprocess.h>
#include <tetwild/Args.h>
#include <tetwild/Common.h>
#include <tetwild/State.h>
#include <tetwild/Logger.h>
#include <tetwild/DistanceQuery.h>
#include <pymesh/MshSaver.h>
#include <igl/fit_plane.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/write_triangle_mesh.h>
#include <igl/unique.h>
#include <igl/unique_simplices.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/writeSTL.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/geometry_nd.h>
#include <unordered_map>

namespace tetwild {

namespace {

void checkBoundary(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const State &state) {
    PyMesh::MshSaver mSaver(state.working_dir+state.postfix+"_boundary.msh", true);
    Eigen::VectorXd oV;
    Eigen::VectorXi oF;
    oV.resize(V.rows() * 3);
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < 3; j++)
            oV(i * 3 + j) = V(i, j);
    }
    oF.resize(F.rows() * 3);
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++)
            oF(i * 3 + j) = F(i, j);
    }

    if (oV.rows() > 0) {
        mSaver.save_mesh(oV, oF, 3, mSaver.TRI);

        Eigen::VectorXd bF(F.rows());
        for (int i = 0; i < bF.rows(); i++)
            bF(i) = 0;
        Eigen::VectorXd bV(V.rows());
        for (int i = 0; i < bV.rows(); i++)
            bV(i) = 0;

        std::vector<std::vector<int>> conn_f4v(V.rows(), std::vector<int>());
        for (int i = 0; i < F.rows(); i++) {
            for (int j = 0; j < 3; j++)
                conn_f4v[F(i, j)].push_back(i);
        }

        for (int i = 0; i < F.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                std::vector<int> tmp;
                std::set_intersection(conn_f4v[F(i, j)].begin(), conn_f4v[F(i, j)].end(),
                                      conn_f4v[F(i, (j + 1) % 3)].begin(), conn_f4v[F(i, (j + 1) % 3)].end(),
                                      std::back_inserter(tmp));
                if (tmp.size() == 1) {
                    bF(tmp[0]) = 1;
//                    logger().debug("boundary tri! {}", tmp[0]);
//                    Triangle_3f tri(Point_3f(V(F(i, 0), 0), V(F(i, 0), 1), V(F(i, 0), 2)),
//                                    Point_3f(V(F(i, 1), 0), V(F(i, 1), 1), V(F(i, 1), 2)),
//                                    Point_3f(V(F(i, 2), 0), V(F(i, 2), 1), V(F(i, 2), 2)));
//                    if(tri.is_degenerate())
//                        logger().debug("degenerate");
                    bV(F(i, j)) = 1;
                    bV(F(i, (j + 1) % 3)) = 1;
                }
            }
        }

        mSaver.save_elem_scalar_field("boundary faces", bF);
        mSaver.save_scalar_field("boundary vertices", bV);
    }

    logger().debug("boundary checked!");
}

} // anonymous namespace

bool Preprocess::init(const Eigen::MatrixXd& V_tmp, const Eigen::MatrixXi& F_tmp,
                      GEO::Mesh& geo_b_mesh, GEO::Mesh& geo_sf_mesh, const Args &args) {

    logger().debug("{} {}", V_tmp.rows(), F_tmp.rows());

    Eigen::VectorXi IV, _;
//        igl::unique_rows(V_tmp, V_in, _, IV);
    igl::remove_duplicate_vertices(V_tmp, F_tmp, 1e-10, V_in, IV, _, F_in);

    if (V_in.rows() == 0 || F_in.rows() == 0)
        return false;

//        for (int i = 0; i < F_in.rows(); i++) {
//            for (int j = 0; j < 3; j++) {
//                F_in(i, j) = IV(F_in(i, j));
//            }
//        }
    logger().debug("#v = {} -> {}", V_tmp.rows(), V_in.rows());
    logger().debug("#f = {} -> {}", F_tmp.rows(), F_in.rows());
//    checkBoundary(V_in, F_in, state);

    ////get GEO meshes
    geo_sf_mesh.vertices.clear();
    geo_sf_mesh.vertices.create_vertices((int) V_in.rows());
    for (int i = 0; i < V_in.rows(); i++) {
        GEO::vec3 &p = geo_sf_mesh.vertices.point(i);
        for (int j = 0; j < 3; j++)
            p[j] = V_in(i, j);
    }
    geo_sf_mesh.facets.clear();
    geo_sf_mesh.facets.create_triangles((int) F_in.rows());
    for (int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            geo_sf_mesh.facets.set_vertex(i, j, F_in(i, j));
    }
    geo_sf_mesh.facets.compute_borders();

    getBoundaryMesh(geo_b_mesh);
    state.is_mesh_closed = (geo_b_mesh.vertices.nb() == 0);

    return true;
}

void Preprocess::getBoundaryMesh(GEO::Mesh& b_mesh) {
    Eigen::MatrixXd& V_sf = V_in;
    Eigen::MatrixXi& F_sf = F_in;

    std::vector<std::vector<int>> conn_f4v(V_sf.rows(), std::vector<int>());
    for (int i = 0; i < F_sf.rows(); i++) {
        for (int j = 0; j < 3; j++)
            conn_f4v[F_sf(i, j)].push_back(i);
    }
    //check isolated vertices
//    for(int i=0;i<conn_f4v.size();i++){
//        if(conn_f4v[i].size()==0)
//            logger().debug("iso");
//    }

    std::vector<std::array<int, 2>> b_edges;
    for (int i = 0; i < F_sf.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            std::vector<int> tmp;
            std::set_intersection(conn_f4v[F_sf(i, j)].begin(), conn_f4v[F_sf(i, j)].end(),
                                  conn_f4v[F_sf(i, (j + 1) % 3)].begin(), conn_f4v[F_sf(i, (j + 1) % 3)].end(),
                                  std::back_inserter(tmp));
            if (tmp.size() == 1)
                b_edges.push_back(std::array<int, 2>({{F_sf(i, j), F_sf(i, (j + 1) % 3)}}));
        }
    }

    if(b_edges.size()==0){
        b_mesh.vertices.clear();
        return;
    }

    std::unordered_set<int> v_ids;
    for (int i = 0; i < b_edges.size(); i++) {
        v_ids.insert(b_edges[i][0]);
        v_ids.insert(b_edges[i][1]);
    }
    std::unordered_map<int, int> v_ids_map;
    int cnt = 0;
    for (int v_id : v_ids) {
        v_ids_map[v_id] = cnt;
        cnt++;
    }

    b_mesh.vertices.clear();
    b_mesh.vertices.create_vertices((int) v_ids.size());
    for (int v_id : v_ids) {
        GEO::vec3 &p = b_mesh.vertices.point(v_ids_map[v_id]);
        for (int j = 0; j < 3; j++)
            p[j] = V_sf(v_id, j);
    }
    b_mesh.facets.clear();
    b_mesh.facets.create_triangles((int) b_edges.size());
    for (int i = 0; i < b_edges.size(); i++) {
        b_mesh.facets.set_vertex(i, 0, v_ids_map[b_edges[i][0]]);
        b_mesh.facets.set_vertex(i, 1, v_ids_map[b_edges[i][1]]);
        b_mesh.facets.set_vertex(i, 2, v_ids_map[b_edges[i][1]]);
    }
}

void Preprocess::process(GEO::Mesh& geo_sf_mesh, std::vector<Point_3>& m_vertices, std::vector<std::array<int, 3>>& m_faces, const Args &args) {
    double eps_scalar = 0.8;
    double eps_scalar_2 = eps_scalar*eps_scalar;

   state.eps *= eps_scalar;
    state.eps_2 *= eps_scalar_2;
//    state.sampling_dist *= eps_scalar*2;

    conn_fs.resize(V_in.size());
    for (int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            conn_fs[F_in(i, j)].insert(i);
    }
    v_is_removed = std::vector<bool>(V_in.rows(), false);
    f_is_removed = std::vector<bool>(F_in.rows(), false);

    // mesh_reorder(geo_sf_mesh, GEO::MESH_ORDER_HILBERT);
    GEO::MeshFacetsAABBWithEps geo_face_tree(geo_sf_mesh);

    std::vector<std::array<int, 2>> edges;
    edges.reserve(F_in.rows()*6);
    for (int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{F_in(i, j), F_in(i, (j + 1) % 3)}};
            if (e[0] > e[1]) e = {{e[1], e[0]}};
            edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    const int edges_size = edges.size();
    for (int i = 0; i < edges_size; i++) {
        double weight = getEdgeLength(edges[i]);
        sm_queue.push(ElementInQueue_sm(edges[i], weight));
        sm_queue.push(ElementInQueue_sm(std::array<int, 2>({{edges[i][1], edges[i][0]}}), weight));
    }

    //simplification
    ts = 0;
    f_tss.resize(F_in.size());
    simplify(geo_sf_mesh, geo_face_tree);

    ////get CGAL surface mesh
    int cnt = 0;
    std::unordered_map<int, int> new_v_ids;

    Eigen::MatrixXd V_out;
    Eigen::MatrixXi F_out;
    V_out.resize(std::count(v_is_removed.begin(), v_is_removed.end(), false), 3);
    F_out.resize(std::count(f_is_removed.begin(), f_is_removed.end(), false), 3);
    for (int i = 0; i < V_in.rows(); i++) {
        if (v_is_removed[i])
            continue;
        new_v_ids[i] = cnt;
        V_out.row(cnt++) = V_in.row(i);
    }

    cnt = 0;
    for (int i = 0; i < F_in.rows(); i++) {
        if (f_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++)
            F_out(cnt, j) = new_v_ids[F_in(i, j)];
        cnt++;
    }
//    igl::writeSTL(state.working_dir+args.postfix+"_simplified.stl", V_out, F_out);
    logger().debug("#v = {}", V_out.rows());
    logger().debug("#f = {}", F_out.rows());

    V_in = V_out;
    F_in = F_out;
    conn_fs.clear();
    conn_fs.resize(V_in.size());
    for (int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            conn_fs[F_in(i, j)].insert(i);
    }
    swap(geo_sf_mesh, geo_face_tree);
    if(args.save_mid_result == 0)
        igl::writeSTL(state.working_dir+state.postfix+"_simplified.stl", V_in, F_in);

//    checkBoundary(V_in, F_in);

    m_vertices.reserve(V_in.rows());
    m_faces.reserve(F_in.rows());
    for (int i = 0; i < V_in.rows(); i++) {
        m_vertices.push_back(Point_3(V_in(i, 0), V_in(i, 1), V_in(i, 2)));
    }
    for (int i = 0; i < F_in.rows(); i++) {
        std::array<int, 3> f = {{F_in(i, 0), F_in(i, 1), F_in(i, 2)}};
        Triangle_3 tr(m_vertices[f[0]], m_vertices[f[1]], m_vertices[f[2]]);
        if (!tr.is_degenerate())//delete all degenerate triangles
            m_faces.push_back(f);
    }
    logger().debug("#v = {}", m_vertices.size());
    logger().debug("#f = {}->{}", F_in.rows(), m_faces.size());

    state.eps /= eps_scalar;
    state.eps_2 /= eps_scalar_2;
//    state.sampling_dist /= eps_scalar*2;

    // igl::write_triangle_mesh("tmp.obj", V_in, F_in);

    // igl::write_triangle_mesh("tmp.obj", V_in, F_in);

    //output colormap
    //    outputSurfaceColormap(geo_face_tree, geo_sf_mesh);
}

void Preprocess::swap(const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree) {
    int cnt = 0;
    for (int i = 0; i < F_in.rows(); i++) {
        bool is_swapped = false;
        for (int j = 0; j < 3; j++) {
            int v_id = F_in(i, j);
            int v1_id = F_in(i, (j + 1) % 3);
            int v2_id = F_in(i, (j + 2) % 3);

            // manifold
            std::vector<int> n12_f_ids;
            setIntersection(conn_fs[v1_id], conn_fs[v2_id], n12_f_ids);
            if (n12_f_ids.size() != 2) {
                continue;
            }
            if (n12_f_ids[1] == i)
                n12_f_ids = {{n12_f_ids[1], n12_f_ids[0]}};
            int v3_id = -1;
            for (int k = 0; k < 3; k++)
                if (F_in(n12_f_ids[1], k) != v1_id && F_in(n12_f_ids[1], k) != v2_id) {
                    v3_id = F_in(n12_f_ids[1], k);
                    break;
                }

            // check quality
            double cos_a = getCosAngle(v_id, v1_id, v2_id);
            double cos_a1 = getCosAngle(v3_id, v1_id, v2_id);
            std::array<GEO::vec3, 2> old_nvs;
            for (int f = 0; f < 2; f++) {
                std::array<GEO::vec3, 3> vs;
                for (int k = 0; k < 3; k++) {
                    vs[k] = GEO::vec3(V_in(F_in(n12_f_ids[f], k), 0), V_in(F_in(n12_f_ids[f], k), 1),
                                      V_in(F_in(n12_f_ids[f], k), 2));
                }
                old_nvs[f] = GEO::Geom::triangle_normal(vs[0], vs[1], vs[2]);
            }
            if (cos_a > -0.999) {
//                continue;
                if (GEO::Geom::cos_angle(old_nvs[0], old_nvs[1]) < 1-1e-6)//not coplanar
                    continue;
            }
            double cos_a_new = getCosAngle(v1_id, v_id, v3_id);
            double cos_a1_new = getCosAngle(v2_id, v_id, v3_id);
            if (std::min(cos_a_new, cos_a1_new) <= std::min(cos_a, cos_a1))
                continue;

            // non flipping
            Eigen::RowVectorXi f1_old = F_in.row(n12_f_ids[0]);
            Eigen::RowVectorXi f2_old = F_in.row(n12_f_ids[1]);
            for (int k = 0; k < 3; k++) {
                if (F_in(n12_f_ids[0], k) == v2_id)
                    F_in(n12_f_ids[0], k) = v3_id;
                if (F_in(n12_f_ids[1], k) == v1_id)
                    F_in(n12_f_ids[1], k) = v_id;
            }
            GEO::vec3 old_nv = cos_a1 < cos_a ? old_nvs[0] : old_nvs[1];
            bool is_filp = false;
            for (int f_id:n12_f_ids) {
                std::array<GEO::vec3, 3> vs;
                for (int k = 0; k < 3; k++) {
                    vs[k] = GEO::vec3(V_in(F_in(f_id, k), 0), V_in(F_in(f_id, k), 1), V_in(F_in(f_id, k), 2));
                }
                GEO::vec3 new_nv = GEO::Geom::triangle_normal(vs[0], vs[1], vs[2]);
                if (GEO::dot(old_nv, new_nv) < 0) {
                    is_filp = true;
                    break;
                }
            }
            if (is_filp) {
                F_in.row(n12_f_ids[0]) = f1_old;
                F_in.row(n12_f_ids[1]) = f2_old;
                continue;
            }

            // non outside envelop
            std::unordered_set<int> new_f_ids;
            new_f_ids.insert(n12_f_ids.begin(), n12_f_ids.end());
            if (isOutEnvelop(new_f_ids, geo_mesh, face_aabb_tree)) {
                F_in.row(n12_f_ids[0]) = f1_old;
                F_in.row(n12_f_ids[1]) = f2_old;
                continue;
            }

            // real update
            conn_fs[v1_id].erase(n12_f_ids[1]);
            conn_fs[v2_id].erase(n12_f_ids[0]);
            conn_fs[v_id].insert(n12_f_ids[1]);
            conn_fs[v3_id].insert(n12_f_ids[0]);
            is_swapped = true;
            break;
        }
        if (is_swapped)
            cnt++;
    }
    logger().debug("{} faces are swapped!!", cnt);
}

double Preprocess::getCosAngle(int v_id, int v1_id, int v2_id) {
    return GEO::Geom::cos_angle(GEO::vec3(V_in(v1_id, 0), V_in(v1_id, 1), V_in(v1_id, 2)) -
                                GEO::vec3(V_in(v_id, 0), V_in(v_id, 1), V_in(v_id, 2)),
                                GEO::vec3(V_in(v2_id, 0), V_in(v2_id, 1), V_in(v2_id, 2)) -
                                GEO::vec3(V_in(v_id, 0), V_in(v_id, 1), V_in(v_id, 2)));
}

void Preprocess::simplify(const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree) {
    int cnt = 0;
//    logger().debug("queue.size() = {}", sm_queue.size());
    while (!sm_queue.empty()) {
        std::array<int, 2> v_ids = sm_queue.top().v_ids;
        double old_weight = sm_queue.top().weight;
        sm_queue.pop();

        if (!isEdgeValid(v_ids, old_weight))
            continue;

        if (!removeAnEdge(v_ids[0], v_ids[1], geo_mesh, face_aabb_tree)) {
            inf_es.push_back(v_ids);
            inf_e_tss.push_back(ts);
        } else {
            cnt++;
            if (cnt % 1000 == 0)
                logger().debug("1000 vertices removed");
        }
    }
    logger().debug("{}", cnt);
    logger().debug("{}", c);

    if (cnt > 0)
        postProcess(geo_mesh, face_aabb_tree);
}

void Preprocess::postProcess(const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree){
    logger().debug("postProcess!");

    std::vector<std::array<int, 2>> tmp_inf_es;
    const int inf_es_size = inf_es.size();
    for (int i = 0; i < inf_es_size; i++) {
        if (!isEdgeValid(inf_es[i]))
            continue;
        bool is_recal = false;
        for (int f_id:conn_fs[inf_es[i][0]]) {
            if (f_tss[f_id] > inf_e_tss[i]) {
                is_recal = true;
                break;
            }
        }
        if (is_recal)
            sm_queue.push(ElementInQueue_sm(inf_es[i], getEdgeLength(inf_es[i])));
        else
            tmp_inf_es.push_back(inf_es[i]);
    }
    std::sort(tmp_inf_es.begin(), tmp_inf_es.end());
    tmp_inf_es.erase(std::unique(tmp_inf_es.begin(), tmp_inf_es.end()), tmp_inf_es.end());
    inf_es = tmp_inf_es;
    ts++;
    inf_e_tss = std::vector<int>(inf_es.size(), ts);

    simplify(geo_mesh, face_aabb_tree);
}

bool Preprocess::removeAnEdge(int v1_id, int v2_id, const GEO::Mesh &geo_mesh, const GEO::MeshFacetsAABBWithEps& face_aabb_tree) {
    if (!isOneRingClean(v1_id) || !isOneRingClean(v2_id))
        return false;

    //check if flip after collapsing
    std::vector<int> n12_f_ids;
    setIntersection(conn_fs[v1_id], conn_fs[v2_id], n12_f_ids);
    if (n12_f_ids.size() != 2) {//!!!
//        logger().debug("error: n12_f_ids.size()!=2");
        return false;
    }

    std::unordered_set<int> new_f_ids;
    for (int f_id:conn_fs[v1_id]) {
        if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1]) {
            new_f_ids.insert(f_id);
        }
    }
    for (int f_id:conn_fs[v2_id]) {
        if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1]) {
            new_f_ids.insert(f_id);
        }
    }

    //check euclidean characteristics (delete degenerate and duplicate elements
    if(!isEuclideanValid(v1_id, v2_id))
        return false;

    //computing normal for checking flipping
    for (int f_id:new_f_ids) {
//    for(int f_id:conn_fs[v1_id]) {
//        if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1]) {
        std::array<GEO::vec3, 3> vs;
        for (int j = 0; j < 3; j++) {
            vs[j] = GEO::vec3(V_in(F_in(f_id, j), 0), V_in(F_in(f_id, j), 1), V_in(F_in(f_id, j), 2));
        }
        GEO::vec3 old_nv = GEO::Geom::triangle_normal(vs[0], vs[1], vs[2]);

        for (int j = 0; j < 3; j++) {
//                if (F_in(f_id, j) == v1_id) {
            if (F_in(f_id, j) == v1_id || F_in(f_id, j) == v2_id) {
                vs[j] = (GEO::vec3(V_in(v1_id, 0), V_in(v1_id, 1), V_in(v1_id, 2)) +
                         GEO::vec3(V_in(v2_id, 0), V_in(v2_id, 1), V_in(v2_id, 2))) / 2;
                break;
            }
        }
        GEO::vec3 new_nv = GEO::Geom::triangle_normal(vs[0], vs[1], vs[2]);
        if (GEO::dot(old_nv, new_nv) < 0)
            return false;
//        }
    }

    //check if go outside of envelop
    Eigen::VectorXd v1_old_p = V_in.row(v1_id);
    Eigen::VectorXd v2_old_p = V_in.row(v2_id);
    V_in.row(v1_id) = (V_in.row(v1_id) + V_in.row(v2_id)) / 2;
    GEO::vec3 mid_p(V_in(v1_id, 0), V_in(v1_id, 1), V_in(v1_id, 2));
    GEO::vec3 nearest_p;
    double _;
    face_aabb_tree.nearest_facet(mid_p, nearest_p, _);//project back to surface
    for(int j=0;j<3;j++)
        V_in(v1_id, j) = nearest_p[j];
    V_in.row(v2_id) = V_in.row(v1_id);
    if (isOutEnvelop(new_f_ids, geo_mesh, face_aabb_tree)) {
        V_in.row(v1_id) = v1_old_p;
        V_in.row(v2_id) = v2_old_p;
        return false;
    }
    c++;

    //real update
    std::unordered_set<int> n_v_ids;//get this info before real update for later usage
    for (int f_id:new_f_ids) {
        for (int j = 0; j < 3; j++) {
            if (F_in(f_id, j) != v1_id && F_in(f_id, j) != v2_id)
                n_v_ids.insert(F_in(f_id, j));
        }
    }

    v_is_removed[v1_id] = true;
    for (int f_id:n12_f_ids) {
        f_is_removed[f_id] = true;
        for (int j = 0; j < 3; j++) {//rm conn_fs
            if (F_in(f_id, j) != v1_id) {
                conn_fs[F_in(f_id, j)].erase(f_id);
//                auto it = std::find(conn_fs[F_in(f_id, j)].begin(), conn_fs[F_in(f_id, j)].end(), f_id);
//                if (it != conn_fs[F_in(f_id, j)].end())
//                    conn_fs[F_in(f_id, j)].erase(it);
            }
        }
    }
    for (int f_id:conn_fs[v1_id]) {//add conn_fs
        if (f_is_removed[f_id])
            continue;
        conn_fs[v2_id].insert(f_id);
        for (int j = 0; j < 3; j++) {
            if (F_in(f_id, j) == v1_id)
                F_in(f_id, j) = v2_id;
        }
    }

    //update timestamps
    ts++;
    for (int f_id:conn_fs[v2_id]) {
        f_tss[f_id] = ts;
    }

    //push new edges into the queue
    for (int v_id:n_v_ids) {
        double weight = getEdgeLength(v2_id, v_id);
        sm_queue.push(ElementInQueue_sm(std::array<int, 2>({{v2_id, v_id}}), weight));
        sm_queue.push(ElementInQueue_sm(std::array<int, 2>({{v_id, v2_id}}), weight));
    }

    return true;
}

bool Preprocess::isEdgeValid(const std::array<int, 2>& v_ids){
    if(v_is_removed[v_ids[0]] || v_is_removed[v_ids[1]])
        return false;

    if(!isHaveCommonEle(conn_fs[v_ids[0]], conn_fs[v_ids[1]]))
        return false;

    return true;
}

bool Preprocess::isEdgeValid(const std::array<int, 2>& v_ids, double old_weight){
    if(v_is_removed[v_ids[0]] || v_is_removed[v_ids[1]])
        return false;

    if(!isHaveCommonEle(conn_fs[v_ids[0]], conn_fs[v_ids[1]]))
        return false;

    if(old_weight!=getEdgeLength(v_ids))//the edge is outdated
        return false;

    return true;
}

double Preprocess::getEdgeLength(const std::array<int, 2>& v_ids){
    return (V_in.row(v_ids[0]) - V_in.row(v_ids[1])).squaredNorm();
}

double Preprocess::getEdgeLength(int v1_id, int v2_id){
    return (V_in.row(v1_id) - V_in.row(v2_id)).squaredNorm();
}

bool Preprocess::isOneRingClean(int v1_id){
    std::vector<std::array<int, 2>> n1_es;
    for (int f_id:conn_fs[v1_id]) {
        for (int j = 0; j < 3; j++) {
            if (F_in(f_id, j) != v1_id) {
                if (F_in(f_id, (j + 1) % 3) == v1_id)
                    n1_es.push_back(std::array<int, 2>({{F_in(f_id, j), 0}}));
                else if (F_in(f_id, (j + 2) % 3) == v1_id)
                    n1_es.push_back(std::array<int, 2>({{F_in(f_id, j), 1}}));
            }
        }
    }
    if (n1_es.size() % 2 != 0)
        return true;
    std::sort(n1_es.begin(), n1_es.end());
    for (int i = 0; i < n1_es.size(); i += 2) {
        if (n1_es[i][0] == n1_es[i + 1][0] && n1_es[i][1] != n1_es[i + 1][1]);
        else
            return false;
    }

    return true;
}


bool Preprocess::isOutEnvelop(const std::unordered_set<int>& new_f_ids,
    const GEO::Mesh &geo_sf_mesh, const GEO::MeshFacetsAABBWithEps& geo_face_tree)
{
    size_t num_querried = 0;
    size_t num_tris = new_f_ids.size();
    size_t num_samples = 0;
    size_t tri_idx = 0;

    static thread_local std::vector<GEO::vec3> ps;
    for (int f_id:new_f_ids) {
        //sample triangles except one-ring of v1v2
        std::array<GEO::vec3, 3> vs = {{
                GEO::vec3(V_in(F_in(f_id, 0), 0), V_in(F_in(f_id, 0), 1), V_in(F_in(f_id, 0), 2)),
                GEO::vec3(V_in(F_in(f_id, 1), 0), V_in(F_in(f_id, 1), 1), V_in(F_in(f_id, 1), 2)),
                GEO::vec3(V_in(F_in(f_id, 2), 0), V_in(F_in(f_id, 2), 1), V_in(F_in(f_id, 2), 2))}};
        ps.clear();
        sampleTriangle(vs, ps, state.sampling_dist);
        ++tri_idx;
        num_samples += ps.size();

//        logger().debug("ps.size = {}", ps.size());
//        logger().debug("is output samples?");
//        int anw = 0;
//        cin >> anw;
//        if (anw != 0) {
////        if (true) {
//            Eigen::MatrixXd V_tmp(ps.size() * 3 + 3, 3);
//            Eigen::MatrixXi F_tmp(ps.size() + 1, 3);
//            for (int i = 0; i < 3; i++) {
//                for (int j = 0; j < 3; j++)
//                    V_tmp(i, j) = vs[i][j];
//                F_tmp(0, i) = i;
//            }
//
//            for (int i = 0; i < ps.size(); i++) {
//                for (int k = 0; k < 3; k++) {
//                    for (int j = 0; j < 3; j++)
//                        V_tmp((1 + i) * 3 + k, j) = ps[i][j];
//                    F_tmp(1 + i, k) = (1 + i) * 3 + k;
//                }
//            }
//            igl::writeSTL(state.working_dir + "_sample.stl", V_tmp, F_tmp);
//        }


//        std::array<double, 3> ls;
//        for (int i = 0; i < 3; i++) {
//            ls[i] = GEO::length(vs[i] - vs[(i + 1) % 3]);
//        }
//        auto min_max = std::minmax_element(ls.begin(), ls.end());
//        int min_i = min_max.first - ls.begin();
//        int max_i = min_max.second - ls.begin();
//
//        double n = ls[max_i] / state.sampling_dist;
//        if (n <= 1) {
//            for (int i = 0; i < 3; i++)
//                ps.push_back(vs[i]);
//        } else {
//            n = int(n) + 1;
//            ps.reserve(n + n + 1);
//            for (int j = 0; j <= n; j++) {
//                ps.push_back(j / n * vs[(min_i + 2) % 3] + (n - j) / n * vs[min_i]);
//                if (j == n)
//                    break;
//                ps.push_back(j / n * vs[(min_i + 2) % 3] + (n - j) / n * vs[(min_i + 1) % 3]);
//            }
//            if (ls[min_i] > state.sampling_dist) {
//                const int ps_size = ps.size();
//                for (int i = 0; i < ps_size - 1; i += 2) {
//                    double m = GEO::length(ps[i] - ps[i + 1]) / state.sampling_dist;
//                    if (m < 1)
//                        break;
//                    m = int(m) + 1;
//                    for (int j = 1; j < m; j++)
//                        ps.push_back(j / m * ps[i] + (m - j) / m * ps[i + 1]);
//                }
//            }
//        }

        //check sampling points
        GEO::vec3 nearest_point;
        double sq_dist = std::numeric_limits<double>::max();
        GEO::index_t prev_facet = GEO::NO_FACET;

        for (const GEO::vec3 &current_point:ps) {
            if (prev_facet != GEO::NO_FACET) {
                get_point_facet_nearest_point(geo_sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
            }
            if (sq_dist > state.eps_2) {
                geo_face_tree.facet_in_envelope_with_hint(
                    current_point, state.eps_2, prev_facet, nearest_point, sq_dist);
            }
            ++num_querried;
            if (sq_dist > state.eps_2) {
                logger().trace("num_triangles {} / {} num_queries {} / {}",
                    tri_idx - 1, num_tris, num_querried, num_samples);
                return true;
            }
        }
    }
    logger().trace("num_triangles {} / {} num_queries {} / {}",
        tri_idx - 1, num_tris, num_querried, num_samples);

    return false;
}

bool Preprocess::isPointOutEnvelop(int v_id, const GEO::MeshFacetsAABBWithEps& geo_face_tree){
    if (geo_face_tree.squared_distance(GEO::vec3(V_in(v_id, 0), V_in(v_id, 1), V_in(v_id, 2))) > state.eps_2)
        return true;
    return false;
}

int calEuclidean(const std::vector<std::array<int, 3>>& fs){
    std::vector<std::array<int, 2>> es;
    es.reserve(fs.size()*3);
    std::unordered_set<int> vs;
    for(int i=0;i<fs.size();i++){
        for(int j=0;j<3;j++){
            vs.insert(fs[i][j]);
            std::array<int, 2> e={{fs[i][j], fs[i][(j+1)%3]}};
            if(e[0]>e[1])
                e={{e[1], e[0]}};
            es.push_back(e);
        }
    }
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());

    return vs.size()-es.size()+fs.size();
}

bool Preprocess::isEuclideanValid(int v1_id, int v2_id){
//    logger().debug("v1:{}", v1_id);
//    for (int f_id:conn_fs[v1_id]) {
//        logger().debug("{}{}{} {}", F_in(f_id, 0), ' ', F_in(f_id, 1), F_in(f_id, 2));
//    }
//    logger().debug("v2:{}", v2_id);
//    for (int f_id:conn_fs[v2_id]) {
//        logger().debug("{}{}{} {}", F_in(f_id, 0), ' ', F_in(f_id, 1), F_in(f_id, 2));
//    }

    std::vector<std::array<int, 3>> fs;
    for(int I=0;I<2;I++) {
        int v_id = I == 0 ? v1_id : v2_id;
        for (int f_id:conn_fs[v_id]) {
            if (F_in(f_id, 0) != F_in(f_id, 1) && F_in(f_id, 1) != F_in(f_id, 2) && F_in(f_id, 0) != F_in(f_id, 2)) {
                std::array<int, 3> f = {{F_in(f_id, 0), F_in(f_id, 1), F_in(f_id, 2)}};
                std::sort(f.begin(), f.end());
                fs.push_back(f);
            }
        }
    }
    std::sort(fs.begin(), fs.end());
    fs.erase(std::unique(fs.begin(), fs.end()), fs.end());
//    logger().debug("fs.size() = {}", fs.size());
    int ec0=calEuclidean(fs);
//    logger().debug("{}", ec0);

    std::vector<std::array<int, 3>> fs1;
    for(int i=0;i<fs.size();i++){
        for(int j=0;j<3;j++){
            if(fs[i][j]==v1_id) {
                fs[i][j] = v2_id;
                break;
            }
        }
//        logger().debug("{} {} {}", fs[i][0], fs[i][1], fs[i][2]);
        if(fs[i][0]!=fs[i][1]&&fs[i][1]!=fs[i][2]&&fs[i][0]!=fs[i][2]){
            std::array<int, 3> f = {{fs[i][0], fs[i][1], fs[i][2]}};
            std::sort(f.begin(), f.end());
            fs1.push_back(f);
        }
    }
    std::sort(fs1.begin(), fs1.end());
    fs1.erase(std::unique(fs1.begin(), fs1.end()), fs1.end());
//    logger().debug("fs1.size() = {}", fs1.size());
    int ec1=calEuclidean(fs1);
//    logger().debug("{}", ec1);

//    pausee();

    if(ec0!=ec1)
        return false;
    return true;
}

void Preprocess::outputSurfaceColormap(const GEO::MeshFacetsAABBWithEps& geo_face_tree, const GEO::Mesh& geo_sf_mesh) {
    Eigen::VectorXd eps_dis(F_in.rows());
    for(int f_id=0;f_id<geo_sf_mesh.facets.nb();f_id++) {
//        if(f_id!=1871)
//            continue;

//    for (int f_id = 0; f_id < F_in.rows(); f_id++) {
        //sample triangles except one-ring of v1v2
        std::vector<GEO::vec3> ps;
        std::array<GEO::vec3, 3> vs = {{
        geo_sf_mesh.vertices.point(geo_sf_mesh.facets.vertex(f_id, 0)),
        geo_sf_mesh.vertices.point(geo_sf_mesh.facets.vertex(f_id, 1)),
        geo_sf_mesh.vertices.point(geo_sf_mesh.facets.vertex(f_id, 2))}};
//                GEO::vec3(V_in(F_in(f_id, 0), 0), V_in(F_in(f_id, 0), 1), V_in(F_in(f_id, 0), 2)),
//                GEO::vec3(V_in(F_in(f_id, 1), 0), V_in(F_in(f_id, 1), 1), V_in(F_in(f_id, 1), 2)),
//                GEO::vec3(V_in(F_in(f_id, 2), 0), V_in(F_in(f_id, 2), 1), V_in(F_in(f_id, 2), 2))};
        std::array<double, 3> ls;
        for (int i = 0; i < 1; i++) {
            ls[i] = GEO::length(vs[i] - vs[(i + 1) % 3]);
            //
            double n = int(ls[i] / state.sampling_dist + 1);
            for (int j = 1; j < n; j++) {
                ps.push_back(double(j) / n * vs[i] + (n - double(j)) / n * vs[(i + 1) % 3]);
            }
            //
        }
//        auto min_max = std::minmax_element(ls.begin(), ls.end());
//        int min_i = min_max.first - ls.begin();
//        int max_i = min_max.second - ls.begin();
//
//        double n = int(ls[max_i] / state.sampling_dist + 1);
//        ps.reserve(2*n);
//        for (int j = 0; j <= n; j++) {
//            ps.push_back(j / n * vs[(min_i + 2) % 3] + (n - j) / n * vs[min_i]);
//            ps.push_back(j / n * vs[(min_i + 2) % 3] + (n - j) / n * vs[(min_i + 1) % 3]);
//        }

//        if(ls[min_i] > state.sampling_dist) {
//            int ps_size = ps.size();
//            for (int i = 0; i < ps_size; i += 2) {
//                double m = int(GEO::length(ps[i] - ps[i + 1]) / state.sampling_dist + 1);
//                if(m==0)
//                    break;
//                for (int j = 1; j < m; j++)
//                    ps.push_back(j / m * ps[i] + (m - j) / m * ps[i + 1]);
//            }
//        }
//        ps.push_back(vs[(min_i + 2) % 3]);
        for(int i=0;i<3;i++) {
            ps.push_back(vs[i]);
        }
//        ps.push_back((vs[0]+vs[1]+vs[2])/3);

        //check sampling points
        GEO::vec3 current_point = ps[0];
        GEO::vec3 nearest_point;
        double sq_dist;
        GEO::index_t prev_facet = geo_face_tree.nearest_facet(current_point, nearest_point, sq_dist);

        double max_dis = 0;
        int cnt=0;
        GEO::vec3 pp;
        std::vector<int> fs;
        for (const GEO::vec3 &current_point:ps) {
            double dis;
            int n_facet = geo_face_tree.nearest_facet(current_point, nearest_point, dis);
//            sq_dist = current_point.distance2(nearest_point);
//            geo_face_tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_dist);
//            double dis = current_point.distance2(nearest_point);
//            if(f_id==2514)
//                logger().debug("{}: {} {} {}", cnt, dis, sq_dist, int(prev_facet));
            if (dis > max_dis) {
                max_dis = dis;
                pp=current_point;
            }
            cnt++;
            fs.push_back(int(n_facet));
        }
        cnt = 0;
        if(f_id==1681) {
            for (const GEO::vec3 &p:ps) {
                logger().debug("{}: {}, {}, {}; {}; {}", cnt, p[0], p[1], p[2], fs[cnt], geo_face_tree.squared_distance(p));
                cnt++;
            }
        }

        eps_dis(f_id) = sqrt(max_dis / state.eps_2);
        if(eps_dis(f_id)>1) {
            logger().debug("ERROR: simplified input goes outside of the envelop");
            logger().debug("{}", f_id);
            logger().debug("{}", eps_dis(f_id));
            logger().debug("{} {}", max_dis, state.eps_2);
            cnt = 0;
            for (const GEO::vec3 &p:ps) {
                logger().debug("{}: {}, {}, {}; {}; {}", cnt, p[0], p[1], p[2], fs[cnt], geo_face_tree.squared_distance(p));
                cnt++;
            }
//            logger().debug("{}", geo_face_tree.squared_distance(pp));
//            double dd;
//            logger().debug("{}", int(geo_face_tree.nearest_facet(pp, nearest_point, dd)));
//            logger().debug("{}", dd);

            std::vector<int> vf={{1681, 1675, 1671, 1666}};
            for(int j=0;j<vf.size();j++) {
                logger().debug("f {}: {}", vf[j], GEO::Geom::triangle_area(
                        geo_sf_mesh.vertices.point(geo_sf_mesh.facets.vertex(vf[j], 0)),
                        geo_sf_mesh.vertices.point(geo_sf_mesh.facets.vertex(vf[j], 1)),
                        geo_sf_mesh.vertices.point(geo_sf_mesh.facets.vertex(vf[j], 2))));

//                logger().debug("{} {}{}{}", geo_sf_mesh.facets.vertex(vf[j], 0), geo_sf_mesh.facets.vertex(vf[j], 1), " "
//, geo_sf_mesh.facets.vertex(vf[j], 2));
                int v1_id = geo_sf_mesh.facets.vertex(vf[j], 0);
                int v2_id = geo_sf_mesh.facets.vertex(vf[j], 1);
                int v3_id = geo_sf_mesh.facets.vertex(vf[j], 2);
                std::array<int, 3> v_ids = {{v1_id ,v2_id, v3_id}};
                for (int k = 0; k < 3; k++) {
                    logger().debug("{}: {} {} {}", v_ids[k], geo_sf_mesh.vertices.point(v_ids[k])[0], geo_sf_mesh.vertices.point(v_ids[k])[1], geo_sf_mesh.vertices.point(v_ids[k])[2]);
                }
            }


            double min_dis=0;
            int vf_id = 0;
            std::vector<double> diss;
            GEO::vec3 nearest_p;
            double _1, _2, _3;
            for(int i=0;i<geo_sf_mesh.facets.nb();i++) {
                double s_dis = GEO::Geom::point_triangle_squared_distance(ps[10],
                                                                          geo_sf_mesh.vertices.point(
                                                                                  geo_sf_mesh.facets.vertex(i, 0)),
                                                                          geo_sf_mesh.vertices.point(
                                                                                  geo_sf_mesh.facets.vertex(i, 1)),
                                                                          geo_sf_mesh.vertices.point(
                                                                                  geo_sf_mesh.facets.vertex(i, 2)),
                                                                          nearest_p, _1, _2, _3);

                diss.push_back(s_dis);
                if (i == 0 || s_dis < min_dis) {
                    min_dis = s_dis;
                    vf_id = i;
                }
            }
            for(int i=0;i<geo_sf_mesh.facets.nb();i++) {
                if(diss[i]==min_dis)
                    logger().debug("vf_id = {}", i);
            }
            logger().debug("something");
            logger().debug("double check 1681 {}", GEO::Geom::point_triangle_squared_distance(ps[10],
                                                                      geo_sf_mesh.vertices.point(
                                                                              geo_sf_mesh.facets.vertex(1681, 0)),
                                                                      geo_sf_mesh.vertices.point(
                                                                              geo_sf_mesh.facets.vertex(1681, 1)),
                                                                      geo_sf_mesh.vertices.point(
                                                                              geo_sf_mesh.facets.vertex(1681, 2)),
                                                                                   nearest_p, _1, _2, _3));

            logger().debug("min_dis = {}", min_dis);
            logger().debug("diag = {}", state.bbox_diag);
        }


    }


    Eigen::VectorXd V_vec(V_in.rows() * 3);
    Eigen::VectorXi F_vec(F_in.rows() * 3);
    for (int i = 0; i < V_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            V_vec(i * 3 + j) = V_in(i, j);
    }
    for (int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            F_vec(i * 3 + j) = F_in(i, j);
    }

    PyMesh::MshSaver mshSaver(state.working_dir + state.postfix + "_sf.msh");
    mshSaver.save_mesh(V_vec, F_vec, 3, mshSaver.TRI);
    mshSaver.save_elem_scalar_field("distance to surface", eps_dis);
}

} // namespace tetwild
