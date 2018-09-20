// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 4/3/17.
//

#include <tetwild/BSPSubdivision.h>
#include <tetwild/MeshConformer.h>
#include <tetwild/Logger.h>

namespace tetwild {

void BSPSubdivision::init() {
    for (int old_n_id = 0; old_n_id < MC.bsp_nodes.size(); old_n_id++) {
        if (MC.bsp_nodes[old_n_id].div_faces.size() == 0) {
            MC.bsp_nodes[old_n_id].is_leaf = true;
            continue;
        }
        processing_n_ids.push(old_n_id);
    }
	float nf=0;
	for (int old_f_id = 0; old_f_id < MC.bsp_faces.size(); old_f_id++) {
		if (MC.bsp_faces[old_f_id].div_faces.size() == 0) {
			continue;
		}
		nf+=1;
	}
	logger().debug("# nodes need subdivision = {}/{}/{}", nf, processing_n_ids.size(), MC.bsp_nodes.size());
}

void BSPSubdivision::subdivideBSPNodes() {
    std::vector<BSPtreeNode> &nodes = MC.bsp_nodes;
    std::vector<BSPFace> &faces = MC.bsp_faces;
    std::vector<BSPEdge> &edges = MC.bsp_edges;
    std::vector<Point_3> &vertices = MC.bsp_vertices;
    const std::vector<Point_3> &div_vertices = MC.m_vertices;
    const std::vector<std::array<int, 3>> &div_faces = MC.m_faces;

    while (!processing_n_ids.empty()) {
        int old_n_id = processing_n_ids.front();
        processing_n_ids.pop();

        ///re-assign divfaces
        int cnt_pos=0, cnt_neg=0;
        Plane_3 pln;
        std::unordered_set<int> v_ids;
        std::unordered_map<int, int> v_sides;
        std::vector<int> rm_df_ids;
        int on_df_id;
        bool is_divide=false;
        for(auto it=nodes[old_n_id].div_faces.begin(); it!=nodes[old_n_id].div_faces.end();it++) {
            pln = Plane_3(div_vertices[div_faces[*it][0]],
                          div_vertices[div_faces[*it][1]],
                          div_vertices[div_faces[*it][2]]);
            ///map sides for vertices
            v_ids.clear();
            for (int i = 0; i < nodes[old_n_id].faces.size(); i++) {
                for (int j = 0; j < faces[nodes[old_n_id].faces[i]].vertices.size(); j++) {
                    v_ids.insert(faces[nodes[old_n_id].faces[i]].vertices[j]);
                }
            }
            v_sides.clear();
            calVertexSides(pln, v_ids, vertices, v_sides);
            cnt_pos=0;
            cnt_neg=0;
            for (auto it = v_sides.begin(); it != v_sides.end(); it++) {
                if (it->second == V_POS)
                    cnt_pos++;
                if (it->second == V_NEG)
                    cnt_neg++;
            }
            if (cnt_pos == 0 || cnt_neg == 0) { //fixed//but how could it happen??
                nodes[old_n_id].is_leaf = true;
                rm_df_ids.push_back(*it);
            }
            else{
                is_divide=true;
                on_df_id=*it;
                break;
            }
        }
        if(!is_divide)
            continue;
        ///from here, the node would definitely be subdivided
        BSPtreeNode pos_node, neg_node;
        BSPFace on_face;
        nodes.push_back(neg_node);
        int new_n_id = nodes.size() - 1;
        std::vector<int> new_n_ids = {old_n_id, new_n_id};

        for (auto it = nodes[old_n_id].div_faces.begin(); it != nodes[old_n_id].div_faces.end(); it++) {
            if(std::find(rm_df_ids.begin(), rm_df_ids.end(), *it)!=rm_df_ids.end())
                continue;
            if (*it == on_df_id) {
                on_face.div_faces.insert(*it);
                continue;
            }

            int side = divfaceSide(pln, div_faces[*it], div_vertices);
            if (side == DIVFACE_POS)
                pos_node.div_faces.insert(*it);
            else if (side == DIVFACE_NEG)
                neg_node.div_faces.insert(*it);
            else if (side == DIVFACE_ON)
                on_face.div_faces.insert(*it);
            else if (side == DIVFACE_CROSS) {
                pos_node.div_faces.insert(*it);
                neg_node.div_faces.insert(*it);
            }
        }

        ///split nodes
        for (int i = 0; i < nodes[old_n_id].faces.size(); i++) {
            int old_f_id = nodes[old_n_id].faces[i];

            ///check if need splitting
            int cnt_pos = 0, cnt_neg = 0, cnt_on = 0;
            for (int j = 0; j < faces[old_f_id].vertices.size(); j++) {
                if (v_sides[faces[old_f_id].vertices[j]] == V_POS)
                    cnt_pos++;
                else if (v_sides[faces[old_f_id].vertices[j]] == V_NEG)
                    cnt_neg++;
                else
                    cnt_on++;
            }
            if (cnt_pos + cnt_on == faces[old_f_id].vertices.size()) {
                pos_node.faces.push_back(old_f_id);
                if (cnt_on == 0) {
                    continue;
                }
            }
            if (cnt_neg + cnt_on == faces[old_f_id].vertices.size()) {
                neg_node.faces.push_back(old_f_id);
                if (cnt_on == 0) {
                    continue;
                }
            }

            ///splitting...
            BSPFace pos_face, neg_face;
            BSPEdge on_edge;
            int new_f_id = faces.size();
            std::vector<int> new_f_ids = {old_f_id, new_f_id};

            bool is_connected=false;
            for (int j = 0; j < faces[old_f_id].edges.size(); j++) {
                int old_e_id = faces[old_f_id].edges[j];

                ///check if need splitting
                std::vector<int> pos_vs, neg_vs, on_vs;
                for (int j = 0; j < edges[old_e_id].vertices.size(); j++) {
                    if (v_sides[edges[old_e_id].vertices[j]] == V_POS)
                        pos_vs.push_back(edges[old_e_id].vertices[j]);
                    else if (v_sides[edges[old_e_id].vertices[j]] == V_NEG)
                        neg_vs.push_back(edges[old_e_id].vertices[j]);
                    else
                        on_vs.push_back(edges[old_e_id].vertices[j]);
                }
                if (on_vs.size() == 2) {
                    if (std::find(on_face.edges.begin(), on_face.edges.end(), old_e_id) == on_face.edges.end())
                        on_face.edges.push_back(old_e_id);
                    is_connected=true;
                }
                if(is_connected)
                    continue;
                if (pos_vs.size() + on_vs.size() == 2) {
                    pos_face.edges.push_back(old_e_id);
                    if (on_vs.size() > 0 &&
                        std::find(on_edge.vertices.begin(), on_edge.vertices.end(), on_vs[0]) == on_edge.vertices.end())
                        on_edge.vertices.push_back(on_vs[0]);
                    continue;
                }
                if (neg_vs.size() + on_vs.size() == 2) {
                    neg_face.edges.push_back(old_e_id);
                    if (on_vs.size() > 0 &&
                        std::find(on_edge.vertices.begin(), on_edge.vertices.end(), on_vs[0]) == on_edge.vertices.end())
                        on_edge.vertices.push_back(on_vs[0]);
                    continue;
                }

                ///splitting...
                BSPEdge pos_edge, neg_edge;
                edges.push_back(neg_edge);
                int new_e_id = edges.size() - 1;
                std::vector<int> new_e_ids = {old_e_id, new_e_id};

                pos_edge.vertices.push_back(pos_vs[0]);
                neg_edge.vertices.push_back(neg_vs[0]);

                int new_v_id = 0;
                Segment_3 seg(vertices[edges[old_e_id].vertices[0]], vertices[edges[old_e_id].vertices[1]]);
                auto result = intersection(seg, pln);
                if (result) {
                    const Point_3 *p = boost::get<Point_3>(&*result);
                    vertices.push_back(*p);

                    new_v_id = vertices.size() - 1;
                    on_edge.vertices.push_back(new_v_id);
                    pos_edge.vertices.push_back(new_v_id);
                    neg_edge.vertices.push_back(new_v_id);

                    v_sides[new_v_id] = V_ON;//fixed
                } else {
                    log_and_throw("error cal p!");
                }

                ///add edges
                pos_edge.conn_faces = edges[old_e_id].conn_faces;
                neg_edge.conn_faces = edges[old_e_id].conn_faces;
                for (auto it = edges[old_e_id].conn_faces.begin(); it != edges[old_e_id].conn_faces.end(); it++) {
                    if (*it == old_f_id)
                        continue;
                    faces[*it].edges.push_back(new_e_id);
                    faces[*it].vertices.push_back(new_v_id);
                }
                edges[new_e_ids[0]] = pos_edge;//if get here, it means that old_edge has been cut into 2
                edges[new_e_ids[1]] = neg_edge;

                ///add edges for faces
                pos_face.edges.push_back(new_e_ids[0]);
                neg_face.edges.push_back(new_e_ids[1]);
            }//split one face end

            if (pos_face.edges.size() == 0 || neg_face.edges.size() == 0)//connected pos/neg
                continue;

            ///from now, the face would definitely be subdivided
            faces.push_back(neg_face);//have to do push_back here!!! Otherwise would producing empty faces!!

            ///clean conn_faces for neg_face's edges//fixed
            for(int j=0;j<neg_face.edges.size();j++){
                auto it = std::find(edges[neg_face.edges[j]].conn_faces.begin(), edges[neg_face.edges[j]].conn_faces.end(), old_f_id);
                if(it!=edges[neg_face.edges[j]].conn_faces.end()) {
                    edges[neg_face.edges[j]].conn_faces.erase(it);
                    edges[neg_face.edges[j]].conn_faces.insert(new_f_id);
                }
            }

            ///add on_edge
            //remove duplicated vertices//fixed
            on_edge.conn_faces = {new_f_ids[0], new_f_ids[1]};
            edges.push_back(on_edge);
            int on_e_id = edges.size() - 1;
            pos_face.edges.push_back(on_e_id);
            neg_face.edges.push_back(on_e_id);
            on_face.edges.push_back(on_e_id);

            ///add faces
            pos_face.conn_nodes = faces[old_f_id].conn_nodes;
            neg_face.conn_nodes = faces[old_f_id].conn_nodes;
            ///cal vertices for faces//fixed
            getVertices(pos_face);
            getVertices(neg_face);
            for (auto it = faces[old_f_id].conn_nodes.begin(); it != faces[old_f_id].conn_nodes.end(); it++) {
                if (*it == old_n_id)
                    continue;
                nodes[*it].faces.push_back(new_f_id);
            }
            ///re-assign divfaces for pos_face & neg_face
            std::unordered_set<int> tmp_df_ids = faces[old_f_id].div_faces;
            for(auto it=tmp_df_ids.begin(); it!=tmp_df_ids.end();it++) {
                int side = divfaceSide(pln, div_faces[*it], div_vertices);
                if (side == DIVFACE_POS)
                    pos_face.div_faces.insert(*it);
                else if (side == DIVFACE_NEG)
                    neg_face.div_faces.insert(*it);
                else if(side==DIVFACE_CROSS){
                    pos_face.div_faces.insert(*it);
                    neg_face.div_faces.insert(*it);
                }
            }
            pos_face.matched_f_id=faces[old_f_id].matched_f_id;
            neg_face.matched_f_id=faces[old_f_id].matched_f_id;

            faces[new_f_ids[0]] = pos_face;
            faces[new_f_ids[1]] = neg_face;

            for(int j=0;j<new_f_ids.size();j++){//fixed
                for(auto it=faces[new_f_ids[j]].edges.begin();it!=faces[new_f_ids[j]].edges.end();it++)///fixed
                    edges[*it].conn_faces.insert(new_f_ids[j]);
            }

            ///add faces for nodes
            pos_node.faces.push_back(new_f_ids[0]);
            neg_node.faces.push_back(new_f_ids[1]);
        }//split one node end


        ///clean conn_nodes for neg_node's faces//fixed //it must be done before adding on_face!! Otherwise, the on_face would be influenced.
        for(int i=0;i<neg_node.faces.size();i++){
            auto it = std::find(faces[neg_node.faces[i]].conn_nodes.begin(), faces[neg_node.faces[i]].conn_nodes.end(), old_n_id);
            if(it!=faces[neg_node.faces[i]].conn_nodes.end()) {
                faces[neg_node.faces[i]].conn_nodes.erase(it);
                faces[neg_node.faces[i]].conn_nodes.insert(new_n_id);
            }
        }

        ///add on_face
        on_face.conn_nodes = {new_n_ids[0], new_n_ids[1]};
        ///cal vertices for faces//fixed
        getVertices(on_face);
        faces.push_back(on_face);
        int on_f_id = faces.size() - 1;
        pos_node.faces.push_back(on_f_id);
        neg_node.faces.push_back(on_f_id);
        for(auto it=faces[on_f_id].edges.begin();it!=faces[on_f_id].edges.end();it++)///fixed
            edges[*it].conn_faces.insert(on_f_id);

        ///add nodes
        nodes[new_n_ids[0]] = pos_node;
        nodes[new_n_ids[1]] = neg_node;

        ///check if divface.size==0
        for (int i = 0; i < new_n_ids.size(); i++) {
            if (nodes[new_n_ids[i]].div_faces.size() == 0)
                nodes[new_n_ids[i]].is_leaf = true;
            else
                processing_n_ids.push(new_n_ids[i]);
        }
    }

}

void BSPSubdivision::calVertexSides(const Plane_3& pln, const std::unordered_set<int>& v_ids, const std::vector<Point_3>& vs,
                                    std::unordered_map<int, int>& v_sides){
    for(auto it=v_ids.begin();it!=v_ids.end();it++){
        CGAL::Oriented_side side=pln.oriented_side(vs[*it]);
        switch (side) {
            case CGAL::ON_ORIENTED_BOUNDARY:
                v_sides[*it]=V_ON;
                break;
            case CGAL::ON_POSITIVE_SIDE:
                v_sides[*it]=V_POS;
                break;
            case CGAL::ON_NEGATIVE_SIDE:
                v_sides[*it]=V_NEG;
        }
    }
}

int BSPSubdivision::divfaceSide(const Plane_3& pln, const std::array<int, 3>& p_ids,
                                    const std::vector<Point_3>& ps) {
    int cnt_pos = 0, cnt_on = 0, cnt_neg = 0;
    for (int i = 0; i < p_ids.size(); i++) {
        switch (pln.oriented_side(ps[p_ids[i]])) {
            case CGAL::ON_ORIENTED_BOUNDARY:
                cnt_on++;
                break;
            case CGAL::ON_POSITIVE_SIDE:
                cnt_pos++;
                break;
            case CGAL::ON_NEGATIVE_SIDE:
                cnt_neg++;
        }
    }

    if (cnt_pos > 0 && cnt_neg > 0)
        return DIVFACE_CROSS;
    if (cnt_on == p_ids.size())
        return DIVFACE_ON;
    if (cnt_neg == 0)
        return DIVFACE_POS;
    if (cnt_pos == 0)
        return DIVFACE_NEG;

    return -1;
}

void BSPSubdivision::getVertices(BSPFace& face){
    std::unordered_set<int> vs;
    for(int i=0;i<face.edges.size();i++){
        vs.insert(MC.bsp_edges[face.edges[i]].vertices[0]);
        vs.insert(MC.bsp_edges[face.edges[i]].vertices[1]);
    }
    for(auto it=vs.begin();it!=vs.end();it++)
        face.vertices.push_back(*it);
}

} // namespace tetwild
