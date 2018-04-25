//
// Created by yihu on 8/22/17.
//

#ifndef NEW_GTET_BSPELEMENTS_H
#define NEW_GTET_BSPELEMENTS_H
#include <vector>
#include <unordered_set>

class BSPEdge{
public:
    std::vector<int> vertices;
    std::unordered_set<int> conn_faces;

    BSPEdge(){}
    BSPEdge(int v1, int v2){
        vertices={v1, v2};
    }
};

class BSPFace{
public:
    std::vector<int> vertices;
    std::vector<int> edges;
    std::unordered_set<int> conn_nodes;
    std::unordered_set<int> div_faces;

    int matched_f_id=-1;
};

class BSPtreeNode{
public:
    bool is_leaf=false;
    std::vector<int> faces;
    std::unordered_set<int> div_faces;
};

#endif //NEW_GTET_BSPELEMENTS_H
