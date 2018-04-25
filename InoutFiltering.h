//
// Created by yihu on 6/13/17.
//

#ifndef NEW_GTET_INOUTFILTERING_H
#define NEW_GTET_INOUTFILTERING_H

#include "heads.h"
#include "TetmeshElements.h"

class InoutFiltering {
public:
    std::vector<TetVertex>& tet_vertices;
    std::vector<std::array<int, 4>>& tets;
    std::vector<std::array<int, 4>>& is_surface_fs;
    std::vector<bool>& v_is_removed;
    std::vector<bool>& t_is_removed;
    std::vector<TetQuality>& tet_qualities;

    std::vector<bool> is_inside;
    InoutFiltering(std::vector<TetVertex>& t_vs, std::vector<std::array<int, 4>>& ts,
                   std::vector<std::array<int, 4>>& is_sf_fs,
                   std::vector<bool>& v_is_rm, std::vector<bool>& t_is_rm, std::vector<TetQuality>& tet_qs):
            tet_vertices(t_vs), tets(ts), is_surface_fs(is_sf_fs), v_is_removed(v_is_rm), t_is_removed(t_is_rm),
            tet_qualities(tet_qs){}

    void getSurface(Eigen::MatrixXd& V_sf, Eigen::MatrixXi& F_sf);
    void filter();

    void outputWindingNumberField(const Eigen::VectorXd& W);
};


#endif //NEW_GTET_INOUTFILTERING_H
