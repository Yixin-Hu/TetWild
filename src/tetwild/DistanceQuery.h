// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Jeremie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Jeremie Dumas on 09/04/18.
//

#pragma once

#include <string>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

namespace tetwild {

inline void get_point_facet_nearest_point(
    const GEO::Mesh& M,
    const GEO::vec3& p,
    GEO::index_t f,
    GEO::vec3& nearest_p,
    double& squared_dist
) {
    using namespace GEO;
    geo_debug_assert(M.facets.nb_vertices(f) == 3);
    index_t c = M.facets.corners_begin(f);
    const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
    ++c;
    const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
    ++c;
    const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
    double lambda1, lambda2, lambda3;  // barycentric coords, not used.
    squared_dist = Geom::point_triangle_squared_distance(
        p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3
    );
}

} // namespace tetwild
