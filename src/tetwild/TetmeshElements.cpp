// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Jeremie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Jeremie Dumas on 09/10/18.
//

#include <tetwild/TetmeshElements.h>
#include <tetwild/Logger.h>
#include <tetwild/Serialization.h>

namespace tetwild {

void TetVertex::printInfo() const {
    logger().debug("is_on_surface = {}", is_on_surface);
    logger().debug("is_on_bbox = {}", is_on_bbox);
    logger().debug("conn_tets = {}", conn_tets);
}

void Stage::serialize(std::string serialize_file) {
    igl::serialize(tet_vertices, "tet_vertices", serialize_file, true);
    igl::serialize(tets, "tets", serialize_file);
    igl::serialize(is_surface_fs, "tets", serialize_file);
    igl::serialize(v_is_removed, "v_is_removed", serialize_file);
    igl::serialize(t_is_removed, "t_is_removed", serialize_file);
    igl::serialize(tet_qualities, "tet_qualities", serialize_file);

    igl::serialize(is_shown, "is_shown", serialize_file);
    igl::serialize(resolution, "resolution", serialize_file);
}

void Stage::deserialize(std::string serialize_file) {
    igl::deserialize(tet_vertices, "tet_vertices", serialize_file);
    igl::deserialize(tets, "tets", serialize_file);
    igl::deserialize(is_surface_fs, "tets", serialize_file);
    igl::deserialize(v_is_removed, "v_is_removed", serialize_file);
    igl::deserialize(t_is_removed, "t_is_removed", serialize_file);
    igl::deserialize(tet_qualities, "tet_qualities", serialize_file);

    igl::deserialize(is_shown, "is_shown", serialize_file);
    igl::deserialize(resolution, "resolution", serialize_file);
}

} // namespace tetwild
