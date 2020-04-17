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

#include <tetwild/State.h>
#include <tetwild/Args.h>
#include <tetwild/Exception.h>
#include <igl/bounding_box_diagonal.h>

namespace tetwild {

State::State(const Args &args, const Eigen::MatrixXd &V)
    : working_dir(args.working_dir)
    , postfix(args.postfix)
    , stat_file(args.csv_file)
    , bbox_diag(igl::bounding_box_diagonal(V))
    , eps_input(bbox_diag * args.eps_rel)
    , eps_delta(args.sampling_dist_rel > 0 ? 0 : eps_input / args.stage / std::sqrt(3))
    , initial_edge_len(args.getAbsoluteEdgeLength(bbox_diag))
{
    if (args.sampling_dist_rel > 0) {
        //for testing only
        sampling_dist = bbox_diag * args.sampling_dist_rel / 100.0;
        eps = bbox_diag * args.eps_rel;
        eps_2 = eps * eps;
        if (args.stage != 1) {
            throw TetWildError("args.stage should be equal to 1.");
        }
    } else {
        // d_err = d/sqrt(3)
        sampling_dist = eps_input / args.stage;
        eps = eps_input - sampling_dist / std::sqrt(3) * (args.stage + 1 - sub_stage);
        eps_2 = eps * eps;
        // eps_delta = sampling_dist / std::sqrt(3);
    }

   // logger().debug("eps = {}", eps);
   // logger().debug("ideal_l = {}", initial_edge_len);
}

} // namespace tetwild
