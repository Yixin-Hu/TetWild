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

namespace tetwild {

struct Args;
struct State;
struct MeshRecord;
class BSPFace;
class MeshConformer;
class EdgeCollapser;
class EdgeSplitter;
class EdgeRemover;
class VertexSmoother;

} // namespace tetwild
