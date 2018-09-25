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

#include <tetwild/Logger.h>
#include <tetwild/DisableWarnings.h>
#include <spdlog/details/registry.h>
#include <tetwild/EnableWarnings.h>
#include <memory>
#include <mutex>
#include <iostream>

namespace tetwild {

std::shared_ptr<spdlog::async_logger> Logger::logger_;

// Some code was copied over from <spdlog/async.h>
void Logger::init(bool use_cout, const std::string &filename, bool truncate) {
	throw TetWildError("not implemented");
}

} // namespace tetwild
