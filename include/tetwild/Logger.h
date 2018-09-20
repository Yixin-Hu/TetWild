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

#include <tetwild/Exception.h>
#include <tetwild/DisableWarnings.h>
#include <spdlog/spdlog.h>
#include <spdlog/async.h>
#include <fmt/ranges.h>
#include <tetwild/EnableWarnings.h>

namespace tetwild {

struct Logger {
	static std::shared_ptr<spdlog::async_logger> logger_;

	// By default, write to stdout, but don't write to any file
	static void init(bool use_cout = true, const std::string &filename = "", bool truncate = true);
};

// Retrieve current logger, or create one if not available
inline spdlog::async_logger & logger() {
	if (!Logger::logger_) {
		Logger::init();
	}
	return *Logger::logger_;
}

template<typename T>
[[noreturn]] void log_and_throw(T x) {
	logger().error(x);
	throw TetWildError(x);
}

} // namespace tetwild
