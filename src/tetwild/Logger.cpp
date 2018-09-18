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
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/details/registry.h>
#include <spdlog/details/thread_pool.h>
#include <tetwild/EnableWarnings.h>
#include <memory>
#include <mutex>
#include <iostream>

namespace tetwild {

std::shared_ptr<spdlog::async_logger> Logger::logger_;

// Some code was copied over from <spdlog/async.h>
void Logger::init(bool use_cout, const std::string &filename, bool truncate) {
	std::vector<spdlog::sink_ptr> sinks;
	if (use_cout) {
		sinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
	}
	if (!filename.empty()) {
		sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, truncate));
	}

	auto &registry_inst = spdlog::details::registry::instance();

	// create global thread pool if not already exists..
	std::lock_guard<std::recursive_mutex> tp_lock(registry_inst.tp_mutex());
	auto tp = registry_inst.get_tp();
	if (tp == nullptr) {
		tp = std::make_shared<spdlog::details::thread_pool>(spdlog::details::default_async_q_size, 1);
		registry_inst.set_tp(tp);
	}

    logger_ = std::make_shared<spdlog::async_logger>("tetwild", sinks.begin(), sinks.end(), std::move(tp), spdlog::async_overflow_policy::block);
    registry_inst.register_and_init(logger_);
}

} // namespace tetwild
