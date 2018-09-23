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
////////////////////////////////////////////////////////////////////////////////
#include "utils.h"
#include <tetwild/Logger.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/basic/progress.h>
////////////////////////////////////////////////////////////////////////////////

namespace tetwild {

void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, GEO::Mesh &M) {
	M.clear();
	// Setup vertices
	M.vertices.create_vertices((int) V.rows());
	for (int i = 0; i < (int) M.vertices.nb(); ++i) {
		GEO::vec3 &p = M.vertices.point(i);
		p[0] = V(i, 0);
		p[1] = V(i, 1);
		p[2] = (V.cols() == 2 ? 0 : V(i, 2));
	}
	// Setup faces
	if (F.cols() == 3) {
		M.facets.create_triangles((int) F.rows());
	} else if (F.cols() == 4) {
		M.facets.create_quads((int) F.rows());
	} else {
		throw std::runtime_error("Mesh faces not supported");
	}
	for (int c = 0; c < (int) M.facets.nb(); ++c) {
		for (int lv = 0; lv < F.cols(); ++lv) {
			M.facets.set_vertex(c, lv, F(c, lv));
		}
	}
    M.facets.connect();
}

// -----------------------------------------------------------------------------

void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &T, GEO::Mesh &M) {
	to_geogram_mesh(V, F, M);
	if (T.cols() == 4) {
		M.cells.create_tets((int) T.rows());
	} else if (T.rows() != 0) {
		throw std::runtime_error("Mesh cells not supported");
	}
	for (int c = 0; c < (int) M.cells.nb(); ++c) {
		for (int lv = 0; lv < T.cols(); ++lv) {
			M.cells.set_vertex(c, lv, T(c, lv));
		}
	}
	M.cells.connect();
}

// -----------------------------------------------------------------------------

namespace {

void create_box_mesh(const Eigen::RowVector3d &pmin, const Eigen::RowVector3d &pmax, double padding, GEO::Mesh &M) {
	Eigen::MatrixXd V(9, 3);
	V << 0, 0, 0,
		0, 0, 1,
		0, 1, 0,
		0, 1, 1,
		1, 0, 0,
		1, 0, 1,
		1, 1, 0,
		1, 1, 1,
		0.494941, 0.652018, 0.319279;
	Eigen::MatrixXi T(12, 4);
	T << 1, 3, 4, 9,
		7, 4, 9, 8,
		2, 9, 4, 8,
		2, 9, 6, 1,
		9, 2, 6, 8,
		2, 9, 1, 4,
		6, 9, 5, 1,
		6, 7, 9, 8,
		6, 7, 5, 9,
		5, 3, 1, 9,
		7, 4, 3, 9,
		3, 5, 7, 9;
	T.array() -= 1;
	V.array() = V.array().rowwise() * (pmax.array() + padding) + (1.0 - V.array()).rowwise() * (pmin.array() - padding);
	Eigen::MatrixXi F(0, 3);
	to_geogram_mesh(V, F, T, M);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

void sample_bbox(const Eigen::MatrixXd &V, int num_samples, double padding,
	Eigen::MatrixXd &P, int num_lloyd, int num_newton)
{
	assert(num_samples > 3);
	// bool was_quiet = GEO::Logger::instance()->is_quiet();
	// GEO::Logger::instance()->set_quiet(true);

	Eigen::RowVector3d pmin = V.colwise().minCoeff();
	Eigen::RowVector3d pmax = V.colwise().maxCoeff();
	GEO::Mesh M;
	create_box_mesh(pmin, pmax, padding, M);
	GEO::CentroidalVoronoiTesselation CVT(&M);
	CVT.set_volumetric(true);

	CVT.compute_initial_sampling(num_samples);
	CVT.resize_points(num_samples + V.rows() + 8);

	// Constrained points
	for (int v = 0; v < V.rows(); ++v) {
		Eigen::RowVector3d p = V.row(v);
		std::copy_n(p.data(), 3, CVT.embedding(num_samples + v));
		CVT.lock_point(num_samples + v);
	}
	// Bounding box
	Eigen::Matrix<double, 8, 3, Eigen::RowMajor> B;
	B << 0, 0, 0,
		0, 0, 1,
		0, 1, 0,
		0, 1, 1,
		1, 0, 0,
		1, 0, 1,
		1, 1, 0,
		1, 1, 1;
	B.array() = B.array().rowwise() * (pmax.array() + padding) + (1.0 - B.array()).rowwise() * (pmin.array() - padding);
	std::copy_n(B.data(), 3*8, CVT.embedding(num_samples + V.rows()));
	for (int v = 0; v < 8; ++v) {
		CVT.lock_point(num_samples + V.rows() + v);
	}
	CVT.set_show_iterations(true);

	// if (num_lloyd > 0) {
	// 	CVT.Lloyd_iterations(num_lloyd);
	// }

	// if (num_newton > 0) {
	// 	CVT.Newton_iterations(num_newton);
	// }

	if (num_lloyd > 0) {
		try {
			GEO::ProgressTask progress("Lloyd", 100);
			CVT.set_progress_logger(&progress);
			CVT.Lloyd_iterations(num_lloyd);
		} catch(const GEO::TaskCanceled&) {
		}
	}

	if (num_newton > 0) {
		try {
			GEO::ProgressTask progress("Newton", 100);
			CVT.set_progress_logger(&progress);
			CVT.Newton_iterations(num_newton);
		} catch(const GEO::TaskCanceled&) {
		}
	}

	P.resize(3, CVT.nb_points());
	std::copy_n(CVT.embedding(0), 3*CVT.nb_points(), P.data());
	P.transposeInPlace();

	// GEO::Logger::instance()->set_quiet(was_quiet);
}

// -----------------------------------------------------------------------------

void resample_surface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int num_samples,
	Eigen::MatrixXd &P, int num_lloyd, int num_newton)
{
	assert(num_samples > 3);
	// bool was_quiet = GEO::Logger::instance()->is_quiet();
	// GEO::Logger::instance()->set_quiet(true);

	GEO::Mesh M;
	to_geogram_mesh(V, F, M);
	GEO::CentroidalVoronoiTesselation CVT(&M);
	CVT.compute_initial_sampling(num_samples);
	CVT.set_show_iterations(true);

	if (num_lloyd > 0) {
		try {
			GEO::ProgressTask progress("Lloyd", 100);
			CVT.set_progress_logger(&progress);
			CVT.Lloyd_iterations(num_lloyd);
		} catch(const GEO::TaskCanceled&) {
		}
	}

	if (num_newton > 0) {
		try {
			GEO::ProgressTask progress("Newton", 100);
			CVT.set_progress_logger(&progress);
			CVT.Newton_iterations(num_newton);
		} catch(const GEO::TaskCanceled&) {
		}
	}

	P.resize(3, num_samples);
	std::copy_n(CVT.embedding(0), 3*num_samples, P.data());
	P.transposeInPlace();

	// GEO::Logger::instance()->set_quiet(was_quiet);
}

// -----------------------------------------------------------------------------

void delaunay_tetrahedralization(const Eigen::MatrixXd &V, Eigen::MatrixXi &T) {
	assert(V.cols() == 3);
	int n = (int) V.rows();

	// Compute tetrahedralization
	GEO::Delaunay_var delaunay = GEO::Delaunay::create(3, "BDEL");
	const Eigen::MatrixXd P = V.transpose();
	delaunay->set_vertices(n, P.data());

	// Extract tetrahedra
	T.resize(delaunay->nb_cells(), 4);
	for (int c = 0; c < (int) delaunay->nb_cells(); ++c) {
		for (int lv = 0; lv < 4; ++lv) {
			GEO::signed_index_t v = delaunay->cell_vertex(c, lv);
			T(c, lv) = v;
		}
	}
}

} // namespace tetwild
