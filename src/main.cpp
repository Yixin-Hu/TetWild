// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <tetwild/tetwild.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>
#include <tetwild/MeshRefinement.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <pymesh/MshSaver.h>
#include <tetwild/DisableWarnings.h>
#include <CLI/CLI.hpp>
#include <tetwild/EnableWarnings.h>

using namespace tetwild;

namespace tetwild {
    void extractFinalTetmesh(MeshRefinement& MR, Eigen::MatrixXd &V_out, Eigen::MatrixXi &T_out, Eigen::VectorXd &A_out, const Args &args);
} // namespace tetwild

void saveFinalTetmesh(const std::string &output_filename,
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXd &A, const Args &args)
{
    logger().debug("Writing mesh to {}...", output_filename);
    std::string output_format = output_filename.substr(output_filename.size() - 4, 4);
    if (output_format == "mesh") {
        std::ofstream f(output_filename);
        f.precision(std::numeric_limits<double>::digits10 + 1);
        f << "MeshVersionFormatted 1" << std::endl;
        f << "Dimension 3" << std::endl;

        f << "Vertices" << std::endl << V.rows() << std::endl;
        for (int i = 0; i < V.rows(); i++)
            f << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << " " << 0 << std::endl;
        f << "Triangles" << std::endl << 0 <<std::endl;
        f << "Tetrahedra" << std::endl;
        f << T.rows() << std::endl;
        for (int i = 0; i < T.rows(); i++) {
            for (int j = 0; j < 4; j++) {
                f << T(i, j) + 1 << " ";
            }
            f << 0 << std::endl;
        }

        f << "End";
        f.close();
    } else {
        PyMesh::MshSaver mSaver(output_filename, true);
        PyMesh::VectorF V_flat(V.size());
        PyMesh::VectorI T_flat(T.size());
        Eigen::MatrixXd VV = V.transpose();
        Eigen::MatrixXi TT = T.transpose();
        std::copy_n(VV.data(), V.size(), V_flat.data());
        std::copy_n(TT.data(), T.size(), T_flat.data());
        mSaver.save_mesh(V_flat, T_flat, 3, mSaver.TET);
        mSaver.save_elem_scalar_field("min_dihedral_angle", A);
    }

    if (args.is_quiet) {
        return;
    }
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    extractSurfaceMesh(V, T, V_sf, F_sf);
    igl::writeOBJ(State::state().working_dir+State::state().postfix+"_sf.obj", V_sf, F_sf);
}

void gtet_new_slz(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, const std::string& slz_file,
                  int max_pass, const std::array<bool, 4>& ops,
                  Eigen::MatrixXd &VO, Eigen::MatrixXi &TO, Eigen::VectorXd &AO,
                  const Args &args = Args())
{
    MeshRefinement MR(args);
    MR.deserialization(VI, FI, slz_file);

//    MR.is_dealing_unrounded = true;
    MR.refine(State::state().ENERGY_AMIPS, ops, false, true);

    extractFinalTetmesh(MR, VO, TO, AO, args); //do winding number and output the tetmesh
}

int main(int argc, char *argv[]) {
    int log_level = 1; // debug
    std::string log_filename = "";
    std::string input_filename;
    std::string output_filename;
    Args args;

    CLI::App app{"RobustTetMeshing"};
    app.add_option("input,--input", input_filename, "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")->required();
    app.add_option("output,--output", output_filename, "Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')");
    app.add_option("--postfix", args.postfix, "Postfix P for output files. (string, optional, default: '_')");
    app.add_option("-l,--ideal-edge-length", args.initial_edge_len_rel, "ideal_edge_length = diag_of_bbox * L / 100. (double, optional, default: 5%)");
    app.add_option("-e,--epsilon", args.eps_rel, "epsilon = diag_of_bbox * EPS / 100. (double, optional, default: 0.1%)");
    app.add_option("--stage", args.stage, "Run pipeline in stage STAGE. (integer, optional, default: 1)");
    app.add_option("--filter-energy", args.filter_energy_thres, "Stop mesh improvement when the maximum energy is smaller than ENERGY. (double, optional, default: 10)");
    app.add_option("--max-pass", args.max_num_passes, "Do PASS mesh improvement passes in maximum. (integer, optional, default: 80)");

    app.add_flag("--is-laplacian", args.smooth_open_boundary, "Do Laplacian smoothing for the surface of output on the holes of input (optional)");
    app.add_option("--targeted-num-v", args.target_num_vertices, "Output tetmesh that contains TV vertices. (integer, optional, tolerance: 5%)");
    app.add_option("--bg-mesh", args.background_mesh, "Background tetmesh BGMESH in .msh format for applying sizing field. (string, optional)");
    app.add_flag("-q,--is-quiet", args.is_quiet, "Mute console output. (optional)");
    app.add_option("--log", log_filename, "Log info to given file.");
    app.add_option("--level", log_level, "Log level (0 = most verbose, 6 = off).");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    Logger::init(!args.is_quiet, log_filename);
    log_level = std::max(0, std::min(6, log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    //initialization
    GEO::initialize();
    State::state().postfix = args.postfix;
    if(args.slz_file != "")
        State::state().working_dir = input_filename.substr(0, args.slz_file.size() - 4);
    else
        State::state().working_dir = input_filename.substr(0, input_filename.size() - 4);

    if(args.csv_file == "")
        State::state().stat_file = State::state().working_dir + State::state().postfix + ".csv";
    else
        State::state().stat_file = args.csv_file;

    if(output_filename.empty()) {
        output_filename = State::state().working_dir + State::state().postfix + ".msh";
    }

    if(args.is_quiet) {
        args.write_csv_file = false;
    }

    //do tetrahedralization
    Eigen::MatrixXd VI, VO;
    Eigen::MatrixXi FI, TO;
    Eigen::VectorXd AO;
    igl::read_triangle_mesh(input_filename, VI, FI);
    if(args.slz_file != "") {
        gtet_new_slz(VI, FI, args.slz_file, args.max_num_passes,
            {{true, false, true, true}}, VO, TO, AO, args);
    } else {
        tetwild::tetrahedralization(VI, FI, VO, TO, AO, args);
    }
    saveFinalTetmesh(output_filename, VO, TO, AO, args);

    spdlog::shutdown();

    return 0;
}
