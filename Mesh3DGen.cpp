#include "tetwild.h"
#include <Mesh3DGen.h>
#include <fstream>

namespace nTopology
{
namespace meshGen
{

  Mesh3D::Mesh3D()
  {}

  Mesh3D::~Mesh3D()
  {}

  void Mesh3D::mesh3DGen(const std::vector<glm::vec3>& vertPositions,
                        const std::vector<uint32_t>& faceIndices,
                        std::vector<glm::vec3>& tetVertPos,
                        std::vector<std::array<int, 4>>& tetIndices,
                        float targetEdgeLength) 
  {

    // Input surface mesh size
    auto numFaces = faceIndices.size() / 3;
    auto numVerts = vertPositions.size();

    // input
    Eigen::MatrixXd v2d(numVerts, 3);
    Eigen::MatrixXi f2d(numFaces, 3);


    // output
    std::vector<std::array<double, 3>> tVerts;
    std::vector<std::array<int, 4>> tTets;

    auto fts = std::vector<std::array<int, 3>>(numFaces);
    auto fvs = std::vector<std::array<double, 3>>(numVerts);

    for (auto i = 0; i != numFaces; ++i)
    {
      std::array<int, 3> fc = { (int)faceIndices[3 * i], (int)faceIndices[3 * i + 1], (int)faceIndices[3 * i + 2] } ;
      fts[i] = fc;
    }
    for (auto i = 0; i != numVerts; ++i)
    {
      auto vp = vertPositions[i];
      std::array<double, 3> vt = { (double)vp.x, (double)vp.y, (double)vp.z };
      fvs[i] = vt;
    }

    // execute
    args.is_quiet = true;
    args.output = "";
    args.max_pass = 80;
    args.targeted_num_v = -1;
    args.filter_energy = 10;
    args.i_ideal_edge_length = targetEdgeLength;
    //args.i_epsilon = eps;

    std::string input = "a";
    tetwild::tetrahedralization(fvs, fts, tVerts, tTets);
    //tetwild::gtet_new(v2d, f2d, tVerts, tTets);

    tetVertPos.resize(tVerts.size());
    tetIndices.resize(tTets.size());

    for (auto i = 0; i != tVerts.size(); ++i)
    {
      tetVertPos[i].x = tVerts[i][0];
      tetVertPos[i].y = tVerts[i][1];
      tetVertPos[i].z = tVerts[i][2];
    }
    for (auto i = 0; i != tTets.size(); ++i)
    {
      for (auto j = 0; j != 4; ++j)
      {
        tetIndices[i][j] = tTets[i][j];
      }
    }
  }

}

}