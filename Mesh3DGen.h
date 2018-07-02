#include "tetwild.h"
#include <glm/glm.hpp>

namespace nTopology
{
namespace meshGen
{

class Mesh3D
{
public:
  Mesh3D();
  virtual ~Mesh3D();

  static void mesh3DGen(const std::vector<glm::vec3>& vertPositions,
                 const std::vector<uint32_t>& faceIndices,
                 std::vector<glm::vec3>& tetVertPos,
                 std::vector<std::array<int, 4>>& tetIndices,
                 float targetEdgeLength);

};

}
}

