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


#include <tetwild/ForwardDecls.h>
#include <tetwild/TetmeshElements.h>
#include <tetwild/Logger.h>
#include <tetwild/Utils.h>
#include <mmg/libmmg.h>
#include <igl/edge_lengths.h>
#include <igl/face_areas.h>
#include <igl/dihedral_angles.h>
#include <Eigen/Dense>

namespace tetwild {

namespace {

#define _MMG5_EPSD      1.e-30
#define _MMG5_EPSD2     1.0e-200
#define _MMG5_EPS       1.e-06
#define _MMG5_EPSOK     1.e-15
#define _MMG5_NULKAL    1.e-30
#define _MMG3D_ALPHAD    20.7846096908265    //0.04811252243247      /* 12*sqrt(3) */

/**
 * \param a pointer toward the coor of the first tetra vertex.
 * \param b pointer toward the coor of the second tetra vertex.
 * \param c pointer toward the coor of the third tetra vertex.
 * \param d pointer toward the coor of the fourth tetra vertex.
 * \return The isotropic quality of the tet.
 *
 * Compute the quality of a tetra given by 4 points a,b,c,d with respect to the
 * isotropic metric \a met.
 *
 */
static
inline double mmg_caltet_iso_4pt(const double *a, const double *b, const double *c, const double *d) {
    double       abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz;
    double       cdx,cdy,cdz;
    double       vol,v1,v2,v3,rap;

    /* volume */
    abx = b[0] - a[0];
    aby = b[1] - a[1];
    abz = b[2] - a[2];
    rap = abx*abx + aby*aby + abz*abz;

    acx = c[0] - a[0];
    acy = c[1] - a[1];
    acz = c[2] - a[2];
    rap += acx*acx + acy*acy + acz*acz;

    adx = d[0] - a[0];
    ady = d[1] - a[1];
    adz = d[2] - a[2];
    rap += adx*adx + ady*ady + adz*adz;

    v1  = acy*adz - acz*ady;
    v2  = acz*adx - acx*adz;
    v3  = acx*ady - acy*adx;
    vol = abx * v1 + aby * v2 + abz * v3;
    if ( vol < _MMG5_EPSD2 )  return(0.0);

    bcx = c[0] - b[0];
    bcy = c[1] - b[1];
    bcz = c[2] - b[2];
    rap += bcx*bcx + bcy*bcy + bcz*bcz;

    bdx = d[0] - b[0];
    bdy = d[1] - b[1];
    bdz = d[2] - b[2];
    rap += bdx*bdx + bdy*bdy + bdz*bdz;

    cdx = d[0] - c[0];
    cdy = d[1] - c[1];
    cdz = d[2] - c[2];
    rap += cdx*cdx + cdy*cdy + cdz*cdz;
    if ( rap < _MMG5_EPSD2 )  return(0.0);

    /* quality = vol / len^3/2 */
    rap = rap * sqrt(rap);
    return(vol / rap);
}

/** Compute 3 * 3 determinant : det(c1-c0,c2-c0,v) */
inline double mmg_det3pt1vec(const double *c0,const double *c1,const double *c2,const double *v) {
    double m00,m10,m20,m01,m11,m21,det;

    m00 = c1[0] - c0[0] ; m01 = c2[0] - c0[0];
    m10 = c1[1] - c0[1] ; m11 = c2[1] - c0[1];
    m20 = c1[2] - c0[2] ; m21 = c2[2] - c0[2];
    det = v[0]*(m10*m21 - m20*m11) -v[1]*(m00*m21-m20*m01) + v[2]*(m00*m11-m10*m01);

    return(det);
}

/** Compute 3 * 3 determinant : det(c1-c0,c2-c0,c3-c0) */
inline double mmg_det4pt(const double *c0,const double *c1,const double *c2,const double *c3) {
  double m[3];

  m[0] = c3[0] - c0[0];
  m[1] = c3[1] - c0[1];
  m[2] = c3[2] - c0[2];

  return( mmg_det3pt1vec(c0,c1,c2,m) );
}

// int _MMG5_minQualCheck ( int iel, double minqual, double alpha )
// {
//   double minqualOnAlpha;

//   minqualOnAlpha = minqual/alpha;

//   if ( minqualOnAlpha < _MMG5_NULKAL ) {
//     fprintf(stderr,"\n  ## Error: %s: too bad quality for the worst element: "
//             "(elt %d -> %15e)\n",__func__,iel,minqual);
//     return(0);
//   }
//   else if ( minqualOnAlpha < _MMG5_EPSOK ) {
//     fprintf(stderr,"\n  ## Warning: %s: very bad quality for the worst element: "
//             "(elt %d -> %15e)\n",__func__,iel,minqual);
//   }

//   return(1);
// }

// int _MMG3D_tetraQual(MMG5_pMesh mesh, MMG5_pSol met,char metRidTyp) {
//   MMG5_pTetra pt;
//   double      minqual;
//   int         k,iel;

//   minqual = 2./_MMG3D_ALPHAD;

//   /*compute tet quality*/
//   for (k=1; k<=mesh->ne; k++) {
//     pt = &mesh->tetra[k];
//      if( !MG_EOK(pt) )   continue;

//      if ( !metRidTyp && met->size == 6 && met->m ) {
//        pt->qual = _MMG5_caltet33_ani(mesh,met,pt);
//      }
//      else
//        pt->qual = _MMG5_orcal(mesh,met,k);

//     if ( pt->qual < minqual ) {
//       minqual = pt->qual;
//       iel     = k;
//     }
//   }

//   return ( _MMG5_minQualCheck(iel,minqual,_MMG3D_ALPHAD) );
// }

std::array<double, 3> get_pos(const TetVertex &p) {
    return {{ p.posf[0], p.posf[1], p.posf[2] }};
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

bool isMeshQualityOk(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed)
{
    double minqual = 2./_MMG3D_ALPHAD;
    for (size_t e = 0; e < tets.size(); ++e) {
        if (tet_is_removed[e]) { continue; }
        auto v1 = get_pos(verts[tets[e][0]]);
        auto v2 = get_pos(verts[tets[e][1]]);
        auto v3 = get_pos(verts[tets[e][2]]);
        auto v4 = get_pos(verts[tets[e][3]]);
        double qual = mmg_caltet_iso_4pt(v1.data(), v2.data(), v3.data(), v4.data());
        if (qual < minqual) {
            minqual = qual;
        }
    }

    double minqualOnAlpha = minqual/_MMG3D_ALPHAD;
    if ( minqualOnAlpha < _MMG5_NULKAL ) {
        return false;
    }
    if ( minqualOnAlpha < _MMG5_EPSOK ) {
        // logger().warn("Worst element quality is very bad, mmg will complain: {}", minqual);
        return false;
    }
    return true;
}

bool isMeshQualityOk(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T) {
    double minqual = 2./_MMG3D_ALPHAD;
    for (int e = 0; e < T.rows(); ++e) {
        Eigen::RowVector3d v1 = V.row(T(e,0));
        Eigen::RowVector3d v2 = V.row(T(e,1));
        Eigen::RowVector3d v3 = V.row(T(e,2));
        Eigen::RowVector3d v4 = V.row(T(e,3));
        double qual = mmg_caltet_iso_4pt(v1.data(), v2.data(), v3.data(), v4.data());
        if (qual < minqual) {
            minqual = qual;
        }
    }

    double minqualOnAlpha = minqual/_MMG3D_ALPHAD;
    if ( minqualOnAlpha < _MMG5_NULKAL ) {
        return false;
    }
    if ( minqualOnAlpha < _MMG5_EPSOK ) {
        // logger().warn("Worst element quality is very bad, mmg will complain: {}", minqual);
        return false;
    }
    return true;
}

// -----------------------------------------------------------------------------

bool checkVolume(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed)
{
    for (size_t e = 0; e < tets.size(); ++e) {
        if (tet_is_removed[e]) { continue; }
        auto v1 = get_pos(verts[tets[e][0]]);
        auto v2 = get_pos(verts[tets[e][1]]);
        auto v3 = get_pos(verts[tets[e][2]]);
        auto v4 = get_pos(verts[tets[e][3]]);
        double vol = mmg_det4pt(v1.data(), v2.data(), v3.data(), v4.data());
        if ( fabs(vol) <= _MMG5_EPSD2 ) {
            // logger().warn("Tet ({},{},{},{}) has zero volume {}", T(e,0), T(e,1), T(e,2), T(e,3), vol);
            return false;
        }
    }
    return true;
}

bool checkVolume(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T) {
    for (int e = 0; e < T.rows(); ++e) {
        Eigen::RowVector3d v1 = V.row(T(e,0));
        Eigen::RowVector3d v2 = V.row(T(e,1));
        Eigen::RowVector3d v3 = V.row(T(e,2));
        Eigen::RowVector3d v4 = V.row(T(e,3));
        double vol = mmg_det4pt(v1.data(), v2.data(), v3.data(), v4.data());
        if ( fabs(vol) <= _MMG5_EPSD2 ) {
            // logger().warn("Tet ({},{},{},{}) has zero volume {}", T(e,0), T(e,1), T(e,2), T(e,3), vol);
            return false;
        }
    }
    return true;
}

// -----------------------------------------------------------------------------

bool hasNoSlivers(const std::vector<TetVertex> &verts,
    const std::vector<std::array<int, 4>> &tets,
    const std::vector<bool> &tet_is_removed,
    double angle_thres)
{
    if (angle_thres == 0.0) { return true; }

    double max_cos = std::cos(angle_thres * M_PI / 180.0);
    double min_cos = -max_cos;

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXd theta, cos_theta;
    extractVolumeMesh(verts, tets, tet_is_removed, V, T);
    igl::dihedral_angles(V, T, theta, cos_theta);

    double c0 = cos_theta.minCoeff();
    double c1 = cos_theta.maxCoeff();
    logger().warn("min_max_dihedral_angles: {} {}", c0, c1);
    return !(c0 < min_cos || c1 > max_cos);
}

} // namespace tetwild

