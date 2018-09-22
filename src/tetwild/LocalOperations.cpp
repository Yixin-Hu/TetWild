// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 5/6/17.
//

#include <tetwild/LocalOperations.h>
#include <tetwild/Common.h>
#include <tetwild/Args.h>
#include <tetwild/Logger.h>
#include <tetwild/DistanceQuery.h>
#include <pymesh/MshSaver.h>
#include <igl/svd3x3.h>
#include <igl/Timer.h>
//#include <igl/face_areas.h>
//#include <igl/dihedral_angles.h>

#define FULL_LOG false
#define CHECK_ENVELOP true

namespace tetwild {

double LocalOperations::comformalAMIPSEnergy_new(const double * T) {
    double helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    double helper_1 = helper_0[2];
    double helper_2 = helper_0[11];
    double helper_3 = helper_0[0];
    double helper_4 = helper_0[3];
    double helper_5 = helper_0[9];
    double helper_6 = 0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
    double helper_7 = helper_0[1];
    double helper_8 = helper_0[4];
    double helper_9 = helper_0[7];
    double helper_10 = helper_0[10];
    double helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 + 0.408248290463863 * helper_8 -
                       1.22474487139159 * helper_9;
    double helper_12 = 0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
    double helper_13 = helper_0[6];
    double helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 + 0.408248290463863 * helper_4 +
                       0.408248290463863 * helper_5;
    double helper_15 = helper_0[5];
    double helper_16 = helper_0[8];
    double helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 - 1.22474487139159 * helper_16 +
                       0.408248290463863 * helper_2;
    double helper_18 = 0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
    double helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
    double helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
    double helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
    return -(helper_1 * (-1.5 * helper_1 + 0.5 * helper_2 + helper_21) +
             helper_10 * (-1.5 * helper_10 + helper_20 + 0.5 * helper_7) +
             helper_13 * (-1.5 * helper_13 + 0.5 * helper_3 + 0.5 * helper_4 + 0.5 * helper_5) +
             helper_15 * (0.5 * helper_1 - 1.5 * helper_15 + 0.5 * helper_16 + 0.5 * helper_2) +
             helper_16 * (0.5 * helper_1 + 0.5 * helper_15 - 1.5 * helper_16 + 0.5 * helper_2) +
             helper_2 * (0.5 * helper_1 - 1.5 * helper_2 + helper_21) +
             helper_3 * (helper_19 - 1.5 * helper_3 + 0.5 * helper_5) +
             helper_4 * (0.5 * helper_13 + 0.5 * helper_3 - 1.5 * helper_4 + 0.5 * helper_5) +
             helper_5 * (helper_19 + 0.5 * helper_3 - 1.5 * helper_5) +
             helper_7 * (0.5 * helper_10 + helper_20 - 1.5 * helper_7) +
             helper_8 * (0.5 * helper_10 + 0.5 * helper_7 - 1.5 * helper_8 + 0.5 * helper_9) +
             helper_9 * (0.5 * helper_10 + 0.5 * helper_7 + 0.5 * helper_8 - 1.5 * helper_9)) *
           pow(pow((helper_1 - helper_2) * (helper_11 * helper_6 - helper_12 * helper_14) -
                   (-helper_10 + helper_7) * (-helper_14 * helper_18 + helper_17 * helper_6) +
                   (helper_3 - helper_5) * (-helper_11 * helper_18 + helper_12 * helper_17), 2), -0.333333333333333);
}

void LocalOperations::comformalAMIPSJacobian_new(const double * T, double *result_0) {
    double helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    double helper_1 = helper_0[1];
    double helper_2 = helper_0[10];
    double helper_3 = helper_1 - helper_2;
    double helper_4 = helper_0[0];
    double helper_5 = helper_0[3];
    double helper_6 = helper_0[9];
    double helper_7 = 0.577350269189626*helper_4 - 1.15470053837925*helper_5 + 0.577350269189626*helper_6;
    double helper_8 = helper_0[2];
    double helper_9 = 0.408248290463863*helper_8;
    double helper_10 = helper_0[5];
    double helper_11 = 0.408248290463863*helper_10;
    double helper_12 = helper_0[8];
    double helper_13 = 1.22474487139159*helper_12;
    double helper_14 = helper_0[11];
    double helper_15 = 0.408248290463863*helper_14;
    double helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
    double helper_17 = 0.577350269189626*helper_8;
    double helper_18 = 1.15470053837925*helper_10;
    double helper_19 = 0.577350269189626*helper_14;
    double helper_20 = helper_17 - helper_18 + helper_19;
    double helper_21 = helper_0[6];
    double helper_22 = -1.22474487139159*helper_21 + 0.408248290463863*helper_4 + 0.408248290463863*helper_5 + 0.408248290463863*helper_6;
    double helper_23 = helper_16*helper_7 - helper_20*helper_22;
    double helper_24 = -helper_14 + helper_8;
    double helper_25 = 0.408248290463863*helper_1;
    double helper_26 = helper_0[4];
    double helper_27 = 0.408248290463863*helper_26;
    double helper_28 = helper_0[7];
    double helper_29 = 1.22474487139159*helper_28;
    double helper_30 = 0.408248290463863*helper_2;
    double helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
    double helper_32 = helper_31*helper_7;
    double helper_33 = 0.577350269189626*helper_1;
    double helper_34 = 1.15470053837925*helper_26;
    double helper_35 = 0.577350269189626*helper_2;
    double helper_36 = helper_33 - helper_34 + helper_35;
    double helper_37 = helper_22*helper_36;
    double helper_38 = helper_4 - helper_6;
    double helper_39 = helper_23*helper_3 - helper_24*(helper_32 - helper_37) - helper_38*(helper_16*helper_36 - helper_20*helper_31);
    double helper_40 = pow(pow(helper_39, 2), -0.333333333333333);
    double helper_41 = 0.707106781186548*helper_10 - 0.707106781186548*helper_12;
    double helper_42 = 0.707106781186548*helper_26 - 0.707106781186548*helper_28;
    double helper_43 = 0.5*helper_21 + 0.5*helper_5;
    double helper_44 = 0.5*helper_26 + 0.5*helper_28;
    double helper_45 = 0.5*helper_10 + 0.5*helper_12;
    double helper_46 = 0.666666666666667*(helper_1*(-1.5*helper_1 + 0.5*helper_2 + helper_44) + helper_10*(-1.5*helper_10 + 0.5*helper_12 + 0.5*helper_14 + 0.5*helper_8) + helper_12*(0.5*helper_10 - 1.5*helper_12 + 0.5*helper_14 + 0.5*helper_8) + helper_14*(-1.5*helper_14 + helper_45 + 0.5*helper_8) + helper_2*(0.5*helper_1 - 1.5*helper_2 + helper_44) + helper_21*(-1.5*helper_21 + 0.5*helper_4 + 0.5*helper_5 + 0.5*helper_6) + helper_26*(0.5*helper_1 + 0.5*helper_2 - 1.5*helper_26 + 0.5*helper_28) + helper_28*(0.5*helper_1 + 0.5*helper_2 + 0.5*helper_26 - 1.5*helper_28) + helper_4*(-1.5*helper_4 + helper_43 + 0.5*helper_6) + helper_5*(0.5*helper_21 + 0.5*helper_4 - 1.5*helper_5 + 0.5*helper_6) + helper_6*(0.5*helper_4 + helper_43 - 1.5*helper_6) + helper_8*(0.5*helper_14 + helper_45 - 1.5*helper_8))/helper_39;
    double helper_47 = -0.707106781186548*helper_21 + 0.707106781186548*helper_5;
    result_0[0] = -helper_40*(1.0*helper_21 - 3.0*helper_4 + helper_46*(helper_41*(-helper_1 + helper_2) - helper_42*(helper_14 - helper_8) - (-helper_17 + helper_18 - helper_19)*(-helper_25 - helper_27 + helper_29 - helper_30) + (-helper_33 + helper_34 - helper_35)*(-helper_11 + helper_13 - helper_15 - helper_9)) + 1.0*helper_5 + 1.0*helper_6);
    result_0[1] = helper_40*(3.0*helper_1 - 1.0*helper_2 - 1.0*helper_26 - 1.0*helper_28 + helper_46*(helper_23 + helper_24*helper_47 - helper_38*helper_41));
    result_0[2] = helper_40*(-1.0*helper_10 - 1.0*helper_12 - 1.0*helper_14 + helper_46*(-helper_3*helper_47 - helper_32 + helper_37 + helper_38*helper_42) + 3.0*helper_8);
}

void LocalOperations::comformalAMIPSHessian_new(const double * T, double *result_0){
    double helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    double helper_1 = helper_0[2];
    double helper_2 = helper_0[11];
    double helper_3 = helper_1 - helper_2;
    double helper_4 = helper_0[0];
    double helper_5 = 0.577350269189626*helper_4;
    double helper_6 = helper_0[3];
    double helper_7 = 1.15470053837925*helper_6;
    double helper_8 = helper_0[9];
    double helper_9 = 0.577350269189626*helper_8;
    double helper_10 = helper_5 - helper_7 + helper_9;
    double helper_11 = helper_0[1];
    double helper_12 = 0.408248290463863*helper_11;
    double helper_13 = helper_0[4];
    double helper_14 = 0.408248290463863*helper_13;
    double helper_15 = helper_0[7];
    double helper_16 = 1.22474487139159*helper_15;
    double helper_17 = helper_0[10];
    double helper_18 = 0.408248290463863*helper_17;
    double helper_19 = helper_12 + helper_14 - helper_16 + helper_18;
    double helper_20 = helper_10*helper_19;
    double helper_21 = 0.577350269189626*helper_11;
    double helper_22 = 1.15470053837925*helper_13;
    double helper_23 = 0.577350269189626*helper_17;
    double helper_24 = helper_21 - helper_22 + helper_23;
    double helper_25 = 0.408248290463863*helper_4;
    double helper_26 = 0.408248290463863*helper_6;
    double helper_27 = helper_0[6];
    double helper_28 = 1.22474487139159*helper_27;
    double helper_29 = 0.408248290463863*helper_8;
    double helper_30 = helper_25 + helper_26 - helper_28 + helper_29;
    double helper_31 = helper_24*helper_30;
    double helper_32 = helper_3*(helper_20 - helper_31);
    double helper_33 = helper_4 - helper_8;
    double helper_34 = 0.408248290463863*helper_1;
    double helper_35 = helper_0[5];
    double helper_36 = 0.408248290463863*helper_35;
    double helper_37 = helper_0[8];
    double helper_38 = 1.22474487139159*helper_37;
    double helper_39 = 0.408248290463863*helper_2;
    double helper_40 = helper_34 + helper_36 - helper_38 + helper_39;
    double helper_41 = helper_24*helper_40;
    double helper_42 = 0.577350269189626*helper_1;
    double helper_43 = 1.15470053837925*helper_35;
    double helper_44 = 0.577350269189626*helper_2;
    double helper_45 = helper_42 - helper_43 + helper_44;
    double helper_46 = helper_19*helper_45;
    double helper_47 = helper_41 - helper_46;
    double helper_48 = helper_33*helper_47;
    double helper_49 = helper_11 - helper_17;
    double helper_50 = helper_10*helper_40;
    double helper_51 = helper_30*helper_45;
    double helper_52 = helper_50 - helper_51;
    double helper_53 = helper_49*helper_52;
    double helper_54 = helper_32 + helper_48 - helper_53;
    double helper_55 = pow(helper_54, 2);
    double helper_56 = pow(helper_55, -0.333333333333333);
    double helper_57 = 1.0*helper_27 - 3.0*helper_4 + 1.0*helper_6 + 1.0*helper_8;
    double helper_58 = 0.707106781186548*helper_13;
    double helper_59 = 0.707106781186548*helper_15;
    double helper_60 = helper_58 - helper_59;
    double helper_61 = helper_3*helper_60;
    double helper_62 = 0.707106781186548*helper_35 - 0.707106781186548*helper_37;
    double helper_63 = helper_49*helper_62;
    double helper_64 = helper_47 + helper_61 - helper_63;
    double helper_65 = 1.33333333333333/helper_54;
    double helper_66 = 1.0/helper_55;
    double helper_67 = 0.5*helper_27 + 0.5*helper_6;
    double helper_68 = -1.5*helper_4 + helper_67 + 0.5*helper_8;
    double helper_69 = 0.5*helper_4 + helper_67 - 1.5*helper_8;
    double helper_70 = -1.5*helper_27 + 0.5*helper_4 + 0.5*helper_6 + 0.5*helper_8;
    double helper_71 = 0.5*helper_27 + 0.5*helper_4 - 1.5*helper_6 + 0.5*helper_8;
    double helper_72 = 0.5*helper_13 + 0.5*helper_15;
    double helper_73 = -1.5*helper_11 + 0.5*helper_17 + helper_72;
    double helper_74 = 0.5*helper_11 - 1.5*helper_17 + helper_72;
    double helper_75 = 0.5*helper_11 + 0.5*helper_13 - 1.5*helper_15 + 0.5*helper_17;
    double helper_76 = 0.5*helper_11 - 1.5*helper_13 + 0.5*helper_15 + 0.5*helper_17;
    double helper_77 = 0.5*helper_35 + 0.5*helper_37;
    double helper_78 = -1.5*helper_1 + 0.5*helper_2 + helper_77;
    double helper_79 = 0.5*helper_1 - 1.5*helper_2 + helper_77;
    double helper_80 = 0.5*helper_1 + 0.5*helper_2 + 0.5*helper_35 - 1.5*helper_37;
    double helper_81 = 0.5*helper_1 + 0.5*helper_2 - 1.5*helper_35 + 0.5*helper_37;
    double helper_82 = helper_1*helper_78 + helper_11*helper_73 + helper_13*helper_76 + helper_15*helper_75 + helper_17*helper_74 + helper_2*helper_79 + helper_27*helper_70 + helper_35*helper_81 + helper_37*helper_80 + helper_4*helper_68 + helper_6*helper_71 + helper_69*helper_8;
    double helper_83 = 0.444444444444444*helper_66*helper_82;
    double helper_84 = helper_66*helper_82;
    double helper_85 = -helper_32 - helper_48 + helper_53;
    double helper_86 = 1.0/helper_85;
    double helper_87 = helper_86*pow(pow(helper_85, 2), -0.333333333333333);
    double helper_88 = 0.707106781186548*helper_6;
    double helper_89 = 0.707106781186548*helper_27;
    double helper_90 = helper_88 - helper_89;
    double helper_91 = 0.666666666666667*helper_10*helper_40 + 0.666666666666667*helper_3*helper_90 - 0.666666666666667*helper_30*helper_45 - 0.666666666666667*helper_33*helper_62;
    double helper_92 = -3.0*helper_11 + 1.0*helper_13 + 1.0*helper_15 + 1.0*helper_17;
    double helper_93 = -helper_11 + helper_17;
    double helper_94 = -helper_1 + helper_2;
    double helper_95 = -helper_21 + helper_22 - helper_23;
    double helper_96 = -helper_34 - helper_36 + helper_38 - helper_39;
    double helper_97 = -helper_42 + helper_43 - helper_44;
    double helper_98 = -helper_12 - helper_14 + helper_16 - helper_18;
    double helper_99 = -0.666666666666667*helper_60*helper_94 + 0.666666666666667*helper_62*helper_93 + 0.666666666666667*helper_95*helper_96 - 0.666666666666667*helper_97*helper_98;
    double helper_100 = helper_3*helper_90;
    double helper_101 = helper_33*helper_62;
    double helper_102 = helper_100 - helper_101 + helper_52;
    double helper_103 = -helper_60*helper_94 + helper_62*helper_93 + helper_95*helper_96 - helper_97*helper_98;
    double helper_104 = 0.444444444444444*helper_102*helper_103*helper_82*helper_86 + helper_57*helper_91 - helper_92*helper_99;
    double helper_105 = 1.85037170770859e-17*helper_1*helper_78 + 1.85037170770859e-17*helper_11*helper_73 + 1.85037170770859e-17*helper_13*helper_76 + 1.85037170770859e-17*helper_15*helper_75 + 1.85037170770859e-17*helper_17*helper_74 + 1.85037170770859e-17*helper_2*helper_79 + 1.85037170770859e-17*helper_27*helper_70 + 1.85037170770859e-17*helper_35*helper_81 + 1.85037170770859e-17*helper_37*helper_80 + 1.85037170770859e-17*helper_4*helper_68 + 1.85037170770859e-17*helper_6*helper_71 + 1.85037170770859e-17*helper_69*helper_8;
    double helper_106 = helper_64*helper_82*helper_86;
    double helper_107 = -0.666666666666667*helper_10*helper_19 + 0.666666666666667*helper_24*helper_30 + 0.666666666666667*helper_33*helper_60 - 0.666666666666667*helper_49*helper_90;
    double helper_108 = -3.0*helper_1 + 1.0*helper_2 + 1.0*helper_35 + 1.0*helper_37;
    double helper_109 = -helper_20 + helper_31 + helper_33*helper_60 - helper_49*helper_90;
    double helper_110 = 0.444444444444444*helper_109*helper_82*helper_86;
    double helper_111 = helper_103*helper_110 + helper_107*helper_57 - helper_108*helper_99;
    double helper_112 = -helper_4 + helper_8;
    double helper_113 = -helper_88 + helper_89;
    double helper_114 = -helper_5 + helper_7 - helper_9;
    double helper_115 = -helper_25 - helper_26 + helper_28 - helper_29;
    double helper_116 = helper_82*helper_86*(helper_112*helper_62 + helper_113*helper_94 + helper_114*helper_96 - helper_115*helper_97);
    double helper_117 = -helper_100 + helper_101 - helper_50 + helper_51;
    double helper_118 = -helper_102*helper_110 + helper_107*helper_92 + helper_108*helper_91;
    double helper_119 = helper_82*helper_86*(helper_112*(-helper_58 + helper_59) - helper_113*helper_93 - helper_114*helper_98 + helper_115*helper_95);
    result_0[0] = helper_56*(helper_57*helper_64*helper_65 - pow(helper_64, 2)*helper_83 + 0.666666666666667*helper_64*helper_84*(-helper_41 + helper_46 - helper_61 + helper_63) + 3.0);
    result_0[1] = helper_87*(helper_104 - helper_105*helper_35 + helper_106*helper_91);
    result_0[2] = helper_87*(helper_106*helper_107 + helper_111);
    result_0[3] = helper_87*(helper_104 + helper_116*helper_99);
    result_0[4] = helper_56*(-pow(helper_117, 2)*helper_83 + helper_117*helper_65*helper_92 + helper_117*helper_84*helper_91 + 3.0);
    result_0[5] = helper_87*(-helper_105*helper_6 - helper_107*helper_116 + helper_118);
    result_0[6] = helper_87*(-helper_105*helper_13 + helper_111 + helper_119*helper_99);
    result_0[7] = helper_87*(helper_118 - helper_119*helper_91);
    result_0[8] = helper_56*(-helper_108*helper_109*helper_65 - 1.11111111111111*pow(helper_109, 2)*helper_84 + 3.0);
}

void LocalOperations::check() {
    ///check correctness
    int n_size=0;
    int d_size=0;
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;

        for (auto it = tet_vertices[i].conn_tets.begin(); it != tet_vertices[i].conn_tets.end(); it++) {
            if (t_is_removed[*it])
                logger().debug("t {} is removed!", *it);
            auto jt=std::find(tets[*it].begin(), tets[*it].end(), i);
            if(jt==tets[*it].end()){
                logger().debug("t {} is not a conn_tet for v {}", *it, i);
            }
        }

        if(tet_vertices[i].conn_tets.size()==0) {
            logger().debug("empty conn_tets: v {}", i);
            assert(tet_vertices[i].conn_tets.size()>0);
        }

//        for(int j=0;j<3;j++) {
//            int tmp_n_size = CGAL::exact(tet_vertices[i].pos[j]).numerator().bit_size();
//            int tmp_d_size = CGAL::exact(tet_vertices[i].pos[j]).denominator().bit_size();
//            if(tmp_n_size>n_size)
//                n_size=tmp_n_size;
//            if(tmp_d_size>d_size)
//                d_size=tmp_d_size;
//        }
    }

    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++) {
            bool is_found = false;
            for (auto it = tet_vertices[tets[i][j]].conn_tets.begin();
                 it != tet_vertices[tets[i][j]].conn_tets.end(); it++) {
                if (*it == i) {
                    is_found = true;
                }
                if (t_is_removed[*it])
                    logger().debug("tet {} is removed!", *it);
            }
            if (!is_found) {
                logger().debug("{} {} {} {}", tets[i][0], tets[i][1], tets[i][2], tets[i][3]);
                logger().debug("tet {} should be conn to v {}", i, tets[i][j]);
            }
        }
    }


}

void LocalOperations::outputInfo(int op_type, double time, bool is_log) {
    logger().debug("outputing info");
    //update min/max dihedral angle infos
    for (int i = 0; i < tets.size(); i++) {
        if (!t_is_removed[i])
            calTetQuality_AD(tets[i], tet_qualities[i]);
    }

    if(args.is_quiet)
        return;

    //some tmp checks for experiments
//    for (int i = 0; i < tets.size(); i++) {
//        if (t_is_removed[i])
//            continue;
//        CGAL::Orientation ori = CGAL::orientation(tet_vertices[tets[i][0]].pos,
//                                                  tet_vertices[tets[i][1]].pos,
//                                                  tet_vertices[tets[i][2]].pos,
//                                                  tet_vertices[tets[i][3]].pos);
//        if (ori != CGAL::POSITIVE) {
//            logger().debug("outputInfo(): tet flipped!!");
//            pausee();
//        }
//
//        //check quality
////        TetQuality tq;
////        calTetQuality_AMIPS(tets[i], tq);
////        if (tq.slim_energy != tet_qualities[i].slim_energy) {
////            logger().debug("quality incorrect!");
////            logger().debug("{} {}", tq.slim_energy, tet_qualities[i].slim_energy);
////            pausee();
////        }
//    }

//    int cnt = 0;
//    for (int i = 0; i < tet_vertices.size(); i++) {
//        if (v_is_removed[i])
//            continue;
//
//        for (int j = 0; j < 3; j++) {
//            if (std::isnan(tet_vertices[i].posf[j])) {
//                logger().debug("v {} is nan", i);
//                pausee();
//            }
//        }
//
//        if (tet_vertices[i].is_on_bbox && tet_vertices[i].is_on_surface) {
//            logger().debug("ERROR: tet_vertices[i].is_on_bbox && tet_vertices[i].is_on_surface");
//            pausee();
//        }
//
//        if (tet_vertices[i].is_on_bbox)
//            cnt++;
//    }
//    logger().debug("on bbox = {}", cnt);

//    std::vector<std::array<int, 4>> fs;
//    for (int i = 0; i < tets.size(); i++) {
//        if (t_is_removed[i])
//            continue;
//        for (int j = 0; j < 4; j++) {
//            if (is_surface_fs[i][j] != state.NOT_SURFACE) {
//                std::array<int, 3> f = {tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4]};
//                std::sort(f.begin(), f.end());
//                fs.push_back(std::array<int, 4>({{f[0], f[1], f[2], is_surface_fs[i][j]}}));
//            }
//        }
//    }
//    if (fs.size() % 2 != 0) {
//        logger().debug("fs.size()%2!=0");
//    }
//    std::sort(fs.begin(), fs.end());
//    for (int i = 0; i < fs.size() - 1; i += 2) {
//        if (fs[i][0] == fs[i + 1][0] && fs[i][1] == fs[i + 1][1] && fs[i][2] == fs[i + 1][2] &&
//            fs[i][3] + fs[i + 1][3] == 0);//good
//        else {
//            logger().debug("{}", i);
//            logger().debug("hehehehe");
//            for (int j = 0; j < 4; j++)
//                logger().debug("{}{}{}{}#vertices outside of envelop = {}", fs[i][j], " ";
//            for (int j = 0; j < 4; j++)
//                std::cout, fs[i + 1][j], " ";
//            pausee();
//        }
//    }

    //check envelop
//    if(op_type != MeshRecord::OpType::OP_OPT_INIT) {
//        cnt = 0;
//        for (int i = 0; i < tet_vertices.size(); i++) {
//            if (!v_is_removed[i] && tet_vertices[i].is_on_surface) {
//                double dis = geo_sf_tree.squared_distance(
//                        GEO::vec3(tet_vertices[i].posf[0], tet_vertices[i].posf[1], tet_vertices[i].posf[2]));
//                if (dis > state.eps_2)
//                    cnt++;
//            }
//        }
//        std::cout, cnt);
//    }

    int cnt = 0;
    int r_cnt = 0;
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (!v_is_removed[i]) {
            cnt++;
            if (tet_vertices[i].is_rounded) {
//                if (tet_vertices[i].pos[0] != tet_vertices[i].posf[0]
//                    || tet_vertices[i].pos[1] != tet_vertices[i].posf[1]
//                    || tet_vertices[i].pos[2] != tet_vertices[i].posf[2]) {
//                    logger().debug("tet_vertices[i].pos!=tet_vertices[i].posf");
//                    logger().debug("{}{}{}{}{}", tet_vertices[i].pos[0] - tet_vertices[i].posf[0], " "
//, tet_vertices[i].pos[1] - tet_vertices[i].posf[1], " "
//, tet_vertices[i].pos[2] - tet_vertices[i].posf[2]);
//
//                }
                r_cnt++;
            }
//            else {
//                if (CGAL::to_double(tet_vertices[i].pos[0]) != tet_vertices[i].posf[0]
//                    || CGAL::to_double(tet_vertices[i].pos[1]) != tet_vertices[i].posf[1]
//                    || CGAL::to_double(tet_vertices[i].pos[2]) != tet_vertices[i].posf[2]) {
//                    logger().debug("CGAL::to_double(tet_vertices[i].pos)!=tet_vertices[i].posf");
//                    logger().debug("{}{}{}{}{}", CGAL::to_double(tet_vertices[i].pos[0]) - tet_vertices[i].posf[0], " "
//, CGAL::to_double(tet_vertices[i].pos[1]) - tet_vertices[i].posf[1], " "
//, CGAL::to_double(tet_vertices[i].pos[2]) - tet_vertices[i].posf[2]);
//                }
//            }
        }
    }


    logger().debug("# vertices = {}({}) {}(r)", cnt, tet_vertices.size(), r_cnt);

    cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (!t_is_removed[i])
            cnt++;
    }
    logger().debug("# tets = {}({})", cnt, tets.size());
    logger().debug("# total operations = {}", counter);
    logger().debug("# accepted operations = {}", suc_counter);


    double min = 10, max = 0;
    double min_avg = 0, max_avg = 0;
    double max_slim_energy = 0, avg_slim_energy = 0;
    std::array<double, 6> cmp_cnt = {{0, 0, 0, 0, 0, 0}};
    cnt = 0;

    for (int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        if (isTetLocked_ui(i))
            continue;

        cnt++;
        if (tet_qualities[i].min_d_angle < min)
            min = tet_qualities[i].min_d_angle;
        if (tet_qualities[i].max_d_angle > max)
            max = tet_qualities[i].max_d_angle;
        if (tet_qualities[i].slim_energy > max_slim_energy)
            max_slim_energy = tet_qualities[i].slim_energy;
        min_avg += tet_qualities[i].min_d_angle;
        max_avg += tet_qualities[i].max_d_angle;
        avg_slim_energy += tet_qualities[i].slim_energy;

        for (int j = 0; j < 3; j++) {
            if (tet_qualities[i].min_d_angle < cmp_d_angles[j])
                cmp_cnt[j]++;
        }
        for (int j = 0; j < 3; j++) {
            if (tet_qualities[i].max_d_angle > cmp_d_angles[j + 3])
                cmp_cnt[j + 3]++;
        }
    }

    logger().debug("min_d_angle = {}, max_d_angle = {}, max_slim_energy = {}", min, max, max_slim_energy);
    logger().debug("avg_min_d_angle = {}, avg_max_d_angle = {}, avg_slim_energy = {}", min_avg / cnt, max_avg / cnt, avg_slim_energy / cnt);
    logger().debug("min_d_angle: <6 {};   <12 {};  <18 {}", cmp_cnt[0] / cnt, cmp_cnt[1] / cnt, cmp_cnt[2] / cnt);
    logger().debug("max_d_angle: >174 {}; >168 {}; >162 {}", cmp_cnt[5] / cnt, cmp_cnt[4] / cnt, cmp_cnt[3] / cnt);

    if(is_log) {
        addRecord(MeshRecord(op_type, time, std::count(v_is_removed.begin(), v_is_removed.end(), false), cnt,
                             min, min_avg / cnt, max, max_avg / cnt, max_slim_energy, avg_slim_energy / cnt), args, state);
    }
}

bool LocalOperations::isTetFlip(const std::array<int, 4>& t) {
    CGAL::Orientation ori;
    bool is_rounded = true;
    for (int j = 0; j < 4; j++)
        if (!tet_vertices[t[j]].is_rounded) {
            is_rounded = false;
            break;
        }
    if (is_rounded)
        ori = CGAL::orientation(tet_vertices[t[0]].posf, tet_vertices[t[1]].posf, tet_vertices[t[2]].posf,
                                tet_vertices[t[3]].posf);
    else
        ori = CGAL::orientation(tet_vertices[t[0]].pos, tet_vertices[t[1]].pos, tet_vertices[t[2]].pos,
                                tet_vertices[t[3]].pos);

    if (ori != CGAL::POSITIVE)
        return true;
    return false;
}

bool LocalOperations::isTetFlip(int t_id){
    return isTetFlip(tets[t_id]);
}

bool LocalOperations::isFlip(const std::vector<std::array<int, 4>>& new_tets) {
    ////check orientation
    for (int i = 0; i < new_tets.size(); i++) {
//        CGAL::Orientation ori = CGAL::orientation(tet_vertices[new_tets[i][0]].pos,
//                                                  tet_vertices[new_tets[i][1]].pos,
//                                                  tet_vertices[new_tets[i][2]].pos,
//                                                  tet_vertices[new_tets[i][3]].pos);
//        if (ori != CGAL::POSITIVE)
//            return true;
        if(isTetFlip(new_tets[i]))
            return true;
    }

    return false;
}

void LocalOperations::getCheckQuality(const std::vector<TetQuality>& tet_qs, TetQuality& tq) {
    double slim_sum = 0, slim_max = 0;
    for (int i = 0; i < tet_qs.size(); i++) {
        if (state.use_energy_max) {
            if (tet_qs[i].slim_energy > slim_max)
                slim_max = tet_qs[i].slim_energy;
        } else
            slim_sum += tet_qs[i].slim_energy * tet_qs[i].volume;
    }
    if (state.use_energy_max)
        tq.slim_energy = slim_max;
    else
        tq.slim_energy = slim_sum;
}

void LocalOperations::getCheckQuality(const std::vector<int>& t_ids, TetQuality& tq){
    double slim_sum = 0, slim_max = 0;
    for (int i = 0; i < t_ids.size(); i++) {
        if (state.use_energy_max) {
            if (tet_qualities[t_ids[i]].slim_energy > slim_max)
                slim_max = tet_qualities[t_ids[i]].slim_energy;
        } else
            slim_sum += tet_qualities[t_ids[i]].slim_energy * tet_qualities[t_ids[i]].volume;
    }
    if (state.use_energy_max)
        tq.slim_energy = slim_max;
    else
        tq.slim_energy = slim_sum;
}

void LocalOperations::getAvgMaxEnergy(double& avg_tq, double& max_tq) {
    avg_tq = 0;
    max_tq = 0;
    int cnt = 0;
    for (unsigned int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        if(isTetLocked_ui(i))
            continue;
        if (tet_qualities[i].slim_energy > max_tq)
            max_tq = tet_qualities[i].slim_energy;
        avg_tq += tet_qualities[i].slim_energy;
        cnt++;
    }
    avg_tq /= cnt;
    if(std::isinf(avg_tq))
        avg_tq = state.MAX_ENERGY;
}

double LocalOperations::getMaxEnergy(){
    double max_tq = 0;
    for (unsigned int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        if(isTetLocked_ui(i))
            continue;
        if (tet_qualities[i].slim_energy > max_tq)
            max_tq = tet_qualities[i].slim_energy;
    }
    return max_tq;
}

double LocalOperations::getSecondMaxEnergy(double max_energy){
    double max_tq = 0;
    for (unsigned int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        if(tet_qualities[i].slim_energy == state.MAX_ENERGY)
            continue;
        if(isTetLocked_ui(i))
            continue;
        if (tet_qualities[i].slim_energy > max_tq)
            max_tq = tet_qualities[i].slim_energy;
    }
    return max_tq;
}

double LocalOperations::getFilterEnergy(bool& is_clean_up) {
    std::array<int, 11> buckets;
    for (int i = 0; i < 11; i++)
        buckets[i] = 0;
    for (unsigned int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        if (tet_qualities[i].slim_energy > args.filter_energy_thres - 1 + 1e10)
            buckets[10]++;
        else {
            for (int j = 0; j < 10; j++) {
                if (tet_qualities[i].slim_energy > args.filter_energy_thres - 1 + pow(10, j)
                    && tet_qualities[i].slim_energy <= args.filter_energy_thres - 1 + pow(10, j + 1)) {
                    buckets[j]++;
                    break;
                }
            }
        }
    }

    std::array<int, 10> tmps1;
    std::array<int, 10> tmps2;
    for (int i = 0; i < 10; i++) {
        tmps1[i] = std::accumulate(buckets.begin(), buckets.begin() + i + 1, 0);
        tmps2[i] = std::accumulate(buckets.begin() + i + 1, buckets.end(), 0);
    }

    if(tmps1[0]>=tmps2[0]) {
        is_clean_up = (tmps2[5] > 0);
        return 8;
    }
    if(tmps1[8]<=tmps2[8])
        return 1e11;

    for (int i = 0; i < 8; i++) {
        if (tmps1[i] < tmps2[i] && tmps1[i + 1] > tmps2[i + 1]){
            return args.filter_energy_thres - 1 + 5 * pow(10, i+1);
        }
    }

    return 8;//would never be execuate, it's fine
}

void LocalOperations::calTetQualities(const std::vector<std::array<int, 4>>& new_tets, std::vector<TetQuality>& tet_qs,
                                      bool all_measure) {
    tet_qs.resize(new_tets.size());
#ifdef TETWILD_WITH_ISPC
    int n = new_tets.size();

    static thread_local std::vector<double> T0;
    static thread_local std::vector<double> T1;
    static thread_local std::vector<double> T2;
    static thread_local std::vector<double> T3;
    static thread_local std::vector<double> T4;
    static thread_local std::vector<double> T5;
    static thread_local std::vector<double> T6;
    static thread_local std::vector<double> T7;
    static thread_local std::vector<double> T8;
    static thread_local std::vector<double> T9;
    static thread_local std::vector<double> T10;
    static thread_local std::vector<double> T11;
    static thread_local std::vector<double> energy;

    if (T0.empty()) {
        // logger().trace("Initial ISPC allocation: n = {}", n);
    } else if (T0.size() != n) {
        // logger().trace("ISPC reallocation: n = {}", n);
    }

    T0.resize(n);
    T1.resize(n);
    T2.resize(n);
    T3.resize(n);
    T4.resize(n);
    T5.resize(n);
    T6.resize(n);
    T7.resize(n);
    T8.resize(n);
    T9.resize(n);
    T10.resize(n);
    T11.resize(n);
    energy.resize(n);

    for (int i = 0; i < n; i++) {
        T0[i] = tet_vertices[new_tets[i][0]].posf[0];
        T1[i] = tet_vertices[new_tets[i][0]].posf[1];
        T2[i] = tet_vertices[new_tets[i][0]].posf[2];
        T3[i] = tet_vertices[new_tets[i][1]].posf[0];
        T4[i] = tet_vertices[new_tets[i][1]].posf[1];
        T5[i] = tet_vertices[new_tets[i][1]].posf[2];
        T6[i] = tet_vertices[new_tets[i][2]].posf[0];
        T7[i] = tet_vertices[new_tets[i][2]].posf[1];
        T8[i] = tet_vertices[new_tets[i][2]].posf[2];
        T9[i] = tet_vertices[new_tets[i][3]].posf[0];
        T10[i] = tet_vertices[new_tets[i][3]].posf[1];
        T11[i] = tet_vertices[new_tets[i][3]].posf[2];
    }

    ispc::energy_ispc(T0.data(), T1.data(), T2.data(), T3.data(), T4.data(),
        T5.data(), T6.data(), T7.data(), T8.data(),
        T9.data(), T10.data(), T11.data(), energy.data(), n);

    for (int i = 0; i < new_tets.size(); i++) {
        CGAL::Orientation ori = CGAL::orientation(tet_vertices[new_tets[i][0]].posf,
                                                  tet_vertices[new_tets[i][1]].posf,
                                                  tet_vertices[new_tets[i][2]].posf,
                                                  tet_vertices[new_tets[i][3]].posf);
        if (ori != CGAL::POSITIVE) {
            tet_qs[i].slim_energy = state.MAX_ENERGY;
            continue;
        } else
            tet_qs[i].slim_energy = energy[i];

        if (std::isinf(energy[i]) || std::isnan(energy[i]))
            tet_qs[i].slim_energy = state.MAX_ENERGY;
    }
#else
    for (int i = 0; i < new_tets.size(); i++) {
        calTetQuality_AMIPS(new_tets[i], tet_qs[i]);
    }
#endif
}

double LocalOperations::calEdgeLength(const std::array<int, 2>& v_ids){
    return CGAL::squared_distance(tet_vertices[v_ids[0]].posf, tet_vertices[v_ids[1]].posf);
}

double LocalOperations::calEdgeLength(int v1_id,int v2_id, bool is_over_refine) {
    return CGAL::squared_distance(tet_vertices[v1_id].posf, tet_vertices[v2_id].posf);
}

void LocalOperations::calTetQuality_AD(const std::array<int, 4>& tet, TetQuality& t_quality) {
    std::array<Vector_3f, 4> nv;
    std::array<double, 4> nv_length;
    std::array<double, 4> heights;
    for (int i = 0; i < 4; i++) {
        Plane_3f pln(tet_vertices[tet[(i + 1) % 4]].posf,
                    tet_vertices[tet[(i + 2) % 4]].posf,
                    tet_vertices[tet[(i + 3) % 4]].posf);
        if(pln.is_degenerate()){
            t_quality.min_d_angle = 0;
            t_quality.max_d_angle = M_PI;
            return;
        }
        Point_3f tmp_p = pln.projection(tet_vertices[tet[i]].posf);
        if(tmp_p == tet_vertices[tet[i]].posf){
            t_quality.min_d_angle = 0;
            t_quality.max_d_angle = M_PI;
            return;
        }
        nv[i] = tet_vertices[tet[i]].posf - tmp_p;
        heights[i] = CGAL::squared_distance(tet_vertices[tet[i]].posf, tmp_p);

//        if(std::isnan(heights[i])){//because pln is degenerate
//            logger().debug("{}", tet_vertices[tet[i]].posf);
//            logger().debug("{}", tet_vertices[tet[(i + 1) % 4]].posf);
//            logger().debug("{}", tet_vertices[tet[(i + 2) % 4]].posf);
//            logger().debug("{}", tet_vertices[tet[(i + 3) % 4]].posf);
//            logger().debug("{}", pln.is_degenerate());
//
//            logger().debug("{}", tmp_p);
//            logger().debug("{}", nv[i]);
//            logger().debug("{}", heights[i]);
//            pausee();
//        }

        //re-scale
        std::array<double, 3> tmp_nv = {{CGAL::abs(nv[i][0]), CGAL::abs(nv[i][1]), CGAL::abs(nv[i][2])}};
        auto tmp = std::max_element(tmp_nv.begin(), tmp_nv.end());
        if(*tmp == 0 || heights[i] == 0){
            t_quality.min_d_angle = 0;
            t_quality.max_d_angle = M_PI;
//            t_quality.asp_ratio_2 = state.MAX_ENERGY;
            return;
        } else if (*tmp < 1e-5) {
            nv[i] = Vector_3f(nv[i][0] / *tmp, nv[i][1] / *tmp, nv[i][2] / *tmp);
            nv_length[i] = sqrt(heights[i] / ((*tmp) * (*tmp)));
        } else {
            nv_length[i] = sqrt(heights[i]);
        }
    }

    std::vector<std::array<int, 2>> opp_edges;
    for (int i = 0; i < 3; i++) {
        opp_edges.push_back(std::array<int, 2>({{0, i + 1}}));
        opp_edges.push_back(std::array<int, 2>({{i + 1, (i + 1) % 3 + 1}}));
    }

    ////compute dihedral angles
    std::array<double, 6> dihedral_angles;
    for (int i = 0; i < (int) opp_edges.size(); i++) {
        double dihedral_angle = -nv[opp_edges[i][0]] * nv[opp_edges[i][1]] /
                                 (nv_length[opp_edges[i][0]] * nv_length[opp_edges[i][1]]);
        if (dihedral_angle > 1)
            dihedral_angles[i] = 0;
        else if (dihedral_angle < -1)
            dihedral_angles[i] = M_PI;
        else
            dihedral_angles[i] = std::acos(dihedral_angle);
    }
//    std::sort(dihedral_angles.begin(), dihedral_angles.end());
    auto it=std::minmax_element(dihedral_angles.begin(), dihedral_angles.end());
    t_quality.min_d_angle = *(it.first);
    t_quality.max_d_angle = *(it.second);

//    std::sort(heights.begin(), heights.end());
//    auto h = std::min_element(heights.begin(), heights.end());
//    t_quality.asp_ratio_2 = max_e_l / *h;
}

void LocalOperations::calTetQuality_AMIPS(const std::array<int, 4>& tet, TetQuality& t_quality) {
    if (energy_type == state.ENERGY_AMIPS) {
        CGAL::Orientation ori = CGAL::orientation(tet_vertices[tet[0]].posf,
                                                  tet_vertices[tet[1]].posf,
                                                  tet_vertices[tet[2]].posf,
                                                  tet_vertices[tet[3]].posf);
        if (ori != CGAL::POSITIVE) {//degenerate in floats
            t_quality.slim_energy = state.MAX_ENERGY;
        } else {
            std::array<double, 12> T;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 3; j++) {
                    T[i*3+j] = tet_vertices[tet[i]].posf[j];
                }
            }
            t_quality.slim_energy = comformalAMIPSEnergy_new(T.data());
            if (std::isinf(t_quality.slim_energy) || std::isnan(t_quality.slim_energy))
                t_quality.slim_energy = state.MAX_ENERGY;
        }
    }
    if(std::isinf(t_quality.slim_energy) || std::isnan(t_quality.slim_energy) || t_quality.slim_energy <= 0)
        t_quality.slim_energy = state.MAX_ENERGY;
}

bool LocalOperations::isEdgeOnSurface(int v1_id, int v2_id) {
    if (!tet_vertices[v1_id].is_on_surface || !tet_vertices[v2_id].is_on_surface)
        return false;

    std::vector<int> t_ids;
    setIntersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, t_ids);
    assert(t_ids.size()!=0);
    return isEdgeOnSurface(v1_id, v2_id, t_ids);
}

bool LocalOperations::isEdgeOnBbox(int v1_id, int v2_id){
    if(!tet_vertices[v1_id].is_on_bbox || !tet_vertices[v2_id].is_on_bbox)
        return false;

    std::vector<int> t_ids;
    setIntersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, t_ids);
    return isEdgeOnBbox(v1_id, v2_id, t_ids);
}

bool LocalOperations::isEdgeOnSurface(int v1_id, int v2_id, const std::vector<int>& t_ids){
    for (int i = 0; i < t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[t_ids[i]][j] != v1_id && tets[t_ids[i]][j] != v2_id) {
                if (is_surface_fs[t_ids[i]][j]!=state.NOT_SURFACE)
                    return true;
            }
        }
    }
    return false;
}

bool LocalOperations::isEdgeOnBbox(int v1_id, int v2_id, const std::vector<int>& t_ids){
    std::unordered_set<int> v_ids;
    for (int i = 0; i < t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[t_ids[i]][j] != v1_id && tets[t_ids[i]][j] != v2_id) {
                v_ids.insert(tets[t_ids[i]][j]);
            }
        }
    }
    if(v_ids.size()!=t_ids.size())
        return true;
    return false;
}

bool LocalOperations::isEdgeOnBoundary(int v1_id, int v2_id) {
//    if (boundary_points.size() == 0)//if it's a closed mesh, then there cannot be any boundary edges.
//        return false;

    if(state.is_mesh_closed)
        return false;

    if (!tet_vertices[v1_id].is_on_boundary || !tet_vertices[v2_id].is_on_boundary)
        return false;

//    return true;

    int cnt = 0;
    for (int t_id: tet_vertices[v1_id].conn_tets) {
        std::array<int, 4> opp_js;
        int ii = 0;
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] == v1_id || tets[t_id][j] == v2_id)
                continue;
            opp_js[ii++] = j;
        }
        if (ii == 2) {
            if (is_surface_fs[t_id][opp_js[0]] != state.NOT_SURFACE)
                cnt++;
            if (is_surface_fs[t_id][opp_js[1]] != state.NOT_SURFACE)
                cnt++;
            if (cnt > 2)
                return false;
        }
    }
    if (cnt == 2) //is boundary edge
        return true;

    return false;
}

bool LocalOperations::isFaceOutEnvelop(const Triangle_3f& tri) {
#if CHECK_ENVELOP
    if(state.use_sampling){
        return isFaceOutEnvelop_sampling(tri);
    }
    return true;
#else
    return false;
#endif

//    EnvelopSide side = getUpperLowerBounds(tri);
//    if (side == EnvelopSide::OUTSIDE)
//        return true;
//    else if (side == EnvelopSide::INSIDE)
//        return false;
//    else {
//        int depth = 0;
//        int cnt = 0, sub_cnt = 0;
//        std::queue<Triangle_3f> tris_queue;
//        tris_queue.push(tri);
//        cnt++;
//        while (!tris_queue.empty()) {
////            logger().debug("depth = {}{}cnt = {}{}sub_cnt = {}", depth, ", "
////, cnt, ", "
////, sub_cnt);
//            if (depth == 6)
//                return true;
////                return false;
//            Triangle_3f &cur_tri = tris_queue.front();
//
//            //subdivide
//            std::array<Point_3f, 3> mps;
//            for (int j = 0; j < 3; j++)
//                mps[j] = CGAL::midpoint(cur_tri[j], cur_tri[(j + 1) % 3]);
//
//            std::array<Triangle_3f, 4> tris;
//            for (int j = 0; j < 3; j++)
//                tris[j] = Triangle_3f(cur_tri[j], mps[(j + 1) % 3], mps[(j - 1 + 3) % 3]);
//            tris[3] = Triangle_3f(mps[0], mps[1], mps[2]);
//
//            for (int j = 0; j < 4; j++) {
////                logger().debug("{}{}depth = {}{}cnt = {}{}sub_cnt = {}", j, ": ";
//                side = getUpperLowerBounds(tris[j]);
//                if (side == EnvelopSide::OUTSIDE)
//                    return true;
//                else if (side == EnvelopSide::UNCERTAIN) {
//                    tris_queue.push(tris[j]);
//                    sub_cnt++;
//                }
//            }
//
//            tris_queue.pop();
//            cnt--;
////            std::cout, depth, ", "
////, cnt, ", "
////, sub_cnt);
//            if (cnt == 0) {
//                cnt = sub_cnt;
//                sub_cnt = 0;
//                depth++;
//            }
////            pausee();
//        }
//        return false;
//    }
}

bool LocalOperations::isPointOutEnvelop(const Point_3f& p) {
#if CHECK_ENVELOP
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    if (geo_sf_tree.squared_distance(geo_p) > state.eps_2)
        return true;

    return false;
#else
    return false;
#endif
}

bool LocalOperations::isFaceOutEnvelop_sampling(const Triangle_3f& tri) {
#if CHECK_ENVELOP
    if (tri.is_degenerate())
        return false;

#if TIMING_BREAKDOWN
    igl_timer0.start();
#endif
    std::array<GEO::vec3, 3> vs = {{GEO::vec3(tri[0][0], tri[0][1], tri[0][2]),
                                    GEO::vec3(tri[1][0], tri[1][1], tri[1][2]),
                                    GEO::vec3(tri[2][0], tri[2][1], tri[2][2])}};
    static thread_local std::vector<GEO::vec3> ps;
    ps.clear();
    sampleTriangle(vs, ps, state.sampling_dist);
#if TIMING_BREAKDOWN
    breakdown_timing0[id_sampling] += igl_timer0.getElapsedTime();
#endif

    size_t num_queries = 0;
    size_t num_samples = ps.size();

    //decide in/out
#if TIMING_BREAKDOWN
    igl_timer0.start();
#endif

    GEO::vec3 current_point = ps[0];
    GEO::vec3 nearest_point;
    double sq_dist = std::numeric_limits<double>::max();
    GEO::index_t prev_facet = GEO::NO_FACET;
    int cnt = 0;
    const unsigned int ps_size = ps.size();
    for (unsigned int i = ps_size / 2; i < ps.size(); i = (i + 1) % ps_size) {//check from the middle
        GEO::vec3 &current_point = ps[i];
        if (prev_facet != GEO::NO_FACET) {
            get_point_facet_nearest_point(geo_sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
        }
        if (sq_dist > state.eps_2) {
            geo_sf_tree.facet_in_envelope_with_hint(
                current_point, state.eps_2, prev_facet, nearest_point, sq_dist);
        }
        ++num_queries;
        if (sq_dist > state.eps_2) {
#if TIMING_BREAKDOWN
            breakdown_timing0[id_aabb] += igl_timer0.getElapsedTime();
#endif
            logger().trace("num_queries {} / {}", num_queries, num_samples);
            return true;
        }
        cnt++;
        if (cnt >= ps_size)
            break;
    }

#if TIMING_BREAKDOWN
    breakdown_timing0[id_aabb] += igl_timer0.getElapsedTime();
#endif

    logger().trace("num_queries {} / {}", num_queries, num_samples);
    return false;
#else
    return false;
#endif
}

bool LocalOperations::isPointOutBoundaryEnvelop(const Point_3f& p) {
#if CHECK_ENVELOP
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    if (geo_b_tree.squared_distance(geo_p) > state.eps_2) {
        return true;
    }
    return false;
#else
    return false;
#endif
}

bool LocalOperations::isBoundarySlide(int v1_id, int v2_id, Point_3f& old_pf){
    return false;

#if CHECK_ENVELOP
    if(state.is_mesh_closed)
        return false;

    std::unordered_set<int> n_v_ids;
    for(int t_id:tet_vertices[v1_id].conn_tets){
        for(int j=0;j<4;j++)
            if(tets[t_id][j]!=v1_id && tets[t_id][j]!=v2_id && tet_vertices[tets[t_id][j]].is_on_boundary)
                n_v_ids.insert(tets[t_id][j]);
    }
    if(n_v_ids.size()==0)
        return false;

#if TIMING_BREAKDOWN
    igl_timer0.start();
#endif
    static thread_local std::vector<GEO::vec3> b_points;
    static thread_local std::vector<GEO::vec3> ps;
    b_points.clear();
    for(int v_id:n_v_ids) {
        if (!isEdgeOnBoundary(v1_id, v_id))
            continue;
        //sample the edge (v1, v) and push the sampling points into vector
        GEO::vec3 p1(tet_vertices[v1_id].posf[0], tet_vertices[v1_id].posf[1], tet_vertices[v1_id].posf[2]);
        GEO::vec3 p2(tet_vertices[v_id].posf[0], tet_vertices[v_id].posf[1], tet_vertices[v_id].posf[2]);
        b_points.push_back(p1);
        b_points.push_back(p2);
        int n = GEO::distance(p1, p2) / state.sampling_dist + 1;
        if (n == 1)
            continue;
        b_points.reserve(b_points.size() + n + 1);
        for (int k = 1; k <= n - 1; k++)
            b_points.push_back(p1 * ((double) k / (double) n) + p2 * ((double) (n - k) / (double) n));
    }

    //sampling faces
    if(v2_id>=0 && tet_vertices[v2_id].is_on_boundary) {
        std::vector<int> n12_t_ids;
        setIntersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, n12_t_ids);
        std::unordered_set<int> n12_v_ids;
        for (int t_id:n12_t_ids) {
            for (int j = 0; j < 4; j++)
                if (tets[t_id][j] != v1_id && tets[t_id][j] != v2_id && tet_vertices[tets[t_id][j]].is_on_boundary)
                    n12_v_ids.insert(tets[t_id][j]);
        }
        bool is_12_on_boundary = false;
        if(n12_v_ids.size()!=0) {
            is_12_on_boundary = isEdgeOnBoundary(v1_id, v2_id);
        }
        for(int v_id:n12_v_ids) {
            if (!isEdgeOnBoundary(v1_id, v_id) || !isEdgeOnBoundary(v2_id, v_id))
                continue;
            if (!is_12_on_boundary) {
                GEO::vec3 p1(tet_vertices[v1_id].posf[0], tet_vertices[v1_id].posf[1], tet_vertices[v1_id].posf[2]);
                GEO::vec3 p2(old_pf[0], old_pf[1], old_pf[2]);
                int n = GEO::distance(p1, p2) / state.sampling_dist + 1;
                b_points.reserve(b_points.size() + n + 1);
                b_points.push_back(p1);
                for (int k = 1; k <= n - 1; k++)
                    b_points.push_back(p1 * ((double) k / (double) n) + p2 * ((double) (n - k) / (double) n));
                b_points.push_back(p2);
            } else {
                Triangle_3f tri(tet_vertices[v_id].posf, tet_vertices[v2_id].posf, old_pf);
                std::array<GEO::vec3, 3> vs = {{GEO::vec3(tri[0][0], tri[0][1], tri[0][2]),
                                                GEO::vec3(tri[1][0], tri[1][1], tri[1][2]),
                                                GEO::vec3(tri[2][0], tri[2][1], tri[2][2])}};
                ps.clear();
                sampleTriangle(vs, ps, state.sampling_dist);

//                sampleTriangle(tri, ps);//CANNOT directly push the sampling points into b_points

                b_points.reserve(b_points.size() + ps.size()); // preallocate memory
                b_points.insert(b_points.end(), ps.begin(), ps.end());
            }
        }
    }
#if TIMING_BREAKDOWN
    breakdown_timing0[id_sampling] += igl_timer0.getElapsedTime();
#endif
    if(b_points.size()==0)
        return false;

#if TIMING_BREAKDOWN
    igl_timer0.start();
#endif
    GEO::vec3 current_point = b_points[0];
    GEO::vec3 nearest_point;
    double sq_dist;
    GEO::index_t prev_facet = geo_b_tree.nearest_facet(current_point, nearest_point, sq_dist);
    int cnt = 0;
    const unsigned int b_points_size = b_points.size();
    for (unsigned int i = b_points_size / 2; ; i = (i + 1) % b_points_size) {
        GEO::vec3 &current_point = b_points[i];
        sq_dist = current_point.distance2(nearest_point);
        geo_b_tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_dist);
        double dis = current_point.distance2(nearest_point);
        if (dis > state.eps_2) {
#if TIMING_BREAKDOWN
            breakdown_timing0[id_aabb] += igl_timer0.getElapsedTime();
#endif
            return true;
        }
        cnt++;
        if (cnt >= b_points.size())
            break;
    }
#if TIMING_BREAKDOWN
    breakdown_timing0[id_aabb] += igl_timer0.getElapsedTime();
#endif

    return false;
#else
    return false;
#endif
}

bool LocalOperations::isTetOnSurface(int t_id){
    for(int i=0;i<4;i++){
        if(is_surface_fs[t_id][i]!=state.NOT_SURFACE)
            return false;
    }
    return true;
}

bool LocalOperations::isTetRounded(int t_id){
    for(int i=0;i<4;i++){
        if(!tet_vertices[tets[t_id][i]].is_rounded)
            return false;
    }
    return true;
}

void LocalOperations::getFaceConnTets(int v1_id, int v2_id, int v3_id, std::vector<int>& t_ids){
    std::vector<int> v1, v2, v3, tmp;
    v1.reserve(tet_vertices[v1_id].conn_tets.size());
    for(auto it=tet_vertices[v1_id].conn_tets.begin();it!=tet_vertices[v1_id].conn_tets.end();it++)
        v1.push_back(*it);
    v2.reserve(tet_vertices[v2_id].conn_tets.size());
    for(auto it=tet_vertices[v2_id].conn_tets.begin();it!=tet_vertices[v2_id].conn_tets.end();it++)
        v2.push_back(*it);
    v3.reserve(tet_vertices[v3_id].conn_tets.size());
    for(auto it=tet_vertices[v3_id].conn_tets.begin();it!=tet_vertices[v3_id].conn_tets.end();it++)
        v3.push_back(*it);

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::sort(v3.begin(), v3.end());

    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(tmp));
    std::set_intersection(v3.begin(), v3.end(), tmp.begin(), tmp.end(), std::back_inserter(t_ids));
}

bool LocalOperations::isIsolated(int v_id) {
    for (auto it = tet_vertices[v_id].conn_tets.begin(); it != tet_vertices[v_id].conn_tets.end(); it++) {
        for (int j = 0; j < 4; j++) {
            if (tets[*it][j] != v_id && is_surface_fs[*it][j] != state.NOT_SURFACE)
                return false;
        }
    }

    return true;
}

bool LocalOperations::isBoundaryPoint(int v_id) {
    if(state.is_mesh_closed)
        return false;
    std::unordered_set<int> n_v_ids;
    for (int t_id:tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++)
            if (tets[t_id][j] != v_id && tet_vertices[tets[t_id][j]].is_on_boundary)
                n_v_ids.insert(tets[t_id][j]);
    }
    for (int n_v_id:n_v_ids) {
        if (isEdgeOnBoundary(n_v_id, v_id))
            return true;
    }
    return false;
}

void LocalOperations::checkUnrounded() {
    bool is_output = false;
    for (unsigned int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        if (!tet_vertices[i].is_rounded) {
            is_output = true;
            break;
        }
    }
    if(!is_output)
        return;

    std::ofstream of;
    of.open(state.working_dir + "unrounded_check.txt");
    int cnt_sf = 0;
    int cnt_b = 0;
    int cnt_all = 0;
    int cnt_sf1 = 0;
    std::vector<double> diss;
    for (unsigned int i = 0; i < tet_vertices.size(); i++) {
        if (v_is_removed[i])
            continue;
        if (tet_vertices[i].is_rounded)
            continue;

        cnt_all++;

        if (tet_vertices[i].is_on_boundary)
            cnt_b++;
        if (tet_vertices[i].is_on_surface) {
            cnt_sf++;
            continue;
        }

        bool is_found = false;
        for (int t_id:tet_vertices[i].conn_tets) {
            for (int j = 0; j < 4; j++) {
                if (tets[t_id][j] == i) {
                    if (is_surface_fs[t_id][j] != state.NOT_SURFACE) {
                        cnt_sf1++;
                        is_found = true;
                    }
                    break;
                }
            }
            if (is_found)
                break;
        }
        if (is_found)
            continue;

        GEO::vec3 geo_p(tet_vertices[i].posf[0], tet_vertices[i].posf[1], tet_vertices[i].posf[2]);
        double dis = sqrt(geo_sf_tree.squared_distance(geo_p));
        diss.push_back(dis);
    }

    of << "# all = " << cnt_all << std::endl;
    of << "# surface = " << cnt_sf << std::endl;
    of << "# boundary = " << cnt_b << std::endl;
//    of<<"Is closed? "<<is_closed<<std::endl;
    of << "# connect to surface = " << cnt_sf1 << std::endl;
    of << "# non-surface = " << diss.size() << std::endl;
    for (double dis:diss) {
        of << dis << std::endl;
    }
}

bool LocalOperations::isLocked_ui(const std::array<int, 2>& e){
    return (tet_vertices[e[0]].is_locked || tet_vertices[e[1]].is_locked);
}

bool LocalOperations::isTetLocked_ui(int tid){
//    return false;

    for(int j=0;j<4;j++)
        if(tet_vertices[tets[tid][j]].is_locked)
            return true;
    return false;
}

void LocalOperations::outputSurfaceColormap(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in, double old_eps) {
    state.sampling_dist /= 2;

    Eigen::VectorXd eps_dis(F_in.rows());
    for (int f_id = 0; f_id < F_in.rows(); f_id++) {
        //sample triangles except one-ring of v1v2
        std::vector<GEO::vec3> ps;
        std::array<GEO::vec3, 3> vs = {{
                GEO::vec3(V_in(F_in(f_id, 0), 0), V_in(F_in(f_id, 0), 1), V_in(F_in(f_id, 0), 2)),
                GEO::vec3(V_in(F_in(f_id, 1), 0), V_in(F_in(f_id, 1), 1), V_in(F_in(f_id, 1), 2)),
                GEO::vec3(V_in(F_in(f_id, 2), 0), V_in(F_in(f_id, 2), 1), V_in(F_in(f_id, 2), 2))}};
//        sampleTriangle(vs, ps);
        double sqrt3_2 = sqrt(3) / 2;

        std::array<double, 3> ls;
        for (int i = 0; i < 3; i++) {
            ls[i] = GEO::length2(vs[i] - vs[(i + 1) % 3]);
        }
        auto min_max = std::minmax_element(ls.begin(), ls.end());
        int min_i = min_max.first - ls.begin();
        int max_i = min_max.second - ls.begin();
        double N = sqrt(ls[max_i]) / state.sampling_dist;
        if (N <= 1) {
            for (int i = 0; i < 3; i++)
                ps.push_back(vs[i]);
            return;
        }
        if (N == int(N))
            N -= 1;

        GEO::vec3 v0 = vs[max_i];
        GEO::vec3 v1 = vs[(max_i + 1) % 3];
        GEO::vec3 v2 = vs[(max_i + 2) % 3];

        GEO::vec3 n_v0v1 = GEO::normalize(v1 - v0);
        for (int n = 0; n <= N; n++) {
            ps.push_back(v0 + n_v0v1 * state.sampling_dist * n);
        }
        ps.push_back(v1);

        double h = GEO::distance(GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
        int M = h / (sqrt3_2 * state.sampling_dist);
        if (M < 1) {
            ps.push_back(v2);
            return;
        }

        GEO::vec3 n_v0v2 = GEO::normalize(v2 - v0);
        GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
        double tan_v0, tan_v1, sin_v0, sin_v1;
        sin_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / (GEO::distance(v0, v2) * GEO::distance(v0, v1));
        tan_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / GEO::dot((v2 - v0), (v1 - v0));
        tan_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / GEO::dot((v2 - v1), (v0 - v1));
        sin_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / (GEO::distance(v1, v2) * GEO::distance(v0, v1));

        for (int m = 1; m <= M; m++) {
            int n = sqrt3_2 / tan_v0 * m + 0.5;
            int n1 = sqrt3_2 / tan_v0 * m;
            if (m % 2 == 0 && n == n1) {
                n += 1;
            }
            GEO::vec3 v0_m = v0 + m * sqrt3_2 * state.sampling_dist / sin_v0 * n_v0v2;
            GEO::vec3 v1_m = v1 + m * sqrt3_2 * state.sampling_dist / sin_v1 * n_v1v2;

            double delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * state.sampling_dist;
            GEO::vec3 v = v0_m + delta_d * n_v0v1;
            int N1 = GEO::distance(v, v1_m) / state.sampling_dist;
            ps.push_back(v0_m);
            for (int i = 0; i <= N1; i++) {
                ps.push_back(v + i * n_v0v1 * state.sampling_dist);
            }
            ps.push_back(v1_m);
        }

        ps.push_back(v2);

//        std::array<double, 3> ls;
//        for (int i = 0; i < 3; i++) {
//            ls[i] = GEO::length(vs[i] - vs[(i + 1) % 3]);
//        }
//        auto min_max = std::minmax_element(ls.begin(), ls.end());
//        int min_i = min_max.first - ls.begin();
//        int max_i = min_max.second - ls.begin();
//
//        double n = int(ls[max_i] / state.sampling_dist + 1);
//        ps.reserve(2 * n);
//        for (int j = 1; j < n; j++) {
//            ps.push_back(j / n * vs[(min_i + 2) % 3] + (n - j) / n * vs[min_i]);
//            ps.push_back(j / n * vs[(min_i + 2) % 3] + (n - j) / n * vs[(min_i + 1) % 3]);
//        }
//        if (ls[min_i] > state.sampling_dist) {
//            int ps_size = ps.size();
//            for (int i = 0; i < ps_size; i += 2) {
//                double m = int(GEO::length(ps[i] - ps[i + 1]) / state.sampling_dist + 1);
//                if (m == 0)
//                    break;
//                for (int j = 1; j < m; j++)
//                    ps.push_back(j / m * ps[i] + (m - j) / m * ps[i + 1]);
//            }
//        }
//        for (int i = 0; i < 3; i++)
//            ps.push_back(vs[i]);

        //check sampling points
        GEO::vec3 current_point = ps[0];
        GEO::vec3 nearest_point;
        double sq_dist;
        GEO::index_t prev_facet = geo_sf_tree.nearest_facet(current_point, nearest_point, sq_dist);

        double max_dis = 0;
        int cnt = 0;
        for (const GEO::vec3 &current_point:ps) {
            sq_dist = current_point.distance2(nearest_point);
            geo_sf_tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_dist);
            double dis = current_point.distance2(nearest_point);
            if (dis > max_dis)
                max_dis = dis;
        }
        eps_dis(f_id) = sqrt(max_dis / (old_eps * old_eps));
    }

    Eigen::VectorXd V_vec(V_in.rows() * 3);
    Eigen::VectorXi F_vec(F_in.rows() * 3);
    for (unsigned int i = 0; i < V_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            V_vec(i * 3 + j) = V_in(i, j);
    }
    for (unsigned int i = 0; i < F_in.rows(); i++) {
        for (int j = 0; j < 3; j++)
            F_vec(i * 3 + j) = F_in(i, j);
    }

    PyMesh::MshSaver mshSaver(state.working_dir + args.postfix + "_sf" + std::to_string(mid_id++) + ".msh");
    mshSaver.save_mesh(V_vec, F_vec, 3, mshSaver.TRI);
    mshSaver.save_elem_scalar_field("distance to surface", eps_dis);

    state.sampling_dist *= 2;
}

} // namespace tetwild
