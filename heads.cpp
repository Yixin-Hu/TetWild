// This file is part of TetWild, a software for generating tetrahedral meshes.
// 
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 3/29/17.
//

#include "heads.h"

//Eigen::MatrixXd G;
std::string g_working_dir = "";
std::string g_stat_file = "";
std::string g_postfix = "";
std::string g_output_file = "";

bool is_using_energy_max = true;

bool is_using_sampling = true;
bool is_use_project = false;
int cnt_geo_aabb = 0;

int NOT_SURFACE = 0;

double g_eps = 0;
double g_eps_2 = 0;
double g_dd = 0;
double g_ideal_l = 0;
double g_diag_l = 0;
bool g_is_close = true;

double g_eps_input = 0;
double g_eps_delta = 0;
int g_cur_stage = 1;

bool is_print_tmp = false;

Args args;
//std::vector<MeshRecord> mesh_records = std::vector<MeshRecord>();
bool is_app_csv=false;
void addRecord(const MeshRecord& record) {
    if (!args.is_output_csv)
        return;
    std::ofstream f;
    if (is_app_csv)
        f.open(g_stat_file, std::ios::app);
    else {
        f.open(g_stat_file);
        is_app_csv = true;
    }
    f << record.op << "," << record.timing << "," << record.n_v << "," << record.n_t << ","
      << record.min_min_d_angle << "," << record.avg_min_d_angle << ","
      << record.max_max_d_angle << "," << record.avg_max_d_angle << ","
      << record.max_energy << "," << record.avg_energy << "\n";
    f.close();
}

void pausee(){
    cout<<"Is pausing... (Enter '0' to exit and other characters to continue.)"<<endl;
    char c;
    cin>>c;
    if(c=='0')
        exit(0);
}

bool isHaveCommonEle(const std::unordered_set<int>& v1, const std::unordered_set<int>& v2) {
    for (auto it = v1.begin(); it != v1.end(); it++)
        if(std::find(v2.begin(), v2.end(), *it)!=v2.end())
            return true;
    return false;
}

void setIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& s) {
    std::unordered_set<int> s_tmp;
    std::vector<int> v1, v2;
    v1.reserve(s1.size());
    for (auto it = s1.begin(); it != s1.end(); it++)
        v1.push_back(*it);
    v2.reserve(s2.size());
    for (auto it = s2.begin(); it != s2.end(); it++)
        v2.push_back(*it);
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(s_tmp, s_tmp.end()));
    s = s_tmp;

//    s.clear();
//    s.reserve(std::min(s1.size(), s2.size()));
//    std::unordered_set<int> s_tmp = s2;
//    int size = s_tmp.size();
//    for(int ele:s1){
//        s_tmp.insert(ele);
//        if(s_tmp.size()>size){
//            size = s_tmp.size();
//        } else
//            s.insert(ele);
//    }
}

void setIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& s) {
//    s.clear();
//    s.reserve(std::min(s1.size(), s2.size()));
//    std::unordered_set<int> s_tmp = s2;
//    int size = s_tmp.size();
//    for(int ele:s1){
//        s_tmp.insert(ele);
//        if(s_tmp.size()>size){
//            size = s_tmp.size();
//        } else
//            s.push_back(ele);
//    }
    std::vector<int> v1, v2;
    v1.reserve(s1.size());
    for(auto it=s1.begin();it!=s1.end();it++)
        v1.push_back(*it);
    v2.reserve(s2.size());
    for(auto it=s2.begin();it!=s2.end();it++)
        v2.push_back(*it);
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(s));
}


void sampleTriangle(const std::array<GEO::vec3, 3>& vs, std::vector<GEO::vec3>& ps) {
    double sqrt3_2 = std::sqrt(3) / 2;

    std::array<double, 3> ls;
    for (int i = 0; i < 3; i++) {
        ls[i] = GEO::length2(vs[i] - vs[(i + 1) % 3]);
    }
    auto min_max = std::minmax_element(ls.begin(), ls.end());
    int min_i = min_max.first - ls.begin();
    int max_i = min_max.second - ls.begin();
    double N = sqrt(ls[max_i]) / g_dd;
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
        ps.push_back(v0 + n_v0v1 * g_dd * n);
    }
    ps.push_back(v1);

    double h = GEO::distance(GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
    int M = h / (sqrt3_2 * g_dd);
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
        GEO::vec3 v0_m = v0 + m * sqrt3_2 * g_dd / sin_v0 * n_v0v2;
        GEO::vec3 v1_m = v1 + m * sqrt3_2 * g_dd / sin_v1 * n_v1v2;
        if (GEO::distance(v0_m, v1_m) <= g_dd)
            break;

        double delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * g_dd;
        GEO::vec3 v = v0_m + delta_d * n_v0v1;
        int N1 = GEO::distance(v, v1_m) / g_dd;
//        ps.push_back(v0_m);
        for (int i = 0; i <= N1; i++) {
            ps.push_back(v + i * n_v0v1 * g_dd);
        }
//        ps.push_back(v1_m);
    }
    ps.push_back(v2);

    //sample edges
    N = sqrt(ls[(max_i + 1) % 3]) / g_dd;
    if (N > 1) {
        if (N == int(N))
            N -= 1;
        GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
        for (int n = 1; n <= N; n++) {
            ps.push_back(v1 + n_v1v2 * g_dd * n);
        }
    }

    N = sqrt(ls[(max_i + 2) % 3]) / g_dd;
    if (N > 1) {
        if (N == int(N))
            N -= 1;
        GEO::vec3 n_v2v0 = GEO::normalize(v0 - v2);
        for (int n = 1; n <= N; n++) {
            ps.push_back(v2 + n_v2v0 * g_dd * n);
        }
    }

//    cout << "ps.size = " << ps.size() << endl;
//    cout << "is output samples?" << endl;
//    int anw = 0;
//    cin >> anw;
//    if (anw != 0) {
////        if (true) {
//        Eigen::MatrixXd V_tmp(ps.size() * 3 + 3, 3);
//        Eigen::MatrixXi F_tmp(ps.size() + 1, 3);
//        for (int i = 0; i < 3; i++) {
//            for (int j = 0; j < 3; j++)
//                V_tmp(i, j) = vs[i][j];
//            F_tmp(0, i) = i;
//        }
//
//        for (int i = 0; i < ps.size(); i++) {
//            for (int k = 0; k < 3; k++) {
//                for (int j = 0; j < 3; j++)
//                    V_tmp((1 + i) * 3 + k, j) = ps[i][j];
//                F_tmp(1 + i, k) = (1 + i) * 3 + k;
//            }
//        }
//        igl::writeSTL(g_working_dir + "_sample.stl", V_tmp, F_tmp);
//    }
}