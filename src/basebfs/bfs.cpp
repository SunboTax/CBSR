#include "basebfs/bfs.h"

#include <bits/types/time_t.h>

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "partition.h"

void biBFS(BiGraph &bg, vector<CCR::queryInfo> &queryInfo, vector<int> &queryRes) {
    LOG("Base BFS query start");
    int i = 0;
    for (auto const &info : queryInfo) {
        auto res = biBFSWorker(bg, info);
        queryRes[i++] = res;
    }
}

bool biBFSWorker(BiGraph &bg, CCR::queryInfo info) {
    auto u = info.u;
    auto w = info.w;
    auto start = info.start;
    auto end = info.end;
    if (start > end) {
        return false;
    }
    if (u == w) {
        return true;
    }

    std::queue<pair<int, time_t>> Q;
    Q.push(make_pair(u, start));
    // 到达某个上部点的最小时间
    vector<int> minT(bg.adj_matrix_u.size(), std::numeric_limits<int>::max());
    vector<int> maxMinT(bg.adj_matrix_l.size(), std::numeric_limits<int>::max());
    minT[u] = 0;

    std::queue<pair<int, time_t>> P;
    P.push(make_pair(w, end));
    vector<int> maxT(bg.adj_matrix_u.size(), 0);
    vector<int> minMaxT(bg.adj_matrix_l.size(), 0);
    maxT[w] = std::numeric_limits<int>::max();
    set<int> s;
    bool firstQ = true;
    while (!Q.empty() and !P.empty()) {
        if (Q.size() <= P.size()) {
            int len = Q.size();
            for (int k = 0; k < len; k++) {
                // while (!Q.empty()) {
                pair<int, time_t> value = Q.front();
                Q.pop();
                int u_ = value.first;
                time_t t_e = value.second;
                // cout<<u_<<" "<<t_e<<"\n";
                auto flag = firstQ;
                firstQ = false;
                for (int i = 0; i < bg.adj_matrix_u[u_].size(); i++) {
                    auto v = bg.adj_matrix_u[u_][i];
                    auto t1 = bg.timeSection_u[u_][i].first;
                    auto t2 = bg.timeSection_u[u_][i].second;
                    if (t1 < t_e) {
                        if (flag and t2 > t_e) {
                        } else {
                            continue;
                        }
                    }
                    vector<long> candi{u_, v, t1, t2};
                    auto hash_res = hash_vector(candi);
                    if (s.count(hash_res)) {
                        continue;
                    }
                    s.insert(hash_res);

                    long t = 0;
                    if (t1 >= maxMinT[v]) {
                        continue;
                    }

                    for (int j = 0; j < bg.adj_matrix_l[v].size(); j++) {
                        auto w_ = bg.adj_matrix_l[v][j];
                        auto t1_ = bg.timeSection_l[v][j].first;
                        auto t2_ = bg.timeSection_l[v][j].second;
                        if (t2_ >= minT[w_]) {
                            continue;
                        }

                        // 判断e和e_是否重叠
                        if (max(t1, t1_) < min(t2, t2_) and t2_ <= end) {
                            // cout<<t1_<<" -- "<<t2_<<"\n";
                            if (w_ == w) {
                                return true;
                            }
                            if (minT[w_] <= maxT[w_]) {
                                return true;
                            }
                            Q.push(make_pair(w_, t2_));
                            minT[w_] = t2_;
                        }

                        if (minT[w_] > t) {
                            t = minT[w_];
                        }
                    }

                    maxMinT[v] = t;
                }
            }
        } else {
            // int len = P.size();
            // for (int k = 0; k < len; k++) {
            while (!P.empty()) {
                pair<int, time_t> value = P.front();
                P.pop();
                int u_ = value.first;
                time_t t_s = value.second;
                for (int i = 0; i < bg.adj_matrix_u[u_].size(); i++) {
                    auto v = bg.adj_matrix_u[u_][i];
                    auto t1 = bg.timeSection_u[u_][i].first;
                    auto t2 = bg.timeSection_u[u_][i].second;
                    if (t2 > t_s) {
                        continue;
                    }
                    vector<long> candi{u_, v, t1, t2};
                    auto hash_res = hash_vector(candi);
                    if (s.count(hash_res)) {
                        continue;
                    }
                    s.insert(hash_res);

                    long t = std::numeric_limits<long>::max();
                    if (t2 <= minMaxT[v]) {
                        continue;
                    }

                    for (int j = 0; j < bg.adj_matrix_l[v].size(); j++) {
                        auto w_ = bg.adj_matrix_l[v][j];
                        auto t1_ = bg.timeSection_l[v][j].first;
                        auto t2_ = bg.timeSection_l[v][j].second;
                        if (t1_ <= maxT[w_]) {
                            continue;
                        }

                        // if overlap
                        if (max(t1, t1_) < min(t2, t2_) and t1_ >= start) {
                            if (w_ == u) {
                                return true;
                            }
                            if (minT[w_] <= maxT[w_]) {
                                return true;
                            }
                            P.push(make_pair(w_, t1_));
                            maxT[w_] = t1_;
                        }

                        if (maxT[w_] < t) {
                            t = maxT[w_];
                        }
                    }
                    minMaxT[v] = t;
                }
            }
        }
    }
    return false;
}