#include "topchain/topChain.h"

#include <assert.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>

#include "partition.h"
namespace topchain {
using namespace std;

bool cmp1(DEG d1, DEG d2)
{
    return d1.d > d2.d;
}

void queryTransform(BiGraph& bg, ptree::Graph& g, int u, int w, time_t start,
    time_t end, int& src, int& dst)
{
    int v1 = -1;
    int v2 = -1;
    // cout << __FILE__ << ":" << __LINE__ << "\n";
    // 找到u在start时刻之后的第一个位置v1
    for (int i = 0; i < bg.adj_matrix_u[u].size(); i++) {
        if (bg.timeSection_u[u][i].second <= start) {
            // cout << __FILE__ << ":" << __LINE__ << "\n";
            if (i == 0) {
                break;
            }
            v1 = bg.adj_matrix_u[u][i - 1];

            int k = i - 2;
            while (bg.timeTable_u[v1].empty() && k >= 0) {
                v1 = bg.adj_matrix_u[u][k];
                k--;
            }

            if (bg.timeTable_u[v1].empty() && k == -1) {
                v1 = -1;
            }
            break;
        }

        if (i == bg.adj_matrix_u[u].size() - 1 && v1 == -1) {
            if (i == 0) {
                break;
            }
            v1 = bg.adj_matrix_u[u][i];
            int k = i - 2;
            while (bg.timeTable_u[v1].empty() && k >= 0) {
                v1 = bg.adj_matrix_u[u][k];
                k--;
            }
            if (bg.timeTable_u[v1].empty() && k == -1) {
                v1 = -1;
                // cout << __FILE__ << ":" << __LINE__ << ":" << i << "---" << k
                // << endl;
            }
            break;
        }
    }

    // 找打w在end时刻之前的最后的位置v2
    int pos2;
    for (int i = 0; i < bg.adj_matrix_u[w].size(); i++) {
        if (bg.timeSection_u[w][i].first < end) {
            if (i == 0) {
                v2 = bg.adj_matrix_u[w][0];
                pos2 = 0;
            } else if (bg.timeSection_u[w][i - 1].first < end) {
                v2 = bg.adj_matrix_u[w][i - 1];
                pos2 = i - 1;
            } else {
                v2 = bg.adj_matrix_u[w][i];
                pos2 = i;
            }

            // v2 = bg.adj_matrix_u[w][i - 1];
            int k = i - 1;
            while (bg.timeTable_u[v2].empty() && k >= 0) {
                v2 = bg.adj_matrix_u[w][k];
                pos2 = k;
                k--;
            }
            if (bg.timeTable_u[v2].empty() && k == -1)
                v2 = -1;
            break;

            if (i == bg.adj_matrix_u[w].size() - 1
                && bg.timeSection_u[w][i].first >= end)
                break;
        }
    }

    // LOG("v1 = {}, v2 = {}", v1, v2);
    if (v1 == -1 || v2 == -1) {
        src = -1;
        dst = -1;
        return;
    }

    // for (auto t : bg.timeTable_u[v1]) cout << "[" << t.first << " " <<
    // t.second
    // << "],"; cout << endl;

    int v1_idx = -1;
    time_t mm = std::numeric_limits<time_t>::max();
    for (int i = 0; i < bg.timeTable_u[v1].size(); i++) {
        if (bg.timeTable_u[v1][i].second > start) {
            if (bg.timeTable_u[v1][i].second < mm) {
                mm = bg.timeTable_u[v1][i].second;
                v1_idx = i;
            }
        }
    }

    if (v1_idx >= 0) {
    } else {
        src = -1;
        dst = -1;
        return;
    }

    int v2_idx = -1;

    mm = 0;

    for (int i = 0; i < bg.timeTable_u[v2].size(); i++) {
        if (bg.timeTable_u[v2][i].first < end) {
            if (bg.timeTable_u[v2][i].first > mm) {
                mm = bg.timeTable_u[v2][i].first;
                v2_idx = i;
            }
        }
    }

    if (v2_idx >= 0) {
        // v2_idx = iter - bg.timeTable_u[v2].begin();
    } else {
        src = -1;
        dst = -1;
        return;
    }

    src = g.cnum(make_pair(v1, bg.timeTable_u[v1][v1_idx].first));
    // cout << __FILE__ << ":" << __LINE__ << "\n";
    dst = g.cnum(make_pair(v2, bg.timeTable_u[v2][v2_idx].first));
}

// void queryByTopChain(BiGraph &bg, ptree::Graph &graph,
//                      std::vector<CCR::queryInfo> &queryInfo, GraphT gt) {
//   int idx = -1;
//   flag = new int[n];
//   for (auto &info : queryInfo) {
//     idx++;
//     auto u = info.u;
//     auto w = info.w;
//     auto start = info.start;
//     auto end = info.end;
//     // 开始时间晚于结束时间，查询结果为假
//     if (start > end) {
//       info.reachable = false;
//       continue;
//     }
//     int src = -1;
//     int dst = -1;
//     queryTransform(bg, graph, u, w, start, end, src, dst);
//     if (src < 0 || dst < 0) {
//       info.reachable = false;
//       continue;
//     }
//     // TODO: 转成topchain需要的src与dst
//     start = min(graph[src].start, start);
//     end = max(graph[dst].start, end);
//     u = gt.getFrom(src, start, end);
//     w = gt.getTo(dst, start, end);

//     bool arrive;
//     if (u < 0 || w < 0) {
//       arrive = 0;
//     } else if (u == w)
//       arrive = 1;
//     else if (reachindex[u].Ldown >= reachindex[w].Ldown)
//       arrive = 0;
//     else if (reachindex[u].flabel1 > reachindex[w].flabel1 ||
//              reachindex[u].flabel2 > reachindex[w].flabel2)
//       arrive = 0;
//     else if (intersect(u, w))
//       arrive = 1;
//     else {
//       mark = idx + 1;
//       arrive = query(u, w);
//     }
//     info.reachable = arrive;
//   }
// }
//----------------------query----------------------------

// int main1(int argc, char* argv[]) {
//   // srand(time(0));
//   string s;
//   k = atoi(argv[1]);
//   // argv[2] for graph
//   s = argv[2];
//   string s2;
//   // argv[3] for query
//   s2 = argv[3];
//   // argv[4] for info
//   string s3 = argv[4];

//   // read
//   // timer.start();
//   readGraph(s);
//   readChainInfo(s3);
//   // timer.stop();
//   //  cout << "read graph time: " << //timer.GetRuntime() << endl;

//   // index
//   // timer.start();
//   findChain();
//   feline();  // feline label
//   handle();
//   // timer.stop();
//   //  cout << "!!!index time : " << //timer.GetRuntime()*1000 << endl;
//   double indexTime = 0;  // timer.GetRuntime() * 1000;
//   copyIndex();

//   // query
//   query(s2);
//   // cout << countQuery << endl;
//   cout << "CONSTR_TIME= " << indexTime << " ms" << endl;

//   long long count1 = 0;
//   for (int i = 0; i < n; ++i) {
//     count1 += Lin[i].size();
//     count1 += Lout[i].size();
//   }
//   // n*k*2
//   cout << "INDEX_SIZE= " << (count1 + n * 1 + n * 2) * 4.0 / 1024 / 1024
//        << " MB" << endl;
// }

} // namespace topchain
