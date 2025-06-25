#include "utils.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <utility>

#include "CCRTree.h"

namespace CCR {
vector<string> split(string line, const char *symbols, int lenth) {
    vector<string> res = vector<string>();
    char *m;
    char *slice = NULL;
    char *buffer_ = (char *)malloc(sizeof(char) * lenth);
    strcpy(buffer_, line.c_str());
    m = strtok_r(buffer_, symbols, &slice);
    while (m != nullptr) {
        res.push_back(m);
        m = strtok_r(nullptr, symbols, &slice);
    }
    return res;
}

Graph transformation(BiGraph &bg, int delta) {
    Graph g;
    // flow[i][j] = <partition.id, <v,startTime>>
    vector<vector<pair<int, pair<time_t, int>>>> flow(bg.adj_matrix_u.size(), vector<pair<int, pair<time_t, int>>>());
    for (int v = 0; v < bg.adj_matrix_l.size(); v++) {
        vector<Edge> edges;
        set<time_t> times;

        for (int u = 0; u < bg.adj_matrix_l[v].size(); u++) {
            edges.emplace_back(
                Edge(bg.adj_matrix_l[v][u], bg.timeSection_l[v][u].first, bg.timeSection_l[v][u].second));
            times.insert(bg.timeSection_l[v][u].first);
            times.insert(bg.timeSection_l[v][u].second);
        }

        sort(edges.begin(), edges.end(), [](Edge &a, Edge &b) { return a.start < b.start; });

        vector<time_t> Tlist;             // 增序排列的时间节点
        unordered_map<time_t, int> Tmap;  // 时间节点映射到时间节点序号的哈希表

        for (auto t : times) {
            Tlist.push_back(t);
            Tmap[t] = Tlist.size() - 1;
        }

        vector<pair<int, int>> partition;  // 连接区间
        for (int i = 0; i < edges.size(); i++) {
            auto edge = edges[i];

            pair<int, int> tmp;
            tmp.first = Tmap[edge.start];
            tmp.second = Tmap[edge.end];
            // 当前区间不与连接区间重合或者当前还没划分区间
            if (partition.size() == 0 || partition.back().second < tmp.first) {
                // 不需要考虑最后一条边：如果最后一条边不与之前的partition相交，那么最后一条边就不会包含可达信息
                if (i < edges.size() - 1) {
                    auto edge_ = edges[i + 1];

                    // 这里只验证后一条边是否与当前边连通，如果不连通，则当前边不包含可达信息
                    if (edge.end > edge_.start) {
                        partition.emplace_back(tmp);
                        flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.first]));
                        flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.second - 1]));
                    }
                }
            } else if (partition.back().second > tmp.first) {
                flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.first]));
                flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.second - 1]));
                partition.back().second = max(partition.back().second, tmp.second);
            }
        }

        for (auto part : partition) {
            g.addVertex(v, Tlist[part.first], Tlist[part.first + 1]);
            bg.timeTable_u[v].push_back(make_pair(Tlist[part.first], Tlist[part.first + 1]));
            for (int k = part.first + 1; k < part.second; k++) {
                g.addVertex(v, Tlist[k], Tlist[k + 1]);
                bg.timeTable_u[v].push_back(make_pair(Tlist[k], Tlist[k + 1]));
                g.addEdge(make_pair(v, Tlist[k - 1]), make_pair(v, Tlist[k]));
            }
        }
    }

    // 降低查询转换过程的平均时间复杂度
    for (int i = 0; i < bg.timeTable_u.size(); i++)
        sort(bg.timeTable_u[i].begin(), bg.timeTable_u[i].end(),
             [](pair<time_t, time_t> &a, pair<time_t, time_t> &b) { return a.first < b.first; });

    // 构建不同下部点之间的连接边
    for (int u = 0; u < bg.adj_matrix_u.size(); u++) {
        sort(flow[u].begin(), flow[u].end(), [](pair<int, pair<time_t, int>> &a, pair<int, pair<time_t, int>> &b) {
            return a.second.second < b.second.second;
        });
        if (flow[u].size() < 3) continue;
        flow[u].erase(flow[u].begin());
        for (int i = 0; i < flow[u].size() / 2; i++) {
            // 单条边时间跨度约束
            if (flow[u][2 * i + 1].second.second - flow[u][2 * i].second.second <= delta &&
                !(flow[u][2 * i + 1].second.first == flow[u][2 * i].second.first &&
                  flow[u][2 * i + 1].first == flow[u][2 * i].first))
                g.addEdge(flow[u][2 * i].second, flow[u][2 * i + 1].second);
        }
    }

    return g;
}

Graph transformation(BiGraph &bg) {
    Graph g;
    // flow[i][j] = <partition.id, <v,startTime>>
    vector<vector<pair<int, pair<time_t, int>>>> flow(bg.adj_matrix_u.size(), vector<pair<int, pair<time_t, int>>>());
    for (int v = 0; v < bg.adj_matrix_l.size(); v++) {
        vector<Edge> edges;
        set<time_t> times;

        for (int u = 0; u < bg.adj_matrix_l[v].size(); u++) {
            edges.emplace_back(
                Edge(bg.adj_matrix_l[v][u], bg.timeSection_l[v][u].first, bg.timeSection_l[v][u].second));
            times.insert(bg.timeSection_l[v][u].first);
            times.insert(bg.timeSection_l[v][u].second);
        }

        sort(edges.begin(), edges.end(), [](Edge &a, Edge &b) { return a.start < b.start; });

        vector<time_t> Tlist;             // 增序排列的时间节点
        unordered_map<time_t, int> Tmap;  // 时间节点映射到时间节点序号的哈希表

        for (auto t : times) {
            Tlist.push_back(t);
            Tmap[t] = Tlist.size() - 1;
        }

        vector<pair<int, int>> partition;  // 连接区间
        for (int i = 0; i < edges.size(); i++) {
            auto edge = edges[i];

            pair<int, int> tmp;
            tmp.first = Tmap[edge.start];
            tmp.second = Tmap[edge.end];
            // 当前区间不与连接区间重合或者当前还没划分区间
            if (partition.size() == 0 || partition.back().second < tmp.first) {
                // 不需要考虑最后一条边：如果最后一条边不与之前的partition相交，那么最后一条边就不会包含可达信息
                if (i < edges.size() - 1) {
                    auto edge_ = edges[i + 1];

                    // 这里只验证后一条边是否与当前边连通，如果不连通，则当前边不包含可达信息
                    if (edge.end > edge_.start) {
                        partition.emplace_back(tmp);
                        flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.first]));
                        flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.second - 1]));
                    }
                }
            } else if (partition.back().second > tmp.first) {
                flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.first]));
                flow[edge.node].emplace_back(partition.size(), make_pair(v, Tlist[tmp.second - 1]));
                partition.back().second = max(partition.back().second, tmp.second);
            }
        }

        for (auto part : partition) {
            g.addVertex(v, Tlist[part.first], Tlist[part.first + 1]);
            bg.timeTable_u[v].push_back(make_pair(Tlist[part.first], Tlist[part.first + 1]));
            for (int k = part.first + 1; k < part.second; k++) {
                g.addVertex(v, Tlist[k], Tlist[k + 1]);
                bg.timeTable_u[v].push_back(make_pair(Tlist[k], Tlist[k + 1]));
                g.addEdge(make_pair(v, Tlist[k - 1]), make_pair(v, Tlist[k]));
            }
        }
    }

    // 降低查询转换过程的平均时间复杂度
    for (int i = 0; i < bg.timeTable_u.size(); i++)
        sort(bg.timeTable_u[i].begin(), bg.timeTable_u[i].end(),
             [](pair<time_t, time_t> &a, pair<time_t, time_t> &b) { return a.first < b.first; });

    // 构建不同下部点之间的连接边
    for (int u = 0; u < bg.adj_matrix_u.size(); u++) {
        sort(flow[u].begin(), flow[u].end(), [](pair<int, pair<time_t, int>> &a, pair<int, pair<time_t, int>> &b) {
            return a.second.second < b.second.second;
        });
        if (flow[u].size() < 3) continue;
        flow[u].erase(flow[u].begin());
        for (int i = 0; i < flow[u].size() / 2; i++) {
            // 单条边时间跨度约束
            if (flow[u][2 * i + 1].second.second - flow[u][2 * i].second.second <= 99999999 &&
                !(flow[u][2 * i + 1].second.first == flow[u][2 * i].second.first &&
                  flow[u][2 * i + 1].first == flow[u][2 * i].first))
                g.addEdge(flow[u][2 * i].second, flow[u][2 * i + 1].second);
        }
    }

    return g;
}

using Quaduple = tuple<int, time_t, time_t, int>;
struct QuadupleHasher {
    std::size_t operator()(const Quaduple &q) const {
        std::size_t h1 = (std::get<0>(q));
        std::size_t h2 = (std::get<1>(q));
        std::size_t h3 = (std::get<2>(q));
        std::size_t h4 = (std::get<3>(q));
        return ((h1 << 32) | h2) + ((h3 << 32) | h4);
    }
};

struct QuadupleEqual {
    bool operator()(const Quaduple &q1, const Quaduple &q2) const {
        return std::get<0>(q1) == std::get<0>(q2) && std::get<1>(q1) == std::get<1>(q2) &&
               std::get<2>(q1) == std::get<2>(q2) && std::get<3>(q1) == std::get<3>(q2);
    }
};

Graph newTransformation(BiGraph &bg, int delta) {
    Graph g;

    // unordered_map<Quaduple, pair<time_t, time_t>, QuadupleHasher, QuadupleEqual> flow;
    // for (int v = 0; v < bg.adj_matrix_l.size(); v++) {
    //     set<time_t> times;
    //     unordered_map<int, set<int>> events;
    //     for (int u = 0; u < bg.adj_matrix_l[v].size(); u++) {
    //         times.insert(bg.timeSection_l[v][u].first);
    //         times.insert(bg.timeSection_l[v][u].second);
    //         int k = u + 1;
    //         events[bg.timeSection_l[v][u].first].insert(k);
    //         events[bg.timeSection_l[v][u].second].insert(-k);
    //     }

    //     vector<time_t> Tlist;
    //     unordered_map<time_t, int> Tmap;

    //     for (auto t : times) {
    //         Tlist.push_back(t);
    //         Tmap[t] = Tlist.size() - 1;
    //     }

    //     for (int i = 0; i < bg.adj_matrix_l[v].size(); i++) {
    //         int u = bg.adj_matrix_l[v][i];
    //         time_t start = bg.timeSection_l[v][i].first;
    //         time_t end = bg.timeSection_l[v][i].second;
    //         time_t start_ = Tlist[Tmap[start] + 1];
    //         time_t end_ = Tlist[Tmap[end] - 1];
    //         flow[make_tuple(u, start, end, v)] = make_pair(start_, end_);
    //     }

    //     set<int> current_edges;
    //     for (int i = 0; i < Tlist.size(); i++) {
    //         auto set_size = current_edges.size();

    //         time_t t = Tlist[i];
    //         for (auto e : events[t]) {
    //             if (e >= 0) {
    //                 current_edges.insert(e);
    //             } else {
    //                 current_edges.erase(-e);
    //             }
    //         }

    //         if (current_edges.empty())
    //             continue;
    //         else {
    //             bg.timeTable_u[v].push_back(make_pair(Tlist[i], Tlist[i + 1]));
    //             if (i == 0) {
    //                 g.addVertex(v, Tlist[i], Tlist[i + 1]);
    //             } else {
    //                 g.addVertex(v, Tlist[i], Tlist[i + 1]);
    //                 if (set_size > 0) {
    //                     g.addEdge(make_tuple(v, Tlist[i - 1], Tlist[i]), make_tuple(v, Tlist[i], Tlist[i + 1]));
    //                 }
    //             }
    //         }
    //     }
    // }

    // for (int u = 0; u < bg.adj_matrix_u.size(); u++) {
    //     map<int, set<int>> events;
    //     for (int v = 0; v < bg.adj_matrix_u[u].size(); v++) {
    //         int k = v + 1;
    //         events[bg.timeSection_u[u][v].first].insert(k);
    //         events[bg.timeSection_u[u][v].second].insert(-k);
    //     }

    //     vector<int> ends;
    //     vector<int> starts;
    //     ends.reserve(50);
    //     starts.reserve(50);

    //     int FIND_END = 1;
    //     int FIND_START = 0;
    //     time_t max_end = 0;
    //     time_t min_start = INT64_MAX;

    //     auto iter = events.begin();
    //     while (iter != events.end()) {
    //         auto kv = *iter;
    //         time_t t = kv.first;
    //         set<int> &es = kv.second;
    //         if (FIND_END) {
    //             for (auto e : es) {
    //                 if (e < 0) {
    //                     int k = (-e) - 1;
    //                     ends.push_back(k);
    //                 }
    //             }
    //             if (!ends.empty()) {
    //                 max_end = max(max_end, t);
    //             }
    //             for (auto e : es) {
    //                 if (e > 0 && !ends.empty()) {
    //                     int k = e - 1;
    //                     starts.push_back(k);
    //                 }
    //             }
    //             if (!starts.empty()) {
    //                 FIND_START = 1;
    //                 FIND_END = 0;
    //                 min_start = min(min_start, t);
    //             }
    //             iter++;

    //         } else if (FIND_START) {
    //             for (auto e : es) {
    //                 if (e < 0) {
    //                     FIND_END = 1;
    //                     FIND_START = 0;
    //                     break;
    //                 }
    //             }

    //             // FIND_END=1时加边，否则继续找开始点
    //             if (FIND_END) {
    //                 if (starts.size() == 1 && ends.size() == 1) {
    //                     int v1 = bg.adj_matrix_u[u][starts[0]];
    //                     int v2 = bg.adj_matrix_u[u][ends[0]];
    //                     time_t t1 = bg.timeSection_u[u][starts[0]].first;
    //                     time_t t2 = bg.timeSection_u[u][starts[0]].second;

    //                     time_t t3 = bg.timeSection_u[u][ends[0]].first;
    //                     time_t t4 = bg.timeSection_u[u][ends[0]].second;

    //                     time_t t5 = flow[make_tuple(u, t1, t2, v1)].first;
    //                     time_t t6 = flow[make_tuple(u, t3, t4, v2)].second;

    //                     g.addEdge(make_tuple(v2, t6, t4), make_tuple(v1, t1, t5));

    //                 } else if (starts.size() == 1 && ends.size() > 1) {
    //                     int v1 = bg.adj_matrix_u[u][starts[0]];
    //                     time_t t1 = bg.timeSection_u[u][starts[0]].first;
    //                     time_t t2 = bg.timeSection_u[u][starts[0]].second;
    //                     time_t t3 = flow[make_tuple(u, t1, t2, v1)].first;

    //                     g.addVertex(v1, min_start, min_start);

    //                     g.addEdge(make_tuple(v1, min_start, min_start), make_tuple(v1, t1, t3));

    //                     for (auto e : ends) {
    //                         int v2 = bg.adj_matrix_u[u][e];
    //                         time_t t4 = bg.timeSection_u[u][e].first;
    //                         time_t t5 = bg.timeSection_u[u][e].second;
    //                         time_t t6 = flow[make_tuple(u, t4, t5, v2)].second;
    //                         g.addEdge(make_tuple(v2, t6, t5), make_tuple(v1, min_start, min_start));
    //                     }

    //                 } else if (starts.size() > 1 && ends.size() == 1) {
    //                     int v1 = bg.adj_matrix_u[u][ends[0]];
    //                     time_t t1 = bg.timeSection_u[u][ends[0]].first;
    //                     time_t t2 = bg.timeSection_u[u][ends[0]].second;
    //                     time_t t3 = flow[make_tuple(u, t1, t2, v1)].second;
    //                     int v_min = INT32_MAX;
    //                     for (auto s : starts) v_min = min(v_min, bg.adj_matrix_u[u][s]);
    //                     g.addVertex(v_min, min_start, min_start);

    //                     g.addEdge(make_tuple(v1, t3, t2), make_tuple(v_min, max_end, max_end));

    //                     for (auto s : starts) {
    //                         int v2 = bg.adj_matrix_u[u][s];
    //                         time_t t4 = bg.timeSection_u[u][s].first;
    //                         time_t t5 = bg.timeSection_u[u][s].second;
    //                         time_t t6 = flow[make_tuple(u, t4, t5, v2)].first;

    //                         g.addEdge(make_tuple(v_min, max_end, max_end), make_tuple(v2, t4, t6));
    //                     }
    //                 } else if (starts.size() > 1 && ends.size() > 1) {
    //                     int v_min = INT32_MAX;
    //                     for (auto s : starts) v_min = min(v_min, bg.adj_matrix_u[u][s]);

    //                     g.addVertex(v_min, min_start, min_start);
    //                     for (auto s : starts) {
    //                         int v1 = bg.adj_matrix_u[u][s];
    //                         time_t t1 = bg.timeSection_u[u][s].first;
    //                         time_t t2 = bg.timeSection_u[u][s].second;
    //                         time_t t3 = flow[make_tuple(u, t1, t2, v1)].first;

    //                         g.addEdge(make_tuple(v_min, min_start, min_start), make_tuple(v1, t1, t3));
    //                     }

    //                     for (auto e : ends) {
    //                         int v2 = bg.adj_matrix_u[u][e];
    //                         time_t t4 = bg.timeSection_u[u][e].first;
    //                         time_t t5 = bg.timeSection_u[u][e].second;
    //                         time_t t6 = flow[make_tuple(u, t4, t5, v2)].second;

    //                         g.addEdge(make_tuple(v2, t6, t5), make_tuple(v_min, min_start, min_start));
    //                     }
    //                 }
    //                 max_end = 0;
    //                 min_start = INT64_MAX;
    //                 starts.clear();
    //                 ends.clear();
    //                 starts.reserve(50);
    //                 ends.reserve(50);
    //             } else {
    //                 for (auto e : es) {
    //                     if (e > 0) {
    //                         int k = e - 1;
    //                         starts.push_back(k);
    //                     }
    //                 }
    //                 min_start = min(min_start, t);
    //                 iter++;
    //             }
    //         }
    //     }
    // }

    // for (int i = 0; i < bg.timeTable_u.size(); i++)
    //     sort(bg.timeTable_u[i].begin(), bg.timeTable_u[i].end(),
    //          [](pair<time_t, time_t> &a, pair<time_t, time_t> &b) { return a.first < b.first; });

    return g;
}

vector<queryInfo> readQuery(string filepath) {
    vector<queryInfo> querys;

    fstream fin(filepath, ios::in);
    int u, w;
    time_t start, end;
    while (fin >> u >> w >> start >> end) {
        queryInfo query;
        query.u = u;
        query.w = w;
        query.start = start;
        query.end = end;
        querys.push_back(query);
    }

    fin.close();
    return querys;
}

vector<queryInfo> readQuery(string filepath, string queryResFilePath, int wantReach) {
    vector<queryInfo> querys;

    fstream fin(filepath, ios::in);
    fstream fin2(queryResFilePath, ios::in);
    int u, w;
    time_t start, end;
    while (fin >> u >> w >> start >> end) {
        int reach;
        fin2 >> reach;
        if (reach == wantReach) {
            queryInfo query;
            query.u = u;
            query.w = w;
            query.start = start;
            query.end = end;
            querys.push_back(query);
        }
        if (querys.size() >= 100000) break;
    }

    fin.close();
    fin2.close();
    return querys;
}

vector<queryInfo> readQueryForTest(string filepath) {
    vector<queryInfo> querys;

    fstream fin(filepath, ios::in);
    int u, w;
    time_t start, end;
    while (fin >> u >> w) {
        queryInfo query;
        query.u = u;
        query.w = w;
        querys.push_back(query);
    }

    fin.close();
    return querys;
}
/*
void GraphInfo(Graph &g, string filename)
{
        fstream ftmp("analyse/" + filename + ".txt", ios::out);
        stringstream stmp;
        for (int u = 0; u < g.n; u++)
        {
                stmp << setw(8) << u << ": ";
                for (int v : g.adj_matrix[u])
                        stmp << setw(8) << v << " ";
                stmp << endl;
        }
        stmp << "}" << endl;
        ftmp << stmp.str();
        ftmp.close();
}

// 绘制前k条路径包含的点所在的子图（包含所有与这些点相关联的边）
void Graphvisual(Graph &g, string filename, int k)
{
        vector<vector<int>> prec_adj(g.n, vector<int>());
        for (int u = 0; u < g.n; u++)
                for (int v : g.adj_matrix[u])
                        prec_adj[v].push_back(u);

        unordered_set<int> nodeSet;

        fstream ftmp("analyse/" + filename + ".dot", ios::out);
        stringstream stmp;
        stmp << "digraph graphe {\ngraph[rankdir =\"LR\"];\nnodesep = 3;\nnode
[fontname=\"Arial\", shape = \"record\", style=filled,size = 5];\nedge
[color=black];" << endl; for (int i = 0; i < k; i++)
        {
                for (int u : g.pathMap[i])
                {
                        stmp << "\t" << u << "[label = \"id=" << u << "$|place-"
<< g.orderToClique[u].node << "|" << g.orderToClique[u].start << "|" <<
g.orderToClique[u].end << "\"];" << endl; nodeSet.insert(u);
                }
                for (int j = 0; j < g.pathMap[i].size() - 1; j++)
                        stmp << "\t" << g.pathMap[i][j] << "->" <<
g.pathMap[i][j + 1] << "[penwidth = 5];" << endl;
        }

        for (int i = 0; i < k; i++)
                for (int j = 0; j < g.pathMap[i].size(); j++)
                {
                        int u = g.pathMap[i][j];
                        for (int v : g.adj_matrix[u])
                        {
                                if (j != g.pathMap.size() - 1 && v ==
g.pathMap[i][j + 1]) continue; if (nodeSet.find(v) == nodeSet.end()) stmp <<
"\t" << v << "[label = \"id=" << v << "$|place-" << g.orderToClique[v].node <<
"|" << g.orderToClique[v].start << "|" << g.orderToClique[v].end << "\"];" <<
endl; stmp << "\t" << u << "->" << v << ";" << endl;
                        }
                        for (int v : prec_adj[u])
                        {
                                if (j != 0 && v == g.pathMap[i][j - 1])
                                        continue;
                                if (nodeSet.find(v) == nodeSet.end())
                                        stmp << "\t" << v << "[label = \"id=" <<
v << "$|place-" << g.orderToClique[v].node << "|" << g.orderToClique[v].start <<
"|" << g.orderToClique[v].end << "\"];" << endl; stmp << "\t" << v << "->" << u
<< ";" << endl;
                        }
                }

        stmp << "}" << endl;
        ftmp << stmp.str();
        ftmp.close();
}
*/

/*
 * 通过将有向图投影为无向图寻找独立子图
 */
vector<subGraph> allSubGraph(Graph &g) {
    vector<vector<int>> connection(g.num_vertices() - 1, vector<int>());
    vector<unordered_set<int>> subDagNodes;
    // // 避开最后添加的虚拟点
    // for (int u = 0; u < g.num_vertices() - 1; u++)
    //     for (int v : g.out_edges(u)) {
    //         if (v == g.num_vertices() - 1) continue;
    //         connection[u].push_back(v);
    //         connection[v].push_back(u);
    //     }

    // vector<int> nodeToSub(g.num_vertices() - 1, -1);
    // unordered_set<int> cliqueID;
    // for (int i = 0; i < g.num_vertices() - 1; i++) cliqueID.insert(i);
    // while (cliqueID.size() > 0) {
    //     subDagNodes.push_back(unordered_set<int>());
    //     int num = subDagNodes.size() - 1;
    //     queue<int> q;
    //     q.push(*cliqueID.begin());
    //     cliqueID.erase(q.front());
    //     while (q.size() > 0) {
    //         int u = q.front();
    //         q.pop();
    //         subDagNodes.back().insert(u);
    //         nodeToSub[u] = num;
    //         for (int v : connection[u])
    //             if (cliqueID.count(v)) {
    //                 q.push(v);
    //                 cliqueID.erase(v);
    //             }
    //     }
    // }

    vector<subGraph> subGraphs(subDagNodes.size(), subGraph());
    // for (int i = 0; i < subDagNodes.size(); i++) {
    //     subGraphs[i].set_gid(i);
    //     for (int u : subDagNodes[i]) {
    //         subGraphs[i].addVertex(g[u]);
    //         for (int v : g.out_edges(u)) {
    //             subGraphs[i].addVertex(g[v]);
    //             subGraphs[i].addEdge(make_tuple(g[u].node, g[u].start, g[u].end),
    //                                  make_tuple(g[v].node, g[v].start, g[v].end));
    //         }
    //     }
    // }

    return subGraphs;
}

void printGraph(Graph &graph, string prefix, int id) {
    fstream ftmp("result/" + prefix + "_" + to_string(id) + ".dot", ios::out);
    stringstream stmp;
    stmp << "digraph graphe {\ngraph[rankdir=\"LR\"];\nnodesep = 3;\nnode "
            "[fontname=\"Arial\", shape = \"record\",style=filled,size = "
            "5];\nedge [color=black];"
         << endl;

    for (int i = 0; i < graph.num_vertices(); i++) {
        stmp << "{rank=\"same\";";

        // for (int j = 0; j < depthMap[i].size(); j++) {
        stmp << graph[i].id << "[label = \"id=" << graph[i].id << "|start=" << graph[i].start
             << "|end=" << graph[i].end;
        stmp << "\"]};" << endl;
        // }
        // for (int node : depthMap[i])
        //   for (int u : graph.in_edges(node))
        //     stmp << u << "->" << node << ";" << endl;
        for (auto edge : graph.out_edges(i)) {
            stmp << i << "->" << edge << ";" << endl;
        }
    }

    stmp << "}" << endl;
    ftmp << stmp.str();
    ftmp.close();
}

void printGraph(Tree &tree, string prefix, int depth, int id) {
    vector<int> nodeDepth(tree.num_vertices(), -1);
    vector<vector<int>> depthMap(depth, vector<int>());
    int max_dep = 0;

    // for (int i = 0; i < tree.num_vertices(); i++)
    // {
    // 	int d = dfs(tree, nodeDepth, i);
    // 	max_dep = std::max(max_dep, d);
    // 	if (d < depth)
    // 		depthMap[d].push_back(i);
    // }
    // max_dep = std::min(max_dep, depth);

    fstream ftmp("result/" + prefix + "_" + to_string(id) + ".dot", ios::out);
    stringstream stmp;
    stmp << "digraph graphe {\ngraph[rankdir=\"LR\"];\nnodesep = 3;\nnode "
            "[fontname=\"Arial\", shape = \"record\",style=filled,size = "
            "5];\nedge [color=black];"
         << endl;

    for (int i = 0; i < max_dep; i++) {
        stmp << "{rank=\"same\";";
        for (int j = 0; j < depthMap[i].size() - 1; j++) stmp << to_string(depthMap[i][j]) << ",";
        stmp << depthMap[i].back() << "}";
        for (int j = 0; j < depthMap[i].size(); j++) {
            stmp << depthMap[i][j] << "[label = \"id=" << depthMap[i][j] << "|depth=" << i
                 << "|out=" << tree.out_degree(depthMap[i][j]);
            if (tree.isStart_[depthMap[i][j]] == 1)
                stmp << ";color=red";
            else if (tree.isStart_[depthMap[i][j]] == 2)
                stmp << ";color=green";
            else if (tree.isStart_[depthMap[i][j]] == 3)
                stmp << ";color=blue";
            stmp << "\"];" << endl;
        }
        for (int node : depthMap[i])
            for (int u : tree.in_edges(node)) stmp << u << "->" << node << ";" << endl;
    }

    stmp << "}" << endl;
    ftmp << stmp.str();
    ftmp.close();
}
void printDAG_metis(Graph &graph, string outputPath) {
    fstream fout(outputPath.c_str(), ios::out);
    stringstream ss;
    ss << graph.num_vertices() << " " << graph.num_edges() << endl;
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (int j : graph.out_edges(i)) {
            ss << (j + 1) << " ";
        }
        for (int j : graph.in_edges(i)) {
            ss << (j + 1) << " ";
        }
        ss << endl;
    }
    fout << ss.str();
    fout.close();
}

void printDAG_1(Graph &graph, string outputPath) {
    fstream fout(outputPath.c_str(), ios::out);
    stringstream ss;
    ss << graph.num_vertices() << " " << graph.num_edges() << endl;
    // std::cout << "edge = " << graph.num_edges() << "\n";
    int size = 0;
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (int j : graph.out_edges(i)) {
            ss << (j + 1) << " ";
            size++;
        }
        ss << endl;
    }
    // std::cout << "edge size = " << size << "\n";
    fout << ss.str();
    fout.close();
}

void printDAG_2(Graph &graph, string outputPath) {
    fstream fout(outputPath.c_str(), ios::out);
    stringstream ss;
    ss << "graph_for_greach" << endl;
    ss << graph.num_vertices() << endl;
    for (int i = 0; i < graph.num_vertices(); i++) {
        ss << i << ": ";
        for (int j : graph.out_edges(i)) ss << j << " ";
        ss << "#" << endl;
    }

    fout << ss.str();
    fout.close();
}

// bool TreeQuery(uint32_t u, uint32_t v, list<res_type>::iterator &res_iter)
// {
//   uint32_t from = u;
//   uint32_t to = v;

//   //
//   正向树上from可以到达to的条件是from的区间包含to的区间；反向树上from可以到达to的条件是to的区间包含from的区间
//   //
//   树上转移点对的first元素都是沿着边相反方向传递的节点，second元素都是沿着边正向传递的节点
//   if (!res_iter->reverse_)
//   {
//     if (CanReach(res_iter->section_table_[from],
//     res_iter->section_table_[to]))
//       return true;
//     else
//     {
//       bool reachable = false;
//       for (auto from : res_iter->black_table_[from].second)
//         for (auto to : res_iter->black_table_[to].first)
//         {
//           reachable =  TreeQuery()
//         }
//     }
//   }
//   else
//   {
//     if (CanReach(res_iter->section_table_[to],
//     res_iter->section_table_[from]))
//       return true;
//     else
//     {
//       from = iter->black_table_[from].second;
//       to = iter->black_table_[to].first;
//     }
//   }

//   return false;
// }

// /*
//  * 查询过程分为三步：
//  * 1. 通过 u 和 start确定此时 u
//  * 连接的下部点v，然后对于v连接的边，对其时间排序后找到start所在的区间，然后找到该区间对应的DAG上的节点
//  * c1
//  * 2. 通过 w 和 end 确定此时 w
//  * 连接的下部点v'，然后对于v'连接的边，对其时间排序后找到end所在的区间，然后找到该区间对应的DAG上的节点
//  * c2
//  * 3. 通过 res 查询 c1 到 c2 的
//  *
//  * 时间复杂度：O( d(u) + d(w) + log(d(v_1))+log(d(v_2)) + k + 2)
//  */
// bool Query(int u, int w, time_t start, time_t end, BiGraph &bg, Graph &g,
//            list<res_type> &res)
// {
//   int i = 0;
//   int j = 0;
//   int v_1 = -1, v_2 = -1, v_3 = -1, v_4 = -1;

//   while (i < bg.adj_matrix_u[u].size())
//   {
//     if (start >= bg.timeSection_u[u][i].first &&
//         start < bg.timeSection_u[u][i].second)
//     {
//       v_1 = bg.adj_matrix_u[u][i];
//       i++;
//     }
//     else if (end >= bg.timeSection_u[u][i].first &&
//              end < bg.timeSection_u[u][i].second)
//     {
//       v_3 = bg.adj_matrix_u[u][i];
//       break;
//     }
//     else
//       i++;
//   }

//   while (j < bg.adj_matrix_u[w].size())
//   {
//     if (start >= bg.timeSection_u[u][j].first &&
//         start < bg.timeSection_u[u][j].second)
//     {
//       v_4 = bg.adj_matrix_u[u][i];
//       j++;
//     }
//     else if (end >= bg.timeSection_u[u][j].first &&
//              end < bg.timeSection_u[u][j].second)
//     {
//       v_2 = bg.adj_matrix_u[u][j];
//       break;
//     }
//     else
//       j++;
//   }

//   if ((v_1 == -1 || v_2 == -1) && (v_3 == -1 || v_4 == -1))
//     return false;

//   auto iter_1 =
//       lower_bound(bg.timeTable[v_1].begin(), bg.timeTable[v_1].end(), start);
//   iter_1--;
//   auto iter_2 =
//       lower_bound(bg.timeTable[v_2].begin(), bg.timeTable[v_2].end(), end);
//   iter_2--;
//   auto iter_3 =
//       lower_bound(bg.timeTable[v_3].begin(), bg.timeTable[v_3].end(), end);
//   iter_3--;
//   auto iter_4 =
//       lower_bound(bg.timeTable[v_4].begin(), bg.timeTable[v_4].end(), start);
//   iter_4--;

//   uint32_t c_1 = g.cnum(make_pair(v_1, *iter_1));
//   uint32_t c_2 = g.cnum(make_pair(v_2, *iter_2));
//   uint32_t c_3 = g.cnum(make_pair(v_3, *iter_3));
//   uint32_t c_4 = g.cnum(make_pair(v_4, *iter_4));

//   return TreeQuery(c_1, c_2, res) || TreeQuery(c_4, c_3, res);
// }

// void ToolPrintGraph(Graph &g, const Tree &tree,
//                     const std::unordered_set<uint64_t> &cross_edges)
// {
//   std::cout << "Source,Target,Color\n";
//   for (auto t : tree.edge_pair_)
//   {
//     if (tree.froot == t.first or tree.froot == t.second)
//       continue;
//     std::cout << t.first << "," << t.second << "," << 3 << "\n";
//   }
//   for (int i = 0; i < g.num_vertices(); i++)
//   {
//     for (auto edge : g.out_edges(i))
//     {
//       auto tmp = MERGE_TO_UINT64(i, edge);
//       if (cross_edges.count(tmp) > 0)
//       {
//         std::cout << i << "," << edge << "," << 2 << "\n";
//       }
//       else
//       {
//         std::cout << i << "," << edge << "," << 1 << "\n";
//       }
//     }
//   }
// }

std::string extractFileName(const std::string &filePath) {
    // 使用 basename 函数提取文件名
    char *baseName = basename(const_cast<char *>(filePath.c_str()));
    return std::string(baseName);
}
}  // namespace CCR
