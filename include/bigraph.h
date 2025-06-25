#ifndef GRAPH_H
#define GRAPH_H

#include <string.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

class BiGraph {
   public:
    vector<vector<int>> adj_matrix_u;
    vector<vector<int>> adj_matrix_l;
    vector<vector<pair<time_t, time_t>>> timeSection_u;
    vector<vector<pair<time_t, time_t>>> timeSection_l;
    vector<vector<time_t>> timeTable;                  // 记录每个下部点的时间戳
    vector<vector<pair<time_t, time_t>>> timeTable_u;  // 记录每个下部点真正用在dag中的时间戳

    BiGraph(){};
    BiGraph(int u, int l) {
        this->adj_matrix_u = vector<vector<int>>(u, vector<int>());
        this->adj_matrix_l = vector<vector<int>>(l, vector<int>());
        this->timeSection_u = vector<vector<pair<time_t, time_t>>>(u, vector<pair<time_t, time_t>>());
        this->timeSection_l = vector<vector<pair<time_t, time_t>>>(l, vector<pair<time_t, time_t>>());
        this->timeTable = vector<vector<time_t>>();
        this->timeTable_u = vector<vector<pair<time_t, time_t>>>();
    }
    BiGraph(string filepath, string mode) {
        ifstream fin(filepath.c_str(), ios::in);
        vector<set<time_t>> timeSet;
        int num_upper, num_lower;

        if (!fin.is_open()) {
            this->adj_matrix_u = vector<vector<int>>();
            this->adj_matrix_l = vector<vector<int>>();
            this->timeSection_u = vector<vector<pair<time_t, time_t>>>();
            this->timeSection_l = vector<vector<pair<time_t, time_t>>>();
            this->timeTable = vector<vector<time_t>>();
            cout << "读取文件失败" << endl;
            exit(1);
        } else if (mode == "g") {
            string buffer;
            fin >> buffer >> num_upper >> num_lower;

            this->adj_matrix_u = vector<vector<int>>(num_upper, vector<int>());
            this->adj_matrix_l = vector<vector<int>>(num_lower, vector<int>());
            this->timeSection_u = vector<vector<pair<time_t, time_t>>>(num_upper, vector<pair<time_t, time_t>>());
            this->timeSection_l = vector<vector<pair<time_t, time_t>>>(num_lower, vector<pair<time_t, time_t>>());
            this->timeTable = vector<vector<time_t>>(num_lower, vector<time_t>());
            this->timeTable_u = vector<vector<pair<time_t, time_t>>>(num_lower, vector<pair<time_t, time_t>>());

            int min_d = 1 * 600;
            int max_d = 8 * 3600;
            float alpha = -2.5;

            timeSet = vector<set<time_t>>(num_lower, set<time_t>());

            int64_t upper, start, lower;
            while (fin >> upper >> start >> lower) {
                time_t end;

                float y = ((float)rand()) / ((float)RAND_MAX);
                double x1 = pow(max_d, (alpha + 1));
                double x2 = pow(min_d, (alpha + 1));
                double x3 = (x2 - x1) * y + x1;
                double x4 = pow(x3, 1 / (alpha + 1));
                time_t duration = (time_t)x4;

                // time_t duration = (time_t)pow(((pow(max_d, (alpha +
                // 1)) - pow(min_d, (alpha + 1))) * y + pow(min_d,
                // (alpha + 1))), (1 / (alpha + 1)));

                end = start + duration;

                this->adj_matrix_u[upper].push_back(lower);
                this->adj_matrix_l[lower].push_back(upper);
                this->timeSection_u[upper].push_back(pair<time_t, time_t>(start, end));
                this->timeSection_l[lower].push_back(pair<time_t, time_t>(start, end));
                timeSet[lower].insert(start);
                timeSet[lower].insert(end);
            }
            for (int i = 0; i < adj_matrix_u.size(); i++) {
                reverse(adj_matrix_u[i].begin(), adj_matrix_u[i].end());
                reverse(timeSection_u[i].begin(), timeSection_u[i].end());
            }
        } else if (mode == "r") {
            string buffer;
            fin >> buffer >> num_upper >> num_lower;

            this->adj_matrix_u = vector<vector<int>>(num_upper, vector<int>());
            this->adj_matrix_l = vector<vector<int>>(num_lower, vector<int>());
            this->timeSection_u = vector<vector<pair<time_t, time_t>>>(num_upper, vector<pair<time_t, time_t>>());
            this->timeSection_l = vector<vector<pair<time_t, time_t>>>(num_lower, vector<pair<time_t, time_t>>());
            this->timeTable = vector<vector<time_t>>(num_lower, vector<time_t>());
            this->timeTable_u = vector<vector<pair<time_t, time_t>>>(num_lower, vector<pair<time_t, time_t>>());

            int upper, lower;
            time_t start, end;

            timeSet = vector<set<time_t>>(num_lower, set<time_t>());

            while (fin >> upper >> lower >> start >> end) {
                this->adj_matrix_u[upper].push_back(lower);
                this->adj_matrix_l[lower].push_back(upper);
                this->timeSection_u[upper].push_back(pair<time_t, time_t>(start, end));
                this->timeSection_l[lower].push_back(pair<time_t, time_t>(start, end));
                timeSet[lower].insert(start);
                timeSet[lower].insert(end);
            }
        }

        // 增序排列时间戳
        for (int i = 0; i < adj_matrix_u.size(); i++) {
            reverse(adj_matrix_u[i].begin(), adj_matrix_u[i].end());
            reverse(timeSection_u[i].begin(), timeSection_u[i].end());
        }

        for (int i = 0; i < num_lower; i++) this->timeTable[i].assign(timeSet[i].begin(), timeSet[i].end());

        fin.close();
    }

    BiGraph(const BiGraph &g) {
        this->adj_matrix_u = *(new vector<vector<int>>(g.adj_matrix_u));
        this->adj_matrix_l = *(new vector<vector<int>>(g.adj_matrix_l));
        this->timeSection_u = *(new vector<vector<pair<time_t, time_t>>>(g.timeSection_u));
        this->timeSection_l = *(new vector<vector<pair<time_t, time_t>>>(g.timeSection_l));
        this->timeTable = *(new vector<vector<time_t>>(g.timeTable));
        this->timeTable_u = *(new vector<vector<pair<time_t, time_t>>>(g.timeTable_u));
    }
    ~BiGraph(){};
};

// class Clique
// {
// public:
//     int node;
//     time_t start;
//     time_t end;
//     int num;

//     Clique() : node(-1), start(-1), end(-1), num(0){};
//     Clique(int node_, time_t start_, time_t end_, int num_) : node(node_),
//     start(start_), end(end_), num(num_){}; bool operator==(const Clique &a)
//     const
//     {
//         return this->node == a.node && this->start == a.start;
//     }

//     bool operator<(const Clique &a) const
//     {
//         return this->start < a.start;
//     }

//     void info()
//     {
//         cout << "Clique" << node << "  " << start << "  " << end << endl;
//     }

//     void info() const
//     {
//         cout << "Clique" << node << "  " << start << "  " << end << endl;
//     }
// };

// class CliqueHash
// {
// public:
//     size_t operator()(const Clique &clique) const
//     {
//         size_t h1 = (uint64_t)clique.node << 32;
//         size_t h2 = clique.start;
//         return h1 | h2;
//     }
// };

// class CliqueEqual
// {
// public:
//     bool operator()(const Clique &a, const Clique &b) const
//     {
//         return a.node == b.node && a.start == b.start && a.end == b.end;
//     }
// };

/*
class Graph
{
public:
    int n;                                                             // 点个数
    vector<unordered_set<int>> adj_matrix;                             // 邻接表
    vector<int> topoOrder;                                             //
拓扑排序 unordered_map<Clique, int, CliqueHash, CliqueEqual> cliqueToOrder; //
从Clique映射为原编号 vector<Clique> orderToClique; vector<vector<int>> pathMap;
    vector<int> nodePath;

    Graph()
    {
        n = 0;
    };

    Graph(const Graph &g)
    {
        this->n = g.n;
        this->adj_matrix = vector<unordered_set<int>>(g.adj_matrix);
        this->topoOrder = vector<int>(g.topoOrder);
        this->cliqueToOrder = unordered_map<Clique, int, CliqueHash,
CliqueEqual>(g.cliqueToOrder); this->orderToClique =
vector<Clique>(g.orderToClique); this->pathMap = vector<vector<int>>(g.pathMap);
        this->nodePath = vector<int>(g.nodePath);
    };

    ~Graph(){};

    void add_node(const Clique &a)
    {
        if (this->cliqueToOrder.find(a) != this->cliqueToOrder.end())
            return;

        this->orderToClique.push_back(a);
        this->cliqueToOrder[a] = this->n;
        this->adj_matrix.push_back(unordered_set<int>());
        this->n++;
    }

    void add_edge(const Clique &a, const Clique &b)
    {
        add_node(a);
        add_node(b);

        this->adj_matrix[this->cliqueToOrder[a]].insert(this->cliqueToOrder[b]);
    }

    void topoSort()
    {
        queue<int> q;
        vector<int> inDegree(this->n, 0);

        for (int u = 0; u < this->n; u++)
            for (int v : this->adj_matrix[u])
                inDegree[v]++;

        for (int i = 0; i < this->n; i++)
            if (inDegree[i] == 0)
                q.push(i);

        while (q.size() > 0)
        {
            int v = q.front();
            q.pop();

            this->topoOrder.push_back(v);
            for (int j : this->adj_matrix[v])
            {
                inDegree[j]--;
                if (inDegree[j] == 0)
                    q.push(j);
            }
        }
    }

    void closure()
    {
        vector<vector<bool>> adj;
        for (int i = 0; i < n; i++)
        {
            vector<bool> tmp(n, false);
            for (auto j : this->adj_matrix[i])
                tmp[j] = true;
            adj.push_back(tmp);
        }

        for (int k = 0; k < n; k++)
            for (int i = 0; i < n; i++)
                if (adj[i][k])
                    for (int j = 0; j < n; j++)
                        adj[i][j] = adj[i][j] | adj[k][j];

        vector<vector<int>> newAdj(n, vector<int>());
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (adj[i][j])
                    newAdj[i].push_back(j);
        this->adj_matrix = newAdj;
    }
};

class subGraph : public Graph
{
public:
    int gID;
};

class cmpGraph
{
public:
    int n;                                 // 点个数
    vector<unordered_set<int>> adj_matrix; // 邻接表
    vector<int> topoOrder;                 // 拓扑排序
    vector<vector<int>> superNodes;        // 超点包含的Clique
    vector<Clique> SntoClique;             // 超点内部序列的第一个节点
    vector<int> parent;                    // Clique所属的超点

    cmpGraph(cmpGraph &g)
    {
        this->n = g.n;
        this->adj_matrix = vector<unordered_set<int>>(g.adj_matrix);
        this->topoOrder = vector<int>(g.topoOrder);
        this->superNodes = vector<vector<int>>(g.superNodes);
        this->SntoClique = vector<Clique>(g.SntoClique);
        this->parent = vector<int>(g.parent);
    }

    cmpGraph(Graph &g)
    {
        parent = vector<int>(g.n, -1);

        // 确定原图每个点的出边和入边数量
        vector<int> in(g.n, 0);
        vector<int> out(g.n, 0);
        for (int i = 0; i < g.n; i++)
        {
            out[i] = g.adj_matrix[i].size();
            for (int j : g.adj_matrix[i])
                in[j]++;
        }

        // 先根据路径压缩Clique
        for (int p = 0; p < g.pathMap.size(); p++)
        {
            vector<int> path = g.pathMap[p];
            vector<int> nodes;
            for (int u : path)
            {
                // 该节点前的一组Clique缩为一个超点
                if (in[u] > 1 && out[u] < 2)
                {
                    if (nodes.size() > 0)
                    {
                        superNodes.push_back(nodes);
                        SntoClique.push_back(Clique(g.orderToClique[nodes.front()].node,
g.orderToClique[nodes.front()].start, g.orderToClique[nodes.back()].end));

                        for (int i : nodes)
                            parent[i] = superNodes.size() - 1;
                    }
                    nodes = {u};
                    if (u == path.back())

                    {
                        superNodes.push_back(nodes);
                        SntoClique.push_back(Clique(superNodes.size() - 1,
g.orderToClique[nodes.front()].start, g.orderToClique[nodes.back()].end)); for
(int i : nodes) parent[i] = superNodes.size() - 1;
                    }
                }
                // 该节点前的一组Clique缩为一个超点，且该节点单独缩为一个超点
                else if (in[u] > 1 && out[u] > 1)
                {
                    if (nodes.size() > 0)
                    {
                        superNodes.push_back(nodes);
                        SntoClique.push_back(Clique(superNodes.size() - 1,
g.orderToClique[nodes.front()].start, g.orderToClique[nodes.back()].end)); for
(int i : nodes) parent[i] = superNodes.size() - 1;
                    }
                    superNodes.push_back({u});
                    SntoClique.push_back(Clique(superNodes.size() - 1,
g.orderToClique[u].start, g.orderToClique[u].end)); parent[u] =
superNodes.size() - 1; nodes.clear();
                }
                // 将包含该节点在内的一组Clique缩为一个超点
                else if (in[u] < 2 && out[u] > 1)
                {
                    nodes.push_back(u);
                    superNodes.push_back(nodes);
                    SntoClique.push_back(Clique(superNodes.size() - 1,
g.orderToClique[nodes.front()].start, g.orderToClique[nodes.back()].end)); for
(int i : nodes) parent[i] = superNodes.size() - 1; nodes.clear();
                }
                // 将该节点添加到当前组内
                else if (in[u] < 2 && out[u] < 2)
                {
                    nodes.push_back(u);
                    if (u == path.back())
                    {
                        superNodes.push_back(nodes);
                        SntoClique.push_back(Clique(superNodes.size() - 1,
g.orderToClique[nodes.front()].start, g.orderToClique[nodes.back()].end)); for
(int i : nodes) parent[i] = superNodes.size() - 1;
                    }
                }
            }
        }

        for (int i = 0; i < g.n; i++)
            if (parent[i] < 0)
                cout << i << endl;

        this->n = superNodes.size();

        // 针对原图邻接表通过parent数组构建超点的邻接表
        this->adj_matrix = vector<unordered_set<int>>(this->n,
unordered_set<int>()); for (int u = 0; u < g.n; u++) for (int v :
g.adj_matrix[u]) if (parent[u] != parent[v])
                    this->adj_matrix[parent[u]].insert(parent[v]);
    }

    void topoSort()
    {
        queue<int> q;
        vector<int> inDegree(this->n, 0);

        for (int u = 0; u < this->n; u++)
            for (int v : this->adj_matrix[u])
                inDegree[v]++;

        for (int i = 0; i < this->n; i++)
            if (inDegree[i] == 0)
                q.push(i);

        while (q.size() > 0)
        {
            int v = q.front();
            q.pop();

            this->topoOrder.push_back(v);
            for (int j : this->adj_matrix[v])
            {
                inDegree[j]--;
                if (inDegree[j] == 0)
                    q.push(j);
            }
        }
    }
};
*/

#endif  // GRAPH_H