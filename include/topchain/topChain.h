#ifndef TOPCHAIN_H
#define TOPCHAIN_H
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

#include "graph.h"
#include "utils.h"
using namespace std;
namespace topchain {

const int infinity = 1e9;

class Edge {
   public:
    void init(int v, int w, int size) {
        this->v = v;
        this->w = w;
        this->interval.resize(size);
    }

   public:
    int v;
    int w;
    vector<int> interval;
};

class Graph {
   public:
    Graph() {}
    Graph(ptree::Graph &g);

   public:
    vector<vector<Edge>> adj;
    int V, static_E, dynamic_E;
};

class GraphT {
   public:
    GraphT() {}
    GraphT(const Graph &g);
    void add_edge(int u, int v);
    void make_unique(vector<int> &data);
    int getVertexID(int u, int offset, int startOrArrival);

    void initial_query(const char *filePath);  // query file
    void initial_query();
    void initial_ds();
    void initial_ds_f();

    void run_earliest_arrival();
    void earliest_arrival(int source);
    void run_fastest();
    void fastest(int source);

    // for testing the correctness
    void run_earliest_arrival(const char *filePath);  // output the result
    void earliest_arrival(int source, FILE *file);
    void run_fastest(const char *filePath);
    void fastest(int source, FILE *file);

    void print_result(const int source, const vector<int> &t_time, FILE *file);
    void print_avg_time();
    void print_avg_time(const char *filePath1, const char *filePath2);
    void show_statistics(const char *filePath1, const char *filePath2);

    bool reachability(int from, int to, int t1, int t2);
    int getFrom(int from, int t1, int t2);
    int getTo(int to, int t1, int t2);
    int output(string s);

   public:
    int V, E, gv;
    vector<vector<int>> adj;
    vector<vector<int>> startT, arrivalT;
    vector<int> curSize, sumSize;
    vector<int> timeId, originalId;
    vector<int> sources;
    int t_start, t_end;
    double time_sum;
    vector<int> a_time, f_time;
    queue<int> myqueue;
    vector<bool> mark_in_Q;  // whether the vertex is in queue
    vector<bool> visited;    // Whether the vertex has been visited

    int originalV;
};

void readGraph(GraphT &g);
void readChainInfo(GraphT &g);
void transform(ptree::Graph &g);
void findChain();
void feline();
void handle();
void copyIndex();
GraphT transToScc(ptree::Graph &g);
void queryByTopChain(BiGraph &bg, ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo, GraphT gt);

class ReachIndex {
   public:
    // uint16_t layerup;
    // uint16_t layerdown;
    int OutLimit;
    int InLimit;
    int HLimit;
    pair<int, int> *Label;
    int Lup;
    int Ldown;
    int flabel1;
    int flabel2;

    ReachIndex() { Label = NULL; }

    ~ReachIndex() {
        if (Label != NULL) delete[] Label;
    }
};
struct DEG {
    int d;
    int index;
};
bool cmp1(DEG d1, DEG d2);

struct TopChainMeta {
    vector<int> inOut;
    vector<int> pNext;
    vector<vector<int>> v;
    vector<vector<int>> vr;
    int n, m;
    vector<int> ind;
    vector<int> outd;
    vector<int> Lup;
    vector<int> Ldown;

    vector<vector<pair<int, int>>> Lin;
    vector<vector<pair<int, int>>> Lout;

    int k = 5;
    vector<DEG> deg;
    vector<int> rank_order;
    vector<int> reverserank_order;
    int chainNum;
    vector<int> vis;
    vector<int> toChain;
    vector<int> toPos;
    vector<vector<int>> toV;
    vector<vector<int>> startT, arrivalT;
    vector<int> curSize, sumSize;
    vector<int> timeId, originalId;
    // copy index
    vector<ReachIndex> reachindex;
    vector<int> flabel1;
    vector<int> flabel2;
    vector<int> topoOrder;
    vector<int> rtopoOrder;
    int maxDown = 0;
    int maxUp = 0;
    int mark = 0;
    vector<int> flagVec;

    void readGraph(GraphT &gt) {
        n = gt.V;
        m = gt.E;
        v.resize(n);
        vr.resize(n);
        ind.resize(n);
        outd.resize(n);

        int u;
        int size;
        int to;
        for (int i = 0; i < n; i++) {
            u = i;
            size = gt.adj[i].size();
            for (int j = 0; j < size; j++) {
                to = gt.adj[i][j];
                v[u].push_back(to);
                vr[to].push_back(u);
                outd[u]++;
                ind[to]++;
            }
        }
        cout << "load done." << endl;
    }
    int getVertexID(int u, int offset, int startOrArrival) {
        //  0 arrival 1 start

        int base = u == 0 ? 0 : sumSize[u - 1];
        const vector<int> &arrival = arrivalT[u];
        const vector<int> &start = startT[u];

        if (startOrArrival == 0) {
            return base + offset;
        } else if (startOrArrival == 1) {
            return base + arrival.size() + offset;
        }
        return -1;
    }

    void readChainInfo(topchain::GraphT &gt) {
        inOut.resize(n);
        pNext.resize(n);

        chainNum = gt.originalV;
        // fstream in;
        // in.open(s.c_str(), ios::in);
        // in >> chainNum;
        startT = vector<vector<int>>(chainNum, vector<int>());
        arrivalT = vector<vector<int>>(chainNum, vector<int>());
        curSize = vector<int>(chainNum, 0);
        sumSize = vector<int>(chainNum, 0);

        originalId.resize(n);

        toV.resize(n);
        toChain.resize(n);
        toPos.resize(n);

        int u;
        int size;
        int tmp;
        for (int i = 0; i < chainNum; ++i) {
            u = i;
            size = gt.arrivalT[i].size();
            for (int j = 0; j < size; ++j) {
                tmp = gt.arrivalT[i][j];
                arrivalT[i].push_back(tmp);
            }
        }
        for (int i = 0; i < chainNum; ++i) {
            u = i;
            size = gt.startT[i].size();

            for (int j = 0; j < size; ++j) {
                tmp = gt.startT[i][j];
                startT[i].push_back(tmp);
            }
        }
        // cout << "chain info load done" << endl;
    }
    void findChain() {
        for (int i = 0; i < chainNum; i++) {
            int cur = startT[i].size() + arrivalT[i].size();
            curSize[i] = cur;
            sumSize[i] = cur + (i > 0 ? sumSize[i - 1] : 0);
        }

        // find chain and pos
        for (int i = 0; i < chainNum; ++i) {
            int p1 = 0;
            int p2 = 0;
            int j = 0;
            while (j < arrivalT[i].size() + startT[i].size()) {
                if (p1 == arrivalT[i].size()) {
                    int u = getVertexID(i, p2, 1);
                    originalId[u] = i;
                    toPos[u] = startT[i][p2];  // toV[i].size();
                    // toV[i].push_back(u);
                    p2++;
                } else if (p2 == startT[i].size()) {
                    int u = getVertexID(i, p1, 0);
                    originalId[u] = i;
                    toPos[u] = arrivalT[i][p1];  // toV[i].size();
                    // toV[i].push_back(u);
                    p1++;
                } else if (arrivalT[i][p1] <= startT[i][p2]) {
                    int u = getVertexID(i, p1, 0);
                    originalId[u] = i;
                    toPos[u] = arrivalT[i][p1];  // toV[i].size();
                    // toV[i].push_back(u);
                    p1++;
                } else {
                    int u = getVertexID(i, p2, 1);
                    originalId[u] = i;
                    toPos[u] = startT[i][p2];  // toV[i].size();
                    // toV[i].push_back(u);
                    p2++;
                }

                j++;
            }
        }

        // assign degree
        deg.resize(chainNum);
        for (int i = 0; i < chainNum; ++i) {
            deg[i].d = startT[i].size() + arrivalT[i].size();
            deg[i].index = i;
        }
    }

    // feline label

    void feline() {
        // cout << "feline begin" << endl;
        vector<int> inds1;
        vector<int> inds2;
        stack<int> stack1;
        stack<int> stack2;
        stack<int> stack3;
        inds1.resize(n);
        inds2.resize(n);
        flabel1.resize(n);
        flabel2.resize(n);
        topoOrder.resize(n);
        rtopoOrder.resize(n);
        Ldown.resize(n);
        Lup.resize(n);
        // copy ind
        for (int i = 0; i < n; ++i) {
            inds1[i] = inds2[i] = ind[i];
        }
        queue<int> q;
        for (int i = n - 1; i >= 0; --i) {
            if (ind[i] == 0) {
                stack1.push(i);
                Ldown[i] = 0;
            }

            if (outd[i] == 0) {
                stack3.push(i);
                Lup[i] = 0;
            }

            // if (ind[i] == 0) q.push(i);
        }

        //   //timer t1;
        //   t1.start();
        int count1 = 0;
        while (!stack1.empty()) {
            // cout << "!!" << endl;
            int tmp = stack1.top();
            stack1.pop();
            topoOrder[count1] = tmp;
            flabel1[tmp] = count1++;
            for (int i = v[tmp].size() - 1; i >= 0; --i) {
                inds1[v[tmp][i]]--;
                if (inds1[v[tmp][i]] == 0) {
                    stack1.push(v[tmp][i]);
                    // Ldown[v[tmp][i]] = Ldown[tmp]+1;
                    // maxDown = max(Ldown[v[tmp][i]], maxDown);
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            int tmp = topoOrder[i];
            if (vr[tmp].size() > 0) {
                int maxi = Ldown[vr[tmp][0]];
                for (int j = 0; j < vr[tmp].size(); ++j) {
                    if (Ldown[vr[tmp][j]] > maxi) maxi = Ldown[vr[tmp][j]];
                }
                Ldown[tmp] = maxi + 1;
                maxDown = max(maxDown, Ldown[tmp]);
            }
        }

        // find reverse topo

        count1 = 0;
        while (!stack3.empty()) {
            int tmp = stack3.top();
            stack3.pop();
            rtopoOrder[count1++] = tmp;
            for (int i = 0; i < vr[tmp].size(); ++i) {
                outd[vr[tmp][i]]--;
                if (outd[vr[tmp][i]] == 0) {
                    stack3.push(vr[tmp][i]);
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            int tmp = rtopoOrder[i];
            if (v[tmp].size() > 0) {
                int maxi = Lup[v[tmp][0]];
                for (int j = 0; j < v[tmp].size(); ++j) {
                    if (Lup[v[tmp][j]] > maxi) maxi = Lup[v[tmp][j]];
                }
                Lup[tmp] = maxi + 1;
                maxUp = max(maxUp, Lup[tmp]);
            }
        }

        count1 = 0;
        for (int i = 0; i < n; ++i) {
            if (ind[i] == 0) stack2.push(i);
        }
        while (!stack2.empty()) {
            int tmp = stack2.top();
            stack2.pop();
            flabel2[tmp] = count1++;
            for (int i = 0; i < v[tmp].size(); ++i) {
                inds2[v[tmp][i]]--;
                if (inds2[v[tmp][i]] == 0) stack2.push(v[tmp][i]);
            }
        }
        // assert(count1 == n);
        // cout << "feline finished" << endl;
    }

    void handle() {
        Lin.resize(n);
        Lout.resize(n);
        queue<int> q;
        queue<int> q1;
        int count1 = 0;

        // find chain

        for (int i = 0; i < n; ++i) toChain[i] = originalId[i];
        cout << "chain Num : " << chainNum << endl;

        /*
        int countt = 0;
        for (int i = 0; i < chainNum; ++ i)
        {
                countt += toV[i].size();
        }
        assert(countt == n);
        */

        //-------------------------------------------------------------
        // assign rank_order

        rank_order.resize(chainNum);

        // rank by degree
        sort(deg.begin(), deg.end(), cmp1);
        for (int i = 0; i < chainNum; ++i) {
            // order by degree
            rank_order[deg[i].index] = i;

            // order randomly
            //	rank_order[i] = i;
        }
        /*
        for (int i = chainNum-1; i > 0; --i)
        {
                int tmp = rand()%(i+1);
                swap(rank_order[i], rank_order[tmp]);
        }
        */

        //------------------

        // topological sort and find Lin

        // cout << "!!!begin time : " << timer.GetRuntime()*1000 << endl;
        cout << "begin calc lin\n";
        cout << " n = " << n << "\n";
        for (int i = 0; i < n; ++i) {
            int tmp = topoOrder[i];
            // IP construct
            vector<pair<int, int>> mid[2];
            Lin[tmp].push_back(make_pair(rank_order[toChain[tmp]], toPos[tmp]));

            mid[0].resize(k);
            mid[1].resize(k);
            int cur = 0;
            int curSize = 0;
            int linTmp = 1;
            int linTmpSize = Lin[tmp].size();
            mid[linTmp] = Lin[tmp];
            mid[linTmp].resize(k);
            for (int j = 0; j < vr[tmp].size(); ++j) {
                int nin = vr[tmp][j];
                // merge tmp nin
                int p1 = 0;
                int p2 = 0;
                curSize = 0;
                // mid[cur].resize(k);
                while (p1 < linTmpSize || p2 < Lin[nin].size()) {
                    // cout << "p1 p2: " << p1 << " " << p2 << endl;
                    if (curSize == k) break;

                    if (p1 == linTmpSize) {
                        mid[cur][curSize++] = Lin[nin][p2];
                        p2++;
                        continue;
                    }
                    if (p2 == Lin[nin].size()) {
                        mid[cur][curSize++] = mid[linTmp][p1];
                        p1++;
                        continue;
                    }

                    if (mid[linTmp][p1].first < Lin[nin][p2].first) {
                        mid[cur][curSize++] = mid[linTmp][p1];
                        p1++;
                    } else if (mid[linTmp][p1].first > Lin[nin][p2].first) {
                        mid[cur][curSize++] = Lin[nin][p2];
                        p2++;
                    } else {
                        mid[cur][curSize++] = max(mid[linTmp][p1], Lin[nin][p2]);
                        p1++;
                        p2++;
                    }
                }
                // Lin[tmp] = mid;
                cur = 1 - cur;
                linTmp = 1 - linTmp;
                // mid[linTmp].resize(curSize);
                linTmpSize = curSize;
            }
            Lin[tmp] = mid[linTmp];
            Lin[tmp].resize(linTmpSize);
        }

        cout << "label in done. " << endl;

        for (int i = 0; i < n; ++i) {
            int tmp = rtopoOrder[i];
            // IP construct
            vector<pair<int, int>> mid[2];
            Lout[tmp].push_back(make_pair(rank_order[toChain[tmp]], toPos[tmp]));

            mid[0].resize(k);
            mid[1].resize(k);
            int cur = 0;
            int curSize = 0;
            int loutTmp = 1;
            int loutTmpSize = Lout[tmp].size();
            mid[loutTmp] = Lout[tmp];
            mid[loutTmp].resize(k);
            for (int j = 0; j < v[tmp].size(); ++j) {
                int nin = v[tmp][j];
                // merge tmp nin
                int p1 = 0;
                int p2 = 0;
                curSize = 0;
                // mid[cur].resize(k);
                while (p1 < loutTmpSize || p2 < Lout[nin].size()) {
                    // cout << "p1 p2: " << p1 << " " << p2 << endl;
                    if (curSize == k) break;

                    if (p1 == loutTmpSize) {
                        mid[cur][curSize++] = Lout[nin][p2];
                        p2++;
                        continue;
                    }
                    if (p2 == Lout[nin].size()) {
                        mid[cur][curSize++] = mid[loutTmp][p1];
                        p1++;
                        continue;
                    }

                    if (mid[loutTmp][p1].first < Lout[nin][p2].first) {
                        mid[cur][curSize++] = mid[loutTmp][p1];
                        p1++;
                    } else if (mid[loutTmp][p1].first > Lout[nin][p2].first) {
                        mid[cur][curSize++] = Lout[nin][p2];
                        p2++;
                    } else {
                        mid[cur][curSize++] = min(mid[loutTmp][p1], Lout[nin][p2]);
                        p1++;
                        p2++;
                    }
                }
                cur = 1 - cur;
                loutTmp = 1 - loutTmp;
                loutTmpSize = curSize;
            }
            Lout[tmp] = mid[loutTmp];
            Lout[tmp].resize(loutTmpSize);
        }
        cout << "label done. " << endl;
        vector<int>().swap(ind);
        vector<int>().swap(outd);
    }

    void copyIndex() {
        reachindex.resize(n);
        for (int i = 0; i < n; ++i) {
            reachindex[i].Lup = Lup[i];
            reachindex[i].Ldown = Ldown[i];
            reachindex[i].flabel1 = flabel1[i];
            reachindex[i].flabel2 = flabel2[i];

            reachindex[i].OutLimit = Lout[i].size();
            reachindex[i].InLimit = Lin[i].size() + reachindex[i].OutLimit;
            // reachindex[i].HLimit=hv.size()+reachindex[i].InLimit;
            reachindex[i].Label = new pair<int, int>[reachindex[i].InLimit];
            // LOG("i = {}, OutLimit = {}, InLimit = {}", i, reachindex[i].OutLimit, reachindex[i].InLimit);
            // LOG("reachindex size = {}", reachindex[i].OutLimit);
            for (int j = 0; j < Lout[i].size(); ++j) reachindex[i].Label[j] = Lout[i][j];
            for (int j = 0; j < Lin[i].size(); ++j) {
                // LOG("j = {}, index = {}, linsize = {}", j, reachindex[i].OutLimit + j, Lin[i].size());
                reachindex[i].Label[reachindex[i].OutLimit + j] = Lin[i][j];
            }
        }
    }
    //--------------------------label------------------------

    bool intersect(const int u1, const int u2) {
        int p1 = 0;
        int p2 = reachindex[u2].OutLimit;

        while (p1 < reachindex[u1].OutLimit && p2 < reachindex[u2].InLimit) {
            // cout << "u1 = " << u1 << " u2 = " << u2 << " p1 first = " << reachindex[u1].Label[p1].first
            //      << " p2 first = " << reachindex[u2].Label[p2].first
            //      << " p1 second = " << reachindex[u1].Label[p1].second
            //      << " p2 second = " << reachindex[u2].Label[p2].second << endl;
            if (reachindex[u1].Label[p1].first == reachindex[u2].Label[p2].first &&
                reachindex[u1].Label[p1].second <= reachindex[u2].Label[p2].second)
                return 1;
            else if (reachindex[u1].Label[p1].first == reachindex[u2].Label[p2].first) {
                p1++;
                p2++;
            } else if (reachindex[u1].Label[p1].first < reachindex[u2].Label[p2].first)
                p1++;
            else
                p2++;
        }
        return 0;
    }

    //----------------------query----------------------------

    bool IPtest(const int u, const int v) {
        int it1 = 0, it2 = 0;
        if (reachindex[u].OutLimit < 1) {
            return false;
        }
        if (reachindex[u].Label[reachindex[u].OutLimit - 1].first > reachindex[v].Label[0].first) {
            while (it1 != reachindex[u].OutLimit && it2 != reachindex[v].OutLimit) {
                if (reachindex[u].Label[it1].first == reachindex[v].Label[it2].first) {
                    if (reachindex[u].Label[it1].second > reachindex[v].Label[it2].second) return true;
                    it1++;
                    it2++;
                } else if (reachindex[u].Label[it1].first < reachindex[v].Label[it2].first) {
                    it1++;
                } else
                    return true;
            }
        }

        if (reachindex[v].InLimit > 0 and Lin[u].size() > 0 and
            reachindex[v].Label[reachindex[v].InLimit - 1].first <= reachindex[u].Label[reachindex[u].OutLimit].first)
            return false;

        it1 = reachindex[v].OutLimit;
        it2 = reachindex[u].OutLimit;

        while (it1 != reachindex[v].InLimit && it2 != reachindex[u].InLimit) {
            if (reachindex[v].Label[it1].first == reachindex[u].Label[it2].first) {
                if (reachindex[v].Label[it1].second < reachindex[u].Label[it2].second) return true;
                it1++;
                it2++;
            } else if (reachindex[v].Label[it1].first < reachindex[u].Label[it2].first) {
                it1++;
            } else
                return true;
        }

        return false;
    }

    int insidetest(int from, int to) {
        // new version
        int l = 0;
        int r = reachindex[from].OutLimit - 1;
        int m;
        while (l <= r) {
            m = (l + r) / 2;
            if (reachindex[from].Label[m].first == rank_order[toChain[to]]) {
                if (reachindex[from].Label[m].second <= toPos[to]) {
                    return 1;
                }
                return 2;
            } else if (reachindex[from].Label[m].first < rank_order[toChain[to]])
                l = m + 1;
            else
                r = m - 1;
        }
        return 0;
    }

    int countQuery = 0;

    bool query(int from, int to) {
        if (from < 0 || to < 0) return 0;
        if (from == to) return 1;
        countQuery++;
        flagVec[from] = mark;

        // cout << "from = " << from << " to = " << to << " from flabel1 = " << reachindex[from].flabel1<<" from flabel2
        // = " << reachindex[from].flabel2
        //      << " to flabel1 = " << reachindex[to].flabel1<<" to flabel2 = " << reachindex[to].flabel2 << endl;
        if (reachindex[from].flabel1 > reachindex[to].flabel1 || reachindex[from].flabel2 > reachindex[to].flabel2) {
            return 0;
        }
        if (IPtest(from, to)) {
            return 0;
        }

        // intersect
        if (intersect(from, to)) {
            return 1;
        }

        for (int i = v[from].size() - 1; i >= 0; --i) {
            int w = v[from][i];
            if (w == to) {
                return 1;
            }

            if (reachindex[w].Ldown < reachindex[to].Ldown && (flagVec[w] != mark))
            //	if (flag[w]!=mark)
            {
                if (query(w, to)) {
                    return 1;
                }
            }
        }
        return 0;
    }
    bool naiveQuery(int from, int to) {
        // flag[from] = mark;
        // if (from == to) return 1;
        // for (int i = 0; i < v[from].size(); ++i) {
        //     int w = v[from][i];
        //     if (flag[w] != mark) {
        //         if (naiveQuery(w, to)) return 1;
        //     }
        // }
        return 0;
    }
    void query(string s) {
        int *flag = new int[n];
        memset(flag, 0, sizeof(int) * n);

        // read query
        vector<int> from;
        vector<int> to;
        fstream in;
        in.open(s.c_str(), ios::in);
        int t1, t2;
        while (in >> t1 >> t2) {
            from.push_back(t1);
            to.push_back(t2);
        }
        in.close();

        cout << "query number: " << from.size() << endl;

        // timer.start();
        countQuery = 0;
        // begin query
        int queryNum = from.size();
        int reachNum = 0;
        for (int i = 0; i < queryNum; ++i) {
            bool arrive;
            if (from[i] == to[i]) arrive = 1;
            // else if (Ldown[from[i]]>=Ldown[to[i]] /*|| Lup[from[i]]<=Lup[to[i]]*/)
            // arrive = 0; else if (flabel1[from[i]]>flabel1[to[i]] ||
            // flabel2[from[i]]>flabel2[to[i]]) arrive = 0; else if
            // (flabel1[from[i]]>flabel1[to[i]]) arrive = 0;
            else if (reachindex[from[i]].Ldown >= reachindex[to[i]].Ldown)
                arrive = 0;
            else if (reachindex[from[i]].flabel1 > reachindex[to[i]].flabel1 ||
                     reachindex[from[i]].flabel2 > reachindex[to[i]].flabel2)
                arrive = 0;
            else if (intersect(from[i], to[i]))
                arrive = 1;
            else {
                mark = i + 1;
                arrive = query(from[i], to[i]);
            }
            reachNum += arrive;
            // if (arrive) cout << from[i] << " " << to[i] << endl;
        }

        // timer.stop();
        cout << "reach number: " << reachNum << endl;
        cout << endl;
        cout << "countQuery= " << countQuery << endl;
        //   cout << "QUERY_TIME = " << //timer.GetRuntime() * 1000 << " ms" <<
        //   endl;

        delete[] flag;
    }
};
}  // namespace topchain
#endif
