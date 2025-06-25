#ifndef _PATH_TREE_H
#define _PATH_TREE_H

#include <cstdint>
#include "DWGraphUtil.h"
#include "DataComp.h"
#include "GraphUtil.h"

// test switch
#define _TEST_
namespace ptree {
class PathTree {
   public:
    Graph &g;
    DWGraph pg;
    Graph ng;  // equivalent weight graph
    DWGraph branch;
    Graph newbranch;
    int maxeid;

    vector<int> nextVertex;
    //		vector<set<int> > out_uncover;
    vector<vector<int> > out_uncover;
    map<int, vector<int> > comp_table;
    vector<vector<int> > pathMap;
    vector<int> grts;  // graph reverse topological sort
    int **labels;
    bool effective;

    map<pair<int, int>, bool> tcm;  // for test
    struct timeval after_time, before_time;
    float run_time;

   public:
    PathTree(Graph &graph);
    PathTree(Graph &graph, vector<int> ts);
    ~PathTree();

    uint64_t getIndexSize() {
        uint64_t res = 0;
        res += g.num_vertices() * 3;
        // out_uncover
        for (auto &t : out_uncover) {
            res += t.size();
        }
        return res;
    }
    void buildWeightPathGraph(int type);
    void buildWeightPathGraph_Pred();  // update weight by pred size
    void buildEquGraph();
    void createLabels(int type, ifstream &cfile, bool compress);
    void displayLabels();
    //	void pathDFS(int vid, int& order, int& first_order);
    void pathDFS(int vid, int &order, int &first_order, vector<bool> &visited);
    void transform(DWGraph dg, Graph &graph);
    void readPathMap(ifstream &cfile);

    void compute_tcm();
    bool reach(int src, int trg);
    bool reach_dc(int src, int trg);
    bool test_reach(int src, int trg);
    bool test_reach_dc(int src, int trg);
    void index_size(int *ind_size);
    void insertSet(set<int> &s1, set<int> &s2);

    void mergeVector(vector<int> &v1, vector<int> &v2);

    void buildEquEdgeset(map<int, set<int> > &pathtopo, Graph &equgraph);

    double cover_ratio();
    double compress_ratio();
};
}  // namespace ptree
#endif
