#ifndef _GRAPH_UTIL_H_
#define _GRAPH_UTIL_H_

#include <sys/time.h>

#include "DWGraph.h"
#include "graph.h"

namespace ptree {
class GraphUtil {
   public:
    static void dfs(Graph &g, int vid, vector<int> &preorder, vector<int> &postorder, vector<bool> &visited);
    static void topological_sort(Graph g, vector<int> &ts);
    static void transitive_closure(Graph g, Graph &tc);
    static void tarjan(Graph &g, int vid, int &index, unordered_map<int, pair<int, int> > &order, vector<int> &sn,
                       multimap<int, int> &sccmap, int &scc);
    static void mergeSCC(Graph &g, int *on, vector<int> &ts);
    static void findTreeCover(Graph g, Graph &tree);
    static void findTreeCover(Graph g, Graph &tree, vector<set<int> > &pred);
    static void findTreeCover(Graph &g, Graph &tree, vector<set<int> > &pred, vector<int> &ts);
    static void compute_pred(Graph g, vector<set<int> > &predMap);
    static void findTreeCoverL(Graph g, Graph &tree);
    static void traverse(Graph &tree, int vid, int &pre_post, vector<bool> &visited);
    static void pre_post_labeling(Graph &tree);
    static void pathDecomposition(Graph &g, vector<vector<int> > &pathMap);
    static void pathDecomposition(Graph &g, vector<vector<int> > &pathMap, vector<int> ts);
    static void treePathDecomposition(Graph tree, Graph &g, vector<vector<int> > &pathMap);

    static void genRandomGraph(int n, double c, char *filename);
};
}  // namespace ptree

#endif
