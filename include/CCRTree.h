#ifndef CCRTREE_H
#define CCRTREE_H
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "graph.h"

#define MERGE_TO_UINT64(high, low) \
    ((static_cast<uint64_t>(high) << 32) | static_cast<uint64_t>(low))

#define GET_HIGH_LOW_32BITS(value, high, low)       \
    do {                                            \
        high = static_cast<int>((value) >> 32);     \
        low = static_cast<int>((value)&0xFFFFFFFF); \
    } while (0)

#define MAX_VAL 10000000
#define MIN_VAL -10000000
namespace CCR {
using namespace ptree;

struct Node {
    int id = -1;
    bool in_tree = false;
};

/**
first = 祖先最近， second = 子孙最近，
value == -1 代表无效
*/
using black_table = std::pair<vector<int>, vector<int>>;

using section_table = std::pair<int32_t, int32_t>;

using NodeList = std::vector<Node>;
struct Tree {
    GRA graph_;
    NodeList vl_;
    int vsize_;
    int root_ = -1;
    int froot { -1 };
    int m = { 1 }; // 用于确定区间
    bool isinverse = false;
    // 每个黑色节点最近的红色祖先节点、红色子孙节点
    std::vector<black_table> black_table_;
    // std::vector<bool> is_cross_node_;
    // TODO：label标签，每个红色节点维护一个区间
    std::vector<section_table> section_table_;

    std::vector<std::pair<int, int>> edge_pair_;

    // 0 黑色节点， 1 CrossEdge的起点， 2 CrossEdge的终点, 3 既是终点又是起点
    std::vector<int> isStart_;

    Tree(int n)
    {
        graph_ = GRA(n, In_OutList());
        vl_ = NodeList(n);
        // color_ = std::vector<bool>(n, false);
        isStart_ = std::vector<int>(n, 0);
        black_table_ = std::vector<black_table>(
            n, std::make_pair(vector<int>(), vector<int>()));
        section_table_ = std::vector<section_table>(n, std::make_pair(-1, -1));
        // is_cross_node_ = std::vector<bool>(n, false);
    }
    Tree() = default;

    inline void
    set_num(int n)
    {
        graph_ = GRA(n, In_OutList());
        vl_ = NodeList(n);
    }
    int
    num_vertices()
    {
        return vl_.size();
    }
    EdgeList
    in_edges(int trg)
    {
        return graph_[trg].inList;
    }
    EdgeList
    out_edges(int trg)
    {
        return graph_[trg].outList;
    }
    int
    in_degree(int trg)
    {
        return graph_[trg].inList.size();
    }
    int
    out_degree(int trg)
    {
        return graph_[trg].outList.size();
    }
    inline bool
    has_node(int num)
    {
        if (num >= vsize_)
            return false;
        if (graph_[num].inList.empty() && graph_[num].outList.empty())
            return false;
        return true;
    }
};

/**
 * 某一层的返回值
 */
struct res_type {
    // 一棵树
    bool reverse_ = false;
    // 每个节点最近的红色祖先节点、红色子孙节点
    std::vector<black_table> black_table_;

    // label标签，每个节点维护一个区间
    std::vector<section_table> section_table_;

    // 记录当前这层的节点编号对应到下一层是多少 (当前层，下一层)
    std::unordered_map<int32_t, int32_t> vertex_table_;
    res_type(res_type&& other) noexcept = default;
    res_type() = default;
};

struct SectionInfo {
    int first;
    int second;
    int id; // 节点id
    SectionInfo(int f, int s, int i)
        : first(f)
        , second(s)
        , id(i)
    {
    }
    SectionInfo() { first = -1, second = -1, id = -1; }
};

// pass by value
std::list<res_type> GetLabels(Graph gGraph, int max_num, double percent);
std::vector<int> FindRootNode(Graph& graph);
std::vector<int> FindRootNode(Graph& graph, bool reverse_);

void BuildTree(Tree& tree, Graph& graph, int node);
void BuildTree(Tree& tree, Graph& graph, int node, bool reverse_);
// void BuildNewGraph(Tree &tree_, Graph &origin_graph_, Graph &graph,
// std::unordered_map<int, int> &vertex_table_, std::unordered_set<uint64_t>
// &covered_edges);
void BuildNewGraph(Tree& tree_, Graph& origin_graph_, Graph& graph,
    std::unordered_map<int, int>& vertex_table_,
    std::unordered_set<uint64_t>& covered_edges,
    bool reverse_);

// 正向树用于判断first能否到达second，反向树用于判断second能否到达first
bool CanReach(const section_table& first_elem,
    const section_table& second_elem);

std::unordered_set<uint64_t> DFSBuild(Tree& tree_, Graph& graph_);
std::unordered_set<uint64_t> DFSBuild(Tree& tree_, Graph& graph_,
    bool reverse_);

void MarkColor(Tree& tree_, Graph& graph_,
    std::unordered_set<uint64_t>& cross_edges);
void MarkColor(Tree& tree_, Graph& graph_,
    std::unordered_set<uint64_t>& cross_edges, bool reverse_);

vector<int> AssignTP(Tree& tree_, int node_, vector<int> parent_,
    Graph& graph_,
    std::unordered_set<uint64_t>& cross_edges);
vector<int> AssignTP(Tree& tree_, int node_, vector<int> parent_,
    Graph& graph_, std::unordered_set<uint64_t>& cross_edges,
    bool reverse_);
} // namespace CCR
#endif // CCRTREE_H