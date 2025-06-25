#ifndef PATITION_H
#define PATITION_H

#include <boost/functional/hash.hpp>
#include <cstdint>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph.h"
#include "pathtree/PathTree.h"
#include "topchain/topChain.h"
#include "utils.h"

template <typename T>
std::size_t hash_vector(const std::vector<T> &vec) {
    std::size_t seed = 0;
    for (const auto &i : vec) {
        boost::hash_combine(seed, i);
    }
    return seed;
}

class Partitioner;
class Block {
   private:
   public:
    Partitioner *global_partition;

    // needs new ds to store extraNodes
    // key = 边界点的lid， value = 等价类编号
    // note: 该点被哪些等价类包含
    unordered_map<int, unordered_set<int>> extraTable;
    // 被单拿出来的node
    // key = 边界图上的lid, value = eqclassGraph节点编号
    unordered_map<int, int> selectedNodes;
    void SplitEqNode(int threshold);
    void ExemptBorderNode(int threshold, ptree::Graph &borderGraph, unordered_map<int, int> &s);

    int getEqGraphID(int gid);
    // key = 等价类编号, value = 虚拟节点的id
    // 第i个等价类 对应的虚拟节点id
    std::unordered_map<int, int> sout;
    std::unordered_map<int, int> sin;
    void allocateID(ptree::Graph &borderGraph);
    vector<unordered_set<int>> bufferO;
    vector<unordered_set<int>> bufferI;
    // 标识等价类是否被遍历过
    std::vector<int> eqClassVisitedO;
    std::vector<int> eqClassVisitedI;
    void setEqVisitedO(int gid, int idx) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderOutIdx;
        this->eqClassVisitedO[eqId] = idx;
    }
    void setEqVisitedI(int gid, int idx) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderInIdx;
        this->eqClassVisitedI[eqId] = idx;
    }

    int getEqVisitedO(int gid) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderOutIdx;
        return this->eqClassVisitedO[eqId];
    }

    int getEqVisitedI(int gid) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderInIdx;
        return this->eqClassVisitedI[eqId];
    }

    uint64_t index_size{0};
    std::unordered_map<int, int> global2local;
    // std::vector<int> global2local;
    ptree::Graph localGraph;
    int crossOutEdge{0};
    int crossInEdge{0};
    int blockID{-1};
    // PathTree pt;
    unique_ptr<ptree::PathTree> pt;

    // 存的是边界节点的id, 不是gid也不是local id
    std::vector<std::vector<int>> borderClassOut;
    std::vector<std::vector<int>> borderClassIn;
    std::vector<std::vector<int>> secondEqNodeRecord;
    void buildBorderLink(ptree::Graph &borderGraph, std::unordered_map<int, int> &borderGlobal2Local, ptree::Graph &g,
                         std::vector<Block *> &blocks);
    void buildBorderLink2(ptree::Graph &borderGraph, ptree::Graph &eqGraph,
                          std::unordered_map<int, int> &borderGlobal2Local, ptree::Graph &g,
                          std::vector<Block *> &blocks);
    void buildBorderLink3(ptree::Graph &eqGraph, ptree::Graph &g, std::vector<Block *> &blocks);
    // 利用另一个方向的等价类来打补丁
    void buildBorderPatch(ptree::Graph &borderGraph, ptree::Graph &eqGraph,
                          std::unordered_map<int, int> &borderGlobal2Local, ptree::Graph &g,
                          std::vector<Block *> &blocks);
    void allocateBaseEqClass(ptree::Graph &eqGraph,ptree::Graph &graph);

    void fetchClass(std::unordered_map<int, int> &borderGlobal2Local);
    void fetchClassByPT(std::unordered_map<int, int> &borderGlobal2Local);
    void fetchClassByBorderBFS(std::unordered_map<int, int> &borderGlobal2Local);
    void fetchClassByBorderBFSOut(ptree::Graph &eqGraph,ptree::Graph &graph);
    void fetchClassByBorderBFSIn(ptree::Graph &eqGraph,ptree::Graph &graph);
    uint64_t getIndexSize() { return pt->getIndexSize(); }
    std::unordered_set<int> borderOut;
    std::unordered_set<int> borderIn;
    int getInBorder(int gid, std::vector<int> &que) {
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderInIdx;
        if (idx < 0) {
            return 0;
        }
        std::copy(borderClassIn[idx].begin(), borderClassIn[idx].end(), que.begin());
        return borderClassIn[idx].size();
    }

    int getOutBorder(int gid, std::vector<int> &que) {
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderOutIdx;
        if (idx < 0) {
            return 0;
        }
        std::copy(borderClassOut[idx].begin(), borderClassOut[idx].end(), que.begin());

        return borderClassOut[idx].size();
    }

    int getInBorder2(int gid) {
        // 先看等价类有没有被访问过
        // 再看每个点
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderInIdx;
        if (idx < 0) {
            return -1;
        }
        return idx;
    }

    int getOutBorder2(int gid) {
        // 先看等价类有没有被访问过
        // 再看每个点
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderOutIdx;

        if (idx < 0) {
            return -1;
        }
        return idx;
    }

    void dumpGraph() {
        for (int i = 0; i < global2local.size(); i++) {
            cout << i << " : " << global2local[i] << "\n";
        }
        std::cout << "edge = \n";
        localGraph.dumpGraph();
    }
    void setBlockID(int id) { blockID = id; }
    int getBlockID() const { return blockID; }
    int getSize() const { return localGraph.num_vertices(); }
    int getCrossOutEdge() const { return crossOutEdge; }
    void setBorderType(int gid, int type) {
        auto lid = this->global2local[gid];
        localGraph.setVertexBorderFlag(lid, type);
        if (type == 1) {
            this->borderOut.insert(lid);
        } else if (type == 2) {
            this->borderIn.insert(lid);
        }
    }
    void addNode(ptree::Vertex v) {
        auto res = localGraph.addVertex(v);
        if (res) {
            // 添加成功
            int len = localGraph.num_vertices();
            global2local[v.id] = len - 1;
        } else {
            // do nothing
        }
    }

    void addEdge(int gfirst, int gsecond) {
        if (global2local[gfirst] < 0 or global2local[gsecond] < 0) {
            ERROR("add edge error, not in the same blk");
            exit(-1);
        }
        int localA = global2local[gfirst];
        int localB = global2local[gsecond];
        localGraph.addEdge(localA, localB);
    }

    void increaseCrossInedge() { crossInEdge++; }
    void increaseCrossOutedge() { crossOutEdge++; }
    void decreaseCrossInedge() { crossInEdge--; }
    void decreaseCrossOutedge() { crossOutEdge--; }
    virtual bool query(int src, int dst) {
        auto s = global2local[src];
        auto d = global2local[dst];
        return pt->reach(s, d);
    }
    virtual void runReachability();
};

class PTBlock : public Block {
    // 父类默认使用pathtree，子类不做任何修改
};

class TopChainBlock : public Block {
    void runReachability() override;
    bool query(int src, int dst) override;
    topchain::TopChainMeta topChainMeta;
    topchain::GraphT gt;
};

class Partitioner {
   public:
    Partitioner() = default;

    Partitioner(int maxs, int mins) : maxBlkSize(maxs), minBlkSize(mins) {}
    Partitioner(int numBlk) : minBlkSize(numBlk) {}

    ~Partitioner() = default;

    void runPartition(ptree::Graph &graph, const string &algo);
    uint64_t getIndexSize() {
        uint64_t res = this->index_size;
        for (auto t : blocks) {
            res += t->getIndexSize();
        }
        return res;
    }

    uint64_t bitElemNumber;

    // std::unordered_map<int, std::vector<int>> borderClassIn;
    // std::unordered_map<int, std::vector<int>> borderClassOut;

    // 看下全局的sout是否可行
    std::unordered_map<size_t, int> GlobalRecordTableOut;
    std::unordered_map<size_t, int> GlobalRecordTableIn;
    std::unordered_map<size_t, int> SingleNode2EqGraphNode;

    // 等价类编号+eqGraph节点编号
    std::unordered_map<int, int> GlobalSout;
    std::unordered_map<int, int> GlobalSin;

    // tools
    void PrintPartitionRatio(ptree::Graph &graph) const;
    void ComputeBorder(ptree::Graph &graph);

    void runLocalReachability();
    void runBorderReachability(ptree::Graph &graph);
    void runQueryWithBfs(BiGraph &bg, ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo,
                         std::vector<int> &queryRes);
    // 处理边界图标签
    void LabelGraph(ptree::Graph &graph);
    int dfsLabel(ptree::Graph &graph, int node, int label);

    std::vector<Block *> blocks;
    // std::vector<int> borderGlobal2Local;
    // ptree::Graph borderGraph;
    ptree::Graph eqClassGraph;
    uint64_t index_size{0};
    int threshold{0};
    topchain::TopChainMeta topChainMeta_partitioner;
    topchain::GraphT gt_partitioner;

   private:
    vector<int> que1;
    vector<int> que2;

    vector<int> visited_outter;
    int maxBlkSize{0};
    int minBlkSize{0};
    unique_ptr<ptree::PathTree> pt;
    int getInBorder(int dst, ptree::Graph &graph);
    int getOutBorder(int src, ptree::Graph &graph);
    bool runBiBFS(int q1hi, int q2hi, int &sum, int idx);
    bool runBiBFSBy2Hop(int q1hi, int q2hi, int &sum, int idx);
    bool runIndexQuery(int src, int dst);
    bool runPTQuery(int src, int dst);
};

// for oreach

#endif
