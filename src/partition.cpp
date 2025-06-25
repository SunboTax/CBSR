#include "partition.h"

#include <bits/types/locale_t.h>
#include <bits/types/time_t.h>

// #include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "graph.h"
#include "utils.h"
#define InnerBFS
// #define PATHTREE
#define BORDERBFS
// #define PathTreeBorder
// #define NOEQCLASS  // 不使用等价类
// #define MTHREAD
#define TWOHOPBORDER
// #define INNODE
// #define NEWSETTING
#define MERGE_TO_UINT64(high, low) ((static_cast<uint64_t>(high) << 32) | static_cast<uint64_t>(low))

void setBitAtPosition(std::vector<uint64_t> &bits, uint64_t pos) {
    // 计算pos在哪个uint64_t元素中
    uint64_t elementIndex = pos / 64;

    // 确保vector有足够的元素
    if (elementIndex >= bits.size()) {
        bits.resize(elementIndex + 1, 0);
    }

    // 计算pos在该元素中的具体位置
    uint64_t bitPosition = pos % 64;

    // 使用位运算将对应位置的bit置为1
    bits[elementIndex] |= (1ULL << bitPosition);
}

struct Element {
    int value;
    size_t vecIndex;   // 当前元素来自哪个vector
    size_t elemIndex;  // 当前元素在其vector中的位置

    bool operator>(const Element &other) const { return value > other.value; }
};

std::vector<int> mergeAndUniqueUsingMinHeap(const std::vector<std::vector<int>> &vectors) {
    // 使用lambda函数来定义比较逻辑，构建一个最小堆
    auto comp = [](const Element &a, const Element &b) { return a.value > b.value; };
    std::priority_queue<Element, std::vector<Element>, decltype(comp)> minHeap(comp);

    // 初始化堆，每个vector的第一个元素加入堆
    for (size_t i = 0; i < vectors.size(); ++i) {
        if (!vectors[i].empty()) {
            minHeap.push({vectors[i][0], i, 0});
        }
    }

    std::vector<int> result;
    int last_added = INT_MIN;  // 用于去重，初始化为一个不可能的值

    while (!minHeap.empty()) {
        Element current = minHeap.top();
        minHeap.pop();

        // 只有当当前值与最后加入的值不同，才添加到结果中
        if (result.empty() || current.value != last_added) {
            result.push_back(current.value);
            last_added = current.value;
        }

        // 如果当前vector还有元素，把下一个元素加入堆
        if (current.elemIndex + 1 < vectors[current.vecIndex].size()) {
            minHeap.push({vectors[current.vecIndex][current.elemIndex + 1], current.vecIndex, current.elemIndex + 1});
        }
    }

    return result;
}
int cnt1 = 0;
int cnt2 = 0;
int cnt3 = 0;
void Partitioner::PrintPartitionRatio(ptree::Graph &graph) const {
    if (this->blocks.empty()) {
        WARN("empty blocks");
        return;
    }
    double ce = 0;
    int partitionMaxIdx = -1;
    auto size = this->blocks.size();
    vector<vector<int>> matrix(size, vector<int>(size, 0));
    vector<int> blk_size(size, 0);
    for (auto node : graph.vertices()) {
        blk_size[node.partition]++;
        partitionMaxIdx = max(partitionMaxIdx, node.partition);
        for (auto edge : graph.out_edges(node.id)) {
            if (graph[edge].partition != node.partition) {
                matrix[node.partition][graph[edge].partition] = 1;
                ce += 1;
            }
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (matrix[i][j] == 1 and matrix[j][i] == 1) {
                WARN("{} and {} has bidirectional edge", i, j);
            }
        }
    }

    LOG("edge cut ratio = {}", ce / graph.num_edges());
    for (auto &blk : this->blocks) {
        LOG("block {} size = {}", blk->getBlockID(), blk->getSize());
    }
}

void Partitioner::runPartition(ptree::Graph &graph, const string &algo) {
    CCR::Timer timer;
    timer.ticker();
    LOG("begin Run Partition");
    auto gsize = graph.num_vertices();
    LOG("graph vertex size = {}, edge size = {}, block size = {}", gsize, graph.num_edges(), this->minBlkSize);
    this->blocks = vector<Block *>(this->minBlkSize);
    int maxVertexSize = (gsize + this->minBlkSize - 1) / this->minBlkSize;
    for (int i = 0; i < minBlkSize; i++) {
        if (algo == "PathTree") {
            this->blocks[i] = new Block();
        } else if (algo == "TopChain") {
            this->blocks[i] = new TopChainBlock();
        } else {
            ERROR("未定义 {}", algo);
            exit(0);
        }
        this->blocks[i]->global_partition = this;
    }

    auto cmp = [](ptree::Vertex &first, ptree::Vertex &second) {
        if (first.start == second.start) {
            return first.end > second.end;
        }
        return first.start > second.start;
    };
    auto checkCondition = [&maxVertexSize, this](Block *blk) {
        if (blk->getSize() >= maxVertexSize) {
            DEBUG("partition {} 超过阈值 {}", blk->getSize(), maxVertexSize);
            return true;
        }
        return false;
    };

    priority_queue<ptree::Vertex, std::vector<ptree::Vertex>, decltype(cmp)> que(cmp);
    for (auto &node : graph.vertices()) {
        if (graph.in_degree(node.id) == 0) {
            graph[node.id].visited = true;
            que.push(node);
        }
    }

    int partitionIdx = 0;
    while (!que.empty()) {
        blocks[partitionIdx]->setBlockID(partitionIdx);
        ptree::Vertex v = que.top();
        que.pop();

        graph[v.id].partition = partitionIdx;
        graph[v.id].gid = v.id;
        blocks[partitionIdx]->addNode(graph[v.id]);
        graph[v.id].lid = blocks[partitionIdx]->localGraph.num_vertices() - 1;
        for (auto t : graph.out_edges(v.id)) {
            if (not graph[t].visited) {
                que.push(graph[t]);
                graph[t].visited = true;
                // tmpBlk.crossOutEdge += graph.out_degree(v.id);
            }
        }
        auto checkRes = checkCondition(blocks[partitionIdx]);
        if (checkRes) {
            partitionIdx++;
        }
    }

    // this->borderGlobal2Local = vector<int>(graph.num_vertices(), -1);
    double cutedge = 0;
    // 处理边界节点
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (auto t : graph.out_edges(i)) {
            int idx = graph[i].partition;
            if (idx == -1 or graph[t].partition == -1) {
                ERROR("自环");
                exit(1);
            }
            if (graph[i].partition == graph[t].partition) {
                blocks[idx]->addEdge(i, t);
            }
        }
    }

    double sum = 0;
    timer.ticker();
    LOG("partition end, time = {}, cut size = {}, cut ratio = {}", timer.get_last_consuming(), cutedge,
        cutedge / graph.num_edges());
}

void Partitioner::runLocalReachability() {
#ifdef MTHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < blocks.size(); i++) {
        blocks[i]->runReachability();
    }
}

void Partitioner::runBorderReachability(ptree::Graph &graph) {
    // 预先分配好节点
    LOG("before allocate, node size = {}", this->eqClassGraph.num_vertices());
    for (int i = 0; i < blocks.size(); i++) {
        // blocks[i]->allocateID(this->eqClassGraph);
    }
    LOG("after allocate, node size = {}", this->eqClassGraph.num_vertices());

    // // 看一眼边界图的映射
    // LOG("mapping");
    // for (auto p : this->borderGlobal2Local) {
    //     cout << p.first << " --> " << p.second << "\n";
    // }
    // 虚拟in节点只能出现一次
    for (size_t i = 0; i < blocks.size(); i++) {
        // blocks[i]->ExemptBorderNode(this->threshold, this->eqClassGraph, s);
    }

    LOG("allocate EQclass Border Node done, #Nodes = {}", this->eqClassGraph.num_vertices());

    // 1. 构建块内部的联通关系
#ifdef MTHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < blocks.size(); i++) {
        // blocks[i]->buildBorderLink2(this->borderGraph, this->eqClassGraph, this->borderGlobal2Local, graph,
        //                             this->blocks);
        // blocks[i]->buildBorderLink(this->borderGraph, this->borderGlobal2Local, graph, this->blocks);
        blocks[i]->buildBorderLink3(this->eqClassGraph, graph, this->blocks);
        // blocks[i]->buildBorderPatch(this->borderGraph, this->eqClassGraph, this->borderGlobal2Local, graph,
        //                             this->blocks);
    }

    // 连接外部关系 out eqclass->in eqclass
    LOG("[AFTER] eq vertex size = {}, edge size = {}", eqClassGraph.num_vertices(), eqClassGraph.num_edges());
    this->index_size += eqClassGraph.num_edges();
    this->index_size += eqClassGraph.num_vertices();
#ifdef PathTreeBorder
    // 2. 用paththee计算边界点的可达性
    LOG("Border #construction begin");
    auto &g = this->eqClassGraph;
    int gsize = g.num_vertices();
    vector<int> reverse_topo_sort;
    ptree::GraphUtil::topological_sort(g, reverse_topo_sort);
    cout << "#DAG vertex size:" << g.num_vertices() << "\t#DAG edges size:" << g.num_edges() << endl;

    pt = make_unique<ptree::PathTree>(g, reverse_topo_sort);
    ifstream cfile;
    pt->createLabels(4, cfile, false);
    LOG("Border #construction done");
#endif
    LOG("cnt1 = {}, cnt2 = {}, cnt3 = {}", cnt1, cnt2, cnt3);

#ifdef TWOHOPBORDER
    // 需要给每个点加上一个时间戳，正常排序即可，无特殊意义
    if (this->eqClassGraph.num_vertices() == 0) {
        return;
    }
    LabelGraph(this->eqClassGraph);
    this->gt_partitioner = topchain::transToScc(this->eqClassGraph);
    this->topChainMeta_partitioner.readGraph(this->gt_partitioner);
    this->topChainMeta_partitioner.readChainInfo(this->gt_partitioner);
    this->topChainMeta_partitioner.findChain();
    this->topChainMeta_partitioner.feline();
    this->topChainMeta_partitioner.handle();
    this->topChainMeta_partitioner.copyIndex();
    this->topChainMeta_partitioner.flag = new int[topChainMeta_partitioner.n];
#endif
}

// 需要使用拓扑排序进行标签赋值
void Partitioner::LabelGraph(ptree::Graph &graph) {
    int start_time = 1;
    auto vec = graph.topological_sort();
    for (auto node : vec) {
        graph[node].start = start_time;
        start_time += 10;
    }
}

void Block::allocateID(ptree::Graph &borderGraph) {
#ifdef INNODE
    int osize = borderClassIn.size();
    int insize = borderClassOut.size();
#else
    int osize = borderClassOut.size();
    int insize = borderClassIn.size();
#endif
    // if (osize == 0) {
    //     return;
    // }
    int borderGraphSize = borderGraph.num_vertices();
    int lastOneID = borderGraphSize + osize - 1;
    borderGraph.addVertex(lastOneID);
    for (int i = 0; i < osize; i++) {
        this->sout.emplace(make_pair(i, i + borderGraphSize));
    }
    borderGraphSize = borderGraph.num_vertices();
    lastOneID = borderGraphSize + insize - 1;
    borderGraph.addVertex(lastOneID);
    for (int i = 0; i < insize; i++) {
        this->sin.emplace(make_pair(i, i + borderGraphSize));
    }
}

int Block::getEqGraphID(int gid) {
    if (gid == -1) {
        exit(-9);
    }
    auto lid = this->global2local[gid];
    int value = -1;
#ifdef INNODE
    value = localGraph[lid].borderInIdx;
#else
    value = localGraph[lid].borderOutIdx;
#endif
    if (this->sout.find(value) != this->sout.end()) {
        return this->sout[value];
    }
    return -1;
}

void Block::buildBorderLink2(ptree::Graph &borderGraph, ptree::Graph &eqClassGraph,
                             std::unordered_map<int, int> &borderGlobal2Local, ptree::Graph &g,
                             std::vector<Block *> &blocks) {
    std::unordered_set<uint64_t> s;
#ifdef INNODE
    auto myPartition = this->blockID;
    for (int i = 0; i < borderClassIn.size(); i++) {
        auto src = this->sout[i];
        for (auto Inode : bufferI[i]) {
            auto dst = sout[Inode];
            eqClassGraph.addEdge(dst, src);
        }
        auto borders = this->borderClassIn[i];
        for (auto bnode : borders) {
            for (auto node : borderGraph.in_edges(bnode)) {
                auto dst_partition = borderGraph[node].partition;
                auto dst = blocks[dst_partition]->getEqGraphID(borderGraph[node].gid);
                if (dst < 0) {
                    continue;
                }
                uint64_t mergeValue = MERGE_TO_UINT64(src, dst);
                if (s.count(mergeValue) > 0) {
                    continue;
                }

                s.insert(mergeValue);
                eqClassGraph.addEdge(dst, src);
            }
        }
    }
#else
    // 所有的out衍生虚拟节点全部连起来
    for (int i = 0; i < borderClassOut.size(); i++) {
        auto &borders = this->borderClassOut[i];
        auto src = this->sout[i];
        for (auto Onode : bufferO[i]) {
            auto dst = sout[Onode];
            eqClassGraph.addEdge(src, dst);
        }
        for (auto bnode : borders) {
            for (auto node : borderGraph.out_edges(bnode)) {
                auto dst_partition = borderGraph[node].partition;

                auto dst = blocks[dst_partition]->getEqGraphID(borderGraph[node].gid);
                if (dst < 0) {
                    continue;
                }
                uint64_t mergeValue = MERGE_TO_UINT64(src, dst);
                if (s.count(mergeValue) > 0) {
                    continue;
                }

                s.insert(mergeValue);
                eqClassGraph.addEdge(src, dst);
            }
        }
    }
#endif
}

void Block::buildBorderPatch(ptree::Graph &borderGraph, ptree::Graph &eqClassGraph,
                             std::unordered_map<int, int> &borderGlobal2Local, ptree::Graph &g,
                             std::vector<Block *> &blocks) {
#ifdef INNODE
    // for (int i = 0; i < borderClassOut.size(); i++) {
    //     auto &borders = this->borderClassOut[i];
    //     // 拿到当前另一个方向的虚拟编号
    //     auto src = this->sin[i];
    //     for (auto bnode : borders) {
    //         for (auto node : borderGraph.out_edges(bnode)) {
    //             auto dst_partition = borderGraph[node].partition;
    //             auto dst = blocks[dst_partition]->getEqGraphID(borderGraph[node].gid);
    //             if (dst < 0) {
    //                 continue;
    //             }
    //             uint64_t mergeValue = MERGE_TO_UINT64(src, dst);
    //             if (s.count(mergeValue) > 0) {
    //                 continue;
    //             }
    //             s.insert(mergeValue);

    //             eqClassGraph.addEdge(src, dst);
    //         }
    //     }
    // }
    for (int i = 0; i < borderClassOut.size(); i++) {
        auto &borders = this->borderClassOut[i];
        std::unordered_set<uint64_t> s;
        vector<int> buffer;
        for (auto bnode : borders) {
            for (auto node : borderGraph.out_edges(bnode)) {
                auto dst_partition = borderGraph[node].partition;
                auto dst = blocks[dst_partition]->getEqGraphID(borderGraph[node].gid);
                if (s.count(dst) > 0) {
                    continue;
                }
                s.insert(dst);
                buffer.push_back(dst);
            }
        }
        this->secondEqNodeRecord.push_back(buffer);
    }
#else

#endif
}

void Block::buildBorderLink3(ptree::Graph &eqClassGraph, ptree::Graph &g, std::vector<Block *> &blocks) {
    auto &table = global_partition->GlobalRecordTableIn;
    for (auto inNode : this->borderIn) {
        auto gid = localGraph[inNode].gid;
        if (global_partition->SingleNode2EqGraphNode.find(gid) != global_partition->SingleNode2EqGraphNode.end()) {
            // in边界点和自己的等价类连起来
            auto src = global_partition->SingleNode2EqGraphNode[gid];
            auto dst = localGraph[inNode].borderOutIdx;
            if (src >= 0 and dst >= 0) {
                eqClassGraph.addEdge(src, dst);
                cnt2++;
            }
        }
    }
}

void Block::ExemptBorderNode(int threshold, ptree::Graph &eqClassGraph, unordered_map<int, int> &s) {
    // 找到需要豁免的边界点，并且分配eqClass节点
    for (auto &it : this->extraTable) {
        // it.first = 边界图lid
        if (it.second.size() > threshold) {
            if (s.find(it.first) != s.end()) {
                // this->selectedNodes.insert({it.first, s[it.first]});
            } else {
                // 分配一个id并记录下来
                auto id = eqClassGraph.num_vertices();
                eqClassGraph.addVertex(id);
                // this->selectedNodes.insert(make_pair(it.first, id));
                s.insert({it.first, id});
                // LOG("exempt code = {}", it.first);
            }
        }
    }
}

void Block::SplitEqNode(int threshold) {}

void Block::buildBorderLink(ptree::Graph &borderGraph, std::unordered_map<int, int> &borderGlobal2Local,
                            ptree::Graph &g, std::vector<Block *> &blocks) {
    // 只构建虚拟节点的连接
    // key = oidx, value = 虚拟节点的id
    // std::unordered_map<int, int> s;
    std::unordered_set<int> s;
    int size = localGraph.num_vertices();
#ifdef INNODE
    for (int i = 0; i < size; i++) {
        // 连接该节点与相关的out节点
        if (localGraph[i].border_out) {
#ifdef NOEQCLASS
            int gid = localGraph[i].gid;
            for (auto &node : this->borderIn) {
                if (this->pt->reach(node, i)) {
                    borderGraph.addEdge(borderGlobal2Local[this->localGraph[node].gid], borderGlobal2Local[gid]);
                    // borderGraph.addEdge(borderGlobal2Local[gid], );
                }
            }
#else

            auto oidx = localGraph[i].borderInIdx;
            if (oidx < 0) {
                continue;
            }
            if (s.count(oidx) == 0) {
                // 加入集合、新建虚拟节点的连接关系
                s.insert(oidx);
                auto idx = this->sout[oidx];

                // 新增一个虚拟节点
                // auto idx = borderGraph.num_vertices();
                s.insert(oidx);
                // borderGraph.addVertex(idx);
                // cout << idx << "  " << borderGraph.num_vertices() << "\n";
                int gid = localGraph[i].gid;
                // 和虚拟节点相连
                borderGraph.addEdge(idx, borderGlobal2Local[gid]);
                // 虚拟节点再和out边界点相连
                auto out_nodes = this->borderClassIn[oidx];
                for (auto dst_bid : out_nodes) {
                    borderGraph.addEdge(dst_bid, idx);
                }
            } else {
                // auto idx = s[oidx];
                auto idx = this->sout[oidx];
                int gid = localGraph[i].gid;
                // 和虚拟节点相连
                borderGraph.addEdge(idx, borderGlobal2Local[gid]);
            }
#endif
        }
    }
#else
    for (int i = 0; i < size; i++) {
        // 连接该节点与相关的out节点
        if (localGraph[i].border_in) {
#ifdef NOEQCLASS
            int gid = this->local2global[i];
            for (auto &node : this->borderOut) {
                if (this->pt->reach(i, node)) {
                    borderGraph.addEdge(borderGlobal2Local[gid], borderGlobal2Local[this->local2global[node]]);
                }
            }
#else

            auto oidx = localGraph[i].borderOutIdx;
            if (oidx < 0) {
                continue;
            }
            if (s.count(oidx) == 0) {
                // 加入集合、新建虚拟节点的连接关系
                s.insert(oidx);
                auto idx = this->sout[oidx];

                // 新增一个虚拟节点
                // auto idx = borderGraph.num_vertices();
                s.insert(oidx);
                // borderGraph.addVertex(idx);
                // cout << idx << "  " << borderGraph.num_vertices() << "\n";
                // int gid = this->local2global[i];
                int gid = localGraph[i].gid;
                // 和虚拟节点相连
                borderGraph.addEdge(borderGlobal2Local[gid], idx);
                // 虚拟节点再和out边界点相连
                auto out_nodes = this->borderClassOut[oidx];
                for (auto dst_bid : out_nodes) {
                    borderGraph.addEdge(idx, dst_bid);
                }
            } else {
                // auto idx = s[oidx];
                auto idx = this->sout[oidx];
                int gid = localGraph[i].gid;

                // 和虚拟节点相连
                borderGraph.addEdge(borderGlobal2Local[gid], idx);
            }
            // auto out_nodes = this->borderClassOut[oidx];
            // int gid = this->local2global[i];
            // for (auto dst_bid : out_nodes) {
            //     borderGraph.addEdge(borderGlobal2Local[gid], dst_bid);
            // }
#endif
        }
    }
#endif
}

void Block::runReachability() {
    CCR::Timer timer;
    auto &g = localGraph;
    int gsize = localGraph.num_vertices();
    vector<int> reverse_topo_sort;
    ptree::GraphUtil::topological_sort(g, reverse_topo_sort);
    cout << "#DAG vertex size:" << g.num_vertices() << "\t#DAG edges size:" << g.num_edges() << endl;

    pt = make_unique<ptree::PathTree>(g, reverse_topo_sort);
    timer.ticker();
    ifstream cfile;
    pt->createLabels(4, cfile, false);
    timer.ticker();
    auto labeling_time = timer.get_last_consuming();
    LOG("#construction time:{} (ms)", labeling_time);
}

void TopChainBlock::runReachability() {
    gt = topchain::transToScc(this->localGraph);
    this->topChainMeta.readGraph(gt);
    this->topChainMeta.readChainInfo(gt);
    this->topChainMeta.findChain();
    this->topChainMeta.feline();
    this->topChainMeta.handle();
    this->topChainMeta.copyIndex();
    this->topChainMeta.flag = new int[topChainMeta.n];
}

bool TopChainBlock::query(int src, int dst) {
    auto s = global2local[src];
    auto d = global2local[dst];
    auto start = localGraph[s].start;
    auto end = localGraph[d].start;

    if (src == dst) {
        return true;
    }
    auto u = gt.getFrom(s, start, end);
    auto w = gt.getTo(d, start, end);
    // cout<<src<<" "<<dst<<" "<<u<<" "<<w<<"\n";

    bool arrive;
    if (u < 0 || w < 0) {
        arrive = 0;
    } else if (u == w)
        arrive = 1;
    else if (topChainMeta.reachindex[u].Ldown >= topChainMeta.reachindex[w].Ldown)
        arrive = 0;
    else if (topChainMeta.reachindex[u].flabel1 > topChainMeta.reachindex[w].flabel1 ||
             topChainMeta.reachindex[u].flabel2 > topChainMeta.reachindex[w].flabel2)
        arrive = 0;
    else if (topChainMeta.intersect(u, w))
        arrive = 1;
    else {
        topChainMeta.mark++;
        arrive = topChainMeta.query(u, w);
    }
    return arrive;
}

bool file_exists(const std::string &filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void queryTransform(BiGraph &bg, ptree::Graph &g, int u, int w, time_t start, time_t end, int &src, int &dst) {
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
                // cout << __FILE__ << ":" << __LINE__ << ":" << i << "---" << k <<
                // endl;
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
            if (bg.timeTable_u[v2].empty() && k == -1) v2 = -1;
            break;

            if (i == bg.adj_matrix_u[w].size() - 1 && bg.timeSection_u[w][i].first >= end) break;
        }
    }

    // LOG("v1 = {}, v2 = {}", v1, v2);
    if (v1 == -1 || v2 == -1) {
        src = -1;
        dst = -1;
        return;
    }

    // for (auto t : bg.timeTable_u[v1]) cout << "[" << t.first << " " << t.second
    // << "],"; cout << endl;

    int v1_idx = -1;
    time_t mm = std::numeric_limits<time_t>::max();
    //   auto iter = std::lower_bound(
    //       bg.timeTable_u[v1].begin(), bg.timeTable_u[v1].end(), start,
    //       [](const auto &a, const auto &b) { return a.second < b; });

    for (int i = 0; i < bg.timeTable_u[v1].size(); i++) {
        if (bg.timeTable_u[v1][i].second > start) {
            if (bg.timeTable_u[v1][i].second < mm) {
                mm = bg.timeTable_u[v1][i].second;
                v1_idx = i;
            }
        }
    }

    //   if (iter != bg.timeTable_u[v1].end()) {
    //     v1_idx = iter - bg.timeTable_u[v1].begin();
    //   }
    if (v1_idx >= 0) {
    } else {
        src = -1;
        dst = -1;
        return;
    }
    // auto iter = std::lower_bound(bg.timeTable_u[v1].begin(),
    // bg.timeTable_u[v1].end(), start); if (iter != bg.timeTable_u[v1].end()) {
    //     v1_idx = iter - bg.timeTable_u[v1].begin();
    // } else {
    //     v1_idx = iter - bg.timeTable_u[v1].begin() - 1;
    // }
    // // cout << __FILE__ << ":" << __LINE__ << "\n";
    int v2_idx = -1;
    //   auto iter = std::lower_bound(
    //       bg.timeTable_u[v2].begin(), bg.timeTable_u[v2].end(), end,
    //       [](const auto &a, const auto &b) { return a.second < b; });
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
        // DEBUG("v1 = {}, v1_idx = {}, start = {}, estart = {}, eend = {}, cstart =
        // {}, cend = {}", v2, v2_idx, start,
        //       bg.timeTable_u[v2][v2_idx].first,
        //       bg.timeTable_u[v2][v2_idx].second, bg.timeSection_u[u][pos2].first,
        //       bg.timeSection_u[u][pos2].second);
        // for (auto t : bg.timeTable_u[v2]) {
        //     cout << "[" << t.first << " " << t.second << "],";
        // }
        // cout << endl;
        // exit(0);
    }
    // iter = std::lower_bound(bg.timeTable_u[v2].begin(),
    // bg.timeTable_u[v2].end(), end); if (iter != bg.timeTable_u[v2].begin()) {
    //     v2_idx = iter - bg.timeTable_u[v2].begin() - 1;
    // } else {
    //     v2_idx = iter - bg.timeTable_u[v2].begin();
    // }

    src = g.cnum(make_pair(v1, bg.timeTable_u[v1][v1_idx].first));
    // cout << __FILE__ << ":" << __LINE__ << "\n";
    dst = g.cnum(make_pair(v2, bg.timeTable_u[v2][v2_idx].first));
    // if (src == 0 or dst == 0) {
    //     cout << v1 << " ==== " << bg.timeTable_u[v1][v1_idx] << "\n";
    //     cout << bg.timeTable_u[v1].size() << "  " << v1_idx << "\n";
    //     cout << v2 << " ==== " << bg.timeTable_u[v2][v2_idx] << "\n";
    //     cout << bg.timeTable_u[v2].size() << "  " << v2_idx << "\n";
    // }
    // cout << __FILE__ << ":" << __LINE__ << "\n";
}

#ifdef NEWSETTING
void queryTransform(BiGraph &bg, ptree::Graph &g, int u, int w, time_t start, time_t end,
                    pair<vector<int>, vector<int>> &querypairs) {
    vector<int> v1s;
    vector<int> v2s;

    for (int i = 0; i < bg.timeSection_u[u].size(); i++)
        if (bg.timeSection_u[u][i].first <= start && bg.timeSection_u[u][i].second > start)
            v1s.push_back(bg.adj_matrix_u[u][i]);

    if (v1s.empty()) {
        time_t max_end = 0;
        int v1 = -1;
        for (int i = 0; i < bg.timeSection_u[u].size(); i++)
            if (bg.timeSection_u[u][i].first > start) {
                max_end = bg.timeSection_u[u][i].second;
                v1 = i;
                break;
            }
        if (v1 >= 0) {
            for (int i = v1; i < bg.timeSection_u[u].size(); i++)
                if (bg.timeSection_u[u][i].first < max_end) v1s.push_back(bg.adj_matrix_u[u][i]);
        }
    }

    // cout<<bg.timeSection_u[u].size()<<"\n";
    for (int i = 0; i < bg.timeSection_u[u].size(); i++)
        if (bg.timeSection_u[u][i].first < end && bg.timeSection_u[u][i].second >= end) {
            // LOG("i = {}, size = {}",i,bg.timeSection_u[u].size());
            v2s.push_back(bg.adj_matrix_u[u][i]);
        }

    if (v2s.empty()) {
        time_t min_start = numeric_limits<time_t>::max();
        int v2 = -1;
        for (int i = bg.timeSection_u[u].size() - 1; i >= 0; i--)
            if (bg.timeSection_u[u][i].second < end) {
                min_start = bg.timeSection_u[u][i].first;
                v2 = i;
                break;
            }
        if (v2 >= 0) {
            for (int i = v2; i >= 0; i--)
                if (bg.timeSection_u[u][i].second > min_start) {
                    // cout<<i<<"  "<<bg.timeSection_u[u].size()<<"=\n";
                    v2s.push_back(bg.adj_matrix_u[u][i]);
                }
        }
    }

    for (auto v1 : v1s) {
        int v1_idx = -1;
        time_t mm = std::numeric_limits<time_t>::max();

        for (int i = 0; i < bg.timeTable_u[v1].size(); i++)
            if (bg.timeTable_u[v1][i].second > start)
                if (bg.timeTable_u[v1][i].second < mm) {
                    mm = bg.timeTable_u[v1][i].second;
                    v1_idx = i;
                }
        if (v1_idx >= 0) {
            int src = g.cnum(make_tuple(v1, bg.timeTable_u[v1][v1_idx].first, bg.timeTable_u[v1][v1_idx].second));
            querypairs.first.push_back(src);
        }
    }

    for (auto v2 : v2s) {
        int v2_idx = -1;
        time_t mm = 0;
        // LOG("v2 = {}, timeTable_u size = {}",v2,bg.timeTable_u.size());
        for (int i = 0; i < bg.timeTable_u[v2].size(); i++)
            if (bg.timeTable_u[v2][i].first < end)
                if (bg.timeTable_u[v2][i].first > mm) {
                    mm = bg.timeTable_u[v2][i].first;
                    v2_idx = i;
                }
        if (v2_idx >= 0) {
            int dst = g.cnum(make_tuple(v2, bg.timeTable_u[v2][v2_idx].first, bg.timeTable_u[v2][v2_idx].second));
            querypairs.second.push_back(v2_idx);
        }
    }
}
#endif

int test_src = 122912;
int test_dst = 365;

unordered_set<int> table;
bool dfs(ptree::Graph &graph, int gid, Partitioner *p) {
    if (table.count(gid) > 0) {
        return false;
    }
    if (gid == test_dst) {
        auto part = graph[gid].partition;
        auto &blk = p->blocks[part];
        auto inidx = blk->localGraph[blk->global2local[gid]].borderInIdx;
        auto xuni = blk->sout[inidx];
        cout << gid << " " << graph[gid].partition << " " << xuni << "\n";
        return true;
    }
    for (auto t : graph.out_edges(gid)) {
        if (dfs(graph, t, p)) {
            auto part = graph[gid].partition;
            auto &blk = p->blocks[part];
            auto inidx = blk->localGraph[blk->global2local[gid]].borderInIdx;
            auto xuni = blk->sout[inidx];
            stringstream ss;
            if (inidx >= 0) {
                auto tmp = blk->borderClassIn[inidx];
                for (auto ttt : tmp) {
                    ss << ttt << ",";
                }
            }
            cout << gid << " " << graph[gid].partition << " " << xuni << " " << ss.str() << "\n";
            return true;
        }
    }
    table.insert(gid);
    return false;
}

void Partitioner::runQueryWithBfs(BiGraph &bg, ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo,
                                  std::vector<int> &queryRes) {
    // Bi BFS
    LOG("query start");
    std::unordered_set<int64_t> s;

    int size = 0;

    int cnt1 = 0;

    for (auto &t : this->blocks) {
        size = max(size, t->getSize());
    }
    int idx = 0;
    int sum_bfs = 0;
    int sum_bfs_reach = 0;
    bool flag = false;
    for (auto const &info : queryInfo) {
        idx++;
        auto u = info.u;
        auto w = info.w;
        auto start = info.start;
        auto end = info.end;
        // 开始时间晚于结束时间，查询结果为假
        if (start > end) {
            queryRes[idx - 1] = false;
            continue;
        }

        int src = u;
        int dst = w;
        queryTransform(bg, graph, u, w, start, end, src, dst);
        if (src < 0 || dst < 0) {
            queryRes[idx - 1] = false;
            continue;
        }

        if (graph[src].partition == graph[dst].partition) {
            // sum_pt_query++;
            if (this->blocks[graph[src].partition]->query(src, dst)) {
                queryRes[idx - 1] = true;
                cnt1++;
            }
        }

        if (queryRes[idx - 1]) {
            continue;
        }

        dst = getInBorder(dst, graph);
        src = getOutBorder(src, graph);
        if (src < 0 || dst < 0) {
            queryRes[idx - 1] = false;
            continue;
        }
        int sum = 0;
#ifdef TWOHOPBORDER
        // queryRes[idx - 1] = runBiBFSBy2Hop(q1hi, q2hi, sum, idx);
        queryRes[idx - 1] = runIndexQuery(src, dst);
        // LOG("res = {}", queryRes[idx - 1]);
#elif defined PathTreeBorder
        queryRes[idx - 1] = runPTQuery(src, dst);
#else
        if (!flag) {
            que1.resize(eqClassGraph.num_vertices());
            que2.resize(eqClassGraph.num_vertices());
            visited_outter = vector<int>(eqClassGraph.num_vertices(), std::numeric_limits<int>::min());
            flag = true;
        }
        queryRes[idx - 1] = runBiBFS(src, dst, sum, idx);
#endif
        // t3.ticker();
        // sum3 += t3.get_last_consuming();
        // LOG("bfs layer = {}, reach = {}",sum,queryRes[idx-1]);
        sum_bfs++;
        if (queryRes[idx - 1]) {
            sum_bfs_reach++;
        }
    }
    LOG("bfs query {}, reach {}, cnt1 = {}", sum_bfs, sum_bfs_reach, cnt1);
}

int Partitioner::getInBorder(int dst, ptree::Graph &graph) {
    auto src_partition = graph[dst].partition;
    return this->blocks[src_partition]->getInBorder2(dst);
}

int Partitioner::getOutBorder(int src, ptree::Graph &graph) {
    auto src_partition = graph[src].partition;
    return this->blocks[src_partition]->getOutBorder2(src);
}

bool Partitioner::runPTQuery(int src, int dst) {
    if (this->pt->reach(src, dst)) {
        return true;
    }
    return false;
}

bool Partitioner::runIndexQuery(int src, int dst) {
    bool arrive;

    if (src == dst) {
        return true;
    }
    auto start = this->eqClassGraph[src].start;
    auto end = this->eqClassGraph[dst].start;
    auto u = this->gt_partitioner.getFrom(src, start, end);
    auto w = this->gt_partitioner.getTo(dst, start, end);
    if (u < 0 || w < 0) {
        arrive = 0;
    } else if (u == w)
        arrive = 1;
    else if (this->topChainMeta_partitioner.reachindex[u].Ldown >= this->topChainMeta_partitioner.reachindex[w].Ldown)
        arrive = 0;
    else if (this->topChainMeta_partitioner.reachindex[u].flabel1 >
                 this->topChainMeta_partitioner.reachindex[w].flabel1 ||
             this->topChainMeta_partitioner.reachindex[u].flabel2 >
                 this->topChainMeta_partitioner.reachindex[w].flabel2)
        arrive = 0;
    else if (this->topChainMeta_partitioner.intersect(u, w))
        arrive = 1;
    else {
        this->topChainMeta_partitioner.mark++;
        arrive = this->topChainMeta_partitioner.query(u, w);
    }
    if (arrive) {
        return true;
    }

    return false;
}

bool Partitioner::runBiBFSBy2Hop(int q1hi, int q2hi, int &sum, int idx) {
    if (q1hi == 0 || q2hi == 0) {
        return false;
    }
    bool arrive = false;
    // 纯单向出发了
#ifdef INNODE
    int dst = que2[q2hi - 1];  // 实际上就一个元素
    int q1lo = 0;
    while (q1hi > q1lo) {
        auto src = que1[q1lo++];
        if (src == dst) {
            return true;
        }
        auto start = this->eqClassGraph[src].start;
        auto end = this->eqClassGraph[dst].start;
        auto u = this->gt_partitioner.getFrom(src, start, end);
        auto w = this->gt_partitioner.getTo(dst, start, end);
        if (u < 0 || w < 0) {
            continue;
        } else if (u == w) {
            return true;
        } else if (this->topChainMeta_partitioner.reachindex[u].Ldown >=
                   this->topChainMeta_partitioner.reachindex[w].Ldown) {
            continue;
        } else if (this->topChainMeta_partitioner.reachindex[u].flabel1 >
                       this->topChainMeta_partitioner.reachindex[w].flabel1 ||
                   this->topChainMeta_partitioner.reachindex[u].flabel2 >
                       this->topChainMeta_partitioner.reachindex[w].flabel2) {
            continue;
        }
        if (this->topChainMeta_partitioner.reachindex[u].flabel1 >
                this->topChainMeta_partitioner.reachindex[w].flabel1 ||
            this->topChainMeta_partitioner.reachindex[u].flabel2 >
                this->topChainMeta_partitioner.reachindex[w].flabel2) {
            continue;
        }
        if (this->topChainMeta_partitioner.IPtest(u, w)) {
            continue;
        }
        if (visited_outter[src] == -idx) {
            return true;
        }
        visited_outter[src] = idx;
        for (auto edge : this->eqClassGraph.out_edges(src)) {
            if (visited_outter[edge] == -idx) {
                return true;
            }
            if (visited_outter[edge] == idx) {
                continue;
            }
            visited_outter[edge] = idx;
            this->que1[q1hi++] = edge;
        }
    }
    this->topChainMeta_partitioner.mark++;
    return false;
#else

#endif
    return false;
}

bool Partitioner::runBiBFS(int src, int dst, int &sum, int idx) {
    // get border graph
    // auto &g = this->borderGraph;
    auto &g = this->eqClassGraph;
    // std::vector<int> visited(borderGraph.num_vertices(), 0);
    int q1lo = 0;
    int q2lo = 0;
    int q1hi = 1;
    int q2hi = 1;

    que1[0] = src;
    que2[0] = dst;
    sum = 0;
    while (q1lo < q1hi and q2lo < q2hi) {
        sum++;
        // DEBUG("in run bfs, src size = {}, dst size = {}", q1hi - q1lo, q2hi -
        // q2lo); unordered_set<int> temp;
        // 正向走，visited 为 1
        if (q1hi - q1lo < q2hi - q2lo) {
            auto origin = q1hi;
            for (int i = q1lo; i < origin; i++) {
                auto cur = que1[q1lo++];
                if (visited_outter[cur] == -idx) {
                    return true;
                }
                visited_outter[cur] = idx;
                for (auto edge : g.out_edges(cur)) {
                    if (visited_outter[edge] == -idx) {
                        return true;
                    } else if (visited_outter[edge] != idx) {
                        visited_outter[edge] = idx;
                        // temp.insert(edge);
                        this->que1[q1hi++] = edge;
                    }
                }
            }
        } else {
            auto origin = q2hi;
            for (int i = q2lo; i < origin; i++) {
                auto cur = que2[q2lo++];
                if (visited_outter[cur] == idx) {
                    return true;
                }
                visited_outter[cur] = -idx;
                for (auto edge : g.in_edges(cur)) {
                    if (visited_outter[edge] == idx) {
                        return true;
                    } else if (visited_outter[edge] != -idx) {
                        // temp.insert(edge);
                        visited_outter[edge] = -idx;
                        this->que2[q2hi++] = edge;
                    }
                }
            }
        }
    }
    return false;
}

string printVector(vector<int> &vec) {
    stringstream ss;
    ss << "[";
    for (auto t : vec) {
        ss << t << ",";
    }
    ss << "]";
    return ss.str();
}

void Block::allocateBaseEqClass(ptree::Graph &eqGraph, ptree::Graph &graph) {
    std::unordered_map<size_t, int> &table = global_partition->GlobalRecordTableOut;
    std::unordered_map<size_t, int> &table2 = global_partition->GlobalRecordTableIn;
    // 分配base单点
    for (int i = 0; i < localGraph.num_vertices(); i++) {
        if (localGraph[i].border_in) {
            if (global_partition->SingleNode2EqGraphNode.find(localGraph[i].gid) ==
                global_partition->SingleNode2EqGraphNode.end()) {
                int idx = eqGraph.num_vertices();
                eqGraph.addVertex(idx);
                global_partition->SingleNode2EqGraphNode[localGraph[i].gid] = idx;
                vector<uint64_t> tmp_vec(this->global_partition->bitElemNumber, 0);
                auto gid = localGraph[i].gid;
                setBitAtPosition(tmp_vec, graph[gid].global_order);
                // graph[gid].bits = tmp_vec;
                auto hash_res = hash_vector(tmp_vec);
                table[hash_res] = idx;
                table2[hash_res] = idx;
            }
        }
    }
}

void Partitioner::ComputeBorder(ptree::Graph &graph) {
    // 把下一个分区的节点也纳入进来
    // 先做正向，反向不一定要做
    unordered_set<int> table;
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (auto edge : graph.out_edges(i)) {
            if (graph[i].partition != graph[edge].partition) {
                int idx = graph[i].partition;
                int toidx = graph[edge].partition;
                // 加点加边
                blocks[idx]->addNode(graph[edge]);
                blocks[idx]->addEdge(i, edge);

                blocks[idx]->setBorderType(edge, 1);
                blocks[toidx]->setBorderType(edge, 2);

                if (table.find(edge) == table.end()) {
                    table.insert(edge);
                    graph[edge].global_order = table.size() - 1;
                }

                // auto &lgraph = blocks[toidx]->localGraph;
                // for (auto edge2 : graph.out_edges(edge)) {
                //     if (graph[edge2].partition == graph[edge].partition and
                //         lgraph[blocks[toidx]->global2local[edge2]].border_in and
                //         graph[edge2].inParitionRecord.find(idx) != graph[edge2].inParitionRecord.end()) {
                //         blocks[idx]->addNode(graph[edge2]);
                //         blocks[idx]->addEdge(edge, edge2);
                //     }
                // }
            }
        }
    }
    LOG("total border edge = {}", table.size());
    this->bitElemNumber = (table.size() + 63) / 64;
    for (int i = 0; i < blocks.size(); i++) {
        blocks[i]->allocateBaseEqClass(eqClassGraph, graph);
    }

#ifdef MTHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < blocks.size(); i++) {
#ifdef PATHTREE
        blocks[i].fetchClassByPT(this->borderGlobal2Local);
#elif defined BORDERBFS
        blocks[i]->fetchClassByBorderBFSOut(this->eqClassGraph, graph);
#else
        blocks[i].fetchClass(this->borderGlobal2Local);
#endif
    }
    for (int i = 0; i < blocks.size(); i++) {
#ifdef PATHTREE
        blocks[i].fetchClassByPT(this->borderGlobal2Local);
#elif defined BORDERBFS
#ifdef MTHREAD
#pragma omp parallel for
#endif
        blocks[i]->fetchClassByBorderBFSIn(this->eqClassGraph, graph);
#else
        blocks[i].fetchClass(this->borderGlobal2Local);
#endif
    }
}

void Block::fetchClass(std::unordered_map<int, int> &borderGlobal2Local) {}

void Block::fetchClassByPT(std::unordered_map<int, int> &borderGlobal2Local) {}

void Block::fetchClassByBorderBFSIn(ptree::Graph &eqGraph, ptree::Graph &graph) {
    // 从入边界点进行处理
    std::unordered_map<size_t, int> &table2 = global_partition->GlobalRecordTableIn;
    std::unordered_map<size_t, int> &table1 = global_partition->GlobalRecordTableOut;

    int size = localGraph.num_vertices();
    LOG("fetchClassByBorderBFS, total node = {}, blk = {}", size, this->blockID);
    vector<int> border_degree(size, 0);
    vector<int> que(size, 0);
    vector<int> visited(size, 0);
    vector<vector<int>> buffer(size);

    int lo = 0;
    int hi = 0;
    for (auto border_node : this->borderIn) {
        // if (localGraph.in_degree(border_node) != 0) {
        //     continue;
        // }
        auto gid = localGraph[border_node].gid;
        graph[gid].bits = vector<uint64_t>(this->global_partition->bitElemNumber, 0);
        setBitAtPosition(graph[gid].bits, graph[gid].global_order);
        que[hi++] = border_node;
        visited[border_node] = -1;
    }
    while (lo < hi) {
        auto node = que[lo++];
        for (auto inEdge : localGraph.out_edges(node)) {
            border_degree[inEdge]++;
            if (visited[inEdge] != -1) {
                visited[inEdge] = -1;
                que[hi++] = inEdge;
            }
        }
    }
    if (not this->borderIn.empty()) {
        for (int i = 0; i < size; i++) {
            if (border_degree[i] > 0) {
                buffer[i].reserve(border_degree[i]);
            }
        }
    }

    lo = 0;
    hi = 0;
    for (auto border_node : this->borderIn) {
        if (border_degree[border_node] == 0) {
            que[hi++] = border_node;
        }
    }

    while (lo < hi) {
        auto node = que[lo++];
        // 收集自己出边传来的class信息
        std::vector<std::vector<int>> allCandidate;
        allCandidate.reserve(1 + border_degree[node]);

        if (localGraph[node].border_in) {
            auto gid = localGraph[node].gid;
            buffer[node].push_back(global_partition->SingleNode2EqGraphNode[gid]);
        }

        // 归并来自子节点的
        // for (auto &kid : buffer[node]) {
        //     allCandidate.push_back(this->global_partition->borderClassIn[kid]);
        // }

        // 归并且去重
        // auto res = mergeAndUniqueUsingMinHeap(allCandidate);
        auto gid = localGraph[node].gid;
        auto res = graph[gid].bits;
        auto hash_res = hash_vector(res);

        if (table2.find(hash_res) == table2.end()) {
            int idx = eqGraph.num_vertices();
            eqGraph.addVertex(idx);
            // global_partition->GlobalSin.insert({idx, idx});
            for (auto &kid : buffer[node]) {
                // int dst = global_partition->GlobalSin[kid];
                int dst = kid;
                eqGraph.addEdge(dst, idx);
            }
            table2[hash_res] = idx;
        }
        // this->global_partition->borderClassIn.insert({table2[hash_res], res});
        localGraph[node].borderInIdx = table2[hash_res];

        // 向前传递信息
        for (auto outEdge : localGraph.out_edges(node)) {
            border_degree[outEdge]--;
            buffer[outEdge].push_back(table2[hash_res]);

            if (border_degree[outEdge] == 0) {
                que[hi++] = outEdge;
            }
            auto parent_gid = localGraph[outEdge].gid;
            if (graph[parent_gid].bits.empty()) {
                graph[parent_gid].bits = res;
            }
            for (int i = 0; i < res.size(); i++) {
                graph[parent_gid].bits[i] |= res[i];
            }
        }
        vector<uint64_t>().swap(graph[gid].bits);
    }
    LOG("blk = {}, after compact, table size = {}", this->blockID, table2.size());
}

void Block::fetchClassByBorderBFSOut(ptree::Graph &eqGraph, ptree::Graph &graph) {
    std::unordered_map<size_t, int> &table = global_partition->GlobalRecordTableOut;
    std::unordered_map<size_t, int> &table2 = global_partition->GlobalRecordTableIn;

    int size = localGraph.num_vertices();
    LOG("fetchClassByBorderBFS, total node = {}", size);
    vector<int> border_degree(size, 0);
    vector<int> que(size, 0);
    vector<int> visited(size, 0);
    vector<vector<int>> buffer(size);

    int lo = 0;
    int hi = 0;
    for (auto border_node : this->borderOut) {
        if (localGraph.out_degree(border_node) != 0) {
            continue;
        }
        que[hi++] = border_node;
        // visited[border_node] = 1;
    }
    // 1.2 收集one hop有哪些节点
    // 1.3 遍历整图收集边界degree
    while (lo < hi) {
        auto node = que[lo++];
        for (auto inEdge : localGraph.in_edges(node)) {
            border_degree[inEdge]++;
            if (visited[inEdge] != 1) {
                visited[inEdge] = 1;
                que[hi++] = inEdge;
            }
        }
    }
    if (not this->borderOut.empty()) {
        for (int i = 0; i < size; i++) {
            if (border_degree[i] > 0) {
                buffer[i].reserve(1 + border_degree[i]);
            }
        }
    }

    // 3. 反向bfs
    // 3.1 先把source找好
    lo = 0;
    hi = 0;
    for (auto border_node : this->borderOut) {
        if (border_degree[border_node] == 0) {
            que[hi++] = border_node;
        }
    }

    while (lo < hi) {
        auto node = que[lo++];
        // 收集自己出边传来的class信息
        std::vector<std::vector<int>> allCandidate;
        allCandidate.reserve(1 + border_degree[node]);
        auto gid = localGraph[node].gid;
        if (localGraph[node].border_out) {
            if (graph[gid].bits.empty()){
                graph[gid].bits = vector<uint64_t>(this->global_partition->bitElemNumber, 0);
            }
            setBitAtPosition(graph[gid].bits, graph[gid].global_order);
            buffer[node].push_back(global_partition->SingleNode2EqGraphNode[localGraph[node].gid]);
        }
        // 归并来自子节点的

        // for (auto &kid : buffer[node]) {
        //     allCandidate.push_back(this->global_partition->borderClassOut[kid]);
        // }
        // 归并且去重
        // auto res = mergeAndUniqueUsingMinHeap(allCandidate);
        auto &res = graph[gid].bits;
        auto hash_res = hash_vector(res);

        if (table.find(hash_res) == table.end()) {
            // 组成了新的等价类
            // 添加节点
            int idx = eqGraph.num_vertices();
            eqGraph.addVertex(idx);
            global_partition->GlobalSout.insert({idx, idx});
            // 前缀树连接边
            // buffer存的等价类编号
            for (auto &kid : buffer[node]) {
                // 获取kid所代表的全局节点
                int dst = kid;
                eqGraph.addEdge(idx, dst);
            }

            // this->global_partition->borderClassOut.emplace_back(res);
            // this->global_partition->borderClassOut.insert({idx, res});
            table[hash_res] = idx;
        }
        localGraph[node].borderOutIdx = table[hash_res];

        // 向前传递信息
        for (auto inEdge : localGraph.in_edges(node)) {
            border_degree[inEdge]--;
            buffer[inEdge].push_back(table[hash_res]);

            if (border_degree[inEdge] == 0) {
                que[hi++] = inEdge;
            }
            auto parent_gid = localGraph[inEdge].gid;
            if (graph[parent_gid].bits.size() == 0) {
                graph[parent_gid].bits = res;
            }
            for (int i = 0; i < res.size(); i++) {
                graph[parent_gid].bits[i] |= res[i];
            }
        }

        vector<uint64_t>().swap(graph[gid].bits);
    }
    LOG("blk = {} ,after compact,  table size = {}", this->blockID, table.size());
    return;
}

void Block::fetchClassByBorderBFS(std::unordered_map<int, int> &borderGlobal2Local) {
    // 1 反向遍历graph,记录每个节点的边界出度
    std::unordered_map<size_t, int> table;
    int size = localGraph.num_vertices();
    LOG("fetchClassByBorderBFS, total node = {}", size);
    vector<int> border_degree(size, 0);
    vector<int> que(size, 0);
    vector<int> visited(size, 0);
    vector<vector<int>> buffer(size);
    bufferO = vector<unordered_set<int>>(size);
    bufferI = vector<unordered_set<int>>(size);
    int lo = 0;
    int hi = 0;
    for (auto border_node : this->borderOut) {
        if (localGraph.out_degree(border_node) != 0) {
            continue;
        }
        que[hi++] = border_node;
        // visited[border_node] = 1;
    }
    // 1.2 收集one hop有哪些节点
    // 1.3 遍历整图收集边界degree
    while (lo < hi) {
        auto node = que[lo++];
        for (auto inEdge : localGraph.in_edges(node)) {
            border_degree[inEdge]++;
            if (visited[inEdge] != 1) {
                visited[inEdge] = 1;
                que[hi++] = inEdge;
            }
        }
    }
    if (not this->borderOut.empty()) {
        for (int i = 0; i < size; i++) {
            if (border_degree[i] > 0) {
                buffer[i].reserve(border_degree[i]);
            }
        }
    }

    // 3. 反向bfs
    // 3.1 先把source找好
    lo = 0;
    hi = 0;
    for (auto border_node : this->borderOut) {
        if (border_degree[border_node] == 0) {
            que[hi++] = border_node;
        }
    }

    while (lo < hi) {
        auto node = que[lo++];
        // 收集自己出边传来的class信息
        std::vector<std::vector<int>> allCandidate;
        allCandidate.reserve(1 + border_degree[node]);
        if (localGraph[node].border_out) {
            allCandidate.push_back({borderGlobal2Local[localGraph[node].gid]});
        }
        // 归并来自子节点的

        for (auto &kid : buffer[node]) {
            allCandidate.push_back(this->borderClassOut[kid]);
        }
        // 归并且去重
        auto res = mergeAndUniqueUsingMinHeap(allCandidate);
        auto hash_res = hash_vector(res);
        if (table.find(hash_res) == table.end()) {
            // 组成了新的等价类
            int idx = this->borderClassOut.size();
            this->bufferO[idx].insert(buffer[node].begin(), buffer[node].end());
            this->borderClassOut.emplace_back(res);
            table[hash_res] = idx;
#ifndef INNODE
            // 遍历该等价类包含的所有边界点lid
            // key = borderGraph lid, value = set of 等价类编号
            for (auto t : res) {
                extraTable[t].insert(idx);
            }
#endif
        }
        localGraph[node].borderOutIdx = table[hash_res];
        // 向前传递信息
        for (auto inEdge : localGraph.in_edges(node)) {
            border_degree[inEdge]--;
            buffer[inEdge].push_back(table[hash_res]);

            if (border_degree[inEdge] == 0) {
                que[hi++] = inEdge;
            }
        }
    }
    this->eqClassVisitedO = vector<int>(table.size());
    LOG("after compact, out node = {}, table size = {}", this->borderClassOut.size(), table.size());

    // 从入边界点进行处理
    table = std::unordered_map<size_t, int>();
    lo = 0;
    hi = 0;
    for (auto border_node : this->borderIn) {
        que[hi++] = border_node;
        visited[border_node] = -1;
    }
    while (lo < hi) {
        auto node = que[lo++];
        for (auto inEdge : localGraph.out_edges(node)) {
            border_degree[inEdge]++;
            if (visited[inEdge] != -1) {
                visited[inEdge] = -1;
                que[hi++] = inEdge;
            }
        }
    }
    if (not this->borderIn.empty()) {
        for (int i = 0; i < size; i++) {
            buffer[i].clear();
            if (border_degree[i] > 0) {
                buffer[i].reserve(border_degree[i]);
            }
        }
    }

    lo = 0;
    hi = 0;
    for (auto border_node : this->borderIn) {
        // cout<<border_node<<" : degree = "<<border_degree[border_node]<<"\n";
        if (border_degree[border_node] == 0) {
            que[hi++] = border_node;
        }
    }

    while (lo < hi) {
        // LOG("lo = {}, hi ={}",lo,hi);
        auto node = que[lo++];
        // 收集自己出边传来的class信息
        std::vector<std::vector<int>> allCandidate;
        allCandidate.reserve(1 + border_degree[node]);
        if (localGraph[node].border_in) {
            allCandidate.push_back({borderGlobal2Local[localGraph[node].gid]});
        }
        // 归并来自子节点的

        for (auto &kid : buffer[node]) {
            allCandidate.push_back(this->borderClassIn[kid]);
        }
        // 归并且去重
        auto res = mergeAndUniqueUsingMinHeap(allCandidate);
        auto hash_res = hash_vector(res);
        // LOG("res size = {}",res.size());
        if (table.find(hash_res) == table.end()) {
            this->index_size += res.size();
            this->index_size += res.size();
            int idx = this->borderClassIn.size();
            this->bufferI[idx].insert(buffer[node].begin(), buffer[node].end());
            this->borderClassIn.emplace_back(res);
            table[hash_res] = idx;
#ifdef INNODE
            // 遍历该等价类包含的所有边界点lid
            // key = borderGraph lid, value = set of 等价类编号
            for (auto t : res) {
                extraTable[t].insert(idx);
            }
#endif
        }
        localGraph[node].borderInIdx = table[hash_res];

        // 向前传递信息
        for (auto outEdge : localGraph.out_edges(node)) {
            border_degree[outEdge]--;
            buffer[outEdge].push_back(table[hash_res]);

            if (border_degree[outEdge] == 0) {
                que[hi++] = outEdge;
            }
        }
    }
    this->eqClassVisitedI = vector<int>(table.size());
    LOG("after compact, in node = {}, table size = {}", this->borderClassIn.size(), table.size());
}
