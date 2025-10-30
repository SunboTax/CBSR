#include "partition.h"
#include "utils.h"

int Partitioner::GetPerfectPartitionNum(ptree::Graph& graph)
{
    int originNum = this->originNum;
    int lo = 1;
    int hi = originNum;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        double rating = GetComuputedBlockNum(graph, mid);
        if (rating > 1) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

// 改成2分的方法
double
Partitioner::GetComuputedBlockNum(ptree::Graph& graph, int num)
{
    //=====================准备工作==========================
    int originNum = num;
    double rating = -1;
    double hyperGraphNum = 0;
    // pre partition
    auto cmp = [](ptree::Vertex& first, ptree::Vertex& second) {
        if (first.start == second.start) {
            return first.end > second.end;
        }
        return first.start > second.start;
    };

    auto gsize = graph.num_vertices();
    int maxVertexSize = (gsize + originNum - 1) / originNum;

    // 分区，只计算hyperGraph的节点数量
    priority_queue<ptree::Vertex, std::vector<ptree::Vertex>, decltype(cmp)>
        que(cmp);
    for (auto& node : graph.vertices()) {
        if (graph.in_degree(node.id) == 0) {
            graph[node.id].visited = true;
            que.push(node);
        }
    }
    this->blocks = vector<Block*>(originNum);
    for (int i = 0; i < originNum; i++) {
        this->blocks[i] = new Block();
        this->blocks[i]->global_partition = this;
    }

    // ===============defer function========================
    defer({
        LOG("clean {} blocks start", originNum);
        for (auto blk : this->blocks) {
            delete blk;
        }
        LOG("clean {} blocks done", originNum);
        this->SingleNode2EqGraphNode.clear();
        LOG("clean SingleNode2EqGraphNode done");
        for (int i = 0; i < graph.num_vertices(); i++) {
            graph[i].visited = false;
        }
        this->GlobalRecordTableOut.clear();
        this->GlobalRecordTableIn.clear();
    });

    //=====================准备工作over==========================

    //=====================分区==========================
    int partitionIdx = 0;
    int nowSum = 0;
    while (!que.empty()) {
        blocks[partitionIdx]->setBlockID(partitionIdx);
        ptree::Vertex v = que.top();
        que.pop();

        graph[v.id].partition = partitionIdx;
        graph[v.id].gid = v.id;
        // blocks[partitionIdx]->addNode(graph[v.id]);
        nowSum++;
        graph[v.id].lid = blocks[partitionIdx]->localGraph.num_vertices() - 1;
        for (auto t : graph.out_edges(v.id)) {
            if (not graph[t].visited) {
                que.push(graph[t]);
                graph[t].visited = true;
                // tmpBlk.crossOutEdge += graph.out_degree(v.id);
            }
        }

        if (nowSum >= maxVertexSize) {
            partitionIdx++;
            nowSum = 0;
        }
    }

    for (int i = 0; i < graph.num_vertices(); i++) {
        for (auto t : graph.out_edges(i)) {
            int idx = graph[i].partition;
            if (idx == -1 or graph[t].partition == -1) {
                ERROR("自环");
                exit(1);
            }
            // if (graph[i].partition == graph[t].partition) {
            //     blocks[idx]->addEdge(i, t);
            // }
        }
    }
    //=====================等价类计算==========================

    //============等价类计算——加点工作===============
    // 将下个分区的节点纳入计算
    unordered_set<int> table;
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (auto edge : graph.out_edges(i)) {
            if (graph[i].partition != graph[edge].partition) {
                int idx = graph[i].partition;
                int toidx = graph[edge].partition;
                // 加点加边
                // blocks[idx]->addNode(graph[edge]);
                // blocks[idx]->addEdge(i, edge);

                // blocks[idx]->setBorderType(edge, 1);
                // blocks[toidx]->setBorderType(edge, 2);

                if (table.find(edge) == table.end()) {
                    table.insert(edge);
                    graph[edge].global_order = table.size() - 1;
                }
            }
        }
    }
    LOG("total border edge = {}", table.size());
    this->bitElemNumber = (table.size() + 63) / 64;

    // //============等价类计算——分配base节点===============
    // CCR::Timer timer;
    // timer.start();
    // for (int i = 0; i < blocks.size(); i++) {
    //     blocks[i]->allocateBaseEqClass4PreCompute(eqClassGraph, graph);
    // }
    // timer.ticker();
    // LOG("allocate base eqClass time = {} (ms)", timer.get_last_consuming());
    // // 获取所有单点数量
    // hyperGraphNum = this->SingleNode2EqGraphNode.size();
    // LOG("before compute eqClass, total hyperGraphNum = {}", hyperGraphNum);
    // // topo逆排序获取总节点数量
    // int sum = 0;
    // for (int i = 0; i < blocks.size(); i++) {
    //     sum += blocks[i]->fetchClassByBorderBFSOut4PreCompute(eqClassGraph,
    //     graph);
    // }

    // hyperGraphNum += sum;
    // sum = 0;
    // for (int i = 0; i < blocks.size(); i++) {
    //     sum += blocks[i]->fetchClassByBorderBFSIn4PreCompute(eqClassGraph,
    //     graph);
    // }
    // LOG("in add hyper {}", sum);
    // hyperGraphNum += sum;
    // timer.ticker();
    // LOG("after compute eqClass, total hyperGraphNum = {}, inner size = {},
    // block num = {}, time = {}", hyperGraphNum,
    //     maxVertexSize, originNum, timer.get_last_consuming());
    // LOG("hyper size = {}, border size = {}, blk = {}", hyperGraphNum,
    // table.size(), originNum);
    hyperGraphNum = table.size() * 8;
    LOG("predict hyper size = {}, border size = {}, blk = {}", hyperGraphNum,
        table.size(), originNum);
    return (maxVertexSize * 1.0) / (hyperGraphNum * 1.0);
}
