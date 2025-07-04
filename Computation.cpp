#include <gflags/gflags.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "basebfs/bfs.h"
#include "compact.h"
#include "graph.h"
#include "nlohmann/json.hpp"
#include "partition.h"
#include "utils.h"
using namespace std;
using json = nlohmann::json;

double getMemoryUsageGB() {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.substr(0, 6) == "VmSize") {
            size_t pos = line.find_first_of("0123456789");
            if (pos != std::string::npos) {
                // 将读取的内存大小从 kB 转换为 GB
                size_t memSizeKB = std::stoull(line.substr(pos));
                return memSizeKB / 1048576.0;  // 1 GB = 1048576 kB
            }
        }
    }
    return 0;  // 如果没有找到 VmSize，返回0
}

double maxMem = 0;
void monitorMemory() {
    while (true) {
        double memoryUsageGB = getMemoryUsageGB();  // 读取内存使用量并转换为 GB
        // std::cout << "Current memory usage: " << memoryUsageGB << " GB" <<
        LOG("Current memory usage = {} GB", memoryUsageGB);
        maxMem = max(maxMem, memoryUsageGB);
        sleep(10);  // 每10秒执行一次
    }
}
ptree::Graph readGraph(string path) {
    ptree::Graph g;
    int num;
    ifstream ifs(path);
    ifs >> num;
    g.addVertex(num - 1);
    for (int i = 0; i < num; i++) {
        g[i].id = i;
        g[i].start = i;
        g[i].node = i;
    }
    int a, b;
    while (ifs >> a >> b) {
        g.addEdge(a, b);
    }

    return g;
}

void prinfGraph(ptree::Graph g) {
    for (int i = 0; i < g.num_vertices(); i++) {
        cout << i << ": ";
        for (auto edge : g.out_edges(i)) {
            cout << edge << " ";
        }
        cout << "\n";
    }
    for (int i = 0; i < g.num_vertices(); i++) {
        cout << i << ": indegree " << g.in_degree(i) << " outdegree " << g.out_degree(i) << "\n";
    }
    cout << "edges = " << g.num_edges() << "\n";
}

void exportQuery(string path, vector<CCR::queryInfo> &pairs) {
    ofstream ofs(path);
    for (auto pair : pairs) {
        ofs << pair.u << " " << pair.w << " " << pair.start << " " << pair.end << "\n";
    }
}

DEFINE_string(inputTestPath, "", "输入文件路径,绝对路径");
DEFINE_string(inputPath, "", "输入文件路径,绝对路径");
DEFINE_string(queryFilePath, "", "查询测试文件路径,绝对路径");
DEFINE_string(outputPath, "", "输出文件路径,绝对路径");
DEFINE_string(algorithm, "no func", "试图使用的算法");

DEFINE_int32(maxBlock, 64, "最大分区数");
DEFINE_int32(reach, 0, "是否可达查询");
DEFINE_int32(minBlock, 4, "最小分区数");
DEFINE_int32(numStep, 100, "循环次数");
DEFINE_int32(delta, std::numeric_limits<int32_t>::max(), "delta 截断");
DEFINE_int32(splitData, 0, "阈值，控制是否进行点切割");

string querypairfilename = "../result/query";

int main(int argc, char *argv[]) {
    if (argc == 1) return 1;
    spdlog::cfg::load_env_levels();
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    auto input_path = FLAGS_inputPath;
    auto output_path = FLAGS_outputPath;
    auto algo = FLAGS_algorithm;
    auto name = CCR::split(input_path, "/.", 100);
    auto fn = CCR::extractFileName(input_path);
    std::thread backThread(monitorMemory);
    backThread.detach();
    // cout << "开始读取 " << fn << endl;
    LOG("开始读取 {}", fn);
    LOG("output path = {}", output_path);
    querypairfilename += "_" + fn;
    BiGraph bg(input_path, "r");

    if (bg.adj_matrix_l.empty()) {
        exit(-1);
    }

    // string queryResFilePath = FLAGS_queryFilePath + "_" + ".queryres";
    string reachQueryPairs = FLAGS_queryFilePath+"_1";
    string unreachQueryPairs = FLAGS_queryFilePath+"_0";
    auto Reachpairs = CCR::readQuery(reachQueryPairs);
    auto UnReachpairs = CCR::readQuery(unreachQueryPairs);
    // 1代表可达查询，0代表不可达查询
    // auto pairs = CCR::readQuery(FLAGS_queryFilePath, queryResFilePath, FLAGS_reach);
    if (Reachpairs.empty()) {
        ERROR("no query");
        exit(-1);
    }
    LOG("origin reach pairs size = {}, unreach size = {}", Reachpairs.size(),UnReachpairs.size());
    while (Reachpairs.size() < 100000) {
        std::vector<CCR::queryInfo> temp;
        for (auto element : Reachpairs) {
            if (Reachpairs.size() + temp.size() >= 100000) {
                break;
            }
            temp.push_back(element);
        }
        Reachpairs.insert(Reachpairs.end(), temp.begin(), temp.end());
    }

    CCR::Timer timer;
    vector<int> queryRes(Reachpairs.size(), 0);
    uint64_t trans_time = 0;
    unsigned bgEdgeNum = 0;
    for (const auto &Neighbors : bg.adj_matrix_l) bgEdgeNum += Neighbors.size();
    timer.start();
    ptree::Graph g = CCR::transformation(bg, FLAGS_delta);
    timer.ticker();
    trans_time = timer.get_last_consuming();
    LOG("transformation time = {}(ms)", timer.get_last_consuming());
   
    Partitioner partitioner(FLAGS_minBlock);
    partitioner.threshold = FLAGS_splitData;
    uint64_t partition_time;
    timer.ticker();
    partitioner.runPartition(g, algo);
    timer.ticker();
    LOG("partition time = {} (ms)", timer.get_last_consuming());
    partition_time = timer.get_last_consuming();

    uint64_t local_algo_time;
    timer.ticker();
    partitioner.runLocalReachability();  // 内部可达算法
    timer.ticker();
    local_algo_time = timer.get_last_consuming();
    LOG("local algorithm {}, time = {} (ms)", algo, timer.get_last_consuming());

    uint64_t eqClassTime;
    timer.ticker();
    partitioner.ComputeBorder(g);  // 找等价类
    timer.ticker();
    eqClassTime = timer.get_last_consuming();
    LOG("ComputeBorder = {} (ms)", timer.get_last_consuming());

    // 边界可达性算法
    uint64_t eqGraphTime;
    timer.ticker();
    partitioner.runBorderReachability(g);  // 构造边界图
    timer.ticker();
    eqGraphTime = timer.get_last_consuming();
    LOG("border Reachability = {} (ms)", timer.get_last_consuming());

    uint64_t queryTime;
    timer.ticker();
    partitioner.runQueryWithBfs(bg, g, Reachpairs, queryRes);
    timer.ticker();
    queryTime = timer.get_last_consuming();

    uint64_t queryTimeNoReach;
    timer.ticker();
    partitioner.runQueryWithBfs(bg, g, UnReachpairs, queryRes);
    timer.ticker();
    queryTimeNoReach = timer.get_last_consuming();

    json j;
    j["Frag_num"] = FLAGS_minBlock;
    j["Trans_time(ms)"] = trans_time;
    j["Partition_time(ms)"] = partition_time;
    j["indexing_method"] = algo;
    j["local_indexing_time(ms)"] = local_algo_time;
    j["EqSet_time(ms)"] = eqClassTime;
    j["ConGraph_indexing_time(ms)"] = eqGraphTime;
    j["Connection_latency_threshold(s)"] =  FLAGS_delta;

    j["ConGraph_node_num"] = std::to_string(partitioner.eqClassGraph.num_vertices());
    j["ConGraph_edge_num"] = std::to_string(partitioner.eqClassGraph.num_edges());
    j["Memory(GB)"] = maxMem;
    j["TotalPosQuery(ms)"] = queryTime;
    j["AvgPosQuery(us)"] = queryTime*1000 / (double)Reachpairs.size();
    j["TotalNegQuery(ms)"] = queryTimeNoReach;
    j["AvgNegQuery(us)"] = queryTimeNoReach*1000 / (double)UnReachpairs.size();
    j["Index_construction_time(ms)"] = trans_time + partition_time + local_algo_time + eqClassTime + eqGraphTime;
    // j["index size"] = partitioner.getIndexSize();
    j["DAG_node_num"] = g.num_vertices();
    j["DAG_edge_num"] = g.num_edges();
    double sum = 0;
    for (int i = 0; i < g.num_vertices(); i++) {
        for (auto edge : g.out_edges(i)) {
            if (g[i].partition != g[edge].partition) {
                sum++;
            }
        }
    }
    LOG("edge cut = {}", sum / g.num_edges());
    j["edeg_cut"] = sum / g.num_edges();
    std::ofstream ofile(FLAGS_outputPath);
    ofile << j << std::endl;
    LOG("run {} query time = {}(ms), total time = {}, build index time = {}, max "
        "memory usage = {}",
        Reachpairs.size(), timer.get_last_consuming(), timer.get_total_consuming(),
        timer.get_total_consuming() - timer.get_last_consuming(), maxMem);
}

