#ifndef UTILS_H
#define UTILS_H

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#define DEBUG(...) \
    SPDLOG_LOGGER_DEBUG(spdlog::default_logger_raw(), __VA_ARGS__)
#define LOG(...) \
    SPDLOG_LOGGER_INFO(spdlog::default_logger_raw(), __VA_ARGS__)
#define WARN(...) \
    SPDLOG_LOGGER_WARN(spdlog::default_logger_raw(), __VA_ARGS__)
#define ERROR(...) \
    SPDLOG_LOGGER_ERROR(spdlog::default_logger_raw(), __VA_ARGS__)
#include <algorithm>
#include <chrono>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph.h"
#include "spdlog/cfg/env.h"
#include "spdlog/spdlog.h"

#define VMRSS_LINE 22
using namespace std;

namespace CCR {

struct Edge {
    int node;
    time_t start;
    time_t end;

    Edge() {};

    Edge(int node, time_t start, time_t end)
    {
        this->node = node;
        this->start = start;
        this->end = end;
    };

    inline bool
    operator<(const Edge& e) const
    {
        if (this->node < e.node)
            return true;
        if (this->node == e.node) {
            if (this->start < e.start)
                return true;
            if (this->start == e.start) {
                if (this->end < e.end)
                    return true;
            }
        }
        return false;
    }

    inline bool
    operator==(const Edge& e) const
    {
        return this->node == e.node && this->start == e.start
            && this->end == e.end;
    }

    inline bool
    operator!=(const Edge& e) const
    {
        return !(*this == e);
    }
};

struct queryInfo {
    int u = -1;
    int w = -1;
    time_t start = 0;
    time_t end = 0;
    bool reachable = false;
};

vector<string> split(string, const char*, int);
ptree::Graph transformation(BiGraph&, int);
ptree::Graph transformation(BiGraph&);
// ptree::Graph newTransformation(BiGraph &);
ptree::Graph newTransformation(BiGraph&, int);
vector<queryInfo> readQuery(string);
vector<queryInfo> readQuery(string, string, int);

vector<queryInfo> readQueryForTest(string);

// Graph dagConstruction(BiGraph &bg);
// void GraphInfo(Graph &g, string filename);
// void Graphvisual(Graph &g, string filename, int k);
vector<ptree::subGraph> allSubGraph(ptree::Graph& g);

// void printGraph(Tree &tree, string prefix, int depth, int id);
// void printGraph(Graph &graph, string prefix, int depth, int id);
void printGraph(ptree::Graph& graph, string prefix, int id);
void printDAG_1(ptree::Graph& graph, string outputPath);
void printDAG_2(ptree::Graph& graph, string outputPath);
void printDAG_metis(ptree::Graph& graph, string outputPath);

// template <typename T>
// int dfs(T &graph, vector<int> &nodeDepth, int node)
// {
//   if (nodeDepth[node] >= 0)
//     return nodeDepth[node];
//   if (graph.in_degree(node) == 0)
//     nodeDepth[node] = 0;
//   else
//   {
//     int d = 0;
//     for (int u : graph.in_edges(node))
//       d = max(d, dfs(graph, nodeDepth, u));
//     nodeDepth[node] = d + 1;
//   }
//   return nodeDepth[node];
// }

class Timer {
private:
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    std::vector<std::chrono::steady_clock::time_point> times;

public:
    inline void
    start()
    {
        start_time = std::chrono::steady_clock::now();
        times.push_back(start_time);
    }
    inline void
    end()
    {
        end_time = std::chrono::steady_clock::now();
        times.push_back(end_time);
    }
    inline void
    ticker()
    {
        times.emplace_back(std::chrono::steady_clock::now());
    }
    inline long long
    get_last_consuming()
    {
        if (times.size() < 2) {
            ticker();
        }
        size_t size = times.size();
        auto first = times[size - 1];
        auto second = times[size - 2];
        return std::chrono::duration_cast<std::chrono::milliseconds>(first
            - second)
            .count();
    }
    inline long long
    get_total_consuming()
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(end_time
            - start_time)
            .count();
    }
    Timer() { times.reserve(1000); }
};

float GetMemoryUsage(int pid);

std::string extractFileName(const std::string& filePath);

} // namespace CCR

// 定义一个Defer类，用于执行延迟操作
class Defer {
private:
    std::function<void()> func; // 存储要延迟执行的函数

public:
    // 构造函数，接受一个函数作为参数
    explicit Defer(std::function<void()> f)
        : func(std::move(f))
    {
    }

    // 析构函数，在对象销毁时执行延迟操作
    ~Defer()
    {
        if (func) {
            func();
        }
    }

    // 禁止拷贝构造和赋值操作
    Defer(const Defer&) = delete;
    Defer& operator=(const Defer&) = delete;

    // 允许移动构造和赋值操作
    Defer(Defer&&) = default;
    Defer& operator=(Defer&&) = default;
};

// 定义一个宏，用于简化defer的使用
#define CONCAT(a, b) a##b
#define DEFER_TOKEN_NAME(line) CONCAT(defer_, line)
#define defer(func) Defer DEFER_TOKEN_NAME(__LINE__)([&]() { func; })

#endif // UTILS_H
