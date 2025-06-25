#include "compact.h"

#include <algorithm>
#include <deque>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph.h"
#include "utils.h"

void Compute2Hop(ptree::Graph &graph, int n) {
    std::unordered_map<string, int> table;
    int size = graph.num_vertices();
    LOG("total node = {}", size);
    for (int i = 0; i < size; i++) {
        // 2 hop
        std::deque<int> deq;
        std::unordered_set<int> s;
        for (auto edge : graph.out_edges(i)) {
            s.insert(edge);
            deq.push_back(edge);
        }
        // for (auto t : deq) {
        //     for (auto edge : graph.out_edges(t)) {
        //         s.insert(edge);
        //     }
        // }
        for (int k = 1; k < n; k++) {
            int len = deq.size();
            for (int j = 0; j < len; j++) {
                auto node = deq.front();
                deq.pop_front();
                s.insert(node);
                for (auto edge : graph.out_edges(node)) {
                    deq.push_back(edge);
                }
            }
        }

        std::vector<int> candidate(s.begin(), s.end());
        std::sort(candidate.begin(), candidate.end());
        std::stringstream ss;
        for (auto t : candidate) {
            ss << t << ",";
        }

        table[ss.str()]++;
    }
    LOG("after compact, node = {}", table.size());
}

void ComputeBorder(ptree::Graph &graph) {
    std::unordered_map<string, int> table;
    int size = graph.num_vertices();
    LOG("total node = {}", size);
    for (int i = 0; i < size; i++) {
        std::deque<int> deq;
        std::set<int> visited;
        auto src_partition = graph[i].partition;
        deq.push_back(i);
        std::unordered_set<int> s;
        while (!deq.empty()) {
            auto node = deq.front();
            deq.pop_front();
            if (graph[node].border_out) {
                s.insert(node);
            }
            for (auto edge : graph.out_edges(node)) {
                if (graph[edge].partition != src_partition) {
                    s.insert(node);
                } else {
                    if (visited.count(edge) == 0) {
                        deq.push_back(edge);
                        visited.insert(edge);
                    }
                }
            }
        }

        std::vector<int> candidate(s.begin(), s.end());
        std::sort(candidate.begin(), candidate.end());
        std::stringstream ss;
        for (auto t : candidate) {
            ss << t << ",";
        }
        table[ss.str()]++;
        if (i % 10000 == 0) {
            LOG("遍历了{}个节点, 新图节点 = {}, 空集大小 = {}", i, table.size(), table[""]);
        }
    }
    LOG("after compact, node = {}", table.size());
    for (auto it : table) {
        cout << it.first << " -- " << it.second << "\n";
    }
}