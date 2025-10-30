#include <bits/types/time_t.h>
#include <gflags/gflags.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CCRTree.h"
#include "graph.h"
#include "utils.h"
double
getMemoryUsageGB()
{
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.substr(0, 6) == "VmSize") {
            size_t pos = line.find_first_of("0123456789");
            if (pos != std::string::npos) {
                size_t memSizeKB = std::stoull(line.substr(pos));
                return memSizeKB / 1048576.0;
            }
        }
    }
    return 0;
}

void monitorMemory()
{
    while (true) {
        double memoryUsageGB = getMemoryUsageGB();

        LOG("Current memory usage = {} GB", memoryUsageGB);

        sleep(10);
    }
}
using namespace std;

DEFINE_string(mode, "", "使用模式");
DEFINE_string(inputPath, "", "输入文件路径,绝对路径");
DEFINE_string(outputPath, "", "输出文件路径,绝对路径");
DEFINE_int32(parts, 4, "分区数量");
DEFINE_int32(queryNum, 4, "生成查询数量");

string outout_file = "";
void calc(ptree::Graph g, string metis_file, int num_part);

std::string
unixTimestampToDateString(const std::time_t timestamp)
{
    std::tm* tm_struct = std::gmtime(&timestamp);
    char buffer[20];

    std::strftime(buffer, sizeof(buffer), "%Y.%m.%d", tm_struct);

    return buffer;
}

struct queryInfo {
    int u = -1;
    int w = -1;
    time_t start = 0;
    time_t end = 0;
    bool reachable = false;
};

void generateRandomQuery(int query_num, ptree::Graph& g, BiGraph& bg,
    string querypairfilename)
{
    int upper_size = bg.adj_matrix_u.size();

    time_t min_time = std::numeric_limits<time_t>::max();
    time_t max_time = std::numeric_limits<time_t>::min();

    for (int i = 0; i < g.num_vertices(); i++) {
        auto value = g[i].start;
        min_time = std::min(value, min_time);
        max_time = std::max(max_time, value);
    }
    ofstream ofs(querypairfilename);
    if (not ofs.is_open()) {
        cout << "error open " << querypairfilename << "\n";
    }
    auto gsize = g.num_vertices();
    for (int i = 0; i < query_num; i++) {
        auto s = lrand48() % upper_size;
        auto t = lrand48() % upper_size;
        time_t start, end;
        start = min_time + lrand48() % (max_time - min_time);
        end = min_time + lrand48() % (max_time - min_time);
        ofs << s << " " << t << " " << start << " " << end << "\n";
    }
}

void generateRandomQuery(int query_num, int gsize, string querypairfilename)
{
    cout << "generating queries..." << endl;
    int left = 0;
    int s, t;
    vector<int> src;
    vector<int> trg;
    while (left < query_num) {
        s = lrand48() % gsize;
        t = lrand48() % gsize;
        if (s == t)
            continue;
        src.push_back(s);
        trg.push_back(t);
        ++left;
    }
    cout << "generate done" << endl;

    ofstream ofs(querypairfilename);
    if (not ofs.is_open()) {
        cout << "error open " << querypairfilename << "\n";
    }

    for (int i = 0; i < src.size(); i++) {
        auto key = MERGE_TO_UINT64(src[i], trg[i]);
        ofs << src[i] << " " << trg[i] << "\n";
    }
}

int main(int argc, char* argv[])
{
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    auto mode = FLAGS_mode;
    if (mode == "PrintMetisInfo") {
        string path = FLAGS_inputPath;
        auto name = CCR::split(path, "/.", 100);
        auto parts = FLAGS_parts;
        auto fn = CCR::extractFileName(path);
        BiGraph bg(path, "r");
        std::list<CCR::res_type> res;
        ptree::Graph g = CCR::transformation(bg);
        std::cout << "graph size = " << g.num_vertices()
                  << " edge size = " << g.num_edges() << "\n";
        string metis_file = "../" + fn + ".metis.part." + to_string(parts);

        calc(g, metis_file, parts);
    } else if (mode == "GenerateQuery") {
        string outpath = FLAGS_outputPath;
        string inpath = FLAGS_inputPath;
        auto queryNum = FLAGS_queryNum;
        std::thread backThread(monitorMemory);
        backThread.detach();
        BiGraph bg(inpath, "r");
        LOG("begin trans");

        ptree::Graph g = CCR::newTransformation(bg, 0);
        generateRandomQuery(queryNum, g, bg, outpath);
    } else if (mode == "FixDataSet") {
        string outpath = FLAGS_outputPath;
        string inpath = FLAGS_inputPath;
        string buffer;
        ifstream fin(inpath);
        ofstream fout(outpath);
        int num_upper, num_lower;
        fin >> buffer >> num_upper >> num_lower;
        stringstream ss;
        ss << buffer << " " << num_upper << " " << num_lower << "\n";
        int64_t upper, start, lower, end;
        int64_t i = 0;
        int min_d = 1 * 600;
        int max_d = 8 * 3600;
        float alpha = -2.5;
        time_t maxDuration = 0;

        while (fin >> upper >> lower >> start >> end) {
            i++;
            time_t end;

            float y = ((float)rand()) / ((float)RAND_MAX);
            double x1 = pow(max_d, (alpha + 1));
            double x2 = pow(min_d, (alpha + 1));
            double x3 = (x2 - x1) * y + x1;
            double x4 = pow(x3, 1 / (alpha + 1));
            time_t duration = (time_t)x4;
            maxDuration = max(maxDuration, duration);
            end = start + duration;
            ss << upper << " " << lower << " " << start << " " << end << "\n";
            if (i % 50000 == 0) {
                fout << ss.str();
                ss.str("");
            }
        }
        fout << ss.str() << "\n";
        ss.str("");

        LOG("generate data done");
    } else if (mode == "GenerateRfile") {
        string outpath = FLAGS_outputPath;
        string inpath = FLAGS_inputPath;
        ifstream fin(inpath);
        ofstream fout(outpath);
        int num_upper, num_lower;

        string buffer;
        int min_d = 1 * 600;
        int max_d = 8 * 3600;
        float alpha = -2.5;
        fin >> buffer >> num_upper >> num_lower;
        int64_t upper, start, lower;
        time_t maxDuration = 0;
        auto timeSection_u = vector<vector<pair<time_t, time_t>>>(
            num_upper, vector<pair<time_t, time_t>>());

        stringstream ss;
        ss << buffer << " " << num_upper << " " << num_lower << "\n";
        long idx = 0;
        while (fin >> upper >> start >> lower) {
            time_t end;

            float y = ((float)rand()) / ((float)RAND_MAX);
            double x1 = pow(max_d, (alpha + 1));
            double x2 = pow(min_d, (alpha + 1));
            double x3 = (x2 - x1) * y + x1;
            double x4 = pow(x3, 1 / (alpha + 1));
            time_t duration = (time_t)x4;

            maxDuration = max(maxDuration, duration);
            if (timeSection_u[upper].size() > 0)
                end = min(start + duration, timeSection_u[upper].back().first);
            else
                end = start + duration;

            timeSection_u[upper].push_back(pair<time_t, time_t>(start, end));

            ss << upper << " " << lower << " " << start << " " << end << "\n";
            idx++;
            if (idx % 500000 == 0) {
                fout << ss.str();
                ss.str("");
            }
        }
        fout << ss.str();
    }
    return 0;
}

struct Elem {
    time_t min_time = numeric_limits<time_t>::max();
    time_t max_time = 0;
    int part_id;
};

void calc(ptree::Graph g, string metis_file, int num_part)
{
    auto cmp = [](Elem& first, Elem& second) {
        if (first.min_time == second.min_time) {
            return first.max_time < second.max_time;
        }
        return first.min_time < second.min_time;
    };
    vector<Elem> elem_vec(num_part);
    for (int i = 0; i < num_part; i++) {
        elem_vec[i].part_id = i;
    }

    ifstream ifs(metis_file);
    if (not ifs.is_open()) {
        std::cout << "文件读取失败\n";
        exit(-1);
    }
    int total_size = g.num_vertices();
    vector<int> part(total_size, 0);
    vector<vector<int>> matrix(num_part, vector<int>(num_part, 0));

    int tmp;
    int cursor = 0;
    while (ifs >> tmp) {
        auto node = g[cursor];
        elem_vec[tmp].min_time = std::min(elem_vec[tmp].min_time, node.start);
        elem_vec[tmp].max_time = std::max(elem_vec[tmp].max_time, node.start);
        part[cursor++] = tmp;
    }
    sort(elem_vec.begin(), elem_vec.end(), cmp);

    uint32_t cross_edge = 0;
    for (int i = 0; i < total_size; i++) {
        auto src_part = part[i];
        for (auto edge : g.out_edges(i)) {
            auto dst_part = part[edge];
            if (src_part != dst_part) {
                matrix[src_part][dst_part]++;
                cross_edge++;
            }
        }
    }
    std::cout << "cross edge = " << cross_edge
              << " cut ratio = " << (double)cross_edge / (double)g.num_edges()
              << "\n";
    for (int i = 0; i < num_part; i++) {
        cout << i << " part min time = " << elem_vec[i].min_time << " = "
             << unixTimestampToDateString(elem_vec[i].min_time) << "\n";
        cout << i << " part max time = " << elem_vec[i].max_time << " = "
             << unixTimestampToDateString(elem_vec[i].max_time) << "\n";
    }
    cout << "\n";
    std::cout << "matrix = \n";
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < num_part; j++) {
            cout << setw(5) << matrix[elem_vec[i].part_id][elem_vec[j].part_id];
        }
        cout << "\n";
    }
    cout << "\n";
}
