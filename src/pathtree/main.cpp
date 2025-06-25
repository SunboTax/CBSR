// #include <stdio.h>
// #include <stdlib.h>
// #include <sys/time.h>

// #include "graph.h"
// #include "pathtree/DWGraph.h"
// #include "pathtree/DWGraphUtil.h"
// #include "pathtree/GraphUtil.h"
// #include "pathtree/PathTree.h"

// //#define TEST

// static void usage() {
//     cout << "\nUsage:\n"
//             "	pathtree [-h] [-n num_query] [-t alg_type]  [-c chainfile] [-r "
//             "resultfile] filename\n"
//             "Description:\n"
//             "	-h	Print the help message.\n"
//             "	-p	Post-processing (Data Compression).\n"
//             "	-n	Set the total number of random queries. The default "
//             "value is 100,000.\n"
//             "	-c	Set the input chain(path) file.\n"
//             "	-r	Set the result filename ('../results.txt' is default "
//             "file).\n"
//             "	-t	Algorithm type(PTree-1 is default algorithm).\n"
//             "		1:	PTree-1\n"
//             "		2:	PTree-2 (chain decomposition).\n"
//             "		3:	PTree-2 (path decomposition).\n"
//             "		4:	PTree-2 .\n"
//             "	-d	debug mode (test reachability).\n"
//             "	filename	The graph file.\n"
//          << endl;
// }

// static void keepResult(char* resultFileName, char* version, char* filename,
//                        float lt, float qt, int isize, double cr) {
//     ofstream out(resultFileName, ios_base::out | ios_base::app);
//     out << version << "\t" << filename << "\t" << lt << "\t" << qt << "\t"
//         << isize << "\t" << cr << endl;
//     out.close();
// }

// template <typename T>
// string to_string(const T& value) {
//     ostringstream oss;
//     oss << value;
//     return oss.str();
// };

// int main4() {
//     Graph g1;
//     for (int i = 0; i < 7; i++) g1.addVertex(i);
//     g1.addEdge(0, 1);
//     g1.addEdge(1, 2);
//     g1.addEdge(2, 3);
//     g1.addEdge(3, 4);
//     g1.addEdge(4, 2);
//     g1.addEdge(4, 5);
//     g1.addEdge(5, 6);
//     g1.addEdge(6, 1);
//     g1.addEdge(6, 4);
//     g1.printGraph();
//     int* on = new int[g1.num_vertices()];
//     vector<int> ts;
//     GraphUtil::mergeSCC(g1, on, ts);
//     g1.printGraph();
//     cout << "********************************" << endl;

//     DWGraph g, branch;
//     for (int i = 0; i < 4; i++) g.addVertex(i);
//     g.addEdge(0, 1, 1, 0);
//     g.addEdge(1, 2, 2, 1);
//     g.addEdge(2, 3, 6, 2);
//     g.addEdge(3, 4, 6, 3);
//     g.addEdge(4, 2, 50, 4);
//     /*
//     g.addEdge(4,5,5,5);
//     g.addEdge(5,6,5,6);
//     g.addEdge(6,1,2,7);
//     g.addEdge(6,4,15,8);
//     */
//     g.printGraph();
//     cout << "================================================" << endl;

//     DWGraphUtil::findMaxBranching(g, branch);
//     branch.printGraph();
// }

// int main3() {
//     int n = 15;
//     double c = 5;
//     string filename;
//     for (int i = 0; i < 20; i++) {
//         filename = "./rand15_5/rand15_" + to_string(i + 1) + ".gra";
//         DWGraphUtil::genRandomGraph(n, c, filename.c_str());
//     }
// }

// #ifdef TEST
// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         usage();
//         return 1;
//     }
//     char* filename = argv[1];
//     ifstream infile(filename);
//     if (!infile) {
//         cout << "Error: Cannot open " << filename << endl;
//         return -1;
//     }
//     DWGraph g(infile);
//     g.printGraph();
//     cout << "================================================" << endl;
//     DWGraph branch;
//     //	DWGraphUtil::processBranch(g, branch);
//     //	DWGraphUtil::findMaxBranching(g, branch);
//     bool res = DWGraphUtil::checkBranching(g, branch);
//     cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
//     branch.printGraph();
//     if (!res) cerr << "\nSome errors found!" << endl;
// }
// #else
// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         usage();
//         return 1;
//     }

//     int i = 1;
//     int query_num = 100000;
//     int alg_type = 1;
//     char* filename = NULL;
//     char* chainfile = NULL;
//     char* resfilename = "../results.txt";
//     bool debug = false;
//     bool compress = false;
//     bool mcc = true;
//     bool allpairs = false;
//     while (i < argc) {
//         if (strcmp("-h", argv[i]) == 0) {
//             usage();
//             return 1;
//         }

//         if (strcmp("-d", argv[i]) == 0) {
//             i++;
//             debug = true;
//         } else if (strcmp("-p", argv[i]) == 0) {
//             i++;
//             compress = true;
//         } else if (strcmp("-m", argv[i]) == 0) {
//             i++;
//             mcc = false;
//         } else if (strcmp("-n", argv[i]) == 0) {
//             i++;
//             query_num = atoi(argv[i++]);
//         } else if (strcmp("-c", argv[i]) == 0) {
//             i++;
//             chainfile = argv[i++];
//         } else if (strcmp("-r", argv[i]) == 0) {
//             i++;
//             resfilename = argv[i++];
//         } else if (strcmp("-t", argv[i]) == 0) {
//             i++;
//             alg_type = atoi(argv[i++]);
//             //	if (alg_type < 1 || alg_type > 3)
//             //		alg_type = 1;
//         } else if (strcmp("-a", argv[i]) == 0) {
//             i++;
//             allpairs = true;
//         } else {
//             filename = argv[i++];
//         }
//     }

//     ifstream infile(filename);
//     if (!infile) {
//         cout << "Error: Cannot open " << filename << endl;
//         return -1;
//     }

//     ifstream cfile;
//     if (chainfile != NULL) {
//         cfile.open(chainfile, ifstream::in);
//         if (!cfile) {
//             cout << "Error: Cannot open " << chainfile << endl;
//             return -1;
//         }
//     }

//     cout << "INPUT: query_num=" << query_num << " alg_type=" << alg_type
//          << " filename=" << filename << endl;

//     Graph g(infile);
//     cout << "#vertex size:" << g.num_vertices()
//          << "\t#edges size:" << g.num_edges() << endl;

//     int s, t;
//     int left = 0;
//     int gsize = g.num_vertices();

//     bool r;
//     struct timeval after_time, before_time;
//     float labeling_time, query_time;
//     int* sccmap = new int[gsize];  // store pair of orignal vertex and
//                                    // corresponding vertex in merged graph
//     vector<int> reverse_topo_sort;

//     // merge strongly connected component
//     if (mcc) {
//         cout << "merging strongly connected component..." << endl;
//         gettimeofday(&before_time, NULL);
//         GraphUtil::mergeSCC(g, sccmap, reverse_topo_sort);
//         gettimeofday(&after_time, NULL);
//         query_time = (after_time.tv_sec - before_time.tv_sec) * 1000.0 +
//                      (after_time.tv_usec - before_time.tv_usec) * 1.0 / 1000.0;
//         cout << "merging time:" << query_time << " (ms)" << endl;
//     } else {
//         GraphUtil::topological_sort(g, reverse_topo_sort);
//         for (int i = 0; i < gsize; i++) sccmap[i] = i;
//     }
//     cout << "#DAG vertex size:" << g.num_vertices()
//          << "\t#DAG edges size:" << g.num_edges() << endl;

//     //	g.printGraph();

//     // generate queries
//     srand48(time(NULL));
//     cout << "generating queries..." << endl;
//     vector<int> src;
//     vector<int> trg;
//     vector<int>::iterator sit, tit;
//     if (!allpairs) {
//         while (left < query_num) {
//             s = lrand48() % gsize;
//             t = lrand48() % gsize;
//             //	if (s == t) continue;
//             src.push_back(s);
//             trg.push_back(t);
//             ++left;
//         }
//     } else {
//         for (int i = 0; i < gsize; i++) {
//             for (int j = 0; j < i; j++) {
//                 src.push_back(i);
//                 trg.push_back(j);
//             }
//         }
//     }

//     PathTree pt(g, reverse_topo_sort);

//     // for test reach
//     if (debug) pt.compute_tcm();

//     // create labels
//     if (alg_type == 1)
//         cout << "PTree-1 ";
//     else
//         cout << "PTree-2 ";
//     cout << "starting creating labels..." << endl;
//     gettimeofday(&before_time, NULL);
//     pt.createLabels(alg_type, cfile, compress);
//     gettimeofday(&after_time, NULL);
//     labeling_time = (after_time.tv_sec - before_time.tv_sec) * 1000.0 +
//                     (after_time.tv_usec - before_time.tv_usec) * 1.0 / 1000.0;
//     cout << "#construction time:" << labeling_time << " (ms)" << endl;

//     // process queries
//     cout << "process queries..." << endl;
//     gettimeofday(&before_time, NULL);
//     int source, target;

//     if (!compress) {
//         if (debug) {
//             cout << "Query and Test Reachability" << endl;
//             for (sit = src.begin(), tit = trg.begin(); sit != src.end();
//                  ++sit, ++tit) {
//                 s = sccmap[*sit];
//                 t = sccmap[*tit];
//                 r = pt.test_reach(s, t);
//                 if (!r) exit(0);
//             }
//         } else {
//             for (sit = src.begin(), tit = trg.begin(); sit != src.end();
//                  ++sit, ++tit) {
//                 s = sccmap[*sit];
//                 t = sccmap[*tit];
//                 r = pt.reach(s, t);
//             }
//         }
//     } else {
//         if (debug) {
//             cout << "Query and Test Reachability" << endl;
//             for (sit = src.begin(), tit = trg.begin(); sit != src.end();
//                  ++sit, ++tit) {
//                 s = sccmap[*sit];
//                 t = sccmap[*tit];
//                 r = pt.test_reach_dc(s, t);
//                 if (!r) exit(0);
//             }
//         } else {
//             for (sit = src.begin(), tit = trg.begin(); sit != src.end();
//                  ++sit, ++tit) {
//                 s = sccmap[*sit];
//                 t = sccmap[*tit];
//                 r = pt.reach_dc(s, t);
//             }
//         }
//     }

//     gettimeofday(&after_time, NULL);
//     query_time = (after_time.tv_sec - before_time.tv_sec) * 1000.0 +
//                  (after_time.tv_usec - before_time.tv_usec) * 1.0 / 1000.0;
//     cout << "#total query running time:" << query_time << " (ms)" << endl;

//     int* ind_size = new int[2];
//     pt.index_size(ind_size);
//     cout << "#transitive closure size = " << ind_size[1] << endl;
//     double cr = pt.cover_ratio();
//     cout << "#Path-tree cover ratio = " << cr << endl;
//     double comp_ratio = pt.compress_ratio();
//     cout << "#TC compression ratio = " << comp_ratio << endl;
//     if (!compress) {
//         if (alg_type == 1)
//             keepResult(resfilename, "PTree-1", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//         else if (alg_type == 2)
//             keepResult(resfilename, "PTree-2(chain)", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//         else if (alg_type == 3)
//             keepResult(resfilename, "PTree-2(path)", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//         else
//             keepResult(resfilename, "PTree-2", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//     } else {
//         if (alg_type == 1)
//             keepResult(resfilename, "PTree-1C", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//         else if (alg_type == 2)
//             keepResult(resfilename, "PTree-2C(chain)", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//         else if (alg_type == 3)
//             keepResult(resfilename, "PTree-2C(path)", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//         else
//             keepResult(resfilename, "PTree-2C", filename, labeling_time,
//                        query_time, ind_size[1], comp_ratio);
//     }
// }
// #endif
