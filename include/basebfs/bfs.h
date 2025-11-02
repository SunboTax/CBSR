#ifndef BFS_H
#define BFS_H
#include <vector>

#include "bigraph.h"
#include "utils.h"

void biBFS(BiGraph& bg, vector<CCR::queryInfo>& queryInfo,
    vector<int>& queryRes);

bool biBFSWorker(BiGraph& bg, CCR::queryInfo info);

#endif