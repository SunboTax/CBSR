
#include "graph.h"

#include <utility>
#include <vector>
namespace ptree {
Graph::Graph()
{
    vsize = 0;
    graph = GRA();
    vl = VertexList();
}

Graph::Graph(int size)
{
    vsize = size;
    vl = VertexList(size);
    graph = GRA(size, In_OutList());
}

Graph::Graph(const Graph& other)
{
    if (this != &other) {
        graph = other.graph;
        vl = other.vl;
        f = other.f;
        vsize = other.vsize;
    }
}

Graph::Graph(GRA& g, VertexList& vlist)
{
    vsize = vlist.size();
    graph = g;
    vl = vlist;
}

Graph::Graph(unordered_map<int, vector<int>>& inlist,
    unordered_map<int, vector<int>>& outlist)
{
    vsize = inlist.size();
    vl = VertexList(vsize);
    graph = GRA(vsize, In_OutList());
    for (int i = 0; i < vsize; i++)
        addVertex(i);
    unordered_map<int, vector<int>>::iterator hit, hit1;
    unordered_map<int, int> indexmap;
    vector<int> vec;
    vector<int>::iterator vit;
    int k;
    for (hit = inlist.begin(), k = 0; hit != inlist.end(); hit++, k++) {
        indexmap[hit->first] = k;
    }
    for (hit = inlist.begin(), hit1 = outlist.begin(), k = 0;
         hit != inlist.end(); hit++, hit1++, k++) {
        vec = hit->second;
        for (vit = vec.begin(); vit != vec.end(); vit++)
            graph[k].inList.push_back(indexmap[*vit]);
        vec = hit1->second;
        for (vit = vec.begin(); vit != vec.end(); vit++)
            graph[k].outList.push_back(indexmap[*vit]);
    }
}

Graph::~Graph() { }

void Graph::writeGraph(ostream& out)
{
    cout << "Graph size = " << graph.size() << endl;
    out << "graph_for_greach" << endl;
    out << vl.size() << endl;

    GRA::iterator git;
    EdgeList el;
    EdgeList::iterator eit;
    VertexList::iterator vit;
    int i = 0;
    for (i = 0; i < vl.size(); i++) {
        out << i << ": ";
        el = graph[i].outList;
        for (eit = el.begin(); eit != el.end(); eit++)
            out << (*eit) << " ";
        out << "#" << endl;
    }
    /*
            cout << "** In List for graph **" << endl;
            for (i = 0; i < vl.size(); i++) {
                    out << i << ": ";
                    el = graph[i].inList;
                    for (eit = el.begin(); eit != el.end(); eit++)
                            out << (*eit) << " ";
                    out << "#" << endl;
            }
    */
}

void Graph::addVRoot()
{
    int vid = vl.size();
    graph.push_back(In_OutList());
    vl.push_back(Vertex());
    vsize = vl.size();

    Vertex v;
    v.node = -1;
    v.start = 0;
    v.end = 0;
    v.visited = false;
    vl[vid] = v;
    vl[vid].id = vid;
    f[make_pair(v.node, v.start)] = vid;
    EdgeList il = EdgeList();
    EdgeList ol = EdgeList();
    In_OutList oil = In_OutList();
    oil.inList = il;
    oil.outList = ol;
    graph[vid] = oil;
}

bool Graph::addVertex(Vertex v)
{
    if (f.find(make_pair(v.node, v.start)) == f.end()) {
        int vid = vl.size();
        graph.push_back(In_OutList());
        vl.push_back(Vertex());
        v.id = vid;
        vl[vid] = v;
        f[make_pair(v.node, v.start)] = vid;
        vsize = vl.size();
        return true;
    }
    return false;
}

void Graph::addVertex(Vertex v, int gid)
{
    if (f.find(make_pair(v.node, v.start)) == f.end()) {
        int vid = vl.size();
        graph.push_back(In_OutList());
        vl.push_back(Vertex());
        v.id = vid;
        v.gid = gid;
        vl[vid] = v;
        f[make_pair(v.node, v.start)] = vid;
        vsize = vl.size();
    }
}

int Graph::getVertexId(Vertex v)
{
    if (f.find(make_pair(v.node, v.start)) != f.end()) {
        return f[make_pair(v.node, v.start)];
    } else {
        return -1;
    }
}

void Graph::addVertex(int node, time_t start, time_t end)
{
    if (f.find(make_pair(node, start)) == f.end()) {
        int vid = vl.size();
        graph.push_back(In_OutList());
        vl.push_back(Vertex());
        f[make_pair(node, start)] = vid;
        vsize = vl.size();

        Vertex v;
        v.visited = false;
        v.node = node;
        v.start = start;
        v.end = end;
        v.id = vid;
        vl[vid] = v;

        // EdgeList il = EdgeList();
        // EdgeList ol = EdgeList();
        // In_OutList oil = In_OutList();
        // oil.inList = il;
        // oil.outList = ol;
        // graph[vid] = oil;
    }
}

void Graph::addVertex(int vid)
{
    if (vid < 0) {
        return;
    }
    if (vid >= vl.size()) {
        int size = vl.size();
        for (int i = 0; i < (vid - size + 1); i++) {
            graph.push_back(In_OutList());
            vl.push_back(Vertex());
        }
        vsize = vl.size();
    }

    Vertex v;
    v.visited = false;
    vl[vid] = v;

    // EdgeList il = EdgeList();
    // EdgeList ol = EdgeList();
    // In_OutList oil = In_OutList();
    // oil.inList = il;
    // oil.outList = ol;
    // graph[vid] = oil;
}

void Graph::extract(unordered_map<int, vector<int>>& inlist,
    unordered_map<int, vector<int>>& outlist)
{
    for (int i = 0; i < vl.size(); i++) {
        inlist[i] = graph[i].inList;
        outlist[i] = graph[i].outList;
    }
    //	printMap(inlist,outlist);
}

void Graph::addVertex(int node, time_t start, time_t end, int num)
{
    if (f.find(make_pair(node, start)) == f.end()) {
        int vid = vl.size();
        graph.push_back(In_OutList());
        vl.push_back(Vertex());
        f[make_pair(node, start)] = vid;
        vsize = vl.size();

        Vertex v;
        v.visited = false;
        v.node = node;
        v.start = start;
        v.end = end;
        v.id = vid;
        v.pnum = num;
        vl[vid] = v;

        EdgeList il = EdgeList();
        EdgeList ol = EdgeList();
        In_OutList oil = In_OutList();
        oil.inList = il;
        oil.outList = ol;
        graph[vid] = oil;
    }
}

// 1: border out
// 2: border in
void Graph::setVertexBorderFlag(int node, int type)
{
    if (type == 1) {
        vl[node].border_out = true;
    } else if (type == 2) {
        vl[node].border_in = true;
    }
}

void Graph::addEdge(int sid, int tid)
{
    // update edge list
    if (find(graph[sid].outList.begin(), graph[sid].outList.end(), tid)
        != graph[sid].outList.end()) {
        return;
    }
    graph[tid].inList.push_back(sid);
    graph[sid].outList.push_back(tid);
}

void Graph::addEdge(pair<int, time_t> a, pair<int, time_t> b)
{
    int sid = cnum(a);
    int tid = cnum(b);
    addEdge(sid, tid);
    // graph[tid].inList.push_back(sid);
    // graph[sid].outList.push_back(tid);
}

// 拓扑排序，对于不存在偏序关系的两个节点，时间增序排列
// 直接对拓扑排序后的节点数组按照start增序重排序即可：对start的排序不会打乱原有拓扑顺序
vector<int>
Graph::topological_sort()
{
    vector<int> order;
    queue<int> q;
    vector<int> ind(vsize, 0);
    for (int i = 0; i < vsize; i++) {
        ind[i] = this->in_degree(i);
        if (ind[i] == 0)
            q.push(i);
    }

    while (q.size() > 0) {
        int v = q.front();
        q.pop();
        order.push_back(v);
        for (int j : this->out_edges(v)) {
            ind[j]--;
            if (ind[j] == 0)
                q.push(j);
        }
    }

    // sort(order.begin(), order.end(), [this](int &a, int &b) { return
    // (*this)[a].start < (*this)[b].start; });

    return order;
}

int Graph::num_vertices() const
{
    return vl.size();
}

int Graph::num_edges()
{
    EdgeList el;
    GRA::iterator git;
    int num = 0;
    for (git = graph.begin(); git != graph.end(); git++) {
        el = git->outList;
        num += el.size();
    }
    return num;
}

// return out edges of specified vertex
EdgeList&
Graph::out_edges(int src)
{
    return graph[src].outList;
}

// return in edges of specified vertex
EdgeList&
Graph::in_edges(int trg)
{
    return graph[trg].inList;
}

int Graph::out_degree(int src)
{
    return graph[src].outList.size();
}

int Graph::in_degree(int trg)
{
    return graph[trg].inList.size();
}

int Graph::cnum(pair<int, time_t> p)
{
    return f[p];
}

// get roots of graph (root is zero in_degree vertex)
vector<int>
Graph::getRoots()
{
    vector<int> roots;
    GRA::iterator git;
    int i = 0;
    for (git = graph.begin(), i = 0; git != graph.end(); git++, i++) {
        if (git->inList.size() == 0)
            roots.push_back(i);
    }

    return roots;
}

// check whether the edge from src to trg is in the graph
bool Graph::hasEdge(int src, int trg)
{
    EdgeList el = graph[src].outList;
    EdgeList::iterator ei;
    for (ei = el.begin(); ei != el.end(); ei++)
        if ((*ei) == trg)
            return true;
    return false;
}

// return vertex list of graph
VertexList&
Graph::vertices()
{
    return vl;
}

Graph&
Graph::operator=(const Graph& g)
{
    if (this != &g) {
        graph = g.graph;
        vl = g.vl;
        vsize = g.vsize;
        f = g.f;
    }
    return *this;
}

Graph&
Graph::operator=(Graph&& other)
{
    if (this != &other) {
        graph = std::move(other.graph);

        // 移动 vl
        vl = std::move(other.vl);

        // 移动 vsize
        vsize = other.vsize;

        // 重置源对象的成员变量
        other.vsize = 0;
    }
    return *this;
}

// get a specified vertex property
Vertex&
Graph::operator[](const int vid)
{
    return vl[vid];
}

// Graph::Graph(unordered_map<int, vector<int>> &inlist, unordered_map<int,
// vector<int>> &outlist)
// {
// 	vsize = inlist.size();
// 	vl = VertexList(vsize);
// 	graph = GRA(vsize, In_OutList());
// 	for (int i = 0; i < vsize; i++)
// 		addVertex(i);
// 	unordered_map<int, vector<int>>::iterator hit, hit1;
// 	unordered_map<int, int> indexmap;
// 	vector<int> vec;
// 	vector<int>::iterator vit;
// 	int k;
// 	for (hit = inlist.begin(), k = 0; hit != inlist.end(); hit++, k++)
// 	{
// 		indexmap[hit->first] = k;
// 	}
// 	for (hit = inlist.begin(), hit1 = outlist.begin(), k = 0; hit !=
// inlist.end(); hit++, hit1++, k++)
// 	{
// 		vec = hit->second;
// 		for (vit = vec.begin(); vit != vec.end(); vit++)
// 			graph[k].inList.push_back(indexmap[*vit]);
// 		vec = hit1->second;
// 		for (vit = vec.begin(); vit != vec.end(); vit++)
// 			graph[k].outList.push_back(indexmap[*vit]);
// 	}
// }

// void Graph::extract(unordered_map<int, vector<int>> &inlist,
// unordered_map<int, vector<int>> &outlist)
// {
// 	for (int i = 0; i < vl.size(); i++)
// 	{
// 		inlist[i] = graph[i].inList;
// 		outlist[i] = graph[i].outList;
// 	}
// 	//	printMap(inlist,outlist);
// }

// // for test
// void Graph::printMap(unordered_map<int, vector<int>> &inlist,
// unordered_map<int, vector<int>> &outlist)
// {
// 	cout << "==============================================" << endl;
// 	unordered_map<int, vector<int>>::iterator hit;
// 	vector<int>::iterator vit;
// 	for (hit = outlist.begin(); hit != outlist.end(); hit++)
// 	{
// 		cout << hit->first << ": ";
// 		vector<int> vec = hit->second;
// 		for (vit = vec.begin(); vit != vec.end(); vit++)
// 			cout << *vit << " ";
// 		cout << "#" << endl;
// 	}
// 	cout << "In List for graph" << endl;
// 	for (hit = inlist.begin(); hit != inlist.end(); hit++)
// 	{
// 		cout << hit->first << ": ";
// 		vector<int> vec = hit->second;
// 		for (vit = vec.begin(); vit != vec.end(); vit++)
// 			cout << *vit << " ";
// 		cout << "#" << endl;
// 	}
// 	cout << "================================================" << endl;
// }

int subGraph::get_gid()
{
    return gid;
}

void subGraph::set_gid(int gid_)
{
    gid = gid_;
}

// for gripp
vector<int>*
Graph::GetChild(int node)
{
    // return (*adj_listp)[node];
    //  vector<unsigned>* ret = new vector<unsigned>((*adj_listp)[node]);
    vector<int>* ret = new vector<int>(this->out_edges(node));

    return ret;
}

bool Graph::CanReach(unsigned nodeA, unsigned nodeB)
{
    std::vector<uint8_t> vis(this->num_vertices(), 0);
    return search(nodeA, nodeB, vis);
}
bool Graph::search(unsigned int nodeA, unsigned int nodeB,
    std::vector<uint8_t>& vis)
{
    if (nodeA == nodeB)
        return true;
    vis[nodeA] = 1;
    // int len = (*adj_listp)[nodeA].size();
    int len = this->out_degree(nodeA);
    EdgeList& edge = this->out_edges(nodeA);
    for (int i = 0; i < len; i++) {
        // if(!vis[(*adj_listp)[nodeA][i]] &&
        // 	search((*adj_listp)[nodeA][i],nodeB))
        // 		return true;
        if (!vis[edge[i]] and search(edge[i], nodeB, vis)) {
            return true;
        }
    }
    return false;
}

#ifdef GRAIL
const double
Graph::actualgap(const int vid)
{
    return vl[vid].mingap;
    //	return vl[vid].mingap - vl[vid].tcs;
}
const double
Graph::tcs(const int vid)
{
    return vl[vid].tcs;
}

bool Graph::incrementalContains(int src, int trg, int cur)
{
    int i;
    for (i = 0; i < cur; i++) {
        if (vl[src].pre->at(i) > vl[trg].pre->at(i))
            return false;
        if (vl[src].post->at(i) < vl[trg].post->at(i))
            return false;
    }
    return true;
}

bool Graph::contains(int src, int trg, int dim)
{
    int i;
    for (i = 0; i < dim; i++) {
        if (vl[src].pre->at(i) > vl[trg].pre->at(i))
            return false;
        if (vl[src].post->at(i) < vl[trg].post->at(i))
            return false;
    }
    return true;
}

#endif

} // namespace ptree