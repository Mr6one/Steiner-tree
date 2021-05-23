#include <iostream>
#include <vector>
#include <deque>
#include <set>
#include <algorithm>
#include <limits>
#include <chrono>
#include<dirent.h>
#include<string.h>
#include <fstream>

using namespace std;

const int MAX_INT = numeric_limits<int>::max();
enum { REQUIRED_POINT, STEINER_POINT };

class DisjointSet {
public:
    DisjointSet(int n) {
        parent.resize(n);
        set_sizes.resize(n, 0);
    }
    void MakeSet(int);
    void UnionSets(int, int);
    int FindSet(int);
private:
    vector<int> parent;
    vector<int> set_sizes;
};

void DisjointSet::MakeSet(int x) {
    parent[x] = x;
    set_sizes[x] += 1;
}

void DisjointSet::UnionSets(int x, int y) {
    int x_parent = FindSet(x);
    int y_parent = FindSet(y);
    if (x_parent != y_parent) {
        if (set_sizes[x_parent] < set_sizes[y_parent]) {
            parent[x_parent] = y_parent;
            set_sizes[y_parent] += set_sizes[x_parent];
        } else {
            parent[y_parent] = x_parent;
            set_sizes[x_parent] += set_sizes[y_parent];
        }
    }
}

int DisjointSet::FindSet(int x) {
    if (parent[x] == x) {
        return x;
    } else {
        parent[x] = FindSet(parent[x]);
        return parent[x];
    }
}

struct Edge {
    Edge(int u, int v, int weight): u(u), v(v), weight(weight){}
    int u;
    int v;
    int weight;

    bool operator < (const Edge& edge) const {
        return weight < edge.weight;
    }
};

class Graph {
public:
    Graph(int n) {
        graph.resize(n);
    }

    void AddVertex(int, int);
    void AddEdge(int, int, int);
    vector<int> GetAdjacentNodes(int);
    deque<Edge> GetEdges();
    size_t GetSize();
    bool IsConnected(int);

private:
    vector<vector<int>> graph;
    deque<Edge> edges;
};

void Graph::AddVertex(int u, int v) {
    graph[u].push_back(v);
    graph[v].push_back(u);
}

vector<int> Graph::GetAdjacentNodes(int node) {
    return graph[node];
}

void Graph::AddEdge(int u, int v, int weight) {
    edges.push_back(Edge(u, v, weight));
    edges.push_back(Edge(v, u, weight));
}

deque<Edge> Graph::GetEdges() {
    return edges;
}

size_t Graph::GetSize() {
    return graph.size();
}

void DFS(Graph& graph, int vertex, vector<bool>& visited) {
    visited[vertex] = true;
    for (int child: graph.GetAdjacentNodes(vertex)) {
        if (!visited[child]) {
            DFS(graph, child, visited);
        }
    }
}

bool Graph::IsConnected(int vertex = 0) {
    vector<bool> visited(graph.size(), false);
    DFS(*this, vertex, visited);
    for (bool visit: visited) {
        if (!visit) {
            return false;
        }
    }
    return true;
}

int FindMST(Graph& graph) {
    int MST_weight = 0;
    DisjointSet disjoint_set(graph.GetSize());
    for (size_t i = 0; i < graph.GetSize(); ++i) {
        disjoint_set.MakeSet(i);
    }
    deque<Edge> remaining_edges = graph.GetEdges();
    sort(remaining_edges.begin(), remaining_edges.end());

    while (!remaining_edges.empty()) {
        Edge current_edge = remaining_edges[0];
        if (disjoint_set.FindSet(current_edge.u) != disjoint_set.FindSet(current_edge.v)) {
            MST_weight += current_edge.weight;
            disjoint_set.UnionSets(current_edge.u, current_edge.v);
        }
        remaining_edges.pop_front();
    }
    return MST_weight;
}

vector<int> FindMinPath(Graph& graph, int vertex) {
    vector<vector<int>> edges(graph.GetSize());
    for (vector<int>& edge: edges) {
        edge.resize(graph.GetSize());
    }

    for (Edge& edge: graph.GetEdges()) {
        edges[edge.u][edge.v] = edge.weight;
    }

    set<pair<int, int>> currentDist;
    vector<int> min_dist(graph.GetSize(), MAX_INT);

    currentDist.insert(make_pair(0, vertex));
    min_dist[vertex] = 0;
    while (!currentDist.empty()) {
        int next_vertex = currentDist.begin()->second;
        currentDist.erase(currentDist.begin());

        for (int adjacent_node: graph.GetAdjacentNodes(next_vertex)) {
            if (min_dist[next_vertex] + edges[next_vertex][adjacent_node] < min_dist[adjacent_node]) {
                currentDist.erase(make_pair(min_dist[adjacent_node], adjacent_node));
                min_dist[adjacent_node] = min_dist[next_vertex] + edges[next_vertex][adjacent_node];
                currentDist.insert(make_pair(min_dist[adjacent_node], adjacent_node));
             }
        }
    }

    return min_dist;
}

vector<vector<int>> FindMinPathes(Graph& graph) {
    vector<vector<int>> distances;
    for (size_t i = 0; i < graph.GetSize(); ++i) {
        distances.push_back(FindMinPath(graph, i));
    }
    return distances;
}

Graph BuildMetricGraph(Graph& graph, vector<int>& points, vector<int>& mapping, int required_points_num) {
    vector<vector<int>> distances = FindMinPathes(graph);

    Graph metric_graph(required_points_num);
    for (size_t u = 0; u < distances.size(); ++u) {
        for (size_t v = u + 1; v < distances.size(); ++v) {
            if (points[u] == REQUIRED_POINT && points[v] == REQUIRED_POINT) {
                metric_graph.AddVertex(mapping[u], mapping[v]);
                metric_graph.AddEdge(mapping[u], mapping[v], distances[u][v]);
            }
        }
    }

    return metric_graph;
}

void FindSubsets(set<set<int>>& subsets, vector<int>& values, size_t pos) {
    if (pos == values.size()) {
        return;
    }
    set<set<int>> new_set;
    for (set<int> subset: subsets) {
        subset.insert(values[pos]);
        new_set.insert(subset);
    }
    new_set.insert({values[pos]});
    new_set.insert(subsets.begin(), subsets.end());
    subsets = new_set;
    FindSubsets(subsets, values, pos + 1);
}

set<set<int>> FindSubsets(vector<int>& values) {
    set<set<int>> subsets;
    FindSubsets(subsets, values, 0);
    return subsets;
}

int FindOpt(Graph& graph, vector<int>& points) {
    vector<vector<int>> edges(graph.GetSize());
    for (vector<int>& edge: edges) {
        edge.resize(graph.GetSize(), -1);
    }

    for (Edge& edge: graph.GetEdges()) {
        edges[edge.u][edge.v] = edge.weight;
    }

    vector<int> steiners_points;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i] == STEINER_POINT) {
            steiners_points.push_back(i);
        }
    }
    set<set<int>> subsets = FindSubsets(steiners_points);

    int min_weight = FindMST(graph);
    for (const set<int>& subset: subsets) {
        vector<int> verticies;
        for (size_t i = 0; i < points.size(); ++i) {
            if (points[i] == REQUIRED_POINT || subset.find(i) == subset.end()) {
                verticies.push_back(i);
            }
        }

        vector<int> mapping(graph.GetSize(), -1);
        for (size_t i = 0; i < verticies.size(); ++i) {
            mapping[verticies[i]] = i;
        }
        Graph new_graph(verticies.size());
        for (size_t u = 0; u < edges.size(); ++u) {
            for (size_t v = u + 1; v < edges[u].size(); ++ v) {
                if (mapping[u] != -1 && mapping[v] != -1 && edges[u][v] != -1) {
                    new_graph.AddVertex(mapping[u], mapping[v]);
                    new_graph.AddEdge(mapping[u], mapping[v], edges[u][v]);
                }
            }
        }
        if (new_graph.IsConnected()) {
            int MST_weight = FindMST(new_graph);
            if (MST_weight < min_weight) {
                min_weight = MST_weight;
            }
        }
    }
    return min_weight;
}

const string input_dir_name = "./Tests/1";
const string output_dir_name = "./Results";

const string output_file_name = "1.csv";

int main()
{
    DIR* pDIR;
    struct dirent* entry;
    if((pDIR = opendir(input_dir_name.c_str()))) {
        while((entry = readdir(pDIR))) {
            if(strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) {
                cout << entry->d_name << endl;
                ifstream input(input_dir_name + "/" + entry->d_name);
                std::cin.rdbuf(input.rdbuf());

                size_t vertex_numbers, edge_numbers;
                cin >> vertex_numbers >> edge_numbers;
                cout << vertex_numbers << ' ' << edge_numbers << endl;

                Graph graph(vertex_numbers);

                for (size_t i = 0; i < edge_numbers; ++i) {
                    char ch;
                    cin >> ch;
                    int u, v;
                    cin >> u >> v;

                    u = u % vertex_numbers;
                    v = v % vertex_numbers;
                    graph.AddVertex(u, v);

                    int weight;
                    cin >> weight;
                    graph.AddEdge(u, v, weight);
                }
                if (!graph.IsConnected()) {
                    cout << "Graph is not connected";
                    return 0;
                }

                vector<int> points(vertex_numbers, STEINER_POINT);
                vector<int> mapping(vertex_numbers, -1);
                int required_points_num;
                cin >> required_points_num;
                for (int i = 0; i < required_points_num; ++i) {
                    char ch;
                    cin >> ch;
                    int u;
                    cin >> u;
                    u = u % vertex_numbers;
                    points[u] = REQUIRED_POINT;
                    mapping[u] = i;
                }

                chrono::steady_clock::time_point time_begin = chrono::steady_clock::now();
                Graph metric_graph = BuildMetricGraph(graph, points, mapping, required_points_num);
                int MST_weight = FindMST(metric_graph);
                chrono::steady_clock::time_point time_end = chrono::steady_clock::now();
                double time_diff = (chrono::duration_cast<chrono::microseconds>(time_end - time_begin).count()) /1000000.0;

                cout << vertex_numbers << " | ";
                cout << edge_numbers << " | ";
                cout << required_points_num << " | ";
                cout << MST_weight << " | ";
                cout << time_diff << " | ";

                ofstream output;
                output.open(output_dir_name + "/" + output_file_name, fstream::app);
                output << to_string(vertex_numbers) + ";" + to_string(edge_numbers) + ";" +
                    to_string(required_points_num) + ";" + to_string(MST_weight) + ";" + to_string(time_diff) + "\n";
                output.close();
            }
        }
        closedir(pDIR);
    }
    /*
    time_begin = chrono::steady_clock::now();
    int min_weight = FindOpt(graph, points);
    time_end = chrono::steady_clock::now();
    double time_diff_2 = (chrono::duration_cast<chrono::microseconds>(time_end - time_begin).count()) /1000000.0;

    cout << min_weight << " | ";
    cout << time_diff_2  << endl;
    */
    return 0;
}
