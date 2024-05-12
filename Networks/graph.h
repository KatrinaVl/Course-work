#include <set>
#include <vector>

class Graph {
    int num_vert;

public:
    Graph(int64_t V);

    void ConnectivityProblem();

    void DFS(size_t now, std::vector<size_t> &res);

    void Dijkstra(int s, bool path_2);

    void Preparation();

    void Start();

    std::vector<std::vector<int>> Z;
    std::vector<std::vector<double>> L;
    std::vector<int64_t> dist;
    std::vector<std::vector<std::pair<double, int>>> table;
    std::vector<bool> visited;
    std::vector<bool> visited_DFS;
    std::set<std::pair<int, int>> que;
    std::vector<int> parents;
    std::vector<int> path_1;
    std::vector<int> path_2;
    double max_lambda;
    double min_lambda = 1e6;
};
