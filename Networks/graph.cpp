#include "graph.h"

Graph::Graph(int64_t n) {
    num_vert = n;
    dist.resize(n, 1e6);
    L.resize(n, std::vector<double>(n, 0));
    Z.resize(n, std::vector<int>(n, 0));
    table.resize(n);
    visited.resize(n, false);
    visited_DFS.resize(n, false);
    parents.resize(n, 0);
}

void Graph::Dijkstra(int s, bool path_2) {
    dist[s] = 0;
    que.insert({0, s});
    parents[s] = -1;
    visited[s] = true;
    while (!que.empty()) {
        int v = que.begin()->second;
        que.erase(que.begin());
        for (auto &i: table[v]) {
            if (dist[i.second] > dist[v] + i.first && !visited[i.second] && i.first != 0) {
                if (path_2 && v == s && i.second == path_1[1]) {
                    continue;
                }
                if (!path_2 || path_2) {
                    que.erase({dist[i.second], i.second});
                    dist[i.second] = dist[v] + i.first;
                    que.insert({dist[i.second], i.second});
                    parents[i.second] = v;
                    visited[i.second] = true;
                }
            }
        }
    }
}

void Graph::DFS(size_t now, std::vector<size_t> &res) {
    visited_DFS[now] = true;
    res.emplace_back(now);
    for (auto neig: table[now]) {
        if (!visited_DFS[neig.second]) {
            DFS(neig.second, res);
        }
    }
}

void Graph::ConnectivityProblem() {
    std::vector<std::vector<size_t>> components;
    for (size_t i = 1; i < visited_DFS.size(); ++i) {
        if (!visited_DFS[i]) {
            std::vector<size_t> connect_comp;
            DFS(i, connect_comp);
            components.emplace_back(connect_comp);
        }
    }
    if (components.size() == 1) {
        return;
    }
    auto vert_1 = components[0][0];
    auto vert_2 = components[1][0];
    table[vert_1][0].first = 0;
    table[vert_1].push_back({1, table[vert_2][0].second});
    table[vert_2][0].first = 0;
    table[vert_2].push_back({1, table[vert_1][0].second});

}

void Graph::Preparation() {
    que.clear();
    for (size_t i = 0; i < num_vert; ++i) {
        parents[i] = 0;
        visited[i] = false;
        dist[i] = 1e6;
    }

}

void Graph::Start() {
    que.clear();
    path_1.clear();
    path_2.clear();
    for (size_t i = 0; i < num_vert; ++i) {
        parents[i] = 0;
        visited[i] = false;
        dist[i] = 1e6;
    }

}
