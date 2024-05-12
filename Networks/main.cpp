#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>

#include "graph.h"

double Random(double a, double b) {
    return rand() * (b - a) / RAND_MAX + a;
}

std::vector<double> FindLambda(std::vector<double> &f, std::vector<double> &v, double &l, double &F) {
    double lambda;
    double step = 0.01;
    double res_max = 1e6;
    std::vector<double> res_vec;
    for (double i = 0.1; i <= 1; i += step) {
        lambda = i;
        std::vector<double> res(f.size(), 0);
        for (size_t j = 0; j < f.size(); ++j) {
            res[j] = (1 - lambda) * f[j] + lambda * v[j];
        }
        double k = 0;
        double f_m = 0;
        for (size_t j = 0; j < res.size(); ++j) {
            k += res[j];
            f_m = std::max(res[j], f_m);
        }
        if (f_m < res_max && f_m != 0) {
            res_vec = res;
            l = lambda;
            res_max = f_m;
        }
    }
    F = res_max;
    return res_vec;
}

int DoDijkstra(Graph &g, int s, int t, std::vector<int> &vec, bool path_2) {
    g.Dijkstra(s, path_2);
    if (g.dist[t] == 1e6) {
        std::cout << "We can't find shortest path" << std::endl;
        return -1;
    }
    for (int v = t; v != s; v = g.parents[v]) {
        vec.push_back(v);
    }
    vec.push_back(s);
    reverse(vec.begin(), vec.end());
    return 0;
}

void DoFlow(std::vector<int> &path, std::vector<double> &f, int n, double &need_lambda) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        double w;
        if (i != path.size() - 2) {
            w = Random(0, need_lambda);
            need_lambda -= w;
        } else {
            w = need_lambda;
        }
        int to = path[i];
        int from = path[i + 1];
        f[to * n + from] = w;

    }
}


int main() {
    int64_t n = 8;
    int num_arcs = 0;
    Graph g(n);
    double SumLambd = 0;
    std::fstream in_Z;
//    example: in_Z.open("/Users/ekaterina/Desktop/Курсач/matrix_name_Z.txt");
    in_Z.open("/Users/ekaterina/Desktop/Курсач/dics_Z.txt");
    if (!in_Z.is_open()) {
        std::cout << "There is a problem with open file test_2.txt" << std::endl;
    }

    std::fstream in_L;
//    example: in_L.open("/Users/ekaterina/Desktop/Курсач/matrix_name_code.txt");
    in_L.open("/Users/ekaterina/Desktop/Курсач/dics_code.txt");
    if (!in_L.is_open()) {
        std::cout << "There is a problem with open file test.txt" << std::endl;
    }
    int x;
    double y;
    // найти матрицу лямбд
    int j = 0;
    while (j < n * n) {
        in_L >> y;
        g.L[j / n][j % n] = y;
        SumLambd += y;
        g.max_lambda = std::max(g.max_lambda, y);
        if (y != 0) {
            g.min_lambda = std::min(g.min_lambda, y);
        }
        ++j;
    }

    int i = 0;
    while (i < n * n) {
        in_Z >> x;
        g.Z[i / n][i % n] = x;
        ++num_arcs;
        if (x != 0) {
            g.table[i / n].push_back({1, i % n});
        }
        ++i;
    }
    // проверить граф на связность
    // если он не связен сделать обмен ссылками
    g.ConnectivityProblem();
    double F_sum = 0;
    double F_max = 0;
    double F_min = 1e6;
    for (int s = 0; s < 8; ++s) {
        for (int t = 0; t < 8; ++t) {
            g.Start();
            double need_lambda = g.L[s][t];
            if (s == t) {
                continue;
            }

            if (DoDijkstra(g, s, t, g.path_1, false) == -1) {
                std::cout << "We can't find first path for " << s << " " << t << std::endl;
                continue;
            }

            // we found first short path
            std::vector<double> f(num_arcs, 0);
            DoFlow(g.path_1, f, n, need_lambda);
            g.Preparation();

            if (DoDijkstra(g, s, t, g.path_2, true) == -1) {
                std::cout << "We can't find second path for " << s << " " << t << std::endl;
                double f_m;
                f_m = g.L[s][t] / (g.path_1.size() - 1);
                F_sum += f_m;
                continue;
            }

            std::vector<double> v(num_arcs, 0);

            need_lambda = g.L[s][t];
            DoFlow(g.path_2, v, n, need_lambda);

            double l;
            double F;
            std::vector<double> weights = FindLambda(f, v, l, F);

            std::cout << "s and t: " << s << " " << t << std::endl;
            std::cout << "path 1: ";
            for (auto i: g.path_1) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
            std::cout << "path 2: ";
            for (auto i: g.path_2) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
            std::cout << "This is F and lambda on this (s, t): " << F << " " << l << std::endl;
            F_sum += F;
            F_max = std::max(F, F_max);
            F_min = std::min(F, F_min);
            std::cout << std::endl;

            // я нашла такую лямбда и такие веса, при которых max_f_ij минимально

        }
    }
    std::cout << "F_sum: " << F_sum << std::endl;

    std::cout << "F_average: " << F_sum / SumLambd << " " << std::endl;
    std::cout << "F_max: " << F_max << std::endl;
    std::cout << "F_min: " << F_min << std::endl;
    in_Z.close();
    in_L.close();
    return 0;
}
