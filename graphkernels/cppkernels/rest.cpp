/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#include "rest.h"

#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>
#include <utility>

using std::vector;
using std::pair;

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

auto order_by_labels(const vector<int>& labels) {
    vector<pair<int, int>> map;
    map.reserve(labels.size());

    auto idx = 0;
    for (const auto label : labels) {
        map.emplace_back(label, idx++);
    }

    sort(map.begin(), map.end());
    return map;
}

auto compute_valid_vertex_pairs(
        const vector<pair<int, int>>& map1,
        const vector<pair<int, int>>& map2) {
    vector<pair<int, int>> pairs;
    pairs.reserve(map1.size() * map2.size());

    const auto comp = [](const auto& p_a, const auto& p_b){
        return p_a.first < p_b.first;
    };

    auto p = map2.cbegin();  // Memoize low limit (see below).
    for (auto i1 = map1.cbegin(); i1 != map1.cend(); ) {
        // Find range of map2 that contains vertices labelled "label1".
        auto [eq_cbegin, eq_cend] = std::equal_range(p, map2.cend(), *i1, comp);

        // Iterate over all equal values in map1.
        const auto label1 = i1->first;
        do {
            const auto num1 = i1->second;

            // Create all pairs between vertex of map1 and range of map2.
            for (auto p = eq_cbegin; p != eq_cend; ++p) {
                pairs.emplace_back(num1, p->second);
            }

            ++i1;
        } while (i1 != map1.cend() && i1->first == label1);

        // All vertices with that label have been exhausted in both maps.
        p = eq_cend;
    }

    sort(pairs.begin(), pairs.end());
    return pairs;
}

auto productAdjacency(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label) {
    // Step 0: Order vertices by labels.
    const auto map1 = order_by_labels(v1_label);
    const auto map2 = order_by_labels(v2_label);

    // Step 1: Compute all valid vertex pairs (vertices of product graph).
    const auto pairs = compute_valid_vertex_pairs(map1, map2);

    // Step 2: Compute new labels for vertices of the product graph.
    Eigen::Matrix<int, -1, -1, Eigen::RowMajor> H(map1.size(), map2.size());

    auto next_label = 0;
    for (const auto& [v1, v2] : pairs) {
        H(v1, v2) = next_label++;
    }

    // Step 3: Compute the adjacency matrix of the direct product graph.
    vector<Eigen::Triplet<double>> v;
    for (auto i = 0; i < e1.rows(); ++i) {
        const auto e1_s = e1(i, 0);
        const auto e1_t = e1(i, 1);
        const auto e1_label = e1(i, 2);

        for (auto j = 0; j < e2.rows(); ++j) {
            if (e1_label == e2(j, 2)) {
                const auto e2_s = e2(j, 0);
                const auto e2_t = e2(j, 1);

                if (v1_label[e1_s] == v2_label[e2_s]
                &&  v1_label[e1_t] == v2_label[e2_t]
                   ) {
                    v.emplace_back(H(e1_s, e2_s), H(e1_t, e2_t), 1.0);
                    v.emplace_back(H(e1_t, e2_t), H(e1_s, e2_s), 1.0);
                }

                if (v1_label[e1_s] == v2_label[e2_t]
                &&  v1_label[e1_t] == v2_label[e2_s]
                   ) {
                    v.emplace_back(H(e1_s, e2_t), H(e1_t, e2_s), 1.0);
                    v.emplace_back(H(e1_t, e2_s), H(e1_s, e2_t), 1.0);
                }
            }
        }
    }

    SparseMatrix<double> Ax(next_label, next_label);
    Ax.setFromTriplets(v.cbegin(), v.cend());

    return Ax;
}

double geometricRandomWalkKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        double lambda,
        int max_iterations,
        double eps) {
    // compute the adjacency matrix Ax of the direct product graph
    const SparseMatrix<double> Lx = lambda * productAdjacency(
            e1, e2, v1_label, v2_label);

    // inverse of I - lambda * Ax by fixed-poInt iterations
    const auto n_rows = Lx.rows();
    const VectorXd ones = VectorXd::Ones(n_rows);
    auto x = ones;
    VectorXd x_pre = VectorXd::Zero(n_rows);

    auto count = 0;
    do {
        x_pre = x;
        x = ones + Lx * x_pre;
        ++count;
    } while (count <= max_iterations && (x - x_pre).squaredNorm() > eps);
    return x.sum();
}

MatrixXd CalculateGeometricRandomWalkKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double lambda,
        int max_iterations,
        double eps) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto j = 0; j < V_label.size(); ++j) {
        for (auto i = 0; i <= j; ++i) {
            K(i, j) = geometricRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], lambda,
                    max_iterations, eps);
        }
    }

    return K.selfadjointView<Eigen::Upper>();
}

double exponentialRandomWalkKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        double beta) {
    // compute the adjacency matrix Ax of the direct product graph
    const MatrixXd Ax = productAdjacency(e1, e2, v1_label, v2_label);

    return Ax.exp().sum();
}

MatrixXd CalculateExponentialRandomWalkKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double beta) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto j = 0; j < V_label.size(); ++j) {
        for (auto i = 0; i <= j; ++i) {
            K(i, j) = exponentialRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], beta);
        }
    }

    return K.selfadjointView<Eigen::Upper>();
}

double kstepRandomWalkKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        const vector<double>& lambda_list) {
    // compute the adjacency matrix Ax of the direct product graph
    const SparseMatrix<double> Ax = productAdjacency(e1, e2, v1_label, v2_label);

    // prepare identity matrix
    const auto n_rows = Ax.rows();
    SparseMatrix<double> I{n_rows, n_rows};
    I.setIdentity();

    auto Sum = SparseMatrix<double>{n_rows, n_rows};
    Sum.setZero();

    // Compute products until k using:
    // https://en.wikipedia.org/wiki/Horner%27s_method
    auto k = lambda_list.size();
    while (k-- > 0) {
        Sum = (Sum * Ax) + lambda_list[k] * I;
    }

    return Sum.sum();
}

MatrixXd CalculateKStepRandomWalkKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        const vector<double>& par) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto j = 0; j < V_label.size(); ++j) {
        for (auto i = 0; i <= j; ++i) {
            K(i, j) = kstepRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], par);
        }
    }

    return K.selfadjointView<Eigen::Upper>();
}
