/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#include "rest.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/Sparse>

#include <utility>

using std::vector;
using std::pair;

using Eigen::FullPivLU;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SelfAdjointEigenSolver;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

// store each valid vertex pair (v_1, v_2) in a vector
auto productMapping(
        const vector<int>& v1_label,
        const vector<int>& v2_label) {
    vector<pair<int, int>> pairs;
    for (auto i = 0; i < v1_label.size(); ++i) {
        for (auto j = 0; j < v2_label.size(); ++j) {
            if (v1_label[i] == v2_label[j]) {
                pairs.emplace_back(i, j);
            }
        }
    }
    return pairs;
}

// compute the adjacency matrix Ax of the direct product graph (sparse)
SparseMatrix<double> productAdjacency(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        const vector<pair<int, int>>& pairs) {
    MatrixXi H = MatrixXi::Zero(v1_label.size(), v2_label.size());

    auto new_label = 0;
    for (const auto& p : pairs) {
        H(p.first, p.second) = new_label++;
    }

    SparseMatrix<double> Ax(new_label, new_label);

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
    const auto pairs = productMapping(v1_label, v2_label);
    const auto n_vx = pairs.size();

    // compute the adjacency matrix Ax of the direct product graph
    SparseMatrix<double> Lx = lambda * productAdjacency(
            e1, e2, v1_label, v2_label, pairs);

    // prepare identity matrix
    SparseMatrix<double> I(n_vx, n_vx);
    I.setIdentity();

    // inverse of I - lambda * Ax by fixed-poInt iterations
    const VectorXd I_vec = VectorXd::Ones(n_vx);
    auto x = I_vec;
    VectorXd x_pre = VectorXd::Zero(n_vx);

    auto count = 0;
    while ((x - x_pre).squaredNorm() > eps) {
        if (count > max_iterations) {
            // cout << "does not converge until " << count - 1 << " iterations" <<
            // endl;
            break;
        }
        x_pre = x;
        x = I_vec + Lx * x_pre;
        ++count;
    }
    return x.sum();
}

MatrixXd CalculateGeometricRandomWalkKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double lambda,
        int max_iterations,
        double eps) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = geometricRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], lambda,
                    max_iterations, eps);
        }
    }

    return K;
}

double exponentialRandomWalkKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        double beta) {
    const auto pairs = productMapping(v1_label, v2_label);
    const auto n_vx = pairs.size();

    // compute the adjacency matrix Ax of the direct product graph
    SparseMatrix<double> Ax = productAdjacency(e1, e2, v1_label, v2_label, pairs);

    // compute e^{beta * Ax}
    SelfAdjointEigenSolver<MatrixXd> es(Ax);
    VectorXd x = (beta * es.eigenvalues()).array().exp();
    MatrixXd D = x.asDiagonal();
    MatrixXd V = es.eigenvectors();

    // prepare identity matrix
    const auto I = MatrixXd::Identity(n_vx, n_vx);

    FullPivLU<MatrixXd> solver(V);
    MatrixXd V_inv = solver.solve(I);
    MatrixXd Res = V * D * V_inv;

    return Res.sum();
}

MatrixXd CalculateExponentialRandomWalkKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double beta) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = exponentialRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], beta);
        }
    }

    return K;
}

double kstepRandomWalkKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        const vector<double>& lambda_list) {
    const auto pairs = productMapping(v1_label, v2_label);
    const auto n_vx = pairs.size();

    // compute the adjacency matrix Ax of the direct product graph
    SparseMatrix<double> Ax = productAdjacency(e1, e2, v1_label, v2_label, pairs);

    // prepare identity matrix
    SparseMatrix<double> I(n_vx, n_vx);
    I.setIdentity();

    auto Sum = SparseMatrix<double>(n_vx, n_vx);
    Sum.setZero();

    // compute products until k using:
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
    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = kstepRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], par);
        }
    }

    return K;
}
