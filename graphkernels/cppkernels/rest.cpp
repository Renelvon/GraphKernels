/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#include "rest.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/Sparse>

using std::vector;

using Eigen::FullPivLU;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SelfAdjointEigenSolver;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

// map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
int productMapping(
        vector<int>& v1_label,
        vector<int>& v2_label,
        MatrixXi& H) {
    auto n_vx = 0;
    for (auto i = 0; i < v1_label.size(); ++i) {
        for (auto j = 0; j < v2_label.size(); ++j) {
            if (v1_label[i] == v2_label[j]) {
                H(i, j) = n_vx++;
            }
        }
    }
    return n_vx;
}

// compute the adjacency matrix Ax of the direct product graph (sparse)
SparseMatrix<double> productAdjacency(
        MatrixXi& e1,
        MatrixXi& e2,
        vector<int>& v1_label,
        vector<int>& v2_label,
        MatrixXi& H) {
  const auto n_vx = v1_label.size() * v2_label.size();

  SparseMatrix<double> Ax(n_vx, n_vx);

  vector<Eigen::Triplet<double>> v;

  for (auto i = 0L; i < e1.rows(); i++) {
    for (auto j = 0L; j < e2.rows(); j++) {
      if (v1_label[e1(i, 0)] == v2_label[e2(j, 0)] &&
          v1_label[e1(i, 1)] == v2_label[e2(j, 1)] && e1(i, 2) == e2(j, 2)) {
        v.emplace_back(H(e1(i, 0), e2(j, 0)), H(e1(i, 1), e2(j, 1)), 1.0);
        v.emplace_back(H(e1(i, 1), e2(j, 1)), H(e1(i, 0), e2(j, 0)), 1.0);
      }
      if (v1_label[e1(i, 0)] == v2_label[e2(j, 1)] &&
          v1_label[e1(i, 1)] == v2_label[e2(j, 0)] && e1(i, 2) == e2(j, 2)) {
        v.emplace_back(H(e1(i, 0), e2(j, 1)), H(e1(i, 1), e2(j, 0)), 1.0);
        v.emplace_back(H(e1(i, 1), e2(j, 0)), H(e1(i, 0), e2(j, 1)), 1.0);
      }
    }
  }
  Ax.setFromTriplets(v.begin(), v.end());
  return Ax;
}

double geometricRandomWalkKernel(
        MatrixXi& e1,
        MatrixXi& e2,
        vector<int>& v1_label,
        vector<int>& v2_label,
        double lambda,
        int max_iterations,
        double eps) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  const auto n_vx = productMapping(v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax = productAdjacency(e1, e2, v1_label, v2_label, H);

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
    x = I_vec + lambda * Ax * x_pre;
    ++count;
  }
  return x.sum();
}

MatrixXd CalculateGeometricRandomWalkKernelPy(
        vector<MatrixXi>& E,
        vector<vector<int>>& V_label,
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

double exponentialRandomWalkKernel(MatrixXi& e1,
                                   MatrixXi& e2,
                                   vector<int>& v1_label,
                                   vector<int>& v2_label,
                                   double beta) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  const auto n_vx = productMapping(v1_label, v2_label, H);

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax = productAdjacency(e1, e2, v1_label, v2_label, H);

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
        vector<MatrixXi>& E,
        vector<vector<int>>& V_label,
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

double kstepRandomWalkKernel(MatrixXi& e1,
                             MatrixXi& e2,
                             vector<int>& v1_label,
                             vector<int>& v2_label,
                             vector<double>& lambda_list) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  const auto n_vx = productMapping(v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax = productAdjacency(e1, e2, v1_label, v2_label, H);

  auto Sum = SparseMatrix<double>{n_vx, n_vx};
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
        vector<MatrixXi>& E,
        vector<vector<int>>& V_label,
        vector<double>& par) {
    MatrixXd K(V_label.size(), V_label.size());
    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = kstepRandomWalkKernel(
                    E[i], E[j], V_label[i], V_label[j], par);
        }
    }

    return K;
}
