/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#include "rest.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/Sparse>

#include <numeric>

using std::iota;

using Eigen::FullPivLU;
using Eigen::SelfAdjointEigenSolver;
using Eigen::SparseMatrix;

using Tuple3_t = Eigen::Triplet<double>;

// map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
int productMapping(vector<int>& v1_label,
                   vector<int>& v2_label,
                   MatrixXi& H) {
  int n_vx = 0;
  for (auto i = 0; i < v1_label.size(); ++i) {
    for (auto j = 0; j < v2_label.size(); ++j) {
      if (v1_label[i] == v2_label[j]) {
        H(i, j) = n_vx;
        n_vx++;
      }
    }
  }
  return n_vx;
}

// compute the adjacency matrix Ax of the direct product graph (sparse)
MatrixXd productAdjacency(MatrixXi& e1,
                          MatrixXi& e2,
                          vector<int>& v1_label,
                          vector<int>& v2_label,
                          MatrixXi& H) {
  auto n_vx = v1_label.size() * v2_label.size();

  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx;

  vector<Tuple3_t> v;

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
  dAx = MatrixXd(Ax);

  return dAx;
}


// geometric random walk karnel
double geometricRandomWalkKernel(MatrixXi& e1,
                                 MatrixXi& e2,
                                 vector<int>& v1_label,
                                 vector<int>& v2_label,
                                 double lambda) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  int n_vx = productMapping(v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx(n_vx, n_vx);

  dAx = productAdjacency(e1, e2, v1_label, v2_label, H);
  Ax = dAx.sparseView();

  // inverse of I - lambda * Ax by fixed-poInt iterations
  VectorXd I_vec(n_vx);
  for (int i = 0; i < n_vx; i++) {
      I_vec[i] = 1;
  }
  VectorXd x = I_vec;
  VectorXd x_pre(n_vx);
  x_pre.setZero();

  double eps = pow(10, -10);
  int count = 0;
  while ((x - x_pre).squaredNorm() > eps) {
    if (count > 100) {
      // cout << "does not converge until " << count - 1 << " iterations" <<
      // endl;
      break;
    }
    x_pre = x;
    x = I_vec + lambda * Ax * x_pre;
    count++;
  }
  return x.sum();
}

// exponential random walk karnel
double exponentialRandomWalkKernel(MatrixXi& e1,
                                   MatrixXi& e2,
                                   vector<int>& v1_label,
                                   vector<int>& v2_label,
                                   double beta) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  int n_vx = productMapping(v1_label, v2_label, H);

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx(n_vx, n_vx);

  dAx = productAdjacency(e1, e2, v1_label, v2_label, H);
  Ax = dAx.sparseView();

  // compute e^{beta * Ax}
  SelfAdjointEigenSolver<MatrixXd> es(Ax);
  VectorXd x = (beta * es.eigenvalues()).array().exp();
  MatrixXd D = x.asDiagonal();
  MatrixXd V = es.eigenvectors();

  MatrixXd I(n_vx, n_vx);
  I.setIdentity();
  FullPivLU<MatrixXd> solver(V);
  MatrixXd V_inv = solver.solve(I);
  MatrixXd Res = V * D * V_inv;

  // compute the total sum
  double K = 0;
  for (auto i = 0L; i < Res.rows(); i++) {
    for (auto j = 0L; j < Res.cols(); j++) {
      K += Res(i, j);
    }
  }

  return K;
}

// k-step product graph karnel
double kstepRandomWalkKernel(MatrixXi& e1,
                             MatrixXi& e2,
                             vector<int>& v1_label,
                             vector<int>& v2_label,
                             vector<double>& lambda_list) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  int n_vx = productMapping(v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx(n_vx, n_vx);

  dAx = productAdjacency(e1, e2, v1_label, v2_label, H);
  Ax = dAx.sparseView();

  // compute products until k
  auto k_max = static_cast<int>(lambda_list.size()) - 1;
  SparseMatrix<double> Ax_pow = I;
  SparseMatrix<double> Sum = lambda_list[0] * I;
  for (int k = 1; k <= k_max; k++) {
    Ax_pow = Ax * Ax_pow;
    Sum += lambda_list[k] * Ax_pow;
  }

  // compute the total sum
  double K = 0;
  for (auto i = 0L; i < Sum.outerSize(); ++i) {
    for (SparseMatrix<double>::InnerIterator it(Sum, i); it; ++it) {
      K += it.value();
    }
  }

  return K;
}

// compute a kernel value of a pair of graphs
double computeKernelValue(MatrixXi& e1,
                          MatrixXi& e2,
                          vector<int>& v1_label,
                          vector<int>& v2_label,
                          vector<double>& par,
                          int kernel_type) {
  double Kval;
  switch (kernel_type) {
    case 8:  // geometric random walk kernel
      Kval = geometricRandomWalkKernel(e1, e2, v1_label, v2_label, par[0]);
      break;
    case 9:  // exponential random walk kernel
      Kval = exponentialRandomWalkKernel(e1, e2, v1_label, v2_label, par[0]);
      break;
    case 10:  // k-step random walk kernel
      Kval = kstepRandomWalkKernel(e1, e2, v1_label, v2_label, par);
      break;
    default:
      Kval = 0;
      break;
  }
  return Kval;
}

MatrixXd CalculateKernelPy(vector<MatrixXi>& E,
                           vector<vector<int>>& V_label,
                           vector<double>& par,
                           int kernel_type) {
  MatrixXd K(V_label.size(), V_label.size());

  vector<int> idx(V_label.size());
  iota(idx.begin(), idx.end(), 0);
  for (auto&& i : idx) {
    for (auto&& j : idx) {
      K(i, j) = computeKernelValue(E[i], E[j], V_label[i], V_label[j], par,
                                   kernel_type);
      K(j, i) = K(i, j);
    }
  }

  return K;
}
