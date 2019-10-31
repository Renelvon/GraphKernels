/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#include "histogram.h"

#include <cmath>

using std::exp;
using std::vector;

using Eigen::MatrixXd;
using Eigen::MatrixXi;

double linear_kernel(vector<int>& h1, vector<int>& h2) {
    auto sum = 0.0;
    for (auto i = 0; i < h1.size(); ++i) {
        sum += static_cast<double>(h1[i]) * h2[i];
    }
    return sum;
}

double rbf_kernel(vector<int>& h1, vector<int>& h2, double gamma) {
    double sum = 0.0;
    for (auto i = 0; i < h1.size(); ++i) {
        const auto diff = static_cast<double>(h1[i]) - h2[i];
        sum += diff * diff;
    }
    return exp(- gamma * sum);
}


// edge histogram karnel
double edgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, double gamma) {
  int e_label_max = 0;
  for (auto i = 0L; i < e1.rows(); i++) {
      if (e1(i, 2) > e_label_max) {
          e_label_max = e1(i, 2);
      }
  }
  for (auto i = 0L; i < e2.rows(); i++) {
      if (e2(i, 2) > e_label_max) {
          e_label_max = e2(i, 2);
      }
  }

  vector<int> h1(e_label_max + 1, 0);
  vector<int> h2(e_label_max + 1, 0);

  for (auto i = 0L; i < e1.rows(); i++) {
    (h1[e1(i, 2)])++;
  }
  for (auto i = 0L; i < e2.rows(); i++) {
    (h2[e2(i, 2)])++;
  }

  if (gamma > 0.0) {
      return rbf_kernel(h1, h2, gamma);
  }

  return linear_kernel(h1, h2);
}

// vertex histogram karnel
double vertexHistogramKernel(vector<int>& v1_label,
                             vector<int>& v2_label,
                             double gamma) {
  int v1_label_max = *max_element(v1_label.begin(), v1_label.end());
  int v2_label_max = *max_element(v2_label.begin(), v2_label.end());
  int v_label_max = v1_label_max > v2_label_max ? v1_label_max : v2_label_max;

  vector<int> h1(v_label_max + 1, 0);
  vector<int> h2(v_label_max + 1, 0);

  for (int i : v1_label) {
    ++h1[i];
  }
  for (int i : v2_label) {
    ++h2[i];
  }

  if (gamma > 0.0) {
      return rbf_kernel(h1, h2, gamma);
  }

  return linear_kernel(h1, h2);
}

// vertex-edge histogram karnel
double vertexEdgeHistogramKernel(MatrixXi& e1,
                                 MatrixXi& e2,
                                 vector<int>& v1_label,
                                 vector<int>& v2_label,
                                 double gamma) {
  int e_label_max = 0;
  for (auto i = 0L; i < e1.rows(); i++) {
      if (e1(i, 2) > e_label_max) {
          e_label_max = e1(i, 2);
      }
  }
  for (auto i = 0L; i < e2.rows(); i++) {
      if (e2(i, 2) > e_label_max) {
          e_label_max = e2(i, 2);
      }
  }
  e_label_max++;

  int v1_label_max = *max_element(v1_label.begin(), v1_label.end());
  int v2_label_max = *max_element(v2_label.begin(), v2_label.end());
  int v_label_max = v1_label_max > v2_label_max ? v1_label_max : v2_label_max;
  v_label_max++;

  vector<int> h1(v_label_max * v_label_max * e_label_max, 0);
  vector<int> h2(v_label_max * v_label_max * e_label_max, 0);

  int v1;
  int v2;
  for (auto i = 0L; i < e1.rows(); i++) {
    v1 = e1(i, 0);
    v2 = e1(i, 1);
    if (v2 > v1) {
      int v_tmp = v1;
      v1 = v2;
      v2 = v_tmp;
    }
    (h1[v1_label[v1] + v1_label[v2] * v_label_max +
        e1(i, 2) * v_label_max * v_label_max])++;
  }
  for (auto i = 0L; i < e2.rows(); i++) {
    v1 = e2(i, 0);
    v2 = e2(i, 1);
    if (v2 > v1) {
      int v_tmp = v1;
      v1 = v2;
      v2 = v_tmp;
    }
    (h2[v2_label[v1] + v2_label[v2] * v_label_max +
        e2(i, 2) * v_label_max * v_label_max])++;
  }

  if (gamma > 0.0) {
      return rbf_kernel(h1, h2, gamma);
  }

  return linear_kernel(h1, h2);
}

// vertex-vertex-edge histogram karnel
double vertexVertexEdgeHistogramKernel(MatrixXi& e1,
                                       MatrixXi& e2,
                                       vector<int>& v1_label,
                                       vector<int>& v2_label,
                                       double lambda) {
  return vertexHistogramKernel(v1_label, v2_label, -1.0) +
         lambda * vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, -1.0);
}

MatrixXd CalculateHistogramKernelPy(
        vector<MatrixXi>& E,
        vector<vector<int>>& V_label,
        double par,
        int kernel_type) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            auto Kval = 0.0;
            switch (kernel_type) {
                // Simple kernels
                case 1:
                    Kval = edgeHistogramKernel(E[i], E[j], -1.0);
                    break;
                case 2:
                    Kval = vertexHistogramKernel(V_label[i], V_label[j], -1.0);
                    break;
                case 3:
                    Kval = vertexEdgeHistogramKernel(
                            E[i], E[j], V_label[i], V_label[j], -1.0);
                    break;
                case 4:
                    Kval = vertexVertexEdgeHistogramKernel(
                            E[i], E[j], V_label[i], V_label[j], par);
                    break;
                // Gaussian kernels
                case 5:
                    Kval = edgeHistogramKernel(E[i], E[j], par);
                    break;
                case 6:
                    Kval = vertexHistogramKernel(
                            V_label[i], V_label[j], par);
                    break;
                case 7:
                    Kval = vertexEdgeHistogramKernel(
                            E[i], E[j], V_label[i], V_label[j], par);
                    break;
                default:
                    Kval = 42.0;  // FIXME: THIS SHOULD NEVER HAPPEN!
            }
            K(j, i) = K(i, j) = Kval;
        }
    }

    return K;
}
