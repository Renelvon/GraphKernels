/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#include "histogram.h"

#include <algorithm>
#include <cmath>
#include <tuple>

using std::exp;
using std::max;
using std::minmax;
using std::vector;

using Eigen::MatrixXd;
using Eigen::MatrixXi;

double linear_kernel(const vector<int>& h1, const vector<int>& h2) {
    auto sum = 0.0;
    for (auto i = 0; i < h1.size(); ++i) {
        sum += static_cast<double>(h1[i]) * h2[i];
    }
    return sum;
}

double rbf_kernel(const vector<int>& h1, const vector<int>& h2, double gamma) {
    double sum = 0.0;
    for (auto i = 0; i < h1.size(); ++i) {
        const auto diff = static_cast<double>(h1[i]) - h2[i];
        sum += diff * diff;
    }
    return exp(- gamma * sum);
}

template< class T>
constexpr T max_plus_one(T a, T b) {return 1 + max(a, b);}

double edgeHistogramKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        double gamma) {
    const auto e12 = e1.col(2);
    const auto e22 = e2.col(2);
    const auto label_high = max_plus_one(e12.maxCoeff(), e22.maxCoeff());

    vector<int> h1(label_high, 0);
    vector<int> h2(label_high, 0);

    for (auto i = 0; i < e12.size(); ++i) {
        ++h1[e12(i)];
    }

    for (auto i = 0; i < e22.size(); ++i) {
        ++h2[e22(i)];
    }

    if (gamma > 0.0) {
        return rbf_kernel(h1, h2, gamma);
    }

    return linear_kernel(h1, h2);
}

double vertexHistogramKernel(
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        double gamma) {
    const auto v1_label_max = *max_element(v1_label.cbegin(), v1_label.cend());
    const auto v2_label_max = *max_element(v2_label.cbegin(), v2_label.cend());
    const auto v_label_high = max_plus_one(v1_label_max, v2_label_max);

    vector<int> h1(v_label_high, 0);
    vector<int> h2(v_label_high, 0);

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

double vertexEdgeHistogramKernel(
        const MatrixXi& e1,
        const MatrixXi& e2,
        const vector<int>& v1_label,
        const vector<int>& v2_label,
        double gamma) {
    const auto e12 = e1.col(2);
    const auto e22 = e2.col(2);
    const auto e_label_high = max_plus_one(e12.maxCoeff(), e22.maxCoeff());

    const auto v1_label_max = *max_element(v1_label.cbegin(), v1_label.cend());
    const auto v2_label_max = *max_element(v2_label.cbegin(), v2_label.cend());
    const auto v_label_high = max_plus_one(v1_label_max, v2_label_max);

    const auto v_size = v_label_high * v_label_high * e_label_high;
    vector<int> h1(v_size, 0);
    vector<int> h2(v_size, 0);

    const auto e10 = e1.col(0);
    const auto e11 = e1.col(1);
    for (auto i = 0; i < e10.size(); ++i) {
        int v_min;
        int v_max;
        std::tie(v_min, v_max) = minmax(e10(i), e11(i));

        ++h1[v1_label[v_max]
           + v1_label[v_min] * v_label_high
           + e12(i)          * v_label_high * v_label_high
        ];
    }

    const auto e20 = e2.col(0);
    const auto e21 = e2.col(1);
    for (auto i = 0; i < e20.size(); ++i) {
        int v_min;
        int v_max;
        std::tie(v_min, v_max) = minmax(e20(i), e21(i));

        ++h2[v2_label[v_max]
           + v2_label[v_min] * v_label_high
           + e22(i)          * v_label_high * v_label_high
        ];
    }

    if (gamma > 0.0) {
        return rbf_kernel(h1, h2, gamma);
    }

    return linear_kernel(h1, h2);
}

MatrixXd CalculateEdgeHistogramKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double par) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = edgeHistogramKernel(E[i], E[j], par);
        }
    }

    return K;
}

MatrixXd CalculateVertexHistogramKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double par) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = vertexHistogramKernel(
                V_label[i], V_label[j], par);
        }
    }

    return K;
}

MatrixXd CalculateVertexEdgeHistogramKernelPy(
        const vector<MatrixXi>& E,
        const vector<vector<int>>& V_label,
        double par) {
    MatrixXd K(V_label.size(), V_label.size());

    for (auto i = 0; i < V_label.size(); ++i) {
        for (auto j = i; j < V_label.size(); ++j) {
            K(j, i) = K(i, j) = vertexEdgeHistogramKernel(
                    E[i], E[j], V_label[i], V_label[j], par);
        }
    }

    return K;
}
