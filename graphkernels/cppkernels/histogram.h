/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
#define GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_

#include <Eigen/Core>

#include <vector>

Eigen::MatrixXd CalculateHistogramKernelPy(
    const std::vector<Eigen::MatrixXi>& E,
    const std::vector<std::vector<int>>& V_label,
    double par,
    int kernel_type);

Eigen::MatrixXd CalculateEdgeHistogramKernelPy(
    const std::vector<Eigen::MatrixXi>& E,
    const std::vector<std::vector<int>>& V_label,
    double par);

Eigen::MatrixXd CalculateVertexHistogramKernelPy(
    const std::vector<Eigen::MatrixXi>& E,
    const std::vector<std::vector<int>>& V_label,
    double par);

Eigen::MatrixXd CalculateVertexEdgeHistogramKernelPy(
    const std::vector<Eigen::MatrixXi>& E,
    const std::vector<std::vector<int>>& V_label,
    double par);

#endif  // GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
