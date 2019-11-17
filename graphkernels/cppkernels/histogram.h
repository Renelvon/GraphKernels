/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
#define GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_

#include <Eigen/Core>

#include <vector>

Eigen::MatrixXd CalculateEdgeHistogramKernelPy(
    const std::vector<Eigen::MatrixXi>& E,
    double par = -1.0);

Eigen::MatrixXd CalculateVertexHistogramKernelPy(
    const std::vector<std::vector<int>>& V_label,
    double par = -1.0);

Eigen::MatrixXd CalculateVertexEdgeHistogramKernelPy(
    const std::vector<Eigen::MatrixXi>& E,
    const std::vector<std::vector<int>>& V_label,
    double par = -1.0);

#endif  // GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
