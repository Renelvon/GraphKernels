/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
#define GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_

#include <Eigen/Core>

#include <vector>

Eigen::MatrixXd CalculateHistogramKernelPy(
    std::vector<Eigen::MatrixXi>& E,
    std::vector<std::vector<int>>& V_label,
    double par,
    int kernel_type);

#endif  // GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
