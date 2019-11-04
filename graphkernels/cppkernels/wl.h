/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_WL_H_
#define GRAPHKERNELS_CPPKERNELS_WL_H_

#include <Eigen/Core>

#include <vector>

Eigen::MatrixXd WLKernelMatrix(
        const std::vector<Eigen::MatrixXi>& E,
        const std::vector<std::vector<int>>& V_label,
        const std::vector<int>& num_v,
        const std::vector<int>& num_e,
        const std::vector<int>& degree_max,
        int h_max);

#endif  // GRAPHKERNELS_CPPKERNELS_WL_H_
