/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_WL_H_
#define GRAPHKERNELS_CPPKERNELS_WL_H_

#include <Eigen/Core>

#include <vector>

Eigen::MatrixXd WLKernelMatrix(std::vector<Eigen::MatrixXi>& E,
                               std::vector<std::vector<int>>& V_label,
                               std::vector<int>& num_v,
                               std::vector<int>& num_e,
                               std::vector<int>& degree_max,
                               int h_max);

#endif  // GRAPHKERNELS_CPPKERNELS_WL_H_
