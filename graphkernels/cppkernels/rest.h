/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_REST_H_
#define GRAPHKERNELS_CPPKERNELS_REST_H_

#include <Eigen/Core>

#include <vector>

// Random Walk kernels
// ===================

Eigen::MatrixXd CalculateKernelPy(std::vector<Eigen::MatrixXi>& E,
                                  std::vector<std::vector<int>>& V_label,
                                  std::vector<double>& par,
                                  int kernel_type);

#endif  // GRAPHKERNELS_CPPKERNELS_REST_H_
