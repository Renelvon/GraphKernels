/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_REST_H_
#define GRAPHKERNELS_CPPKERNELS_REST_H_

#include <Eigen/Core>

#include <vector>

// Random Walk kernels
// ===================

Eigen::MatrixXd CalculateGeometricRandomWalkKernelPy(
        std::vector<Eigen::MatrixXi>& E,
        std::vector<std::vector<int>>& V_label,
        double par,
        int max_iterations,
        double eps);

Eigen::MatrixXd CalculateExponentialRandomWalkKernelPy(
        std::vector<Eigen::MatrixXi>& E,
        std::vector<std::vector<int>>& V_label,
        double par);

Eigen::MatrixXd CalculateKStepRandomWalkKernelPy(
        std::vector<Eigen::MatrixXi>& E,
        std::vector<std::vector<int>>& V_label,
        std::vector<double>& par);

#endif  // GRAPHKERNELS_CPPKERNELS_REST_H_
